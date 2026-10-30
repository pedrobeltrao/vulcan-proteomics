#!/usr/bin/env python3
"""G3 (DESIGN §6.4, A1.4): retrieval symmetry of the K1 top-K pruning stage (production K = 500) on null spectra,
binned by score.  For the first two G2 mixtures (one per D0 host; scripts/m3e/g2mix.py) the null-provenance spectra
(offset + species) are searched exhaustively (K2 over every candidate in the precursor windows,
vulcan.gpu.engine.search_exhaustive) and with the production top-K search; per spectrum the exhaustive best target
and best decoy (hyperscore over all K2-rescored candidates, >= min_matched_peaks) and the pruned search's best target /
decoy.  Recall_T = P(pruned best target == exhaustive best target), same for decoys, in deciles of the exhaustive
best score (vulcan.gates.stats.g3_gate: |recall_T - recall_D| <= 0.01 in every bin).  Reported for all null spectra
and for the A1.3-clean species subset (exhaustive best target/decoy pairs not near the source proteome).  The
learned scorer does not enter K1 / top-K, so G3 is identical with and without the RT features.
    python scripts/m3e/g3.py --out results/.../g3.json
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
sys.path.insert(0, os.path.join(HERE, "..", "m2"))
sys.path.insert(0, HERE)
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

from vulcan.gates import stats as GS  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
MIX = os.path.join(W, "null_mix")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--mixtures", nargs="+", default=["g2mix00_h0", "g2mix01_h1"])
    a = ap.parse_args()
    import torch
    from e2_correctness import _best_per_label, _rescored
    from gates_a1 import homolog_filters
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    import m3_data as MD
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    lab = MD.labels(idx)
    isdec = lab["is_decoy"]
    p0 = E.SearchParams()
    gi = E.GpuIndex(l1, p0, C.STATIC, C.VAR)
    hfs = homolog_filters()
    acc = {k: [] for k in ("bt_x", "bd_x", "bt_k", "bd_k", "prov", "clean")}
    for name in a.mixtures:
        sp = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        n = len(sp)
        nul = np.flatnonzero(prov >= 1)
        _, _, dx = E.search_exhaustive(sp, gi, l1, p0, spec_subset=nul)
        sx, hx, rx = _rescored(dx, l1, p0)
        btx, bdx = _best_per_label(sx, hx, isdec[rx], n)
        _, _, dk = E.search(sp, gi, l1, p0)
        sk, hk, rk = _rescored(dk, l1, p0)
        btk, bdk = _best_per_label(sk, hk, isdec[rk], n)
        # rows of the exhaustive best target / decoy (for the homolog filter)
        clean = np.ones(n, bool)
        for lab_dec, best in ((False, btx), (True, bdx)):
            m = isdec[rx] == lab_dec
            hit = m & (hx == best[sx])
            brow = np.full(n, -1, np.int64)
            brow[sx[hit]] = rx[hit]
            for src, hf in hfs.items():
                sel = np.flatnonzero((prov == src) & (brow >= 0))
                if sel.size == 0:
                    continue
                rr = np.unique(np.r_[brow[sel], l1.partner[brow[sel]]])
                rr = rr[rr >= 0]
                bad = rr[hf.near([l1.sequence(int(r)) for r in rr])]
                clean[sel[np.isin(brow[sel], bad) | np.isin(l1.partner[brow[sel]], bad)]] = False
        for k, v in (("bt_x", btx), ("bd_x", bdx), ("bt_k", btk), ("bd_k", bdk), ("prov", prov), ("clean", clean)):
            acc[k].append(v[nul])
        print(name, "nulls", nul.size, "with exhaustive candidates", int(np.isfinite(np.maximum(btx, bdx))[nul].sum()), flush=True)
        torch.cuda.empty_cache()
    A = {k: np.concatenate(v) for k, v in acc.items()}
    out = {"K": p0.topk, "mixtures": a.mixtures,
           "all_nulls": GS.g3_gate(A["bt_x"], A["bd_x"], A["bt_k"], A["bd_k"]),
           "clean_species_nulls": GS.g3_gate(*(np.where((A["prov"] >= 2) & A["clean"], A[k], -np.inf)
                                               for k in ("bt_x", "bd_x", "bt_k", "bd_k"))),
           "offset_nulls_with_candidates": int(np.isfinite(np.maximum(A["bt_x"], A["bd_x"])[A["prov"] == 1]).sum())}
    out["pass"] = out["all_nulls"]["pass"] and out["clean_species_nulls"]["pass"]
    out["note"] = "scorer-independent (K1/top-K precede L3): identical with and without RT features"
    for k in ("all_nulls", "clean_species_nulls"):
        print(k, out[k]["pass"], [(round(b["recall_target"], 4), round(b["recall_decoy"], 4)) for b in out[k]["bins"]])
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
