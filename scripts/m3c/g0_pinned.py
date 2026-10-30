#!/usr/bin/env python3
"""3c follow-up, item 1: G0 (label blindness) for the frozen IT-CID scorers with every prediction read from the PINNED
prediction store (vulcan.predict.pinned), on the same two cases as 3c's G0 (scripts/m3c/gates_itcid.py g0):
  LF_D10cf  D10cf fraction 10 vs the E2 Pfu E2001 index and E2's swap index (FASTA shuffled, pairs permuted, T/D
            roles swapped)
  TMT_D4    D4 fraction A1 vs TMT-static versions of the same two indexes (built by 3c, same seeds)
Procedure: (1) search both sides; (2) fill pass: the store is filled with every (model, peptide, charge, CE) key the
two sides need (Koina queried at most once per key; label-free); (3) frozen pass: for EACH side a NEW PinnedStore
instance in mode "frozen" (no Koina access possible; shards re-read from disk and sha256-verified) computes every L3
feature and all 4 frozen scores; PSM tables, features and scores must be bit-identical (strict, as 3c).
    python scripts/m3c/g0_pinned.py  -> results/2026-10-27_3c_itcid/followup/g0_pinned.json
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import gates_itcid as G  # noqa: E402
import pinned_cache as PC  # noqa: E402

import numpy as np  # noqa: E402


def main():
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.predict.pinned import PinnedPredCache, PinnedStore
    from vulcan.rescore.itcid import ItcidFeatureBuilder
    ms = G.models()
    out = {"store": PC.PINNED}
    t_all = time.perf_counter()
    for tag, store, tmt in (("LF_D10cf", C.d10_stores()[10], False), ("TMT_D4", C.d4_stores()[0], True)):
        if tmt:
            base = os.path.join(C.WORK, "index", "g0_pfu_tmt_E2001")
            alt = os.path.join(C.WORK, "index", "g0_pfu_tmt_swapped")
            st, var = C.mods("human_tmt")
        else:
            base = "/tmp/vulcan_e2/index/pfu_contam_E2001"
            alt = "/tmp/vulcan_e2/g0/index_swapped"
            st, var = C.mods("human")
        kind = "tmt" if tmt else "lf"
        p = E.SearchParams(**C.SEARCH)
        sp = E.load_store(store)
        sides = []
        for idx in (base, alt):
            l1 = L1Index(idx)
            if tmt:
                l1.mass32 = (l1.mass32.astype(np.float32) + np.float32(C.TMT)).astype(np.float32)
                l1.mask = l1.mask | np.uint64(1 << 63)
            gi = E.GpuIndex(l1, p, st, var)
            tbl, _, _ = E.search(sp, gi, l1, p)
            del gi
            torch.cuda.empty_cache()
            sides.append((l1, tbl, G.psm_dict(tbl)))
        # (2) fill pass (label-free; one Koina query per missing key)
        fill = PinnedStore(PC.PINNED, "fill")
        for l1, _, psm in sides:
            fb = ItcidFeatureBuilder(lambda nce: PinnedPredCache(fill, kind, st, var, ce=nce), st, var)
            fb.run(l1, psm, sp)
        snap = PinnedStore(PC.PINNED, "frozen").snapshot_id()
        # (3) frozen pass, a fresh store instance per side
        res = []
        for l1, tbl, psm in sides:
            fr = PinnedStore(PC.PINNED, "frozen")
            fb = ItcidFeatureBuilder(lambda nce, fr=fr: PinnedPredCache(fr, kind, st, var, ce=nce), st, var)
            feat, _ = fb.run(l1, psm, sp)
            res.append((tbl, feat, {s: m.apply(psm, feat) for s, m in ms.items()}))
        (t0, f0, s0), (t1, f1, s1) = res
        r = {"store_spectra": store, "n_psms": t0.num_rows, "snapshot_id": snap,
             "koina_keys_filled_for_this_case": fill.stats["filled"],
             "psm_tables_bit_identical": bool(t0.equals(t1)),
             "features_bit_identical": {k: bool(np.array_equal(f0[k], f1[k], equal_nan=True)) for k in f0},
             "features_max_abs_diff": {k: float(np.nanmax(np.abs(f0[k].astype(float) - f1[k].astype(float))))
                                       if f0[k].size else 0.0 for k in f0},
             "static_scores_bit_identical": {k: bool(np.array_equal(s0[k], s1[k])) for k in s0}}
        r["pass"] = r["psm_tables_bit_identical"] and all(r["features_bit_identical"].values()) and \
            all(r["static_scores_bit_identical"].values())
        out[tag] = r
        print(tag, "G0 pinned", r["pass"], flush=True)
    out["pass"] = all(v["pass"] for v in out.values() if isinstance(v, dict))
    out["wall_s"] = time.perf_counter() - t_all
    os.makedirs(PC.FU_RES, exist_ok=True)
    json.dump(out, open(os.path.join(PC.FU_RES, "g0_pinned.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
