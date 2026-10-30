#!/usr/bin/env python3
"""D0 mass-offset / open search IDs with CONGA filtering + stem-level competition + Group-walk (DESIGN L4), ours vs
Sage 0.14.7 open search through the SAME pipeline (Sage's PSMs mapped to our index rows; same score = hyperscore).

    python scripts/m4/open_d0.py --ours closed=DIR offset=DIR open=DIR --sage-open DIR --out RES/open_d0.json

Each input dir holds <run>.psms.parquet with the top `report` PSMs per spectrum (we use top 5).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

from vulcan.gpu import sage_bridge as SB  # noqa: E402
from vulcan.index.build import L1Index, _read  # noqa: E402
from vulcan.ptm import conga as CG  # noqa: E402
from vulcan.ptm.localize import RESIDUE  # noqa: E402

TOPN = 5


def ion_set_cache(l1, rows, bin_da=0.02):
    cache = {}
    out = []
    for r in rows:
        r = int(r)
        if r not in cache:
            s = l1.sequence(r)
            m = int(l1.mask[r])
            rm = np.array([RESIDUE[c] + C.STATIC.get(c, 0.0) + (C.VAR_CLOSED.get(c, 0.0) if m >> k & 1 else 0.0)
                           for k, c in enumerate(s)])
            cache[r] = CG.by_ion_bins(s, rm, 0.0, bin_da)
        out.append(cache[r])
    return out


def run_conga(spec_key, rows, score, delta, lab, l1, label):
    t0 = time.perf_counter()
    ions = ion_set_cache(l1, rows)
    keep = CG.neighbour_filter(spec_key, score, ions, 0.05)
    t_nf = time.perf_counter() - t0
    rows_k = rows[keep]
    stem = np.array([l1.sequence(int(r)).replace("I", "L") for r in rows_k])
    pair = np.array([s.replace("I", "L") for s in lab["fwd_seq"][rows_k]])
    dec = lab["is_decoy"][rows_k]
    psm_id = np.flatnonzero(keep)
    res = CG.conga(psm_id, stem, pair, score[keep], dec, delta[keep])
    w = res["row"]
    out = {k: v for k, v in res.items() if k not in ("row", "group", "q", "q_pooled", "is_decoy", "score")}
    acc = (~res["is_decoy"]) & (res["q"] <= 0.01)
    shifts = np.asarray(delta[keep])[w][acc]
    b = CG.shift_bin(shifts)
    vals, cnt = np.unique(b[b != 0], return_counts=True)
    out.update(label=label, n_psms_in=int(len(rows)), n_after_neighbour_filter=int(keep.sum()),
               accepted_stems_1pct_shifted=int((b != 0).sum()), accepted_stems_1pct_narrow=int((b == 0).sum()),
               n_shift_bins_1Da_with_ge1=int(len(vals)), n_shift_bins_1Da_with_ge5=int((cnt >= 5).sum()),
               top_shift_bins=[(int(v), int(c)) for v, c in sorted(zip(vals, cnt), key=lambda x: -x[1])[:15]],
               time_neighbour_filter_s=t_nf, time_total_s=time.perf_counter() - t0)
    return out


def ours(d, l1, lab, label):
    t = pa.concat_tables([pq.read_table(f) for f in sorted(glob.glob(os.path.join(d, "*.psms.parquet")))])
    t = t.filter(pc.less_equal(t.column("rank"), TOPN))
    sk = np.array([f"{a}:{b}" for a, b in zip(t.column("run").to_pylist(), t.column("spec_idx").to_pylist())])
    rows = t.column("pf_row").to_numpy()
    if "mass_shift" in t.column_names:
        delta = t.column("mass_shift").to_numpy().astype(np.float64)
    else:
        delta = (t.column("expmass").to_numpy() - t.column("isotope_error").to_numpy() * 1.00335
                 - t.column("calcmass").to_numpy()).astype(np.float64)
    return run_conga(sk, rows, t.column("hyperscore").to_numpy(), delta, lab, l1, label)


def sage(d, l1, lab, key):
    s = SB.read_results(d).to_pandas()
    s = s[s["rank"] <= TOPN]
    pm = [SB.parse_peptide(x, C.VAR_CLOSED) for x in s.peptide]
    rows = np.array([key[k] for k in pm])
    assert np.array_equal(lab["is_decoy"][rows], (s.label.values == -1)), "label mapping mismatch"
    sk = (s.filename + ":" + s.scannr.astype(str)).values
    delta = (s.expmass.values - s.calcmass.values).astype(np.float64)
    out = run_conga(sk, rows, s.hyperscore.values.astype(np.float64), delta, lab, l1, "sage_open_same_pipeline")
    r1 = s[s["rank"] == 1]
    out["sage_own"] = {"peptides_q01_rank1_targets": int(r1[(r1.label == 1) & (r1.peptide_q <= 0.01)].peptide.nunique()),
                       "psms_spectrum_q01": int(((r1.label == 1) & (r1.spectrum_q <= 0.01)).sum()),
                       "note": "Sage's own LDA + its q-values (no CONGA, no stem competition): reference only"}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", nargs="+", required=True, help="label=DIR")
    ap.add_argument("--sage-open", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    l1 = L1Index(C.D0_INDEX)
    labt = _read(os.path.join(C.D0_INDEX, "labels.arrow"))
    lab = {"is_decoy": labt.column("is_decoy").to_numpy(zero_copy_only=False),
           "fwd_seq": np.array(labt.column("fwd_seq").to_pylist(), dtype=object)}
    res = {}
    for x in a.ours:
        k, d = x.split("=", 1)
        res[f"ours_{k}"] = ours(d, l1, lab, f"ours_{k}")
        print(k, json.dumps({kk: v for kk, v in res[f"ours_{k}"].items() if kk != "group_pair_counts"}), flush=True)
    if a.sage_open:
        key = {}
        s = SB.read_results(a.sage_open).to_pandas()
        need = {SB.parse_peptide(x, C.VAR_CLOSED) for x in s[s["rank"] <= TOPN].peptide}
        for i in range(len(l1)):
            k = (l1.sequence(i), int(l1.mask[i]))
            if k in need:
                key[k] = i
        res["sage_open"] = sage(a.sage_open, l1, lab, key)
        print("sage", json.dumps({kk: v for kk, v in res["sage_open"].items() if kk != "group_pair_counts"}))
    C.dump(res, a.out)


if __name__ == "__main__":
    main()
