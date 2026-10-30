#!/usr/bin/env python3
"""M6d: precursor-level epsilon / CV of M6 vs M6d on the SAME precursors (those with >= 1 value per condition in the M6
version), so that added precursors (charge siblings) do not change the denominator.  Inputs: --save-long tables of
scripts/m6d/lfq_eval.py.     scripts/m2b/vpy.sh scripts/m6d/same_set.py DATASET SUBSET OUT.json"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as DS  # noqa: E402
import evalq as E  # noqa: E402

import pandas as pd  # noqa: E402

ds, sub, out = sys.argv[1:4]
cfg = DS.get(ds)
runs, A, B = cfg.subset(sub)
res = {}
for mbr in ("off", "on"):
    Ls = {v: pd.read_parquet(f"/tmp/vulcan_m6d/long/{ds}_{sub}_{v}_{mbr}.parquet") for v in ("m6", "m6d")}
    mats = {v: E.matrix(L, "precursor") for v, L in Ls.items()}
    M6, _ = mats["m6"]
    M6 = M6.reindex(columns=A + B)
    keys = sorted(M6.index[(M6[A].notna().sum(axis=1) >= 1) & (M6[B].notna().sum(axis=1) >= 1)])
    r = {}
    for v, (M, P) in mats.items():
        o, _ = E.score(M, P, A, B, cfg.expected, rows=keys)
        r[v] = o
        # MS2-linked precursors only (sibling values dropped) - full depth
        L2 = Ls[v][Ls[v]["provenance"] != "ms2_sibling"]
        M2, P2 = E.matrix(L2, "precursor")
        r[v + "_no_sibling_values"], _ = E.score(M2, P2, A, B, cfg.expected)
    res[f"mbr_{mbr}"] = {"n_keys": len(keys), **r}
    print(mbr, {v: (round(x["eps"], 4), round(x.get("cv", x.get("cv2")), 4), x["n"]) for v, x in r.items()})
json.dump(res, open(out, "w"), indent=1, default=float)
