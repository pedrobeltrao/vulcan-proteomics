#!/usr/bin/env python3
"""Koina repeatability diagnostic (for G0): predict the same peptidoforms twice, in the same order and batches, with
the Prosit TMT-CID and CID models; count peptidoforms whose returned intensities differ (server-side float
nondeterminism, independent of any label).  -> results/2026-10-27_3c_itcid/koina_repeat.json"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import numpy as np  # noqa: E402
from vulcan.predict.koina import KoinaProsit  # noqa: E402

out = {}
for fam in ("human_tmt", "human"):
    l1 = C.load_l1(fam, "prod")
    rows = np.random.default_rng(11).choice(len(l1), 20000, replace=False)
    st, var = C.mods(fam)
    kp = KoinaProsit("tmt" if C.is_tmt(fam) else "lf", st, var)
    res = [kp.predict(l1.residues, l1.seq_off[rows], l1.length[rows], l1.mask[rows].view(np.int64),
                      np.full(len(rows), 2), 35.0) for _ in range(2)]
    a, b = res
    fo = a["ms2_off"]
    diff = [not np.array_equal(a["ms2"][fo[i]:fo[i + 1]], b["ms2"][fo[i]:fo[i + 1]]) for i in range(len(rows))]
    out[fam] = {"n": len(rows), "n_ms2_differ": int(sum(diff)),
                "max_abs_diff": float(np.abs(a["ms2"].astype(float) - b["ms2"].astype(float)).max()),
                "rt_identical": bool(np.array_equal(a["rt"], b["rt"], equal_nan=True))}
    print(fam, out[fam], flush=True)
json.dump(out, open(os.path.join(C.RES, "koina_repeat.json"), "w"), indent=1)
