#!/usr/bin/env python3
"""Evaluate the Sage / Sage+mokapot comparator runs of 3c (sage_comp.py) with the M1e analysis code
(scripts/baselines/analyze.py: tool's own peptide q, I/L-collapsed sequences; picked protein groups; arm P (3 seeds)
and arm C entrapment with the §6.3 rule), pointed at the 3c output directory.

    python scripts/m3c/sage_eval.py  -> results/2026-10-27_3c_itcid/sage_eval.json
"""
from __future__ import annotations

import functools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import sage_comp as SC  # noqa: E402  (registers D10cf in run_searches.DATASETS)

sys.path.insert(0, os.path.join(C.ROOT, "scripts", "baselines"))
import analyze as A  # noqa: E402

A.load = functools.partial(A.load, cache=SC.CACHE)
PIPES = ["sage", "sage+mokapot"]


def main():
    out = {}
    for ds in ("D4", "D10cf"):
        for pipe in PIPES:
            res = A.load(pipe, ds, "prod")
            if res is None:
                continue
            P = A.eval_P(pipe, ds)
            Cc = A.eval_C(pipe, ds)
            rec = {"ids": A.id_counts(res, with_tdc=True), "proteins": A.protein_ids_production(res),
                   "verdicts": A.verdicts_for(P, Cc), "arm_P": P and P["alpha"], "arm_C": Cc and Cc["alpha"],
                   "units": "peptides = I/L-collapsed stripped sequences at the tool's own peptide q (M1e units)"}
            try:
                rec["protein_entrapment"] = A.eval_protein(ds, pipe)
            except Exception as e:  # noqa: BLE001
                rec["protein_entrapment"] = f"failed: {e}"
            out.setdefault(ds, {})[pipe] = rec
            v = rec["verdicts"]["0.01"]
            print(ds, pipe, {a: rec["ids"][a]["peptides"] for a in ("0.001", "0.01", "0.05")}, rec["proteins"]["0.01"],
                  v, flush=True)
    json.dump(out, open(os.path.join(C.RES, "sage_eval.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
