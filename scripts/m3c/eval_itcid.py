#!/usr/bin/env python3
"""3c dev evaluation (D4, D10cf): IDs at q <= {0.1, 1, 5}% (peptidoforms, paired double competition pooled over the
fractions of a set; plus I/L-collapsed unique sequences = M1e/Sage units) and picked protein groups, with §6.3 dev
entrapment verdicts (arm P seeds 2001-2003, arm C Arabidopsis; CI convention (a)) for hyperscore and the frozen IT-CID
static scorers.  Same estimator code as M3 (scripts/m3/evaluate.py functions, with the 3c work dirs).

    python scripts/m3c/eval_itcid.py [--datasets D4 D10cf] [--scorers hyperscore itcid_gbm_pred ...]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m3"))
import common as C  # noqa: E402
import evaluate as M3E  # noqa: E402  (M3 estimator wiring)

import numpy as np  # noqa: E402

from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402

M3E.C = C                                           # index_dir(family, kind) / psm_dir(tag) of 3c
M3E.FAM = {k: v["family"] for k, v in C.DEV.items()}
STATIC = ["itcid_linear_base", "itcid_gbm_base", "itcid_linear_pred", "itcid_gbm_pred"]
M3E.STATIC = STATIC
ALPHAS = (0.001, 0.01, 0.05)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["D4", "D10cf"])
    ap.add_argument("--scorers", nargs="+", default=["hyperscore"] + STATIC)
    ap.add_argument("--models", default=os.path.join(C.RES, "models"))
    ap.add_argument("--out", default=os.path.join(C.RES, "eval.json"))
    a = ap.parse_args()
    models = {s: StaticModel.load(os.path.join(a.models, s + ".pkl")) for s in a.scorers if s in STATIC}
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    r_c = C.c_arm_r()
    for ds in a.datasets:
        t0 = time.perf_counter()
        searches = {k: M3E.Search(ds, k) for k in C.KINDS}
        print(ds, "loaded", round(time.perf_counter() - t0, 1), "s", flush=True)
        res_ds = out.setdefault(ds, {})
        for sc in a.scorers:
            t1 = time.perf_counter()
            mclass = "fixed" if sc == "hyperscore" else "static_low_capacity"
            ws = {k: M3E.scored(s, sc, models) for k, s in searches.items()}
            w = ws["prod"]
            sp = searches["prod"]
            rec = {"model_class": mclass, "ids": {}, "ids_unique_seq_IL": {}, "proteins": {}, "entrapment": {},
                   "protein_entrapment": {}}
            tseq = sp.lab["target_seq"][sp.psm["pf_row"][w["row"]]]
            tseq_il = np.array([s.replace("I", "L") for s in tseq], dtype=object)
            for al in ALPHAS:
                acc = (w["q"] <= al) & ~w["is_decoy"]
                rec["ids"][str(al)] = int(acc.sum())
                rec["ids_unique_seq_IL"][str(al)] = int(np.unique(tseq_il[acc].astype(str)).size)
            pl = M3E.protein_level(sp, w, ALPHAS)
            rec["proteins"] = {str(al): v["n_groups"] for al, v in pl.items()}
            pe = {k: M3E.protein_level(searches[k], ws[k], ALPHAS, entrapment=True) for k in ("E2001", "E2002", "E2003")}
            for al in ALPHAS:
                reps = [ENT.Replicate(M3E.peptide_entrapment(searches[k], ws[k], al, "paired"), int(k[1:]), 0)
                        for k in ("E2001", "E2002", "E2003")]
                armP = ENT.evaluate_arm(reps, al, "paired")
                armC = ENT.evaluate_arm([ENT.Replicate(M3E.peptide_entrapment(searches["C"], ws["C"], al, "combined"),
                                                       -1, 0)], al, "combined", r=r_c)
                armP["consistent_with_control"] = bool(armP["ci95"][0] <= al)
                vs = {"P": armP["verdict"], "C": armC["verdict"]}
                rec["entrapment"][str(al)] = {"P": armP, "C": armC, "overall": ENT.combine_arms(vs, mclass)}
                preps = [ENT.Replicate(pe[k][al], int(k[1:]), 0) for k in ("E2001", "E2002", "E2003")]
                rec["protein_entrapment"][str(al)] = ENT.evaluate_arm(preps, al, "paired")
            rec["eval_wall_s"] = time.perf_counter() - t1
            res_ds[sc] = rec
            e1 = rec["entrapment"]["0.01"]
            print(ds, sc, rec["ids"], rec["ids_unique_seq_IL"]["0.01"], rec["proteins"], "P",
                  round(e1["P"]["mean_fdp"], 4), e1["P"]["verdict"], "C", round(e1["C"]["mean_fdp"], 4),
                  e1["C"]["verdict"], "->", e1["overall"], flush=True)
            json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
