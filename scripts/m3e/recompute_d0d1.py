#!/usr/bin/env python3
"""A1.6: recompute the M3 D0/D1 entrapment verdicts under both CI conventions ((a) bootstrap + between-seed,
(b) bootstrap only), with the M3 evaluation code path unchanged otherwise (scripts/m3/evaluate.py: same searches,
features, scorers, estimators, n_boot 2000, seed 11).  Arm C has one replicate, so (a) == (b) there.  Peptide level
(arm P, 3 entrapment seeds) and protein level (arm P).  The verdict under the convention chosen by
scripts/m3e/ci_convention_sim.py is the A1.6 verdict; "consistent with control" is reported under both but never ranks.
    conda_envs/pred/bin/python scripts/m3e/recompute_d0d1.py --out results/.../d0d1_ci_conventions.json
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
import evaluate as EV  # noqa: E402
import m3_common as C  # noqa: E402

from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402

KEEP = ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t", "ci95_a", "ci95_b", "var_boot",
        "var_between_over_J", "fdp_per_replicate")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    scorers = ["hyperscore"] + EV.STATIC + list(EV.RESETS)
    models = {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in EV.STATIC}
    out = {}
    for ds in ("D0", "D1"):
        searches = {k: EV.Search(ds, k) for k in ("prod", "E2001", "E2002", "E2003", "C")}
        r_c = C.c_arm_r(EV.FAM[ds])
        for sc in scorers:
            mclass = "fixed" if sc == "hyperscore" else "static_low_capacity"
            ws = {k: EV.scored(s, sc, models) for k, s in searches.items()}
            rec = {}
            if sc not in EV.RESETS:
                pe = {k: EV.protein_level(searches[k], ws[k], EV.ALPHAS, entrapment=True) for k in ("E2001", "E2002", "E2003")}
            for al in EV.ALPHAS:
                reps = [ENT.Replicate(EV.peptide_entrapment(searches[k], ws[k], al, "paired"), int(k[1:]), 0)
                        for k in ("E2001", "E2002", "E2003")]
                armC = ENT.evaluate_arm([ENT.Replicate(EV.peptide_entrapment(searches["C"], ws["C"], al, "combined"), -1, 0)],
                                        al, "combined", r=r_c)
                r = {"C": {k: armC[k] for k in KEEP if k in armC}}
                for cv in ("a", "b"):
                    armP = ENT.evaluate_arm(reps, al, "paired", ci_convention=cv)
                    d = {k: armP[k] for k in KEEP if k in armP}
                    d["consistent_with_control"] = bool(armP["ci95"][0] <= al)
                    d["overall_P_and_C"] = ENT.combine_arms({"P": armP["verdict"], "C": armC["verdict"]}, mclass)
                    r[f"P_{cv}"] = d
                    if sc not in EV.RESETS:
                        preps = [ENT.Replicate(pe[k][al], int(k[1:]), 0) for k in ("E2001", "E2002", "E2003")]
                        pp = ENT.evaluate_arm(preps, al, "paired", ci_convention=cv)
                        r[f"protein_P_{cv}"] = {k: pp[k] for k in KEEP if k in pp}
                rec[str(al)] = r
                print(ds, sc, al, "P(a)", r["P_a"]["verdict"], round(r["P_a"]["ci95"][0], 5), "P(b)", r["P_b"]["verdict"],
                      round(r["P_b"]["ci95"][0], 5), "C", r["C"]["verdict"], flush=True)
            out.setdefault(ds, {})[sc] = rec
        del searches
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
