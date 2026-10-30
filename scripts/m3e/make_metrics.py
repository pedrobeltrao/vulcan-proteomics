#!/usr/bin/env python3
"""Assemble results/2026-10-27_M3e_hardening/metrics.json from the component outputs (gates_a1.json, g3.json,
g4.json, armM.json, ci_convention_sim.json, d0d1_ci_conventions.json, d8_confirm.json, bench_prod.json)."""
import json
import os

import numpy as np

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "2026-10-27_M3e_hardening")


def j(n):
    p = os.path.join(R, n)
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    m = {"frozen_scorer": {"file": "results/2026-10-27_M3_scoring/models/static_gbm_pred.pkl"}}
    meta = json.load(open(os.path.join(R, "..", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl.json")))
    m["frozen_scorer"]["sha256"] = meta["sha256"]
    sim = j("ci_convention_sim.json")
    conv = sim["chosen"] if sim else None
    m["ci_convention"] = None if sim is None else {
        "chosen": conv, "eligible": sim["eligible"], "mean_abs_dev_from_095": sim["mean_abs_dev_from_095"],
        "coverage": {k: {"a": v["coverage_a"], "b": v["coverage_b"], "R": v["R"]} for k, v in sim["scenarios"].items()}}
    g = j("gates_a1.json")
    g3, g4 = j("g3.json"), j("g4.json")
    gates = {}
    if g:
        m["nulls"] = g["null_counts"]
        for p, lab in (("gbm", "with_RT"), ("gbm_noRT", "without_RT")):
            L = g["learned"][p]
            t = L["tail_G2_1pct"]
            gates[lab] = {
                "G1": {"pass": L["G1"]["pass"], "rho_top_decile": L["G1"]["rho"], "ci99": L["G1"]["ci"], "n": L["G1"]["n"]},
                "G2": {"pass": L["G2"]["pass"], "overall_diff": L["G2"]["strata"]["0"]["overall"]["diff"] if "0" in L["G2"]["strata"] else
                       list(L["G2"]["strata"].values())[0]["overall"]["diff"],
                       "overall_ucb": list(L["G2"]["strata"].values())[0]["overall"]["ucb"],
                       "top_ucb": list(L["G2"]["strata"].values())[0]["top"]["ucb"], "frac_winner_changed": L["frac_winner_changed"]},
                "tail_G2": {"verdict": t["verdict"], "n_tail": t["n_tail"], "frac_target": t["frac_target"],
                            "frac_target_ucb": t["frac_target_ucb"], "diff_vs_hyperscore": t["diff"], "diff_ucb": t["diff_ucb"],
                            "n_tail_hyperscore": t["n_tail_fixed"], "frac_target_hyperscore": t["frac_target_fixed"]},
                "G5": {"pass": L["G5"]["pass"], "auc": L["G5"]["all"]["auc"], "ci95": L["G5"]["all"].get("ci"),
                       "auc_top_decile": L["G5"]["top_decile"]["auc"]},
            }
    if g3:
        for lab in ("with_RT", "without_RT"):
            gates.setdefault(lab, {})["G3"] = {"pass": g3["pass"], "max_abs_diff_all": max(abs(b["diff"]) for b in g3["all_nulls"]["bins"]),
                                               "max_abs_diff_clean": max(abs(b["diff"]) for b in g3["clean_species_nulls"]["bins"]),
                                               "note": g3["note"]}
    if g4:
        for p, lab in (("gbm", "with_RT"), ("gbm_noRT", "without_RT")):
            gates.setdefault(lab, {})["G4"] = {"pass": g4[p]["pass"], "tau": g4[p]["tau"],
                                               "ks_vs_rest": {k: v["ks_vs_rest"] for k, v in g4[p]["strata"].items()}}
    m["gates"] = gates
    a = j("armM.json")
    if a:
        m["arm_M"] = {sc: {al: {"verdict_" + cv: e["verdict"] for cv, e in r["entrapment"][al].items()} |
                                {"mean_fdp_C": r["entrapment"][al]["a"]["mean_fdp"], "mean_lb": r["entrapment"][al]["a"]["mean_lb"],
                                 "ci95_a": r["entrapment"][al]["a"]["ci95"], "ci95_b": r["entrapment"][al]["b"]["ci95"],
                                 "N_T_mean": float(np.mean([v["N_T"] for v in r["ids_target_side"][al].values()])),
                                 "N_E_mean": float(np.mean([v["N_E"] for v in r["ids_target_side"][al].values()]))}
                           for al in ("0.01", "0.05")} for sc, r in a["scorers"].items()}
    d = j("d0d1_ci_conventions.json")
    if d:
        m["d0d1_recomputed"] = {ds: {sc: {al: {"P_a": v["P_a"]["verdict"], "P_b": v["P_b"]["verdict"], "C": v["C"]["verdict"],
                                               "P_mean": v["P_a"]["mean_fdp"], "P_lo_a": v["P_a"]["ci95"][0], "P_lo_b": v["P_b"]["ci95"][0],
                                               "overall_a": v["P_a"]["overall_P_and_C"], "overall_b": v["P_b"]["overall_P_and_C"]}
                                          for al, v in r.items() if al != "0.001"} for sc, r in dd.items()} for ds, dd in d.items()}
    d8 = j("d8_confirm.json")
    if d8:
        runs = d8["runs"]
        s = {}
        for al in ("0.001", "0.01", "0.05"):
            h = np.array([d8["per_run"][r]["hyperscore"]["ids"][al] for r in runs])
            gg = np.array([d8["per_run"][r]["static_gbm_pred"]["ids"][al] for r in runs])
            ok = h > 0
            rat = gg[ok] / h[ok]
            rng = np.random.default_rng(11)
            bs = [np.exp(np.mean(np.log(rng.choice(rat, rat.size)))) for _ in range(10000)]
            s[al] = {"hyperscore_peptidoforms": int(h.sum()), "static_gbm_pred_peptidoforms": int(gg.sum()),
                     "hyperscore_proteins": int(sum(d8["per_run"][r]["hyperscore"]["proteins"][al] for r in runs)),
                     "static_gbm_pred_proteins": int(sum(d8["per_run"][r]["static_gbm_pred"]["proteins"][al] for r in runs)),
                     "per_run_ratio_geomean": float(np.exp(np.mean(np.log(rat)))),
                     "per_run_ratio_ci95_run_bootstrap": np.percentile(bs, [2.5, 97.5]).tolist(),
                     "n_runs_ratio": int(ok.sum()), "n_runs_gbm_gt_hyperscore": int((gg > h).sum())}
            for sc in ("hyperscore", "static_gbm_pred"):
                p = d8["pooled"][sc][al]
                s[al][f"{sc}_pooled"] = {"P_mean": p["P"]["a"]["mean_fdp"], "P_verdict_a": p["P"]["a"]["verdict"],
                                         "P_verdict_b": p["P"]["b"]["verdict"], "P_ci95_a": p["P"]["a"]["ci95"],
                                         "P_ci95_b": p["P"]["b"]["ci95"], "C_mean": p["C"]["mean_fdp"],
                                         "C_verdict": p["C"]["verdict"], "LB": p["P"]["a"]["mean_lb"]}
        m["d8_confirmation"] = {"runs": runs, "summary": s,
                                "per_run_ids": {r: {sc: d8["per_run"][r][sc]["ids"] for sc in ("hyperscore", "static_gbm_pred")}
                                                for r in runs}}
    b = j("bench_prod.json")
    if b:
        m["speed_production"] = {"summary": b["summary_measured_runs"], "n_runs": b["n_measured"], "one_time": b["one_time"],
                                 "nce_per_run": {k: v["nce"] for k, v in b["runs"].items()}, "loadavg_end": b["loadavg_end"]}
    json.dump(m, open(os.path.join(R, "metrics.json"), "w"), indent=1, default=float)
    print(json.dumps(m.get("gates"), indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
