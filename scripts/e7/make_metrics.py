#!/usr/bin/env python3
"""Collect results/2026-10-28_E7_chimeric_v2/metrics.json from the driver outputs (t1/<ds>/metrics.json), the stratum
read-outs (strata_<ds>_<method>.json), gates_v2.json, runtime.json and models/training_log.json."""
from __future__ import annotations

import json
import os
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2")
NUMS = ["vulcan_v2_sep", "vulcan_v2", "vulcan_chim_v1_sep", "vulcan_chim_v1"]


def load(p):
    p = os.path.join(RES, p)
    return json.load(open(p)) if os.path.exists(p) else None


def ratios(m, num):
    ni = (m.get("noninferiority_by_numerator") or {}).get(num) or (m["noninferiority"] if m.get("noninferiority_numerator") == num else None)
    if not ni:
        return None
    out = {}
    for comp, lv in ni.items():
        out[comp] = {}
        for level, d in lv.items():
            out[comp][level] = {
                "matched_fdp_0.01": {k: d["matched_fdp"]["0.01"][k] for k in ("ratio", "ci95", "ours_ids", "comparator_ids", "lower_ci_ge_0.97")},
                "matched_fdp_0.05": {k: d["matched_fdp"]["0.05"][k] for k in ("ratio", "ci95")},
                "own_q_0.01": {k: d["own_q"]["0.01"][k] for k in ("ratio", "ci95", "comparator_verdict_deciding", "counts_for_noninferiority")}}
    return out


def ids_verdicts(m):
    out = {}
    for meth, r in m["methods"].items():
        if "ids_dataset_level" not in r:
            out[meth] = r
            continue
        out[meth] = {lv: {"ids_q0.01_mean": r["ids_dataset_level"][lv]["0.01"]["mean"],
                          "ids_q0.001_mean": r["ids_dataset_level"][lv]["0.001"]["mean"],
                          "verdict_0.01": {k: r["verdicts"][lv]["0.01"][k] for k in ("deciding", "deciding_conservative")},
                          "P_mean_fdp_0.01": r["verdicts"][lv]["0.01"]["P"]["a"]["mean_fdp"],
                          "P_ci95_0.01_a": r["verdicts"][lv]["0.01"]["P"]["a"]["ci95"],
                          "C_mean_fdp_0.01": (r["verdicts"][lv]["0.01"].get("C") or {}).get("a", {}).get("mean_fdp")}
                     for lv in r["ids_dataset_level"]}
    return out


def strata_summary(s):
    if not s:
        return None
    out = {}
    for k, v in s["by_stratum"].items():
        out[k] = {"ids_prod": v["ids_prod_dataset_level_mean_over_decoy_seeds"],
                  "decoy_estimated_fdr_in_stratum": {al: x["D_over_T"] for al, x in v["decoy_estimated_fdr_in_stratum_prod"].items()},
                  "verdicts": {al: {"P_mean_fdp": x["P"]["a"]["mean_fdp"], "P_ci95": x["P"]["a"]["ci95"],
                                    "P_verdict_a": x["P"]["a"]["verdict"], "P_verdict_b": x["P"]["b"]["verdict"],
                                    "P_lb_ci95": x["P"]["a"]["lb_ci95"], "C_mean_fdp": x["C"]["a"]["mean_fdp"],
                                    "C_verdict_a": x["C"]["a"]["verdict"], "overall_a": x["overall_a"], "overall_b": x["overall_b"],
                                    "mean_n_t_P": x["P"]["a"]["mean_n_t"]}
                               for al, x in v["verdicts"].items()}}
    return {"q_at_matched_fdp": s["q_at_matched_fdp"], "closed_v1_ids_prod": s["closed_v1_ids_prod_mean"], "strata": out,
            "neighbour_filter_prod_1001": s["competition_diag_prod"].get("prod_1001")}


def main():
    out = {"experiment": "E7 chimeric search + scorer v2 (dev only; D0/D1 selection set, D8h confirmation)",
           "git_sha": subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
           "datasets": {}}
    for ds in ("D0", "D1", "D8h"):
        m = load(f"t1/{ds}/metrics.json")
        if not m:
            continue
        out["datasets"][ds] = {"runs": m["runs"], "ids_and_verdicts": ids_verdicts(m),
                               "ratios": {n: ratios(m, n) for n in NUMS},
                               "strata": {meth: strata_summary(load(f"strata_{ds}_{meth}.json"))
                                          for meth in ("vulcan_chim_v1", "vulcan_chim_v1_sep", "vulcan_v2", "vulcan_v2_sep")
                                          if load(f"strata_{ds}_{meth}.json")}}
    for ds in ("D0", "D1", "D8h"):
        for v in ("e7", "e7sep"):
            m = load(f"t1/{ds}_{v}/metrics.json")
            if m:
                out["datasets"].setdefault(ds, {}).setdefault("e7_runs", {})[v] = {"numerator": m.get("noninferiority_numerator"),
                                                                                  "ratios_vs_closed_v1": ratios(m, m.get("noninferiority_numerator"))}
    g = load("gates_v2.json")
    if g:
        x = g["G2_original_G5"]
        out["gates"] = {"G0": {k: g["G0"][k] for k in ("pass", "chimeric_psm_tables_bit_identical", "columns_not_identical",
                                                       "arrow_equals_raw", "nce_same", "scores_bit_identical", "n_psms",
                                                       "n_secondary_psms")},
                        "G2_original": {s: {"pass": v["g2"]["pass"], "overall": v["g2"]["strata"]["0"]["overall"],
                                            "top": v["g2"]["strata"]["0"]["top"], "by_subset": v["g2_by_subset"],
                                            "frac_winner_changed": v["frac_winner_changed"]}
                                        for s, v in x["scorers"].items()},
                        "G5": {s: v["g5"] for s, v in x["scorers"].items()},
                        "G5_fixed_reference": x["scorers"]["v2"]["g5_fixed_reference"],
                        "null_design": {k: x[k] for k in ("n_null_competitions", "n_offset", "n_species_pfu",
                                                          "n_species_yeast", "n_secondary_units", "p_target_win_fixed")},
                        "pending": g.get("pending")}
    rt = load("runtime.json")
    if rt:
        out["runtime_D0_per_run_median"] = {r: v["median"] for r, v in rt["runs"].items()}
        out["runtime_meta"] = {k: rt[k] for k in ("hw_threads", "loadavg_start", "loadavg_end")}
    tl = load("models/training_log.json")
    if tl:
        out["training"] = {k: tl[k] for k in ("n_runs", "n_psms", "excluded_runs", "chosen", "best_cv", "within_tolerance",
                                              "hyperscore_ids_all_runs")}
        out["training"]["cv"] = {k: v["heldout_ids"] for k, v in tl["cv"].items()}
    json.dump(out, open(os.path.join(RES, "metrics.json"), "w"), indent=1, default=float)
    print("wrote", os.path.join(RES, "metrics.json"))


if __name__ == "__main__":
    main()
