#!/usr/bin/env python3
"""Collect results/2026-10-28_E7b_fair_speed/metrics.json from the driver outputs, equality records, benches and
cohort stats.  Missing inputs are recorded as null."""
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_E7b_fair_speed")
OURS = "vulcan_v2_sep"
LEVELS = ("peptidoform", "sequence", "protein")


def jl(p):
    p = os.path.join(RES, p)
    return json.load(open(p)) if os.path.exists(p) else None


def comparator(ds):
    m = jl(f"t1/{ds}/metrics.json")
    if m is None:
        return None
    out = {"runs": m["runs"], "methods": {}, "ratios_v2_sep_over": {}}
    for meth, rec in m["methods"].items():
        if "ids_dataset_level" not in rec:
            out["methods"][meth] = rec
            continue
        r = {"ids": {lv: {a: round(v["mean"], 1) for a, v in rec["ids_dataset_level"][lv].items()} for lv in LEVELS},
             "entrapment": {}}
        for lv in LEVELS:
            r["entrapment"][lv] = {}
            for a, v in rec["verdicts"][lv].items():
                e = {"deciding": v["deciding"], "P_mean_fdp": v["P"]["a"]["mean_fdp"], "P_ci95_a": v["P"]["a"]["ci95"],
                     "P_verdict_a": v["P"]["a"]["verdict"]}
                if "C" in v:
                    e.update(C_mean_fdp=v["C"]["a"]["mean_fdp"], C_ci95_a=v["C"]["a"]["ci95"],
                             C_verdict_a=v["C"]["a"]["verdict"])
                r["entrapment"][lv][a] = e
        out["methods"][meth] = r
    for comp, rec in m.get("noninferiority", {}).items():
        out["ratios_v2_sep_over"][comp] = {
            lv: {"own_q": {a: {k: v[k] for k in ("ratio", "ci95", "ours_ids", "comparator_ids",
                                                 "comparator_verdict_deciding", "counts_for_noninferiority")}
                           for a, v in d["own_q"].items()},
                 "matched_fdp": {f: {k: v[k] for k in ("ratio", "ci95", "ours_ids", "comparator_ids",
                                                       "ours_q_at_match", "comparator_q_at_match")}
                                 for f, v in d["matched_fdp"].items()}}
            for lv, d in rec.items()}
    return out


def r_rule(comp):
    """Pre-registered: R=5 also on D8h iff sage_chim5+pred+mokapot > sage_chim2+pred+mokapot at matched FDP 1%
    (peptidoform, point) on BOTH D0 and D1."""
    res = {}
    for ds in ("D0", "D1"):
        c = (comp.get(ds) or {}).get("ratios_v2_sep_over", {})
        try:
            n2 = c["sage_chim2+pred+mokapot"]["peptidoform"]["matched_fdp"]["0.01"]["comparator_ids"]
            n5 = c["sage_chim5+pred+mokapot"]["peptidoform"]["matched_fdp"]["0.01"]["comparator_ids"]
            res[ds] = {"R2": n2, "R5": n5, "R5_more": n5 > n2}
        except KeyError:
            res[ds] = None
    ok = all(v is not None for v in res.values())
    res["run_R5_on_D8h"] = bool(ok and all(v["R5_more"] for k, v in res.items() if k in ("D0", "D1")))
    return res


def bench_summary(b):
    if b is None:
        return None
    return {"median_over_13_runs": b["summary_median_over_runs"], "sum_over_13_runs": b["summary_sum_over_runs"],
            "threads": b["threads"], "loadavg_start": b["loadavg_start"], "loadavg_end": b["loadavg_end"],
            "per_dataset_median": {ds: {k: {m: float(sorted(r["median"][k][m] for r in d["runs"].values())[
                len(d["runs"]) // 2]) for m in ("wall_s", "core_s")} for k in ("ms1", "search", "features", "model",
                                                                              "compete", "proteins", "total")}
                                    for ds, d in b["datasets"].items()}}


def cohort(cfg):
    c = jl(f"cohort_{cfg}.json")
    if c is None:
        return None
    keep = {k: c.get(k) for k in ("n_runs", "workers", "makespan_s", "runs_per_hour", "runs_per_hour_steady",
                                  "bottleneck", "stage_busy_fraction", "core_s", "gpu_s", "ms1_stage", "chimeric")}
    if "scoring" in c:
        s = c["scoring"]
        keep["scoring"] = {k: s.get(k) for k in ("n_scored", "gpu_stage_wall_s_median", "post_wall_s_median",
                                                 "post_core_s_sum", "ids_1pct_sum", "proteins_1pct_sum")}
    keep["hw"] = c.get("hw")
    return keep


def main():
    comp = {ds: comparator(ds) for ds in ("D0", "D1", "D8h")}
    eq_s, eq_m = jl("equality_scoring.json"), jl("equality_ms1.json")
    out = {
        "experiment": "E7b fair chimeric comparator + v2_sep speed (dev only)",
        "preregistration": "config.yaml (commit 7a2d6ee, before any Sage chimera search / comparator evaluation)",
        "comparator": comp,
        "R_rule": r_rule(comp),
        "sage_rank_strata": {ds: jl(f"sage_rank_strata_{ds}.json") for ds in ("D0", "D1", "D8h")},
        "speed": {
            "bench_before": bench_summary(jl("bench_before.json")),
            "bench_after": bench_summary(jl("bench_after.json")),
            "cohorts": {c: cohort(c) for c in ("chim", "chim_ms1ref", "closed", "base12")},
        },
        "equivalence": {
            "ms1": None if eq_m is None else {"all_ok": eq_m["all_ok"], "n_runs": eq_m["n_runs"]},
            "scoring": None if eq_s is None else {
                ds: {k: v for k, v in r.items() if k.startswith("all_") or k == "lib_misses_total"}
                for ds, r in eq_s["datasets"].items()},
            "rt_history_check": jl("check_rt_history_E9_0008.json"),
            "history_check_all_D8h": (lambda j: None if j is None else {k: v for k, v in j.items() if k != "pairs"})(
                jl("check_history_D8h.json")),
            "cohort_reference": (lambda j: None if j is None else {"all_ok": j["all_ok"], "n_runs": len(j["runs"])})(
                jl("equality_cohort_reference.json")),
            "cohort_copies": {c: (lambda j: None if j is None else {k: j[k] for k in ("all_ok", "n_runs", "n_bases")})(
                jl(f"equality_cohort_copies_{c}.json")) for c in ("chim", "chim_ms1ref")},
        },
    }
    json.dump(out, open(os.path.join(RES, "metrics.json"), "w"), indent=1, default=float)
    print(json.dumps(out["R_rule"]))


if __name__ == "__main__":
    sys.exit(main())
