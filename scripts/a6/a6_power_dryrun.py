#!/usr/bin/env python3
"""A6.8.2 step 2-3: a6_power re-run with the DRY-RUN yields (pilot condition 6) and the throughput tag criterion.

LABEL-FREE: inputs are the dry run's label-free counts (A, A_cand, tail proxy per DEV run) and metadata.
reviews/checks/a6_power.py is imported and run UNCHANGED, except that its yield input file is the dry-run yields
(written in a6_cert_yield.json's format to results/2026-10-28_A6_hn1/power/) and its output path is redirected there.
The one-engine 'cert_hs' sensitivity cell keeps the drafting-time counts (reviews/checks/a6_cert_yield.json); it is a
sensitivity cell only and is not used by the tag rule.

Tag rule (A6.8.2 step 3): N = r x (16:30 - tag time).  Tag iff P(v1 minimum within N runs) >= 0.8 at x1 yields AND
>= 0.5 at x0.5 yields, with the same 400 walks per cell and seeds as a6_power (walk seeds 9001+i, rng seed 131).

    scripts/m2b/vpy.sh scripts/a6/a6_power_dryrun.py --r RUNS_PER_HOUR [--tag-time 15:15]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_A6_hn1")
PW = os.path.join(RES, "power")
sys.path.insert(0, os.path.join(ROOT, "reviews", "checks"))


def _counts(path):
    """{acc: {run: counts}} from a dry_run.json or from an event log (JSON lines; 'admit' events)."""
    if path.endswith(".json"):
        return {acc: rec["counts"] for acc, rec in json.load(open(path))["projects"].items()}
    out = {}
    for x in open(path):
        if x.startswith("{"):
            e = json.loads(x)
            if e.get("ev") == "admit":
                out[e["acc"]] = e["counts"]
    return out


def write_yields(dry_paths, out_path):
    ref = json.load(open(os.path.join(ROOT, "reviews", "checks", "a6_cert_yield.json")))["per_run"]
    per = {}
    allc = {}
    for pth in dry_paths:
        allc.update(_counts(pth))
    for acc, counts in allc.items():
        for run, c in counts.items():
            per[run] = {"project": acc, "run": run,
                        "cert_agree": {"n_certified": c["A"], "A_cand": c["A_cand"], "tail_proxy_own_thr": c["tail_proxy"]},
                        "cert_hs": ref[run]["cert_hs"] if run in ref else
                        {"n_certified": c["n_certified_hs"], "A_cand": c["A_cand"], "tail_proxy_own_thr": c["tail_proxy"]}}
    tot = {k: int(sum(v["cert_agree"][k] for v in per.values())) for k in ("n_certified", "A_cand", "tail_proxy_own_thr")}
    json.dump({"source": list(dry_paths), "per_run": per, "totals_cert_agree": tot}, open(out_path, "w"), indent=1)
    return tot, len(per)


def p_min_grid(AP6, accs, keys, usable, y, budgets, nwalk=400):
    rng = np.random.default_rng(AP6.SEED)
    ws = [AP6.walk(rng, 9001 + i, accs, keys, usable, y) for i in range(nwalk)]
    return {int(b): float(np.mean([w["snap_min"] is not None and w["snap_min"]["runs"] <= b for w in ws]))
            for b in budgets}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", nargs="+", default=[os.path.join(RES, "dry_run", "DEV12_events_part1.log"),
                                                 os.path.join(RES, "dry_run", "DEV12_part2", "dry_run.json")])
    ap.add_argument("--r", type=float, required=True, help="measured end-to-end admitted runs per hour")
    ap.add_argument("--tag-time", default="15:15")
    ap.add_argument("--walk-end", default="16:30")
    ap.add_argument("--skip-full", action="store_true", help="only the walk grid (no gate re-simulation)")
    a = ap.parse_args()
    os.makedirs(PW, exist_ok=True)
    t0 = time.perf_counter()
    tot, nrun = write_yields(a.dry, os.path.join(PW, "a6_cert_yield.json"))
    import a6_power as AP6
    AP6.HERE = PW                                   # yields() reads HERE/a6_cert_yield.json -> the dry-run yields
    AP6.OUT = os.path.join(PW, "a6_power_dryrun.json")
    hh = lambda s: int(s.split(":")[0]) + int(s.split(":")[1]) / 60  # noqa: E731
    window_h = hh(a.walk_end) - hh(a.tag_time)
    N = a.r * window_h
    accs, keys, usable, pinfo = AP6.pool()
    y = AP6.yields("cert_agree")
    budgets = sorted(set(list(range(20, 401, 10)) + [int(math.floor(N))]))
    grid = {}
    for ys in (1.0, 0.5):
        yy = dict(y, log_A_mean_project=y["log_A_mean_project"] + math.log(ys))
        grid[f"x{ys}"] = p_min_grid(AP6, accs, keys, usable, yy, budgets)
    nN = int(math.floor(N))
    p1, p05 = grid["x1.0"][nN], grid["x0.5"][nN]
    decision = "GO-to-tag" if (p1 >= 0.8 and p05 >= 0.5) else "not tagged: infeasible"
    res = {"dry_run_yields_totals": tot, "n_dev_runs": nrun, "yields_fit": y, "pool": pinfo, "r_runs_per_h": a.r,
           "tag_time": a.tag_time, "walk_end": a.walk_end, "window_h": window_h, "N_runs": N,
           "P_v1_min_within_N": {"x1": p1, "x0.5": p05}, "rule": "P(x1) >= 0.8 and P(x0.5) >= 0.5",
           "decision": decision, "P_v1_min_grid": grid,
           "r_needed": {"x1_0.8": next((b for b in budgets if grid["x1.0"][b] >= 0.8), None),
                        "x0.5_0.5": next((b for b in budgets if grid["x0.5"][b] >= 0.5), None)}}
    need = max(x for x in res["r_needed"].values() if x is not None) if all(res["r_needed"].values()) else None
    res["min_runs_needed_both"] = need
    res["min_r_needed_runs_per_h"] = need / window_h if need else None
    json.dump(res, open(os.path.join(PW, "tag_criterion.json"), "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k not in ("P_v1_min_grid", "yields_fit", "pool")}, default=float),
          flush=True)
    if not a.skip_full:
        AP6.main()                                  # full a6_power re-run (walk cells + gate behaviour) on dry-run yields
    print("wall", round(time.perf_counter() - t0), "s")


if __name__ == "__main__":
    main()
