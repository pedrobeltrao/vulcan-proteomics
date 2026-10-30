#!/usr/bin/env python3
"""Collect the C9 pilot outputs (label-free only) into results/2026-10-28_C9_pilot/metrics.json."""
from __future__ import annotations

import glob
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_C9_pilot")


def main():
    runs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(RES, "per_run", "*.json")))
            if not os.path.basename(f).startswith("_")]
    meta = json.load(open(os.path.join(RES, "per_run", "_meta.json")))
    keys = ("n_ms2", "n_certified", "certified_single_ms1_feature", "A_cand", "A_cand_single_ms1_feature",
            "tail_proxy_own_thr", "tail_proxy_D0max_thr", "tail_proxy_own_thr_single_ms1_feature")
    tot = {k: int(sum(r[k] for r in runs)) for k in keys}
    for t in ("0.15", "0.1", "0.05"):
        for k in ("A_cand", "tail_proxy_own_thr"):
            tot[f"ladder_f<{t}_{k}"] = int(sum(r["diag_ladder_threshold"][t][k] for r in runs))
    wa = [sum(r["diag_winner_f_hist_A_cand"][i] for r in runs) for i in range(4)]
    wt = [sum(r["diag_winner_f_hist_tail"][i] for r in runs) for i in range(4)]
    co = {k: sum(r["diag_coisolation"][k] for r in runs) for k in runs[0]["diag_coisolation"]}
    per_run = []
    for r in runs:
        t = r["timing_wall_cpu_s"]
        per_run.append({"project": r["project"], "run": r["run"], **{k: r[k] for k in keys},
                        "frac_certified_with_candidate": r["frac_certified_with_candidate"],
                        "host_thr_own_run": r["host_thr_own_run"],
                        "A_cand_f<0.05": r["diag_ladder_threshold"]["0.05"]["A_cand"],
                        "tail_proxy_own_thr_f<0.05": r["diag_ladder_threshold"]["0.05"]["tail_proxy_own_thr"],
                        "wall_s": round(sum(v[0] for v in t.values()), 1),
                        "process_cpu_s": round(sum(v[1] for v in t.values()), 1),
                        "stage_wall_cpu_s": t})
    n = len(runs)
    cost = {"runs": n, "mean_wall_s_per_run": sum(p["wall_s"] for p in per_run) / n,
            "mean_process_cpu_s_per_run": sum(p["process_cpu_s"] for p in per_run) / n,
            "median_wall_s_per_run": sorted(p["wall_s"] for p in per_run)[n // 2],
            "stage_mean_wall_s": {k: sum(r["timing_wall_cpu_s"].get(k, [0, 0])[0] for r in runs) / n
                                  for k in runs[0]["timing_wall_cpu_s"]},
            "setup": meta, "note": "single process, warm GPU index (setup once); download not included; "
                                   "process CPU excludes GPU kernel time; shared node, load average in setup"}
    out = {"set": "C9 removed-peptide human null pilot (DEV only: D0 PXD001468 x2, D8h 10 runs / 10 projects)",
           "label_discipline": "no per-label quantity of any null candidate computed; full-index labels read only "
                               "for certification, host thresholds and the identified set",
           "totals": tot,
           "rates": {"certified_per_ms2": tot["n_certified"] / tot["n_ms2"],
                     "A_cand_per_certified": tot["A_cand"] / tot["n_certified"],
                     "tail_proxy_own_per_A_cand": tot["tail_proxy_own_thr"] / tot["A_cand"],
                     "tail_proxy_D0max_per_A_cand": tot["tail_proxy_D0max_thr"] / tot["A_cand"],
                     "tail_proxy_own_per_A_cand_f<0.05": tot["ladder_f<0.05_tail_proxy_own_thr"] /
                     tot["ladder_f<0.05_A_cand"],
                     "DN2_reference_tail_proxy_per_A_cand": 0.000895 / 0.391,
                     "per_run_mean": {k: tot[k] / n for k in keys}},
           "diag_winner_ladder_f": {"bins": [0, 0.05, 0.10, 0.15, 0.20], "A_cand_winners": wa, "tail_winners": wt,
                                    "tail_rate_by_bin": [wt[i] / wa[i] for i in range(4)]},
           "diag_coisolation": co | {"tail_rate_single": co["single_tail"] / co["single_A_cand"],
                                     "tail_rate_multi": co["multi_tail"] / co["multi_A_cand"]},
           "per_run": per_run, "cost": cost}
    for name in ("pool_count.json", "certify_strict.json"):
        f = os.path.join(RES, name)
        if os.path.exists(f):
            d = json.load(open(f))
            d.pop("e9_usable_raw_per_project", None) if False else None
            if name == "pool_count.json":
                d["e9_usable_raw_per_project"] = {k: v for k, v in d["e9_usable_raw_per_project"].items()
                                                  if k != "values"}
            out[name[:-5]] = d
    f = os.path.join(RES, "c9_power.json")
    if os.path.exists(f):
        d = json.load(open(f))
        out["power"] = {"summary": d["summary"], "gates_at_v1_minimum": d["gates_at_v1_minimum"],
                        "walk_cells": d["walk_cells"], "criteria": d["criteria"], "seed": d["seed"]}
    json.dump(out, open(os.path.join(RES, "metrics.json"), "w"), indent=1, default=float)
    print(json.dumps({"totals": tot, "rates": out["rates"], "winner_f": out["diag_winner_ladder_f"],
                      "coiso": out["diag_coisolation"], "cost": {k: v for k, v in cost.items() if k != "setup"}},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
