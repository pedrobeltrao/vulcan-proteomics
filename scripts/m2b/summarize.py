#!/usr/bin/env python3
"""Aggregate M2b cohort runs (results/2026-10-27_M2b_cohort/runs/*.json) + side measurements into metrics.json."""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "2026-10-27_M2b_cohort")


def stage_table(d: dict, lustre_GBps: float | None) -> dict:
    """Per-stage capacity in runs/hour if that stage alone were the limit (mean per-run cost, stage parallelism)."""
    runs = d["runs"]
    n = len(runs)
    w = d["workers"]
    ing = np.array([r["wall_s"] for r in runs])
    ms1 = np.array([r.get("ms1_wall_s", 0.0) for r in runs])
    core = np.array([r["core_s"] for r in runs])
    ms1c = np.array([r.get("ms1_core_s", 0.0) for r in runs])
    load = np.array([r["load_t"][1] - r["load_t"][0] for r in runs])
    srch = np.array([r["search_wall_s"] for r in runs])
    gms = np.array([sum(r["gpu_ms"].values()) for r in runs]) / 1e3
    wr = np.array([r["write_t"][1] - r["write_t"][0] for r in runs])
    rawb = np.array([r["raw_bytes"] for r in runs])
    t = {
        "decode_workers (wall incl. RAW copy)": {"parallel": w, "mean_s_per_run": float(ing.mean()),
                                                 "core_s_per_run": float(core.mean()),
                                                 "capacity_runs_per_h": float(w * 3600 / ing.mean())},
        "loader threads (store -> MS2 arrays)": {"parallel": 2, "mean_s_per_run": float(load.mean()),
                                                 "capacity_runs_per_h": float(2 * 3600 / load.mean())},
        "GPU stage (K0..features, serial)": {"parallel": 1, "mean_s_per_run": float(srch.mean()),
                                             "capacity_runs_per_h": float(3600 / srch.mean())},
        "GPU busy (CUDA events)": {"parallel": 1, "mean_s_per_run": float(gms.mean()),
                                   "capacity_runs_per_h": float(3600 / gms.mean())},
        "writer threads (parquet)": {"parallel": 2, "mean_s_per_run": float(wr.mean()),
                                     "capacity_runs_per_h": float(2 * 3600 / max(wr.mean(), 1e-9))},
    }
    if ms1.sum() > 0:
        t["MS1 stage (in decode workers)"] = {"parallel": w, "mean_s_per_run": float(ms1.mean()),
                                              "core_s_per_run": float(ms1c.mean()),
                                              "capacity_runs_per_h": float(w * 3600 / ms1.mean())}
        t["decode_workers (wall incl. RAW copy)"]["capacity_runs_per_h"] = float(w * 3600 / (ing + ms1).mean())
    if lustre_GBps:
        t["Lustre RAW read (measured O_DIRECT aggregate)"] = {
            "GBps": lustre_GBps, "mean_raw_GB_per_run": float(rawb.mean() / 1e9),
            "capacity_runs_per_h": float(lustre_GBps * 3600 / (rawb.mean() / 1e9))}
    b = min(t, key=lambda k: t[k]["capacity_runs_per_h"])
    return {"stages": t, "binding_stage": b, "n_runs": n}


def main():
    lustre = None
    lf = os.path.join(RES, "lustre_read.json")
    if os.path.exists(lf):
        lustre = json.load(open(lf)).get("aggregate_GBps_best")
    out = {"configs": {}, "scaling_c50": [], "notes": "dev data only; repeated files are page-cache warm (see NOTES)"}
    for f in sorted(glob.glob(os.path.join(RES, "runs", "*.json"))):
        d = json.load(open(f))
        tag = os.path.basename(f)[:-5]
        keep = {k: d[k] for k in ("n_runs", "workers", "makespan_s", "runs_per_hour", "runs_per_hour_steady",
                                  "ms2_per_s", "n_ms2", "time_to_first_result_s", "setup", "core_s", "gpu_s",
                                  "host_busy_cores_mean", "stage_busy_fraction", "bottleneck", "per_run_median",
                                  "per_run_sum", "ms1_stage", "stage_raw", "family_top2", "gpu_max_mem_bytes")
                if k in d}
        if "ms1_share_of_worker_busy" in d:
            keep["ms1_share_of_worker_busy"] = d["ms1_share_of_worker_busy"]
        keep["process"] = d.get("process")
        keep["stage_capacity"] = stage_table(d, lustre)
        out["configs"][tag] = keep
        if tag.startswith("c50_"):
            out["scaling_c50"].append({"workers": d["workers"], "makespan_s": d["makespan_s"],
                                       "runs_per_hour": d["runs_per_hour"],
                                       "decode_core_s_per_run": d["core_s"]["decode_jobs_sum"] / d["n_runs"],
                                       "total_core_s_per_run": d["core_s"]["per_run"],
                                       "core_s_per_1e5_ms2": d["core_s"]["per_1e5_ms2"],
                                       "gpu_util_mean_pct": d["gpu_s"]["util_mean_pct_after_setup"],
                                       "bottleneck": d["bottleneck"]})
    out["scaling_c50"].sort(key=lambda r: r["workers"])
    if out["configs"]:
        out["hw"] = json.load(open(sorted(glob.glob(os.path.join(RES, "runs", "*.json")))[0])).get("hw")
    for extra in ("startup.json", "ingest_profile.json", "h2d.json", "family_validation_D1.json",
                  "family_validation_D0.json", "lustre_read.json"):
        p = os.path.join(RES, extra)
        if os.path.exists(p):
            out[extra[:-5]] = json.load(open(p))
    with open(os.path.join(RES, "metrics.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    for r in out["scaling_c50"]:
        print(r)
    for k, v in out["configs"].items():
        print(k, v["runs_per_hour"], v["stage_capacity"]["binding_stage"])


if __name__ == "__main__":
    sys.exit(main())
