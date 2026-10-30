#!/usr/bin/env python3
"""A6.8.2 step 2: end-to-end admission throughput r from the dry-run event log (label-free; DEV data).

Two estimates, both reported; the tag rule uses the SMALLER N:
  * direct: admitted runs / (time of the last admission - start), on the pipelined D8h dry run (10 runs, PRIDE
    download included).  Includes the pipeline fill latency once, as the real walk would.
  * steady state: r_ss = min over resources of capacity / per-run use, with the measured per-run costs:
      - download: per-stream rate measured here; aggregate ceiling 0.197 GB/s (E9, 64 streams, sustained);
        files downloaded per admitted run from a6_power (p50 113 files at 102 runs);
      - CPU: decode + TRFP + Sage core-s per run (measured) + admission pool core-s, against the cores left free by
        the other node users (32 - mean load of non-A6 processes at the start, floored at 0);
      - the serial GPU admission process: wall per run of certification + removal + reduced-index search; the swap
        test runs in a second process (bounded separately by its own wall per run);
    N_ss = r_ss x (window - fill latency), fill latency = time to the first admission in the direct run.

    scripts/m2b/vpy.sh scripts/a6/throughput.py [--dry results/2026-10-28_A6_hn1/dry_run/D8h/dry_run.json]
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_A6_hn1")
NCORES = 32
DL_CEILING_GBS = 0.197          # E9 sustained aggregate (results/2026-10-27_E9_scale NOTES §1)
FILES_PER_RUN = 113 / 102       # a6_power p50 (files downloaded / runs at the v1 minimum)
MEAN_GB_PER_FILE = 1.08         # E9 D8 mean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", default=os.path.join(RES, "dry_run", "D8h", "dry_run.json"))
    ap.add_argument("--window-h", type=float, default=1.25, help="tag 15:15 -> walk end 16:30")
    ap.add_argument("--load-other", type=float, default=None, help="mean load of non-A6 processes (default: start)")
    ap.add_argument("--log", help="event log (JSON lines) instead of dry_run.json (a run that did not finish)")
    ap.add_argument("--pool-core-s", type=float, default=None,
                    help="measured near-test pool CPU (s) of the run, admission + swap tests (from /proc)")
    a = ap.parse_args()
    if a.log:
        E = [json.loads(x) for x in open(a.log) if x.startswith("{")]
        d = {"events": E}
    else:
        d = json.load(open(a.dry))
        E = d["events"]
    by = lambda k: [e for e in E if e["ev"] == k]  # noqa: E731
    start = next(e for e in E if e["ev"] == "start")
    adm = by("admitted")
    n_runs = sum(e["n_runs"] for e in adm)
    t_first, t_last = min(e["t"] for e in adm), max(e["t"] for e in adm)
    r_direct = n_runs / (t_last / 3600)
    dl = [e for e in by("download") if not e.get("reused")]
    dl_rate = float(np.mean([e["bytes"] / e["wall_s"] for e in dl])) / 1e6 if dl else None
    core = {"decode": [e["core_s"] for e in by("decode") if not e.get("reused") and e.get("rc") == 0 and e.get("core_s")],
            "trfp": [e["core_s"] for e in by("trfp") if not e.get("reused") and e.get("rc") == 0 and e.get("core_s")],
            "sage": [(e.get("core_s") or 0) / max(e["n_runs"], 1) for e in by("sage")]}
    admit = by("admit")
    serial_wall = sum(e["wall_s"] for e in admit) / n_runs
    swaps = by("swap_test")
    swap_wall = sum(e["wall_total_s"] for e in swaps) / n_runs if swaps else 0.0
    # admission pool CPU is not in the GPU process's cpu_s; bound it by (pool size x near-test wall)
    # pool workers are reaped at process exit: children CPU minus the timed subprocesses (decode, TRFP, Sage)
    sub = sum(sum(v) for k, v in core.items() if k != "sage") + sum((e.get("core_s") or 0) for e in by("sage"))
    pool_core = (a.pool_core_s if a.pool_core_s is not None else max(d.get("children_cpu_s", 0.0) - sub, 0.0)) / n_runs
    if swaps:
        pool_core /= 2          # the pool CPU covers admission AND swap-test near tests; re-added below
    core_per_run = sum(float(np.mean(v)) if v else 0.0 for v in core.values()) + pool_core * (2 if swaps else 1) + \
        sum(e["cpu_s"] for e in admit) / n_runs
    load_other = a.load_other if a.load_other is not None else start["loadavg"][1]
    free = max(NCORES - load_other, 0.0)
    r_cpu = free * 3600 / core_per_run if core_per_run else float("inf")
    r_dl = DL_CEILING_GBS * 3600 / (MEAN_GB_PER_FILE * FILES_PER_RUN)
    r_serial = 3600 / serial_wall
    r_swap = 3600 / swap_wall if swap_wall else float("inf")
    r_ss = min(r_cpu, r_dl, r_serial, r_swap)
    fill_h = t_first / 3600
    # swap tests run in the GPU process in the dry run: add their main-process CPU (wall as an upper bound)
    out = {"n_runs_admitted": n_runs, "wall_to_last_admission_s": t_last, "fill_latency_s": t_first,
           "r_direct_runs_per_h": r_direct, "N_direct": r_direct * a.window_h,
           "per_run": {"download_MBps_single_stream": dl_rate,
                       "core_s_mean": {k: float(np.mean(v)) if v else None for k, v in core.items()},
                       "admission_serial_wall_s": serial_wall, "swap_test_wall_s": swap_wall,
                       "admission_pool_core_s_bound": pool_core, "core_s_total": core_per_run},
           "load_other": load_other, "free_cores": free,
           "r_bounds_runs_per_h": {"cpu": r_cpu, "download_aggregate": r_dl, "serial_admission": r_serial,
                                   "swap_process": r_swap},
           "r_steady_state": r_ss, "N_steady_state": r_ss * max(a.window_h - fill_h, 0.0),
           "window_h": a.window_h}
    # CPU bound under other node loads (the walk shares the node with the A5 freeze/G0 and the speed track)
    out["r_cpu_by_load_other"] = {str(L): max(NCORES - L, 0.0) * 3600 / core_per_run for L in (0, 10, 20, 30, 40)}
    out["N_steady_state_by_load_other"] = {
        L: min(r, r_dl, r_serial, r_swap) * max(a.window_h - fill_h, 0.0) for L, r in out["r_cpu_by_load_other"].items()}
    out["N_used"] = min(out["N_direct"], out["N_steady_state"])
    out["r_used"] = out["N_used"] / a.window_h
    json.dump(out, open(os.path.join(RES, "throughput.json"), "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
