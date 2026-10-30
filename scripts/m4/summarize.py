#!/usr/bin/env python3
"""Collect M4 outputs into results/2026-10-27_M4_ptm/metrics.json (+ copies of stats/e8/open json)."""
import glob
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

W = C.WORK
R = C.RES
m = {"note": "dev only (D0, D1-free, D2 dev pools); development curves, not claims (DESIGN §6.6). T2 untouched."}
speed = {}
for d in sorted(glob.glob(os.path.join(W, "bench", "*", "stats.json"))):
    name = d.split(os.sep)[-2]
    s = json.load(open(d))
    shutil.copy(d, os.path.join(R, f"stats_{name}.json"))
    warm = s["passes"][-1]
    speed[name] = {
        "mode": s["mode"], "n_peptidoforms": s["setup"]["n_peptidoforms"], "n_postings": s["setup"]["n_postings"],
        "gpu_index_bytes": s["setup"]["gpu_index_bytes"], "gpu_index_build_s": s["setup"]["gpu_index_build_s"],
        "l1_load_s": s["setup"]["l1_load_s"], "gpu_max_mem_bytes": s["gpu_max_mem_bytes"],
        "process_wall_s_2reps": s["process_wall_s"], "process_core_s_2reps": s["process_cpu_core_s"],
        "warm_runs": [{"run": r["run"], "n_spectra": r["counts"]["n_spectra"], "wall_s": r["wall_s_total"],
                       "search_s": r["search_wall_s"]["total"], "core_s": r["cpu_core_s"],
                       "gpu_ms": r["gpu_ms"], "spectra_per_s_search": r["spectra_per_s_search"],
                       "k1_window_candidates": r["counts"].get("k1_window_candidates"),
                       "n_rescored": r["counts"].get("n_rescored")} for r in warm],
    }
    speed[name]["warm_total_wall_s"] = sum(r["wall_s"] for r in speed[name]["warm_runs"])
    speed[name]["warm_total_core_s"] = sum(r["core_s"] for r in speed[name]["warm_runs"])
    hw = s.get("hw")
m["speed_ours"] = speed
m["hw"] = hw
tim = [json.loads(l) for l in open(os.path.join(R, "timings.jsonl"))]
m["speed_sage"] = [{k: t.get(k) for k in ("label", "wall_s", "core_s", "user_s", "sys_s", "threads_env", "start",
                                          "loadavg_before")} for t in tim]
idx_meta = json.load(open(os.path.join(C.index_dir("d2_phospho_prod"), "meta.json")))
m["phospho_index"] = {k: idx_meta[k] for k in ("n_targets", "n_peptidoforms", "n_peptidoforms_by_family", "timing_s",
                                               "content_hash")}
m["phospho_index"]["params"] = {k: idx_meta["params"][k] for k in ("variable_mods", "max_variable_mods",
                                                                   "static_mods", "enzyme")}
m["phospho_index"]["build_note"] = "CPU build, unlocked, node loadavg ~15 (other agents)"
for f in ("open_d0.json", "e8.json"):
    p = os.path.join(R, f)
    if os.path.exists(p):
        m[f[:-5]] = json.load(open(p))
C.dump(m, os.path.join(R, "metrics.json"))
print("ok")
