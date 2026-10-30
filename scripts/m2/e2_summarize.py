#!/usr/bin/env python3
"""Collect E2 outputs into results/2026-10-27_E2_gpu_search/metrics.json (+ git_sha).

Inputs: timings.jsonl (timed.py records), ours stats.json per bench rep, Sage logs, correctness.json,
topk_sweep.json, ids.json.
"""
from __future__ import annotations

import json
import os
import re
import statistics as S
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_common as C  # noqa: E402


def med(xs):
    xs = [x for x in xs if x is not None]
    return S.median(xs) if xs else None


def sage_log(path):
    t = open(path).read()
    g = lambda pat, f=float: (f(m.group(1)) if (m := re.search(pat, t)) else None)  # noqa: E731
    return {"index_build_ms": g(r"fragments, \d+ peptides in (\d+)ms"), "n_fragments": g(r"generated (\d+) fragments", int),
            "n_peptides": g(r"fragments, (\d+) peptides", int), "file_io_ms": g(r"file IO:\s+(\d+) ms"),
            "search_ms": g(r"- search:\s+(\d+) ms"), "search_spectra_per_s": g(r"ms \((\d+) spectra/s\)"),
            "finished_s": g(r"finished in (\d+)s")}


def main():
    recs = [json.loads(l) for l in open(os.path.join(C.RES, "timings.jsonl"))]
    out = {"experiment": "E2", "dev_only": True}
    hw = None
    speed = {}
    for ds in C.DATASETS:
        n_ms2 = 0
        for raw in C.DATASETS[ds]["raw"]:
            n_ms2 += json.load(open(C.store_path(raw) + ".json"))["n_ms2"]
        d = {"n_ms2_spectra": n_ms2, "n_files": len(C.DATASETS[ds]["raw"])}
        for tool, lab in (("ours", f"E2_ours_{ds}_prod"), ("sage_t32", f"E2_sage_t32_{ds}_prod"),
                          ("sage_t4", f"E2_sage_t4_{ds}_prod")):
            rr = [r for r in recs if r["label"] == lab and r["rc"] == 0]
            if not rr:
                continue
            hw = hw or rr[0]["hw"]
            e = {"n_reps": len(rr), "process_wall_s_median": med([r["wall_s"] for r in rr]),
                 "process_core_s_median": med([r["core_s"] for r in rr]), "max_rss_mb": max(r["max_rss_mb"] for r in rr),
                 "lock_wait_s": [r["lock_wait_s"] for r in rr], "loadavg_before": [r["loadavg_before"][0] for r in rr],
                 "threads_env": rr[0]["threads_env"]}
            e["end_to_end_spectra_per_s"] = n_ms2 / e["process_wall_s_median"]
            e["core_s_per_1e5_ms2"] = e["process_core_s_median"] * 1e5 / n_ms2
            if tool == "ours":
                st = [json.load(open(r["stats"])) for r in rr if os.path.exists(r["stats"])]
                first = [sum(x["wall_s_total"] for x in s["passes"][0]) for s in st]
                warm = [sum(x["wall_s_total"] for x in p) for s in st for p in s["passes"][1:]]
                srch = [sum(x["search_wall_s"]["total"] for x in p) for s in st for p in s["passes"][1:]]
                gpu = [sum(x["gpu_s"] for x in p) for s in st for p in s["passes"][1:]]
                core = [sum(x["cpu_core_s"] for x in p) for s in st for p in s["passes"][1:]]
                ph = {}
                for s in st:
                    for p in s["passes"][1:]:
                        for x in p:
                            for k, v in x["gpu_ms"].items():
                                ph.setdefault(k, []).append(v)
                # the bench process runs `--repeat` passes; the single-pass process cost removes the extra warm passes
                sp_wall = [r["wall_s"] - sum(x["wall_s_total"] for p in s["passes"][1:] for x in p) for r, s in zip(rr, st)]
                sp_core = [r["core_s"] - sum(x["cpu_core_s"] for p in s["passes"][1:] for x in p) for r, s in zip(rr, st)]
                e["process_wall_s_single_pass_median"] = med(sp_wall)
                e["process_core_s_single_pass_median"] = med(sp_core)
                e["end_to_end_spectra_per_s"] = n_ms2 / e["process_wall_s_single_pass_median"]
                e["core_s_per_1e5_ms2"] = e["process_core_s_single_pass_median"] * 1e5 / n_ms2
                e.update({
                    "setup_median": {k: med([s["setup"][k] for s in st]) for k in
                                     ("python_import_s", "cuda_init_s", "l1_load_s", "gpu_index_build_s")},
                    "gpu_index_bytes": st[0]["setup"]["gpu_index_bytes"], "n_postings": st[0]["setup"]["n_postings"],
                    "gpu_max_mem_gb": max(s["gpu_max_mem_bytes"] for s in st) / 1e9,
                    "runs_first_pass_wall_s_median": med(first),
                    "runs_warm_wall_s_median (store load + K0..K2 + features + parquet)": med(warm),
                    "runs_warm_search_wall_s_median (K0..K2 + features)": med(srch),
                    "runs_warm_gpu_event_s_median": med(gpu),
                    "runs_warm_cpu_core_s_median": med(core),
                    "h2d_bytes_per_dataset": sum(x["h2d_bytes"] for x in st[0]["passes"][0]),
                    "gpu_index_h2d_bytes": st[0]["setup"]["gpu_index_h2d_bytes"],
                    "warm_spectra_per_s (incl. store load)": n_ms2 / med(warm),
                    "warm_spectra_per_s (search only)": n_ms2 / med(srch),
                    "warm_phase_gpu_event_ms_median_per_file": {k: med(v) for k, v in ph.items()},
                })
            else:
                logs = [sage_log(os.path.join(r["outdir"], "sage.log")) for r in rr]
                e.update({k: med([l[k] for l in logs]) for k in logs[0]})
                e["search_phase_spectra_per_s"] = e.pop("search_spectra_per_s")
            d[tool] = e
        if "ours" in d and "sage_t32" in d:
            d["ratio_process_wall_sage32_over_ours(single pass)"] = d["sage_t32"]["process_wall_s_median"] / d["ours"]["process_wall_s_single_pass_median"]
            d["ratio_core_s_sage32_over_ours(single pass)"] = d["sage_t32"]["process_core_s_median"] / d["ours"]["process_core_s_single_pass_median"]
            d["ratio_search_phase_sage32_over_ours_warm"] = (d["sage_t32"]["search_ms"] / 1e3) / \
                d["ours"]["runs_warm_search_wall_s_median (K0..K2 + features)"]
        if "ours" in d and "sage_t4" in d:
            d["ratio_process_wall_sage4_over_ours(single pass)"] = d["sage_t4"]["process_wall_s_median"] / d["ours"]["process_wall_s_single_pass_median"]
            d["ratio_core_s_sage4_over_ours(single pass)"] = d["sage_t4"]["process_core_s_median"] / d["ours"]["process_core_s_single_pass_median"]
        speed[ds] = d
    # pipelined multi-run steady state (vulcan.gpu.batch): D0 stores x5 = 10 runs in one process
    batch = {}
    for r in recs:
        if r["label"].startswith("E2_ours_batch") and r["rc"] == 0:
            tag = r["label"].replace("E2_ours_batch_", "")
            bdir = os.path.join(C.WORK, "bench", ("batch_io1" if "io1" in tag else "batch_io2") + ("_b" if "rerun" in tag else ""))
            f = os.path.join(bdir, "batch_stats.json")
            if os.path.exists(f):
                b = json.load(open(f))
                per = [x["search_wall_s"] for x in b["runs"]]
                batch[tag] = {"process_wall_s": r["wall_s"], "process_core_s": r["core_s"], "setup_s": b["setup_s"],
                              "steady_wall_s_10_runs": b["wall_s"], "steady_cpu_core_s": b["cpu_core_s"],
                              "n_spectra": b["n_spectra"], "steady_spectra_per_s": b["spectra_per_s_steady"],
                              "per_run_search_wall_s": per,
                              "valid": "io1" not in tag,
                              "note": ("EXCLUDED: ~190 s setup + ~6 s first-search stall in cl_sync_io_wait (Lustre) while other "
                                       "agents ran I/O-heavy jobs; an unlocked rerun right after was normal (NOTES.md)")
                              if "io1" in tag else ""}
    out["pipelined_batch_D0x5"] = batch
    out["hw"] = hw
    out["speed"] = speed
    for name in ("correctness", "topk_sweep", "ids"):
        f = os.path.join(C.RES, f"{name}.json")
        if os.path.exists(f):
            out[name] = json.load(open(f))
    with open(os.path.join(C.RES, "metrics.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=C.ROOT).stdout.strip()
    with open(os.path.join(C.RES, "git_sha"), "w") as fh:
        fh.write(sha + "\n")
    print(json.dumps(speed, indent=1, default=str))


if __name__ == "__main__":
    main()
