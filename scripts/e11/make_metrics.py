#!/usr/bin/env python3
"""E11 metrics.json from the result files in results/2026-10-28_E11_gpu_ms1 (equality, per-run benches, cohorts)."""
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_E11_gpu_ms1")


def load(name):
    p = os.path.join(RES, name)
    return json.load(open(p)) if os.path.exists(p) else None


def per_run(b, stage="ms1"):
    out = {}
    for ds, d in b["datasets"].items():
        for run, r in d["runs"].items():
            m = r["median"]
            out[run] = {"dataset": ds, **{k: m[stage][k] for k in ("wall_s", "core_s", "gpu_s") if k in m[stage]},
                        "sub_wall_s": m[stage].get("sub_wall_s"), "total_core_s": m["total"]["core_s"],
                        "total_wall_s": m["total"]["wall_s"]}
    return out


def bench_block(cpu, gpu):
    if not cpu or not gpu:
        return None
    sc, sg = cpu["summary_sum_over_runs"], gpu["summary_sum_over_runs"]
    stages = {}
    for k in sc:
        stages[k] = {"cpu": {m: v for m, v in sc[k].items() if m != "sub_wall_s"},
                     "gpu": {m: v for m, v in sg[k].items() if m != "sub_wall_s"}}
        if sc[k].get("core_s"):
            stages[k]["core_s_ratio_cpu_over_gpu"] = sc[k]["core_s"] / max(1e-9, sg[k]["core_s"])
            stages[k]["wall_s_ratio_cpu_over_gpu"] = sc[k]["wall_s"] / max(1e-9, sg[k]["wall_s"])
    pc, pg = per_run(cpu), per_run(gpu)
    d0 = [r for r in pc if pc[r]["dataset"] == "D0"]
    mean = lambda d, k: sum(d[r][k] for r in d0) / len(d0)  # noqa: E731
    return {"threads": [cpu["threads"], gpu["threads"]], "sum_over_13_runs": stages,
            "ms1_sub_wall_s_sum": {"cpu": sc["ms1"].get("sub_wall_s"), "gpu": sg["ms1"].get("sub_wall_s")},
            "per_D0_run_mean": {"ms1_wall_s": [mean(pc, "wall_s"), mean(pg, "wall_s")],
                                "ms1_core_s": [mean(pc, "core_s"), mean(pg, "core_s")],
                                "ms1_gpu_s": mean(pg, "gpu_s"),
                                "total_wall_s": [mean(pc, "total_wall_s"), mean(pg, "total_wall_s")],
                                "total_core_s": [mean(pc, "total_core_s"), mean(pg, "total_core_s")]},
            "per_run": {"cpu": pc, "gpu": pg},
            "counts_equal": all(cpu["datasets"][ds]["runs"][r]["counts"] == gpu["datasets"][ds]["runs"][r]["counts"]
                                for ds in cpu["datasets"] for r in cpu["datasets"][ds]["runs"]),
            "loadavg": {"cpu": [cpu["loadavg_start"], cpu["loadavg_end"]],
                        "gpu": [gpu["loadavg_start"], gpu["loadavg_end"]]},
            "hw": gpu["hw"]}


def cohort_block(name):
    c = load(name)
    if not c:
        return None
    keys = ("n_runs", "workers", "makespan_s", "runs_per_hour", "runs_per_hour_steady", "bottleneck",
            "stage_busy_fraction", "ms1_share_of_worker_busy")
    r = {k: c.get(k) for k in keys}
    r["core_s"] = c["core_s"]
    r["gpu_s"] = {k: v for k, v in c["gpu_s"].items() if k != "util_samples"}
    r["per_run_median"] = c["per_run_median"]
    r["per_run_sum"] = c["per_run_sum"]
    r["ms1_backend"] = c.get("ms1_backend", "cpu")
    r["ms1_gpu"] = c.get("ms1_gpu")
    r["gpu_max_mem_bytes_main_process"] = c.get("gpu_max_mem_bytes")
    r["failed"] = c.get("failed")
    r["decode_pool_events"] = c.get("decode_pool_events")
    r["n_ms1_core_s_per_run"] = c["per_run_sum"]["ms1_core_s"] / max(1, c["n_runs"])
    sc = c.get("scoring") or {}
    r["scoring_summary"] = {k: v for k, v in sc.items() if not isinstance(v, (list, dict))} if sc else None
    return r


def main():
    eq = load("equality_gpu.json")
    m = {
        "experiment": "E11 GPU MS1 feature finder (opt-in ms1_backend=gpu), dev only (D0, D1, D8h)",
        "equality_per_run_tables": None if not eq else {
            "n_runs": eq["n_runs"], "all_ok": eq["all_ok"], "n_runs_vs_E7_tables": eq["n_runs_vs_E7"],
            "checked": "Hills arrays (12), feature table (all columns + schema metadata), peak_feature, peak_iso, "
                       "params, hypothesis table (all columns + metadata) GPU vs CPU; feature + hypothesis tables "
                       "GPU vs E7 stored tables",
            "id_set_jaccard": 1.0 if eq["all_ok"] else None,
            "v2_sep_id_ratio_change": 0.0 if eq["all_ok"] else None},
        "bench_t8_full_chimeric": bench_block(load("bench_cpu_t8.json"), load("bench_gpu_t8.json")),
        "bench_t1_ms1_only": bench_block(load("bench_cpu_t1_ms1.json"), load("bench_gpu_t1_ms1.json")),
        "cohort_chim_gpu": cohort_block("cohort_chim_gpu.json"),
        "cohort_chim_cpu": cohort_block("cohort_chim_cpu.json"),
        "equality_cohort_gpu_vs_E7b_base12": None,
    }
    ec = load("equality_cohort_gpu.json")
    if ec:
        m["equality_cohort_gpu_vs_E7b_base12"] = {k: ec.get(k) for k in (
            "n_runs", "n_bases", "all_ok", "n_base_named_runs_bitwise_equal", "n_copies_identity_only",
            "all_ok_up_to_run_identity", "differing_columns_any_run")}
    m["gpu_mem_worker"] = load("gpu_mem_worker.json")
    cg, cc = m["cohort_chim_gpu"], m["cohort_chim_cpu"]
    if cg and cc:
        m["cohort_gpu_over_cpu"] = {"runs_per_hour": cg["runs_per_hour"] / cc["runs_per_hour"],
                                    "runs_per_hour_steady": cg["runs_per_hour_steady"] / cc["runs_per_hour_steady"],
                                    "core_s_per_run_cpu_over_gpu": cc["core_s"]["per_run"] / cg["core_s"]["per_run"],
                                    "ms1_core_s_per_run": [cc["n_ms1_core_s_per_run"], cg["n_ms1_core_s_per_run"]]}
    runs = [json.loads(x) for x in open(os.path.join(RES, "cohort_runs.jsonl"))] \
        if os.path.exists(os.path.join(RES, "cohort_runs.jsonl")) else []
    m["cohort_wrapper_records"] = [{k: v for k, v in r.items() if k != "samples_t_gpumem_tmpfree_load1"}
                                   for r in runs]
    json.dump(m, open(os.path.join(RES, "metrics.json"), "w"), indent=1, default=float)
    print(json.dumps({k: (v if not isinstance(v, dict) else list(v)) for k, v in m.items()}, default=str)[:2000])


if __name__ == "__main__":
    main()
