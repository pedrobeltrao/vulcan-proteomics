#!/usr/bin/env python3
"""Assemble results/2026-10-27_scoring_speed/metrics.json from the measured JSONs (no new measurement)."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_common as S  # noqa: E402

R = S.RES


def J(n):
    p = os.path.join(R, n)
    return json.load(open(p)) if os.path.exists(p) else None


def med(runs, f):
    x = [f(r) for r in runs]
    x = [v for v in x if v is not None]
    return float(np.median(x)) if x else None


def g(r, k, m="wall_s"):
    v = r["sub_stages"].get(k)
    return None if v is None else v.get(m)


def before_stages(runs):
    def sub(r, m):
        s = r["sub_stages"]
        z = lambda k: (s.get(k) or {}).get(m, 0.0)  # noqa: E731
        return {
            "predictor_inference_rt_ms2_incl_nce_grid": z("nce.predict_ms2") + z("main.predict_ms2") + z("main.predict_rt"),
            "nce_selection_total": z("nce.total"),
            "fragment_matching_features_chosen_nce": z("main.lookup_ms2_incl_predict") + z("main.fragment_mz")
            + z("main.observed") + z("main.fragment_features"),
            "rt_lookup_and_calibration": z("main.lookup_rt_incl_predict") + z("main.rt_calibrate"),
            "peaks": z("peaks"),
            "features_total": z("features_total"),
            "gbm_apply": z("gbm_apply_total"),
            "double_competition": z("double_competition"),
            "picked_protein": z("protein_total"),
            "library_load_io": z("library_load"),
            "load_store_io": z("load_store"),
            "search": z("search"),
        }
    return sub


def after_stages(r, m):
    s = r["sub_stages"]
    z = lambda k: (s.get(k) or {}).get(m, (s.get(k) or {}).get("wall_s", 0.0))  # noqa: E731
    return {
        "predictor_inference_rt_ms2_incl_nce_grid": z("features.of_which_predict"),
        "nce_selection_total": z("features.nce_select"),
        "fragment_matching_features_chosen_nce": z("features.fragment_features"),
        "rt_lookup_and_calibration": z("features.rt"),
        "peaks": z("features.peaks"),
        "features_total": z("features_total"),
        "gbm_apply": z("gbm_apply_total"),
        "double_competition": z("double_competition"),
        "picked_protein": z("protein_total"),
        "library_load_io": 0.0,
        "load_store_io": z("load_store"),
        "search": z("search"),
    }


def table(d, fn):
    runs = [r for r in d["runs"].values() if r["tag"] != "warmup"]
    out = {}
    keys = fn(runs[0], "wall_s").keys()
    for k in keys:
        out[k] = {"median_wall_s": med(runs, lambda r: fn(r, "wall_s")[k]),
                  "sum_wall_s": float(sum(fn(r, "wall_s")[k] for r in runs))}
    for k in ("scoring_total", "end_to_end_total"):
        out[k] = {"median_wall_s": med(runs, lambda r: r[k]["wall_s"]),
                  "median_core_s": med(runs, lambda r: r[k]["cpu_core_s"]),
                  "max_wall_s": float(max(r[k]["wall_s"] for r in runs)),
                  "sum_wall_s": float(sum(r[k]["wall_s"] for r in runs)),
                  "sum_core_s": float(sum(r[k]["cpu_core_s"] for r in runs))}
    out["n_measured_runs"] = len(runs)
    out["per_run"] = {n: {"tag": r["tag"], "n_ms2": r["n_ms2"], "n_psms": r["n_psms"], "nce": r["nce"],
                          "scoring_wall_s": r["scoring_total"]["wall_s"], "scoring_core_s": r["scoring_total"]["cpu_core_s"],
                          "ids_1pct": r["ids_1pct"], "proteins_1pct": r["proteins_1pct"],
                          "ms2_predicted": r["cache_stats"]["n_predicted_ms2"]} for n, r in d["runs"].items()}
    return out


def main():
    m = {"scope": "production scoring of the frozen static_gbm_pred (features + NCE choice + GBM + double competition "
                  "+ picked protein), persistent process, dev data (D0, D1, D8 holdout 10 runs)"}
    bef = J("profile_before.json")
    aft = J("bench_after_library.json")
    afc = J("bench_after_firstcontact.json")
    tb = table(bef, before_stages(None))
    ta = table(aft, after_stages)
    m["before_frozen_firstcontact"] = tb
    m["after_library"] = ta
    if afc:
        m["after_new_code_firstcontact_cache"] = table(afc, after_stages)
    ah = J("bench_after_library_host.json")
    if ah:
        m["after_library_host_resident_placement"] = table(ah, after_stages)
    m["after_library"]["placement"] = aft.get("library_placement", "gpu")
    sp = {}
    for k in ta:
        if isinstance(ta[k], dict) and "median_wall_s" in ta[k] and tb.get(k, {}).get("median_wall_s"):
            a, b = ta[k]["median_wall_s"], tb[k]["median_wall_s"]
            sp[k] = {"before_median_wall_s": b, "after_median_wall_s": a, "speedup_median": b / a if a else None}
            if "sum_wall_s" in ta[k]:
                sp[k]["speedup_sum"] = tb[k]["sum_wall_s"] / ta[k]["sum_wall_s"] if ta[k]["sum_wall_s"] else None
    sp["scoring_total"]["core_s_before_median"] = tb["scoring_total"]["median_core_s"]
    sp["scoring_total"]["core_s_after_median"] = ta["scoring_total"]["median_core_s"]
    sp["scoring_total"]["core_speedup_median"] = tb["scoring_total"]["median_core_s"] / ta["scoring_total"]["median_core_s"]
    m["speedup"] = sp
    m["one_time"] = {"before": bef["one_time"], "after": aft["one_time"],
                     "after_library_resident_bytes": aft.get("library_resident_bytes"),
                     "after_max_rss_bytes": aft.get("max_rss_bytes"), "after_gpu_max_mem_bytes": aft.get("gpu_max_mem_bytes")}
    gp = {}
    for impl in ("frozen", "fast"):
        d = J(f"gpuprof_{impl}.json")
        if d:
            rr = [r for r in d["runs"].values() if r["tag"] != "warmup"]
            gp[impl] = {k: {"median": med(rr, lambda r: r[k]), "sum": float(sum(r[k] for r in rr))}
                        for k in ("search_gpu_s", "features_gpu_s", "gbm_gpu_s", "scoring_gpu_s")}
            gp[impl]["cache"] = "firstcontact (M3e state)" if impl == "frozen" else "library"
    m["gpu_s_cuda_kernels"] = gp
    eq = J("equality.json")
    if eq:
        runs = eq["runs"]
        m["equality"] = {
            "library_root": eq["library_root"], "n_runs": len(runs),
            "all_features_bitwise_identical_to_reference": all(v["bitwise_eq_reference"] for r in runs.values()
                                                               for v in r["features"].values()),
            "all_features_bitwise_identical_old_code_same_library": all(
                v["bitwise_eq_old_code_same_library"] for r in runs.values() for v in r["features"].values()),
            "all_scores_bitwise_identical": all(r["score_bitwise_eq_reference"] for r in runs.values()),
            "gbm_gpu_vs_sklearn_bitwise": all(r["gbm_gpu_vs_sklearn_same_features_bitwise"] for r in runs.values()),
            "nce_choice_identical": all(r["nce_new"] == r["nce_reference"] for r in runs.values()),
            "winners_and_q_identical": all(r["winners_identical"] for r in runs.values()),
            "ids_q1pct_sets_identical": {ds: all(r["ids_q0.01"]["identical"] for r in runs.values() if r["ds"] == ds)
                                         for ds in ("D0", "D1", "D8")},
            "ids_q0.1pct_and_5pct_sets_identical": all(r["ids_q0.001"]["identical"] and r["ids_q0.05"]["identical"]
                                                       for r in runs.values()),
            "proteins_1pct_identical": all(r["proteins_1pct"]["reference"] == r["proteins_1pct"]["new"]
                                           for r in runs.values()),
            "protein_glue_identical": all(r["protein_glue_identical_n_groups_1pct"] for r in runs.values()),
            "search_reproduces_reference_psms": all(r["search_reproduces_reference_psms"] for r in runs.values()),
            "ids_1pct": {n: r["ids_q0.01"]["new"] for n, r in runs.items()},
            "proteins_1pct": {n: r["proteins_1pct"]["new"] for n, r in runs.items()},
            "float_tolerance_needed": 0.0}
    if afc:
        ref_ids = {n: r["ids_q0.01"]["new"] for n, r in eq["runs"].items()} if eq else {}
        m["firstcontact_new_code_ids_1pct_equal_reference_counts"] = all(
            r["ids_1pct"] == ref_ids.get(n) for n, r in afc["runs"].items())
    m["predictor_precision"] = J("predictor_precision.json")
    lb = J("library_build.json")
    if lb:
        m["library_build"] = {k: lb.get(k) for k in ("root", "source", "charges", "copy_s", "consolidate_s", "total_s",
                                                     "bytes", "keys")}
        m["library_build"]["digest_sha256"] = lb["digest"]["sha256"]
        m["library_build"]["precompute_steps"] = lb["steps"]
        m["library_build"]["note"] = "one-time per proteome x instrument class; GPU shared (unlocked), loadavg recorded"
    co = J("cohort_summary.json")
    if co:
        cs = {}
        for cfg, v in co.items():
            st = J(f"cohort_{cfg}.json")
            if st is None:
                continue
            runs = st["runs"]
            rec = {"n_runs": v["n_runs"], "makespan_s": v["makespan_s"], "runs_per_hour": v["runs_per_hour"],
                   "runs_per_hour_steady_10_90pct_ingest": v["runs_per_hour_steady"],
                   "runs_per_hour_incl_process_startup": v["runs_per_hour_incl_startup"],
                   "bottleneck": v["bottleneck"], "stage_busy_fraction": v["stage_busy_fraction"],
                   "core_s_per_run": v["core_s"]["per_run"], "gpu_s_events_per_run": v["gpu_s"]["per_run"],
                   "loadavg_end": v["loadavg_end"], "gpu_max_mem_bytes": st.get("gpu_max_mem_bytes"),
                   "decode_timeouts": st["decode_pool_events"].get("timeouts"),
                   "search_threads": st.get("search_threads")}
            if "scoring" in st:
                g = np.array([r["score_gpu_stage_wall_s"] for r in runs])
                sw = np.array([r["search_wall_s"] for r in runs])
                rec.update({"scorer_setup": st["setup"].get("scorer"),
                            "score_gpu_thread_wall_s_median": float(np.median(g)),
                            "gpu_thread_capacity_runs_per_hour": float(3600 / np.mean(g + sw)),
                            "post_wall_s_median": float(np.median([r["post"]["post_wall_s"] for r in runs])),
                            "post_core_s_per_run": float(np.mean([r["post"]["post_core_s"] for r in runs])),
                            "vs_reference_all_identical": all(
                                x.get("psm_table_identical_to_reference") and x.get("ids_1pct_set_identical_to_reference",
                                                                                  cfg.startswith("frozen"))
                                and x.get("score_bitwise_equal_frozen_reference", cfg.startswith("frozen"))
                                for x in v.get("vs_reference", {}).values()),
                            "ids_1pct_by_run": {r["run"]: r["post"]["ids_1pct"] for r in runs if "_r" != r["run"][-4:-2]},
                            "proteins_1pct_by_run": {r["run"]: r["post"]["proteins_1pct"] for r in runs
                                                     if "_r" != r["run"][-4:-2]}})
            cs[cfg] = rec
        m["cohort"] = cs
    m["hw"] = aft["hw"]
    m["loadavg_end"] = {"before": bef.get("loadavg_end"), "after": aft.get("loadavg_end")}
    json.dump(m, open(os.path.join(R, "metrics.json"), "w"), indent=1, default=float)
    for k, v in sp.items():
        print(f"{k:45s} {v['before_median_wall_s']:.3f} -> {v['after_median_wall_s']:.3f}  x{(v['speedup_median'] or 0):.1f}")


if __name__ == "__main__":
    main()
