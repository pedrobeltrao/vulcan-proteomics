#!/usr/bin/env python3
"""E9 metrics.json: download throughput, end-to-end runs/hour, per-stage costs, E9 table (FDP-hat vs N), failures, hw.

    python3 scripts/e9/summarize.py --analysis /tmp/vulcan_e9/e9.json [--analysis-hs ...] --out results/2026-10-27_E9_scale
"""
import argparse
import glob
import json
import os
import statistics as S
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def deep_merge(a, b):
    """Recursively merge b into a (later analysis files add variants/keys without dropping earlier ones)."""
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict):
            deep_merge(a[k], v)
        else:
            a[k] = v
    return a


def load(p):
    return json.load(open(p)) if p and os.path.exists(p) else {}


def timeline(work):
    return [json.loads(x) for x in open(os.path.join(work, "timeline.jsonl"))]


def sessions(ev):
    """Split the timeline into stream sessions (stream_start .. next stream_start)."""
    out, cur = [], None
    for e in ev:
        if e["ev"] == "stream_start":
            cur = {"start": e, "events": []}
            out.append(cur)
        elif cur is not None:
            cur["events"].append(e)
    return out


def stage_costs(work):
    ing_core, ing_wall, gpu_s, srch, mk, nr, ms2 = [], [], [], [], [], 0, 0
    c_wall, c_gpu, c_core = [], [], []
    failed_ingest = []
    for f in sorted(glob.glob(os.path.join(work, "batches", "b*", "cohort_stats.json"))):
        d = json.load(open(f))
        failed_ingest += d.get("failed", [])
        for r in d["runs"]:
            ing_core.append(r["core_s"])
            ing_wall.append(r["wall_s"])
            gpu_s.append(sum(r["gpu_ms"].values()) / 1e3)
            srch.append(r["search_wall_s"])
            ms2 += r.get("n_ms2", 0)
        nr += d["n_runs"]
        mk.append({"n": d["n_runs"], "makespan_s": d["makespan_s"], "core_s_total": d["core_s"]["total"],
                   "setup_s": d["setup"]["setup_total_s"]})
    for f in sorted(glob.glob(os.path.join(work, "batches", "b*", "C", "stats.json"))):
        d = json.load(open(f))
        for r in d["passes"][0]:
            c_wall.append(r["wall_s_total"])
            c_gpu.append(r["gpu_s"])
            c_core.append(r["cpu_core_s"])
    sc_tot, sc_feat, sc_pred, sc_apply, sc_err = [], [], [], [], []
    for f in sorted(glob.glob(os.path.join(work, "queue", "b*.done"))):
        d = json.load(open(f))
        for r in d["runs"]:
            if "error" in r:
                sc_err.append({"run": r["run"], "error": r["error"]})
            for arm in ("P", "C"):
                x = r.get(arm)
                if isinstance(x, dict):
                    sc_tot.append(x["total_s"])
                    sc_feat.append(x["features_s"])
                    sc_pred.append(x["predict_s"])
                    sc_apply.append(x["apply_s"])
    med = lambda x: S.median(x) if x else None  # noqa: E731
    mean = lambda x: S.fmean(x) if x else None  # noqa: E731
    return {
        "n_runs_searched": nr, "n_ms2": ms2,
        "decode_core_s_per_run": {"mean": mean(ing_core), "median": med(ing_core)},
        "decode_wall_s_per_run": {"mean": mean(ing_wall), "median": med(ing_wall)},
        "search_P_gpu_s_per_run": {"mean": mean(gpu_s), "median": med(gpu_s)},
        "search_P_wall_s_per_run": {"mean": mean(srch), "median": med(srch)},
        "search_C_wall_s_per_run_incl_store_load": {"mean": mean(c_wall), "median": med(c_wall)},
        "search_C_gpu_s_per_run": {"mean": mean(c_gpu), "median": med(c_gpu)},
        "score_s_per_run_arm": {"mean": mean(sc_tot), "median": med(sc_tot),
                                "features_mean": mean(sc_feat), "predict_mean": mean(sc_pred), "gbm_apply_mean": mean(sc_apply),
                                "n": len(sc_tot)},
        "cohort_batches": mk,
        "failed_ingest": failed_ingest, "failed_scoring": sc_err,
    }


def throughput(ev, work):
    out = []
    for s in sessions(ev):
        e = s["events"]
        dl = [x for x in e if x["ev"] == "download"]
        if not dl:
            continue
        t0 = s["start"]["t"]
        t_dl_end = max(x["t"] for x in dl)
        b = sum(x["bytes"] for x in dl)
        scored = [x for x in e if x["ev"] == "batch_scored"]
        rec = {"range": [s["start"]["start"], s["start"]["end"]], "streams": s["start"]["streams"],
               "n_downloaded": len(dl), "bytes": b, "download_wall_s": t_dl_end - t0,
               "download_gb_per_s": b / max(t_dl_end - t0, 1e-9) / 1e9,
               "download_runs_per_h": len(dl) / max(t_dl_end - t0, 1e-9) * 3600,
               "pride_sha1_mismatch": sum(not x.get("pride_sha1_match", True) for x in dl),
               "download_failed": [x for x in e if x["ev"] == "download_failed"],
               "cohort_retries": [x for x in e if x["ev"] == "cohort_retry"]}
        end = [x for x in e if x["ev"] == "stream_end"]
        if end:
            rec["stream_wall_s"] = end[0]["wall_s"]
        if scored:
            rec["last_scored_after_s"] = max(x["t"] for x in scored) - t0
        out.append(rec)
    return out


def table(an):
    rows = []
    for n in sorted((k for k in an if k.startswith("N")), key=lambda k: int(k[1:])):
        for arm in ("P", "C"):
            if arm not in an[n]:
                continue
            blocks = [(v, an[n][arm][v]) for v in ("score_raw", "score_cal", "score_bulk", "hyperscore")
                      if v in an[n][arm]]
            for v, r in blocks:
                for a in ("0.001", "0.01", "0.05"):
                    p = r["peptide"][a]
                    row = {"N": int(n[1:]), "runs_present": an[n][arm].get("n_runs_present"),
                           "projects": an[n][arm].get("n_projects"), "arm": arm, "score": v, "alpha": float(a),
                           "level": "peptidoform"}
                    for k in ("global", "naive"):
                        row[k] = {x: p[k][x] for x in ("mean_fdp", "ci95", "mean_lb", "lb_ci95", "verdict", "n_targets",
                                                       "n_entrapment")}
                    rows.append(row)
                    if "protein" in r and a in r["protein"]:
                        pr = r["protein"][a]
                        prow = dict(row, level="protein_group")
                        for k in ("global", "naive"):
                            prow[k] = {x: pr[k].get(x) for x in ("mean_fdp", "ci95", "mean_lb", "lb_ci95", "verdict",
                                                                  "n_target_groups")}
                        rows.append(prow)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/vulcan_e9")
    ap.add_argument("--analysis", nargs="+", required=True, help="analyze.py outputs, deep-merged in order")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "2026-10-27_E9_scale"))
    a = ap.parse_args()
    an = {}
    for p in a.analysis:
        deep_merge(an, load(p))
    ev = timeline(a.work)
    sys.path.insert(0, ROOT)
    try:
        hw = json.loads(subprocess.run([os.path.join(ROOT, "conda_envs", "pe", "bin", "python"), "-c",
                                        "import json;from vulcan.hw import fingerprint;print(json.dumps(fingerprint(),default=str))"],
                                       capture_output=True, text=True, cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT}).stdout)
    except Exception:  # noqa: BLE001
        hw = None
    res = os.path.join(ROOT, "results", "2026-10-27_E9_scale")
    m = {
        "experiment": "E9 (DESIGN §7) + L6 global FDR, dev-grade, dataset D8",
        "hw": hw,
        "download_probe": [load(os.path.join(res, f)) for f in ("download_probe.json", "download_probe2.json")],
        "stream_sessions": throughput(ev, a.work),
        "stage_costs": stage_costs(a.work),
        "e9_table": table(an),
        "fair_coin_family_counts": {n: {arm: {v: an[n][arm][v].get("family_counts") for v in ("score_raw", "hyperscore")
                                              if v in an[n][arm]} for arm in an[n] if arm in ("P", "C")}
                                    for n in an if n.startswith("N")},
        "comparability_G4_style": {n: {arm: {"ks_excess_median": an[n][arm].get("ks_excess_median"),
                                             "raw_per_run": an[n][arm]["comparability"]["score_raw"]["per_run"].get("summary"),
                                             "raw_per_run_charge": an[n][arm]["comparability"]["score_raw"]["per_run_charge"].get("summary"),
                                             "raw_per_instrument": an[n][arm]["comparability"]["score_raw"]["per_instrument"]}
                                       for arm in an[n] if "comparability" in an[n][arm]}
                                   for n in an if n.startswith("N")},
        "l6_outputs": {n: {arm: {v: {"expected_false": an[n][arm][v]["expected_false"],
                                     "run_cutoff_1pct": an[n][arm][v]["run_cutoff_1pct"],
                                     "global_cutoff_1pct": an[n][arm][v]["global_cutoff_1pct"],
                                     "timing_s": an[n][arm][v]["timing_s"]}
                                 for v in ("score_raw",) if v in an[n][arm]} for arm in an[n]}
                       for n in an if n.startswith("N")},
        "A1_1_confirmation_first10": an.get("A1_1_confirmation_first10"),
        "analysis_meta": an.get("meta"),
    }
    json.dump(m, open(os.path.join(a.out, "metrics.json"), "w"), indent=1, default=float)
    print(json.dumps({"sessions": m["stream_sessions"], "stage": {k: v for k, v in m["stage_costs"].items()
                                                                  if k not in ("cohort_batches", "failed_ingest")}},
                     indent=1, default=float)[:4000])


if __name__ == "__main__":
    sys.exit(main())
