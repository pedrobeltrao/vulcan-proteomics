#!/usr/bin/env python3
"""E11 MS1 equality record (dev only: D0, D1, D8h): GPU backend vs CPU reference, bitwise.

For every run (<work>/<ds>/stage/<run>.vspec.arrow):
  (a) find_features(backend="gpu") vs find_features() (CPU default, itself bitwise equal to the original kernels,
      E7b): every Hills array, the feature table (values + schema metadata), peak_feature / peak_iso, params;
  (b) build_hypotheses on each result: every column of the hypothesis table (values + validity) + metadata;
  (c) both vs the tables E7 wrote with the original code (<work>/<ds>/ms1/<run>.{features,hypotheses}.arrow), the
      same check as E7b's equality_ms1.json.
Writes --out JSON.  Same equality helpers as scripts/e7b/equal_ms1.py.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "e7b"))

import numpy as np  # noqa: E402
from equal_ms1 import read_arrow, table_equal  # noqa: E402

HILL_FIELDS = ("off", "p_scan", "p_mz", "p_int", "p_src", "h_mz", "h_start", "h_end", "h_apex_scan", "h_apex_int",
               "h_area", "h_len")


def same(a, b) -> bool:
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1", "D8h"])
    ap.add_argument("--work-root", default=f"/tmp/vulcan_t1_{os.getuid()}")
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args()
    from vulcan import hw
    n_thr = hw.pin_thread_pools(a.threads)
    from vulcan.ms1.features import FeatureParams, find_features, load_ms1
    from vulcan.ms1.hypotheses import build_hypotheses
    out = {"hw": hw.fingerprint(), "threads": n_thr, "loadavg_start": os.getloadavg(), "runs": {}}
    for ds in a.datasets:
        W = os.path.join(a.work_root, ds)
        for f in sorted(os.listdir(os.path.join(W, "stage"))):
            if not f.endswith(".vspec.arrow"):
                continue
            run = f[:-len(".vspec.arrow")]
            store = os.path.join(W, "stage", f)
            t0 = time.time()
            ms1 = load_ms1(store)
            fc = find_features(ms1, FeatureParams())
            fg = find_features(ms1, FeatureParams(), backend="gpu")
            hc, _ = build_hypotheses(store, fc)
            hg, _ = build_hypotheses(store, fg)
            rec = {"dataset": ds, "backend_ran": fg.timings.get("backend"),
                   "hills_gpu_vs_cpu": {k: same(getattr(fc.hills, k), getattr(fg.hills, k)) for k in HILL_FIELDS},
                   "features_gpu_vs_cpu": table_equal(fg.features, fc.features),
                   "peak_feature_equal": same(fc.peak_feature, fg.peak_feature),
                   "peak_iso_equal": same(fc.peak_iso, fg.peak_iso),
                   "params_equal": fc.params == fg.params,
                   "hypotheses_gpu_vs_cpu": table_equal(hg, hc)}
            e7f = os.path.join(W, "ms1", run + ".features.arrow")
            e7h = os.path.join(W, "ms1", run + ".hypotheses.arrow")
            if os.path.exists(e7f) and os.path.exists(e7h):
                rec["features_gpu_vs_E7"] = table_equal(fg.features, read_arrow(e7f))
                rec["hypotheses_gpu_vs_E7"] = table_equal(hg, read_arrow(e7h))
            rec["n_peaks"] = int(len(ms1.mz))
            rec["n_hills"] = int(fg.hills.n)
            rec["n_candidates"] = fg.params["n_candidates"]
            rec["n_features"] = fg.features.num_rows
            rec["n_hypotheses"] = hg.num_rows
            rec["wall_s"] = time.time() - t0
            tb = [rec["features_gpu_vs_cpu"], rec["hypotheses_gpu_vs_cpu"]] + \
                ([rec["features_gpu_vs_E7"], rec["hypotheses_gpu_vs_E7"]] if "features_gpu_vs_E7" in rec else [])
            rec["ok"] = (rec["backend_ran"] == "gpu" and all(rec["hills_gpu_vs_cpu"].values())
                         and rec["peak_feature_equal"] and rec["peak_iso_equal"] and rec["params_equal"]
                         and all(t["columns_equal"] and t["metadata_equal"] for t in tb))
            out["runs"][run] = rec
            print(ds, run[:34], rec["ok"], rec["n_hills"], rec["n_features"], rec["n_hypotheses"],
                  round(rec["wall_s"], 1), flush=True)
            del fc, fg, hc, hg, ms1
    out["all_ok"] = all(r["ok"] for r in out["runs"].values())
    out["n_runs"] = len(out["runs"])
    out["n_runs_vs_E7"] = sum("features_gpu_vs_E7" in r for r in out["runs"].values())
    out["loadavg_end"] = os.getloadavg()
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print("ALL_OK", out["all_ok"], out["n_runs"], "vs_E7:", out["n_runs_vs_E7"])


if __name__ == "__main__":
    main()
