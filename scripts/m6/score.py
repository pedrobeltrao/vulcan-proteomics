#!/usr/bin/env python3
"""L3 predicted features + frozen M3 static scorer (D7 default: static_gbm_pred) for an M6 PSM directory.

    conda_envs/pred/bin/python scripts/m6/score.py --psms DIR --index INDEX [--model static_gbm_pred]

Writes DIR/<run>.feat.parquet and DIR/<run>.score.parquet (column `score`, same row order as the PSM table) and
DIR/score_info.json.  The model is the frozen OT-HCD model of M3 (never refit here).  D9 is Orbitrap Astral HCD:
the model is applied out of its training class (stated in NOTES); the predictor runs with instrument 'QE'.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

t_start = time.perf_counter()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--model", default="static_gbm_pred")
    ap.add_argument("--dataset", default="D9")
    a = ap.parse_args()
    from vulcan.predict.apd import APD, check_allowed
    check_allowed(a.dataset)
    import torch
    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    from vulcan.rescore.static import StaticModel
    t_import = time.perf_counter() - t_start
    model = StaticModel.load(os.path.join(C.ROOT, "results", "2026-10-27_M3_scoring", "models", a.model + ".pkl"))
    t0 = time.perf_counter()
    pred = APD()
    t_model = time.perf_counter() - t0
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PredCache(os.path.join(C.WORK, "predcache"), pred, nce, "QE", C.STATIC, C.VAR)
        return caches[nce]

    fb = FeatureBuilder(factory, C.STATIC, C.VAR)
    l1 = L1Index(a.index)
    infos = {"python_import_s": t_import, "model_load_s": t_model, "predictor": pred.version, "model": a.model,
             "model_meta": {k: model.meta.get(k) for k in ("name", "instrument_class", "sha256") if k in model.meta},
             "runs": {}}
    for f in sorted(glob.glob(os.path.join(a.psms, "*.psms.parquet"))):
        run = os.path.basename(f)[:-len(".psms.parquet")]
        t = pq.read_table(f)
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
        t1 = time.perf_counter()
        sp = load_store(C.store(run))
        t_load = time.perf_counter() - t1
        t1 = time.perf_counter()
        feats, info = fb.run(l1, psm, sp, nce=None)
        t_feat = time.perf_counter() - t1
        t1 = time.perf_counter()
        s = model.apply(psm, feats)
        t_apply = time.perf_counter() - t1
        pq.write_table(pa.table(feats), os.path.join(a.psms, run + ".feat.parquet"))
        pq.write_table(pa.table({"score": s.astype(np.float64)}), os.path.join(a.psms, run + ".score.parquet"))
        info["timing_s"].update({"load_store": t_load, "features_total": t_feat, "model_apply": t_apply})
        info["n_psms"] = int(len(psm["pf_row"]))
        infos["runs"][run] = info
        print(run, json.dumps({k: info[k] for k in ("nce", "rt_mad_min", "cache_stats")}, default=float),
              round(t_feat, 2), round(t_apply, 2), flush=True)
    infos["gpu_max_mem_bytes"] = torch.cuda.max_memory_allocated()
    with open(os.path.join(a.psms, "score_info.json"), "w") as fh:
        json.dump(infos, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
