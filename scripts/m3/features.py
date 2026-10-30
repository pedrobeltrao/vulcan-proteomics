#!/usr/bin/env python3
"""Predicted-feature assembly (L3, label-blind) for PSM directories written by vulcan.gpu.run.

    conda_envs/pred/bin/python scripts/m3/features.py --psms DIR [DIR ...] --index-of DIR=INDEX ... [--dataset D0]

For every <run>.psms.parquet in DIR writes <run>.feat.parquet (same row order) and feat_info.json (NCE choice,
RT calibration, cache hit rate, per-stage timing).  The store of each run is found in the E2 / M3 store dirs.
Decision D4: refuses D2/T2 (vulcan.predict.apd.check_allowed).
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
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402


def find_store(run):
    for d in (os.path.join(C.E2.WORK, "store"), os.path.join(C.WORK, "store_train"), os.path.join(C.WORK, "store_null")):
        p = os.path.join(d, run + ".vspec.arrow")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(run)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", nargs="+", required=True)
    ap.add_argument("--index", nargs="+", required=True, help="one index dir per --psms dir")
    ap.add_argument("--dataset", default="")
    ap.add_argument("--nce", type=float, default=None, help="fixed NCE (default: per-run label-blind grid choice)")
    ap.add_argument("--spectra", default=None, help="optional pickled Spectra override dir (G2 nulls)")
    a = ap.parse_args()
    from vulcan.predict.apd import APD, check_allowed
    check_allowed(a.dataset)
    import torch
    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    t_import = time.perf_counter() - t_start
    t0 = time.perf_counter()
    pred = APD()
    t_model = time.perf_counter() - t0
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PredCache(os.path.join(C.WORK, "predcache"), pred, nce, "QE", C.STATIC, C.VAR)
        return caches[nce]

    fb = FeatureBuilder(factory, C.STATIC, C.VAR)
    for d, idx in zip(a.psms, a.index):
        l1 = L1Index(idx)
        infos = {"python_import_s": t_import, "model_load_s": t_model, "predictor": pred.version, "runs": {}}
        for f in sorted(glob.glob(os.path.join(d, "*.psms.parquet"))):
            run = os.path.basename(f)[:-len(".psms.parquet")]
            t = pq.read_table(f)
            psm = {k: t.column(k).to_numpy() for k in ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt")}
            if a.spectra:
                import pickle
                sp = pickle.load(open(os.path.join(a.spectra, run + ".spectra.pkl"), "rb"))
            else:
                t1 = time.perf_counter()
                sp = load_store(find_store(run))
                t_load = time.perf_counter() - t1
            feats, info = fb.run(l1, psm, sp, nce=a.nce)
            info["timing_s"]["load_store"] = t_load if not a.spectra else 0.0
            info["n_psms"] = int(len(psm["pf_row"]))
            info["n_spectra"] = int(len(sp))
            pq.write_table(pa.table(feats), os.path.join(d, run + ".feat.parquet"))
            infos["runs"][run] = info
            print(run, json.dumps({k: info[k] for k in ("nce", "rt_mad_min", "cache_stats")}, default=float), flush=True)
        infos["gpu_max_mem_bytes"] = torch.cuda.max_memory_allocated()
        with open(os.path.join(d, "feat_info.json"), "w") as fh:
            json.dump(infos, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
