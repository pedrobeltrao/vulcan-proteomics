#!/usr/bin/env python3
"""Runtime impact of chimeric search + v2 scoring (dev D0, prod_1001 index), one persistent process, warm caches,
timed sections under `flock /tmp/vulcan_bench.lock` (run this script under the lock).

Per run, median of 3 repeats after one warm-up:
  closed   engine.search (report 5) + v1 features (FeatureBuilder, PredCache) + frozen v1 GBM
  chimeric MS1 features + hypotheses (vulcan.ms1, in-process, 8 threads) + chimeric search + v2 features + v2 model
Wall and process core-s (RUSAGE_SELF), GPU-s (CUDA events of the engine).

    flock /tmp/vulcan_bench.lock /tmp/vulcan_env_<uid>/pred/bin/python scripts/e7/bench_runtime.py
"""
from __future__ import annotations

import json
import os
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

W = f"/tmp/vulcan_t1_{os.getuid()}/D0"
RES = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2")
RUNS = ["b1906_293T_proteinID_01A_QE3_122212", "b1937_293T_proteinID_01B_QE3_122212"]
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}


def cs():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main():
    from vulcan import hw
    n_thr = hw.pin_thread_pools(8)
    import torch

    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ms1.features import FeatureParams, find_features
    from vulcan.ms1.hypotheses import build_hypotheses
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    from vulcan.rescore.features_v2 import FeatureBuilderV2, model_psm
    from vulcan.rescore.static import StaticModel
    v1 = StaticModel.load(os.path.join(ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl"))
    v2 = StaticModel.load(os.path.join(RES, "models", "static_v2.pkl"))
    pred = APD()
    caches = {}

    def fac(nce):
        if nce not in caches:
            caches[nce] = PredCache(os.path.join(W, "predcache"), pred, nce, "QE", STATIC, VAR)
        return caches[nce]
    fb1, fb2 = FeatureBuilder(fac, STATIC, VAR), FeatureBuilderV2(fac, STATIC, VAR)
    p = E.SearchParams(report_psms=5)
    l1 = L1Index(os.path.join(W, "index", "prod_1001"))
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    out = {"hw_threads": n_thr, "loadavg_start": os.getloadavg(), "runs": {}}
    for run in RUNS:
        store = os.path.join(W, "stage", run + ".vspec.arrow")
        sp = E.load_store(store)
        rec = {"closed": [], "chimeric": []}
        for rep in range(4):
            # closed
            torch.cuda.synchronize()
            t0, c0 = time.perf_counter(), cs()
            tbl, st, _ = E.search(sp, gi, l1, p)
            t1 = time.perf_counter()
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
            f, _ = fb1.run(l1, psm, sp)
            s = v1.apply(psm, f)
            torch.cuda.synchronize()
            t2, c2 = time.perf_counter(), cs()
            closed = {"search_s": t1 - t0, "score_s": t2 - t1, "total_s": t2 - t0, "core_s": c2 - c0,
                      "gpu_s_search": sum(st.gpu_ms.values()) / 1e3, "n_psms": len(s)}
            # chimeric
            t0, c0 = time.perf_counter(), cs()
            fr = find_features(store, FeatureParams())
            hyp, _ = build_hypotheses(store, fr)
            t1 = time.perf_counter()
            groups = CS.precursor_groups(hyp, sp)
            tbl2, st2 = CS.search_chimeric(sp, gi, l1, p, groups)
            torch.cuda.synchronize()
            t2 = time.perf_counter()
            psm2 = {c: tbl2.column(c).to_numpy(zero_copy_only=False) for c in tbl2.column_names}
            ff = fr.features
            mf = {c: ff.column(c).to_numpy() for c in ("trace_cos", "n_isotopes", "rt_apex", "rt_start", "rt_end",
                                                       "intensity_apex")}
            f2, _ = fb2.run_v2(l1, psm2, sp, mf)
            s2 = v2.apply(model_psm(psm2), f2)
            torch.cuda.synchronize()
            t3, c3 = time.perf_counter(), cs()
            chim = {"ms1_s": t1 - t0, "search_s": t2 - t1, "score_s": t3 - t2, "total_s": t3 - t0, "core_s": c3 - c0,
                    "gpu_s_search": st2["primary"]["gpu_s"] + st2["secondary"]["gpu_s"],
                    "n_hypotheses": st2["primary"]["n_hypotheses"] + st2["secondary"]["n_hypotheses"],
                    "n_hypotheses_closed": st2["primary"]["n_hypotheses"], "n_psms": len(s2)}
            if rep > 0:
                rec["closed"].append(closed)
                rec["chimeric"].append(chim)
            print(run, rep, closed, chim, flush=True)
        med = {k: {x: float(np.median([r[x] for r in v])) for x in v[0]} for k, v in rec.items()}
        out["runs"][run] = {"median": med, "repeats": rec}
    out["loadavg_end"] = os.getloadavg()
    json.dump(out, open(os.path.join(RES, "runtime.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
