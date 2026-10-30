#!/usr/bin/env python3
"""E11 per-run cost of the chimeric v2_sep production path with the MS1 stage on the CPU (E7b default) or the GPU.

Same protocol and stages as scripts/e7b/bench_run.py (impl "after"), so the numbers line up with E7b bench_after.json:
one persistent process, per run one warm-up + --reps measured repeats under `flock /tmp/vulcan_bench.lock` (taken
per run), dataset's T1 prod_1001 index, consolidated prediction cache, stages ms1 / search / features / model /
compete / proteins with wall, process core-s (RUSAGE_SELF, all threads) and GPU-s.  Only the MS1 stage differs:
  --ms1-backend cpu  find_features() (E7b fast exact numba path)
  --ms1-backend gpu  find_features(backend="gpu") (E11; bitwise-equal tables, equality_gpu.json)
MS1 sub-stage wall times (find_features timings) and the MS1 kernels' GPU time (CUDA events) are recorded too.
Threads: --threads (default vulcan.hw.n_threads()); E7b used 8.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import resource
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}


def cs():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


class R:
    pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms1-backend", choices=["cpu", "gpu"], required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1", "D8h"])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--work-root", default=f"/tmp/vulcan_t1_{os.getuid()}")
    ap.add_argument("--lib-root", default=os.path.join(ROOT, "data", "cache", "e7b_predlib"))
    ap.add_argument("--lock", default=os.environ.get("VULCAN_BENCH_LOCK", "/tmp/vulcan_bench.lock"))
    ap.add_argument("--ms1-only", action="store_true", help="time only the MS1 stage (no search/scoring)")
    a = ap.parse_args()
    from vulcan import hw
    n_thr = hw.pin_thread_pools(a.threads)
    import torch

    from bench.t1.ours_chimeric import MODEL_V2
    from vulcan.chimeric import search as CS
    from vulcan.chimeric.compete import IonSets
    from vulcan.chimeric.fast_compete import chimeric_competition_fast
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ms1.features import FeatureParams, find_features
    from vulcan.ms1.hypotheses import build_hypotheses
    from vulcan.predict.apd import APD
    from vulcan.predict.library import LibCache
    from vulcan.rescore import post as P
    from vulcan.rescore.fast import FastFeatureBuilder, FastStatic
    from vulcan.rescore.features_v2 import model_psm
    from vulcan.rescore.static import StaticModel
    backend = a.ms1_backend
    if backend == "gpu":
        from vulcan.ms1 import gpu as G
        G.lib()
    v2 = StaticModel.load(MODEL_V2)
    pred = APD()
    p = E.SearchParams(report_psms=5)
    out = {"ms1_backend": backend, "threads": n_thr, "hw": hw.fingerprint(), "reps": a.reps, "datasets": {},
           "ms1_only": a.ms1_only, "loadavg_start": os.getloadavg()}
    lockf = open(a.lock, "a")
    for ds in a.datasets:
        W = os.path.join(a.work_root, ds)
        lib = os.path.join(a.lib_root, ds)
        caches = {}

        def factory(nce):
            if nce not in caches:
                caches[nce] = LibCache(lib, pred, nce, "QE", STATIC, VAR, device="cuda")
            return caches[nce]
        fb = FastFeatureBuilder(factory, STATIC, VAR)
        model = FastStatic(v2)
        t0 = time.time()
        if not a.ms1_only:
            index_dir = os.path.join(W, "index", "prod_1001")
            l1 = L1Index(index_dir)
            gi = E.GpuIndex(l1, p, STATIC, VAR)
            ctx = P.IndexCtx(index_dir)
            ions = IonSets(l1)
        setup_s = time.time() - t0
        runs = sorted(f[:-len(".vspec.arrow")] for f in os.listdir(os.path.join(W, "stage"))
                      if f.endswith(".vspec.arrow"))
        drec = {"setup_s": setup_s, "runs": {}}
        for run in runs:
            store = os.path.join(W, "stage", run + ".vspec.arrow")
            sp = None if a.ms1_only else E.load_store(store)
            reps = []
            fcntl.flock(lockf, fcntl.LOCK_EX)
            try:
                la0 = os.getloadavg()
                for k in range(a.reps + 1):
                    st = {}

                    def tick(name, t, c):
                        st[name] = {"wall_s": time.perf_counter() - t, "core_s": cs() - c}
                        return time.perf_counter(), cs()
                    torch.cuda.synchronize()
                    t, c = time.perf_counter(), cs()
                    fr = find_features(store, FeatureParams(), backend=backend, hill_peaks=False)
                    hyp, _ = build_hypotheses(store, fr)
                    t, c = tick("ms1", t, c)
                    st["ms1"]["sub_wall_s"] = {kk: v for kk, v in fr.timings.items() if isinstance(v, float)}
                    if backend == "gpu":
                        st["ms1"]["gpu_s"] = sum(fr.timings["gpu_ms"].values()) / 1e3
                        st["ms1"]["gpu_ms"] = fr.timings["gpu_ms"]
                    counts = {"n_features": fr.features.num_rows, "n_hypotheses_ms1": hyp.num_rows,
                              "n_hills": int(fr.hills.n), "n_candidates": fr.params["n_candidates"]}
                    if not a.ms1_only:
                        groups = CS.precursor_groups(hyp, sp)
                        tbl, sst = CS.search_chimeric(sp, gi, l1, p, groups)
                        torch.cuda.synchronize()
                        t, c = tick("search", t, c)
                        st["search"]["gpu_s"] = sst["primary"]["gpu_s"] + sst.get("secondary", {}).get("gpu_s", 0.0)
                        psm = {cn: tbl.column(cn).to_numpy(zero_copy_only=False) for cn in tbl.column_names}
                        e0, e1 = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
                        e0.record()
                        feat, info = fb.run(l1, psm, sp, primary=psm["group"] == 0)
                        torch.cuda.synchronize()
                        t, c = tick("features", t, c)
                        s = np.asarray(model.apply(model_psm(psm), feat), np.float64)
                        e1.record()
                        torch.cuda.synchronize()
                        t, c = tick("model", t, c)
                        st["features"]["gpu_stream_s_features_model"] = e0.elapsed_time(e1) / 1e3
                        row = psm["pf_row"]
                        lab = ctx.lab
                        cc = chimeric_competition_fast(psm["spec_idx"].astype(np.int64), psm["group"],
                                                       psm["chim_rank"], psm["pf_id"], lab["tdc_pair"][row], s,
                                                       lab["is_decoy"][row], row, ions, seed=11)
                        t, c = tick("compete", t, c)
                        rr = R()
                        rr.row, rr.score, rr.is_decoy, rr.q = cc["win"], cc["res"].score, cc["res"].is_decoy, cc["q_sep"]
                        res, prot = P.protein_groups(ctx, row, psm["pf_id"], rr)
                        t, c = tick("proteins", t, c)
                        counts.update({"n_psms": int(len(s)), "nce": info["nce"],
                                       "ids_q_sep_1pct": int(((cc["q_sep"] <= 0.01) & ~cc["res"].is_decoy).sum()),
                                       "protein_groups_1pct": prot[0.01]["n_groups"]})
                    stages = [kk for kk in st]
                    st["total"] = {"wall_s": sum(st[kk]["wall_s"] for kk in stages),
                                   "core_s": sum(st[kk]["core_s"] for kk in stages)}
                    st["counts"] = counts
                    print(ds, run[:28], k, {kk: round(v["wall_s"], 3) for kk, v in st.items() if "wall_s" in v},
                          round(st["total"]["core_s"], 2), flush=True)
                    if k > 0:
                        reps.append(st)
                la1 = os.getloadavg()
            finally:
                fcntl.flock(lockf, fcntl.LOCK_UN)
            stages = [kk for kk in reps[0] if kk != "counts"]
            med = {kk: {m: float(np.median([r[kk][m] for r in reps])) for m in reps[0][kk]
                        if isinstance(reps[0][kk][m], float)} for kk in stages}
            med["ms1"]["sub_wall_s"] = {m: float(np.median([r["ms1"]["sub_wall_s"].get(m, 0.0) for r in reps]))
                                        for m in reps[0]["ms1"]["sub_wall_s"]}
            drec["runs"][run] = {"median": med, "counts": reps[0]["counts"], "repeats": reps,
                                 "loadavg_locked": [la0, la1]}
        out["datasets"][ds] = drec
        if not a.ms1_only:
            del gi, l1
        del caches, fb
        torch.cuda.empty_cache()
    allr = [r for d in out["datasets"].values() for r in d["runs"].values()]
    stages = [kk for kk in allr[0]["median"]]
    out["summary_median_over_runs"] = {kk: {m: float(np.median([r["median"][kk][m] for r in allr]))
                                            for m in allr[0]["median"][kk] if m != "sub_wall_s"} for kk in stages}
    out["summary_sum_over_runs"] = {kk: {m: float(np.sum([r["median"][kk][m] for r in allr]))
                                         for m in allr[0]["median"][kk] if m != "sub_wall_s"} for kk in stages}
    out["summary_sum_over_runs"]["ms1"]["sub_wall_s"] = {
        m: float(np.sum([r["median"]["ms1"]["sub_wall_s"].get(m, 0.0) for r in allr]))
        for m in allr[0]["median"]["ms1"]["sub_wall_s"]}
    out["loadavg_end"] = os.getloadavg()
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(out["summary_sum_over_runs"], indent=1))


if __name__ == "__main__":
    main()
