#!/usr/bin/env python3
"""E7b per-run cost of the chimeric v2_sep production path, before vs after (dev only).

One persistent process per invocation; per run: one warm-up, then --reps measured repeats; each run's repeats are
timed while holding `flock /tmp/vulcan_bench.lock` (taken per run, so other agents are blocked only briefly).
Stages (per run):
  ms1       vulcan.ms1 features + precursor-hypothesis table (in-process, --threads numba threads)
            before: exact_reference isotope kernels; after: fast exact path (bitwise equal, equality_ms1.json)
  search    chimeric search (closed primary + secondary groups; unchanged code in both)
  features  v1 feature set with primary-only calibrations
            before: FeatureBuilderV2.run_v2 + PredCache (the E7 path; also computes the 19 unused v2 extras)
            after:  FastFeatureBuilder.run(primary=) + LibCache mirrored on the GPU (equality_scoring.json)
  model     frozen v2 linear model (same in both)
  compete   per-run chimeric competition + per-stratum TDC (q_sep)
            before: chimeric_competition (incl. Group-walk); after: chimeric_competition_fast (no Group-walk)
  proteins  picked protein groups (vulcan.rescore.post.protein_groups, unchanged estimator)
Wall and process core-s (RUSAGE_SELF, all threads) per stage; GPU-s of the search kernels (engine CUDA events);
CUDA-event stream time for features+model (upper bound on their GPU-s).  Index: the dataset's T1 prod_1001.
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
    ap.add_argument("--impl", choices=["before", "after"], required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1", "D8h"])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--work-root", default=f"/tmp/vulcan_t1_{os.getuid()}")
    ap.add_argument("--lock", default="/tmp/vulcan_bench.lock")
    a = ap.parse_args()
    from vulcan import hw
    n_thr = hw.pin_thread_pools(a.threads)
    import torch

    from bench.t1.ours_chimeric import MODEL_V2
    from vulcan.chimeric import search as CS
    from vulcan.chimeric.compete import IonSets, chimeric_competition
    from vulcan.chimeric.fast_compete import chimeric_competition_fast
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ms1.features import FeatureParams, find_features
    from vulcan.ms1.hypotheses import build_hypotheses
    from vulcan.predict.apd import APD
    from vulcan.rescore import post as P
    from vulcan.rescore.features_v2 import FeatureBuilderV2, model_psm
    from vulcan.rescore.static import StaticModel
    before = a.impl == "before"
    v2 = StaticModel.load(MODEL_V2)
    pred = APD()
    p = E.SearchParams(report_psms=5)
    out = {"impl": a.impl, "threads": n_thr, "hw": hw.fingerprint(), "reps": a.reps, "datasets": {},
           "loadavg_start": os.getloadavg()}
    lockf = open(a.lock, "a")
    for ds in a.datasets:
        W = os.path.join(a.work_root, ds)
        lib = os.path.join(ROOT, "data", "cache", "e7b_predlib", ds)
        caches = {}
        if before:
            from vulcan.predict.cache import PredCache

            def factory(nce):
                if nce not in caches:
                    caches[nce] = PredCache(lib, pred, nce, "QE", STATIC, VAR)
                return caches[nce]
            fb = FeatureBuilderV2(factory, STATIC, VAR)
            model = v2
        else:
            from vulcan.predict.library import LibCache
            from vulcan.rescore.fast import FastFeatureBuilder, FastStatic

            def factory(nce):
                if nce not in caches:
                    caches[nce] = LibCache(lib, pred, nce, "QE", STATIC, VAR, device="cuda")
                return caches[nce]
            fb = FastFeatureBuilder(factory, STATIC, VAR)
            model = FastStatic(v2)
        t0 = time.time()
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
            sp = E.load_store(store)
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
                    fr = find_features(store, FeatureParams(), exact_reference=before)
                    hyp, _ = build_hypotheses(store, fr)
                    t, c = tick("ms1", t, c)
                    groups = CS.precursor_groups(hyp, sp)
                    tbl, sst = CS.search_chimeric(sp, gi, l1, p, groups)
                    torch.cuda.synchronize()
                    t, c = tick("search", t, c)
                    st["search"]["gpu_s"] = sst["primary"]["gpu_s"] + sst.get("secondary", {}).get("gpu_s", 0.0)
                    psm = {cn: tbl.column(cn).to_numpy(zero_copy_only=False) for cn in tbl.column_names}
                    e0, e1 = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
                    e0.record()
                    if before:
                        ff = fr.features
                        mf = {cn: ff.column(cn).to_numpy() for cn in ("trace_cos", "n_isotopes", "rt_apex",
                                                                      "rt_start", "rt_end", "intensity_apex")}
                        feat, info = fb.run_v2(l1, psm, sp, mf)
                    else:
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
                    args = (psm["spec_idx"].astype(np.int64), psm["group"], psm["chim_rank"], psm["pf_id"],
                            lab["tdc_pair"][row], s, lab["is_decoy"][row], row)
                    if before:
                        cc = chimeric_competition(*args, IonSets(l1), seed=11)
                    else:
                        cc = chimeric_competition_fast(*args, ions, seed=11)
                    t, c = tick("compete", t, c)
                    rr = R()
                    rr.row, rr.score, rr.is_decoy, rr.q = cc["win"], cc["res"].score, cc["res"].is_decoy, cc["q_sep"]
                    res, prot = P.protein_groups(ctx, row, psm["pf_id"], rr)
                    t, c = tick("proteins", t, c)
                    tot = {"wall_s": sum(v["wall_s"] for v in st.values()),
                           "core_s": sum(v["core_s"] for v in st.values())}
                    st["total"] = tot
                    st["counts"] = {"n_features": fr.features.num_rows, "n_hypotheses_ms1": hyp.num_rows,
                                    "n_psms": int(len(s)), "nce": info["nce"],
                                    "ids_q_sep_1pct": int(((cc["q_sep"] <= 0.01) & ~cc["res"].is_decoy).sum()),
                                    "protein_groups_1pct": prot[0.01]["n_groups"]}
                    print(ds, run[:28], k, {kk: round(v["wall_s"], 3) for kk, v in st.items() if "wall_s" in v},
                          round(tot["core_s"], 1), flush=True)
                    if k > 0:
                        reps.append(st)
                la1 = os.getloadavg()
            finally:
                fcntl.flock(lockf, fcntl.LOCK_UN)
            stages = [k for k in reps[0] if k != "counts"]
            med = {k: {m: float(np.median([r[k][m] for r in reps])) for m in reps[0][k]} for k in stages}
            drec["runs"][run] = {"median": med, "counts": reps[0]["counts"], "repeats": reps,
                                 "loadavg_locked": [la0, la1]}
        out["datasets"][ds] = drec
        del gi, l1, caches, fb
        torch.cuda.empty_cache()
    allr = [r for d in out["datasets"].values() for r in d["runs"].values()]
    stages = list(allr[0]["median"])
    out["summary_median_over_runs"] = {k: {m: float(np.median([r["median"][k][m] for r in allr]))
                                           for m in allr[0]["median"][k]} for k in stages}
    out["summary_sum_over_runs"] = {k: {m: float(np.sum([r["median"][k][m] for r in allr]))
                                        for m in allr[0]["median"][k]} for k in stages}
    out["loadavg_end"] = os.getloadavg()
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(out["summary_median_over_runs"], indent=1))


if __name__ == "__main__":
    main()
