#!/usr/bin/env python3
"""Production-configuration speed of the default scorer (review M6; A1 hardening).  MEASURED, not derived.

One persistent process (pred env), the production path for every run:
  load store -> GPU closed search (prod index, report_psms 5) -> L3 features with the per-run label-blind NCE choice
  (grid {20..35}, the production setting; predictions come from the cache at the CHOSEN NCE, misses predicted lazily
  and appended) -> frozen static_gbm_pred (CPU sklearn) -> paired double competition + q -> picked protein groups.
Cache state = production "first contact": a fresh cache root holding only the shipped human tryptic library
(NCE 29, QE; the M3 library dir copied), every other NCE starts empty and fills lazily as runs are processed (as a
persistent production worker would).  Runs: D0 run 1 as warm-up (excluded; first-call kernel set-up), then the 10 D8
holdout runs and D0 run 2, each timed (wall via perf_counter with cuda.synchronize, CPU core-s via getrusage SELF =
all threads).  One-time costs reported separately: Python imports, model load, GPU index upload, library load.
Index build (L1, once per proteome) is not per-run; its recorded build time is quoted from the index meta.
Run under `flock /tmp/vulcan_bench.lock`.
    conda_envs/pred/bin/python scripts/m3e/bench_prod.py --out results/.../bench_prod.json
"""
from __future__ import annotations

import time

T_START = time.perf_counter()
import argparse  # noqa: E402
import glob  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import resource  # noqa: E402
import shutil  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
LIB_KEY = "a9d5af1427616250"      # QE, NCE 29 human tryptic library (results/2026-10-27_M3_scoring/predictor_library.json)


def cpu():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


class Clock:
    def __init__(self, torch):
        self.torch, self.rec = torch, {}

    def __call__(self, name):
        clock = self

        class _C:
            def __enter__(s):
                clock.torch.cuda.synchronize()
                s.w, s.c = time.perf_counter(), cpu()

            def __exit__(s, *a):
                clock.torch.cuda.synchronize()
                clock.rec[name] = {"wall_s": time.perf_counter() - s.w, "cpu_core_s": cpu() - s.c}
        return _C()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    c0 = cpu()
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    from vulcan.rescore.static import StaticModel
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    import evaluate as EV
    import m3_data as MD
    one = {"python_imports": {"wall_s": time.perf_counter() - T_START, "cpu_core_s": cpu() - c0}}
    clk = Clock(torch)
    root = os.path.join(W, "predcache_bench")
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root)
    shutil.copytree(os.path.join(C.WORK, "predcache", LIB_KEY), os.path.join(root, LIB_KEY))
    with clk("model_load"):
        pred = APD()
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PredCache(root, pred, nce, "QE", C.STATIC, C.VAR)
        return caches[nce]
    idx = C.index_dir("human_contam", "prod")
    p = E.SearchParams(report_psms=C.REPORT)
    with clk("index_load_upload"):
        l1 = L1Index(idx)
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        lab = MD.labels(idx)
    with clk("library_load_nce29"):
        cache29 = factory(29.0)
        for kind in ("rt", "ms2_z2", "ms2_z3"):
            cache29._tab(kind)
    one.update(clk.rec)
    fb = FeatureBuilder(factory, C.STATIC, C.VAR)
    model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
    d0 = [C.E2.store_path(r) for r in C.DEV["D0"]["raw"]]
    d8 = sorted(glob.glob(os.path.join(W, "store_d8", "*.vspec.arrow")))
    plan = [("warmup", d0[0])] + [("D8", f) for f in d8] + [("D0", d0[1])]
    runs = {}
    for tag, path in plan:
        clk.rec = {}
        t0, c_0 = time.perf_counter(), cpu()
        with clk("load_store"):
            sp = E.load_store(path)
        with clk("search"):
            tbl, st, _ = E.search(sp, gi, l1, p)
        with clk("psm_to_numpy"):
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        with clk("features_incl_nce_choice"):
            feat, info = fb.run(l1, psm, sp)
        with clk("gbm_apply"):
            s = model.apply(psm, feat)
        with clk("double_competition"):
            spec = psm["spec_idx"].astype(np.int64)
            r = double_competition(spec, psm["pf_id"], lab["tdc_pair"][psm["pf_row"]], s, lab["is_decoy"][psm["pf_row"]],
                                   seed=11, tiebreak=hash_tiebreak(psm["pf_id"], 11))
        with clk("picked_protein"):
            srch = type("S", (), {})()
            srch.psm, srch.lab, srch.index, srch.unit = psm, lab, idx, psm["pf_id"]
            prot = EV.protein_level(srch, {"row": r.row, "score": r.score, "is_decoy": r.is_decoy, "q": r.q}, (0.01,))
        tot = {"wall_s": time.perf_counter() - t0, "cpu_core_s": cpu() - c_0}
        scoring = {k: sum(clk.rec[x][k] for x in ("features_incl_nce_choice", "gbm_apply", "double_competition",
                                                    "picked_protein")) for k in ("wall_s", "cpu_core_s")}
        name = os.path.basename(path)[:-len(".vspec.arrow")]
        runs[name] = {"tag": tag, "n_ms2": len(sp), "n_psms": tbl.num_rows, "nce": info["nce"],
                      "cache_stats": info.get("cache_stats"), "feature_timing_s": info.get("timing_s"),
                      "stages": dict(clk.rec), "scoring_total": scoring, "end_to_end_total": tot,
                      "ids_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum()), "proteins_1pct": prot[0.01]["n_groups"]}
        print(tag, name, len(sp), "NCE", info["nce"], "scoring", {k: round(v, 2) for k, v in scoring.items()},
              "e2e", {k: round(v, 2) for k, v in tot.items()}, flush=True)
    meas = [v for v in runs.values() if v["tag"] != "warmup"]
    summ = {}
    for key in ("scoring_total", "end_to_end_total"):
        for m in ("wall_s", "cpu_core_s"):
            x = np.array([v[key][m] for v in meas])
            summ[f"{key}.{m}"] = {"median": float(np.median(x)), "min": float(x.min()), "max": float(x.max()),
                                  "mean": float(x.mean())}
            x1e5 = x / np.array([v["n_ms2"] for v in meas]) * 1e5
            summ[f"{key}.{m}_per_1e5_ms2"] = {"median": float(np.median(x1e5)), "min": float(x1e5.min()),
                                              "max": float(x1e5.max())}
    from vulcan.hw import fingerprint
    meta = json.load(open(os.path.join(idx, "meta.json")))
    out = {"one_time": one, "runs": runs, "summary_measured_runs": summ, "n_measured": len(meas),
           "index_build_recorded_s": meta.get("timing", meta.get("build_s")), "hw": fingerprint(),
           "loadavg_end": os.getloadavg(), "cache_root": root, "library_key": LIB_KEY}
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
