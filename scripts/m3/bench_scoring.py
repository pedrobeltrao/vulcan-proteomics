#!/usr/bin/env python3
"""Timed M3 scoring stages on D0 (2 runs, prod index, top-5 PSMs), one process, warm-up excluded and reported.

Stages per run: store load, NCE choice, predictor inference (cold = empty cache; warm = precomputed tryptic library),
fragment feature assembly (GPU), RT calibration, static model apply (linear on GPU and CPU; GBM on CPU), RESET.
Wall + CPU core-s per stage; hwinfo fingerprint.

    flock /tmp/vulcan_bench.lock conda_envs/pred/bin/python scripts/m3/bench_scoring.py
"""
from __future__ import annotations

import json
import os
import resource
import shutil
import sys
import time

t_start = time.perf_counter()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402


def cs():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


class T:
    def __init__(self, rec, name):
        self.rec, self.name = rec, name

    def __enter__(self):
        import torch
        torch.cuda.synchronize()
        self.t, self.c = time.perf_counter(), cs()

    def __exit__(self, *a):
        import torch
        torch.cuda.synchronize()
        self.rec[self.name] = {"wall_s": time.perf_counter() - self.t, "core_s": cs() - self.c}


def main():
    import torch
    from vulcan.gpu.engine import load_store
    from vulcan.hw import fingerprint
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder, SpecPeaks
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, StaticModel, feature_matrix
    from vulcan.stats.reset import reset
    t_import = time.perf_counter() - t_start
    out = {"hw": fingerprint(), "python_import_s": t_import, "runs": {}}
    t0 = time.perf_counter()
    pred = APD()
    out["model_load_s"] = time.perf_counter() - t0
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    runs = MD.load_runs(C.psm_dir("D0_prod"))
    lab = MD.labels(C.index_dir("human_contam", "prod"))
    cold_root = os.path.join(C.WORK, "bench_cold_cache")
    lib_root = os.path.join(C.WORK, "predcache")
    models = {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in
              ("static_linear_pred", "static_gbm_pred")}
    # warm-up: kernels for every peptide length, both models (excluded)
    w = np.random.default_rng(0).choice(len(l1), 20000, replace=False)
    t0 = time.perf_counter()
    pred.predict(l1.residues, l1.seq_off[w], l1.length[w], l1.mask[w].view(np.int64), np.full(len(w), 2), 29.0)
    torch.cuda.synchronize()
    out["predictor_warmup_s"] = time.perf_counter() - t0
    for run, psm, _ in runs:
        rec = {}
        with T(rec, "load_store"):
            sp = load_store(C.E2.store_path(next(r for r in C.DEV["D0"]["raw"] if run in r)))
        lib_caches = {}
        for mode, root in (("cold", cold_root), ("warm_library_first", lib_root), ("warm_library_resident", lib_root)):
            if mode == "cold":
                shutil.rmtree(cold_root, ignore_errors=True)
            caches = lib_caches if root == lib_root else {}

            def fac(nce, root=root, caches=caches):
                if nce not in caches:
                    caches[nce] = PredCache(root, pred, nce, "QE", C.STATIC, C.VAR)
                return caches[nce]
            fb = FeatureBuilder(fac, C.STATIC, C.VAR)
            r2 = {}
            with T(r2, "features_total"):
                feat, info = fb.run(l1, psm, sp, nce=29.0)       # NCE fixed = the library NCE
            r2["stages_s"] = info["timing_s"]
            r2["cache_stats"] = info["cache_stats"]
            rec[mode] = r2
        with T(rec, "nce_select_grid"):
            fb.choose_nce(l1, psm, sp, SpecPeaks(sp, fb.dev))
        with T(rec, "feature_matrix"):
            X = feature_matrix(psm, feat, BASE_FEATURES + PRED_FEATURES)
        m = models["static_linear_pred"]
        with T(rec, "apply_linear_cpu"):
            m.score(X)
        with T(rec, "apply_linear_gpu"):
            Xt = torch.from_numpy(X).cuda()
            s = ((Xt - torch.from_numpy(m.mu).cuda()) / torch.from_numpy(m.sd).cuda()) @ torch.from_numpy(m.model["w"]).cuda()
            s.cpu()
        with T(rec, "apply_gbm_cpu"):
            models["static_gbm_pred"].score(X)
        # L4: paired double competition on hyperscore, then RESET (LDA) on its winners
        rows = psm["pf_row"]
        with T(rec, "double_competition"):
            dc = double_competition(psm["spec_idx"], psm["pf_id"], lab["tdc_pair"][rows], psm["hyperscore"],
                                    lab["is_decoy"][rows], seed=11, tiebreak=hash_tiebreak(psm["pf_id"], 11))
        with T(rec, "reset_lda"):
            reset(X[dc.row], dc.is_decoy, init_score=dc.score, seed=11)
        rec["n_psms"] = int(len(psm["rank"]))
        rec["n_spectra"] = int(len(sp))
        out["runs"][run] = rec
        print(run, json.dumps({k: (v if not isinstance(v, dict) else v.get("wall_s", v)) for k, v in rec.items()},
                              default=float)[:1500], flush=True)
    out["process_wall_s"] = time.perf_counter() - t_start
    out["process_core_s"] = cs()
    json.dump(out, open(os.path.join(C.RES, "bench_scoring.json"), "w"), indent=1, default=float)
    sys.stdout.flush()
    os._exit(0)          # interpreter shutdown hung once (torch/peptdeep teardown) while holding the bench lock


if __name__ == "__main__":
    main()
