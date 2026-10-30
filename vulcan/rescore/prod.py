"""Production scorer for a persistent process (cohort GPU process): frozen static_gbm_pred over L3 predicted
features, then paired double competition + picked protein groups.  Same outputs as the frozen M3/M3e path
(scripts/m3e/bench_prod.py); faster components in vulcan.rescore.fast / vulcan.predict.library / vulcan.rescore.post.

    sc = ProdScorer(lib_root, model_path)            # one-time: predictor + model load, library mmap
    s, feat, info = sc.score(l1, psm, sp)            # GPU process: features (+ per-run NCE choice) + GBM
    r, prot, summ = post.post_stage(ctx, ...)        # host (any process): competition + proteins
"""
from __future__ import annotations

import time

import numpy as np

from vulcan.rescore import post as P

STATIC = {"C": 57.021464}
VAR = {"M": 15.9949}


class ProdScorer:
    def __init__(self, lib_root: str, model_path: str, instrument: str = "QE", static_mods=None, var_mods=None,
                 dev: str = "cuda", log=print, placement: str = "gpu"):
        """placement of the predicted library: 'gpu' (mirrored in device memory, lookups on the GPU; host keeps
        only the memory map), 'host' (read into process memory), 'mmap' (page cache only)."""
        if placement not in ("gpu", "host", "mmap"):
            raise ValueError(placement)
        self.placement = placement
        import torch
        from vulcan.predict.apd import APD
        from vulcan.predict.library import LibCache
        from vulcan.rescore.fast import FastFeatureBuilder, FastStatic
        from vulcan.rescore.static import StaticModel
        self.torch = torch
        self.static, self.var = static_mods or STATIC, var_mods or VAR
        self.instrument = instrument
        t0 = time.perf_counter()
        self.pred = APD(device=dev)
        self.t_predictor_load = time.perf_counter() - t0
        self.caches = {}
        self.lib_root = lib_root

        def factory(nce):
            if nce not in self.caches:
                self.caches[nce] = LibCache(lib_root, self.pred, nce, instrument, self.static, self.var,
                                            resident=placement == "host",
                                            device=dev if placement == "gpu" else None)
            return self.caches[nce]
        self.factory = factory
        self.fb = FastFeatureBuilder(factory, self.static, self.var, dev=dev)
        t1 = time.perf_counter()
        self.model_static = StaticModel.load(model_path)
        self.model = FastStatic(self.model_static, dev)
        self.t_model_load = time.perf_counter() - t1
        self.ctx = {}
        self.log = log

    def warm(self, nces=None, kinds=("rt", "ms2_z2", "ms2_z3", "ms2_z4")):
        """Open the library tables of every grid NCE (resident: read into memory); returns seconds."""
        from vulcan.rescore.features import NCE_GRID
        t0 = time.perf_counter()
        self.library_bytes = sum(self.factory(nce).warm(kinds) for nce in (nces or NCE_GRID))
        return time.perf_counter() - t0

    def index_ctx(self, index_dir) -> P.IndexCtx:
        if index_dir not in self.ctx:
            self.ctx[index_dir] = P.IndexCtx(index_dir)
        return self.ctx[index_dir]

    def score(self, l1, psm: dict, sp, nce=None):
        """-> (score float64 [n], features dict, info).  info['timing_s'] has the sub-stages."""
        torch = self.torch
        t0 = time.perf_counter()
        feat, info = self.fb.run(l1, psm, sp, nce=nce)
        t1 = time.perf_counter()
        s = self.model.apply(psm, feat)
        torch.cuda.current_stream().synchronize()          # this thread's stream only (cohort GPU threads)
        info["timing_s"]["gbm"] = time.perf_counter() - t1
        info["timing_s"]["features_total"] = t1 - t0
        return np.asarray(s, np.float64), feat, info
