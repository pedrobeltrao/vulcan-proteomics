"""Production scoring inside the persistent cohort process (vulcan.cohort --score-lib ...).

GPU process (cohort GPU-stage thread, right after the search of a run): L3 features with the per-run label-blind
NCE choice from the resident predicted library + frozen static_gbm_pred on the GPU (vulcan.rescore.prod).
Host stage (paired double competition + picked protein groups, vulcan.rescore.post) runs in a small pool of
torch-free worker processes, so it overlaps the next runs' search/scoring and the decode workers.

Per run, OUT/<run>.psms.parquet gains a `score` column (frozen static_gbm_pred), and the post workers write
OUT/<run>.peptides.parquet (pair-competition winners: pf_id, spec_idx, score, is_decoy, q) and
OUT/<run>.proteins.parquet (picked protein groups: proteins, score, is_decoy, q, is_obsolete).
"""
from __future__ import annotations

import os
import resource
import time

import numpy as np

PSM_NEED = ("spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt")

_CTX = {}
REPORT_PSMS = 5     # the frozen scorer is defined on top-5 PSM tables (M3: scripts/m3/m3_common.py REPORT)


def _core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


CHIM_NEED = ("group", "chim_rank")


def _post_init(index_dir, chimeric=False):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from vulcan.rescore.post import IndexCtx
    _CTX["ctx"] = IndexCtx(index_dir)
    if chimeric:
        from vulcan.chimeric.compete import IonSets
        from vulcan.index.build import L1Index
        _CTX["ions"] = IonSets(L1Index(index_dir))


class _Res:
    pass


def _post_job_chim(run, pf_row, pf_id, spec_idx, group, chim_rank, score, out_dir, write=True):
    """E7b chimeric post stage: chimeric competition (stage-1 group winners, CONGA neighbour filter, paired double
    competition, rank strata, per-stratum TDC = q_sep; vulcan.chimeric.fast_compete) + picked protein groups."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    from vulcan.chimeric.fast_compete import chimeric_competition_fast
    from vulcan.rescore import post as P
    t0, c0 = time.time(), _core_s()
    ctx = _CTX["ctx"]
    lab = ctx.lab
    cc = chimeric_competition_fast(np.asarray(spec_idx, np.int64), group, chim_rank, pf_id, lab["tdc_pair"][pf_row],
                                   score, lab["is_decoy"][pf_row], pf_row, _CTX["ions"], seed=11)
    r = _Res()
    r.row, r.score, r.is_decoy, r.q = cc["win"], cc["res"].score, cc["res"].is_decoy, cc["q_sep"]
    t1 = time.time()
    res, prot = P.protein_groups(ctx, pf_row, pf_id, r)
    t2 = time.time()
    if write:
        pq.write_table(pa.table({"pf_id": pf_id[r.row], "spec_idx": np.asarray(spec_idx)[r.row], "score": r.score,
                                 "is_decoy": r.is_decoy, "q": r.q, "stratum": cc["stratum"],
                                 "chim_rank": np.asarray(chim_rank)[r.row]}),
                       os.path.join(out_dir, f"{run}.peptides.parquet"))
        pq.write_table(pa.table({"proteins": [";".join(g) for g in res.groups], "score": res.score,
                                 "is_decoy": res.is_decoy, "q": res.q, "is_obsolete": res.is_obsolete}),
                       os.path.join(out_dir, f"{run}.proteins.parquet"))
    return {"run": run, "ids_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum()), "proteins_1pct": prot[0.01]["n_groups"],
            "compete_s": t1 - t0, "protein_s": t2 - t1, "post_wall_s": time.time() - t0, "post_core_s": _core_s() - c0,
            "t_end": time.time()}


def _post_job(run, pf_row, pf_id, spec_idx, score, out_dir, write=True):
    import pyarrow as pa
    import pyarrow.parquet as pq

    from vulcan.rescore import post as P
    t0, c0 = time.time(), _core_s()
    ctx = _CTX["ctx"]
    r = P.compete(ctx, pf_row, pf_id, spec_idx, score)
    t1 = time.time()
    res, prot = P.protein_groups(ctx, pf_row, pf_id, r)
    t2 = time.time()
    if write:
        pq.write_table(pa.table({"pf_id": pf_id[r.row], "spec_idx": np.asarray(spec_idx)[r.row], "score": r.score,
                                 "is_decoy": r.is_decoy, "q": r.q}), os.path.join(out_dir, f"{run}.peptides.parquet"))
        pq.write_table(pa.table({"proteins": [";".join(g) for g in res.groups], "score": res.score,
                                 "is_decoy": res.is_decoy, "q": res.q, "is_obsolete": res.is_obsolete}),
                       os.path.join(out_dir, f"{run}.proteins.parquet"))
    return {"run": run, "ids_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum()), "proteins_1pct": prot[0.01]["n_groups"],
            "compete_s": t1 - t0, "protein_s": t2 - t1, "post_wall_s": time.time() - t0, "post_core_s": _core_s() - c0,
            "t_end": time.time()}


class CohortScorer:
    """Callable hook for vulcan.cohort.run_cohort(post_search=...): (rec, tbl, sp) -> (tbl + score, future)."""

    def __init__(self, l1, index_dir, out_dir, lib_root, model_path, post_workers=2, instrument="QE", log=print,
                 write=True, impl="fast", placement="gpu", chimeric=False):
        """impl='fast' (default): vulcan.rescore.prod + post workers.  impl='frozen' (comparison only): the frozen
        M3/M3e components unchanged (FeatureBuilder + PredCache on `lib_root`, sklearn GBM on the host, competition
        and proteins inline in the GPU-stage thread), i.e. the pre-optimization cost inside the cohort."""
        from concurrent.futures import ProcessPoolExecutor
        import multiprocessing as mp
        t0 = time.perf_counter()
        # start the torch-free post workers first (they load labels + protein map while the predictor loads)
        self.pool = None
        if impl != "frozen":
            self.pool = ProcessPoolExecutor(max(1, post_workers), mp_context=mp.get_context("spawn"),
                                            initializer=_post_init, initargs=(index_dir, chimeric))
            self._prime = [self.pool.submit(os.getpid) for _ in range(max(1, post_workers))]
        self.impl = impl
        self.chimeric = chimeric
        if chimeric and impl != "fast":
            raise ValueError("chimeric scoring is implemented for impl='fast' only")
        if impl == "frozen":
            self.sc = _FrozenScorer(lib_root, model_path, instrument, index_dir)
        else:
            from vulcan.rescore.prod import ProdScorer
            self.sc = ProdScorer(lib_root, model_path, instrument, log=log, placement=placement)
        t1 = time.perf_counter()
        self.warm_s = self.sc.warm()
        self.l1, self.out_dir, self.write = l1, out_dir, write
        self.setup = {"predictor_load_s": self.sc.t_predictor_load, "model_load_s": self.sc.t_model_load,
                      "library_warm_s": self.warm_s, "library_bytes": getattr(self.sc, "library_bytes", None), "scorer_setup_s": time.perf_counter() - t0,
                      "post_workers": max(1, post_workers), "lib_root": lib_root, "model": model_path, "impl": impl,
                      "placement": placement if impl != "frozen" else "PredCache"}
        self.setup["scorer_setup_excl_warm_s"] = t1 - t0
        self.results = {}
        log(f"[cohort] scorer ready in {self.setup['scorer_setup_s']:.1f}s (library warm {self.warm_s:.1f}s)")

    def __call__(self, rec, tbl, sp):
        import pyarrow as pa
        t0, c0 = time.time(), _core_s()
        if self.chimeric:
            from vulcan.rescore.features_v2 import model_psm
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in PSM_NEED + CHIM_NEED}
            t1 = time.perf_counter()
            feat, info = self.sc.fb.run(self.l1, psm, sp, primary=psm["group"] == 0)
            t2 = time.perf_counter()
            s = np.asarray(self.sc.model.apply(model_psm(psm), feat), np.float64)
            self.sc.torch.cuda.current_stream().synchronize()
            info["timing_s"]["gbm"] = time.perf_counter() - t2
            info["timing_s"]["features_total"] = t2 - t1
            tbl = tbl.append_column("score", pa.array(s))
            fut = self.pool.submit(_post_job_chim, rec["run"], psm["pf_row"], psm["pf_id"], psm["spec_idx"],
                                   psm["group"], psm["chim_rank"], s, self.out_dir, self.write)
        else:
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in PSM_NEED}
            s, feat, info = self.sc.score(self.l1, psm, sp)
            tbl = tbl.append_column("score", pa.array(s))
        if self.chimeric:
            pass                                        # post job submitted above
        elif self.impl == "frozen":                     # inline, as the frozen bench path
            fut = _Done(self.sc.post(rec["run"], psm, s))
        else:
            fut = self.pool.submit(_post_job, rec["run"], psm["pf_row"], psm["pf_id"], psm["spec_idx"], s,
                                   self.out_dir, self.write)
        rec["score_gpu_stage_wall_s"] = time.time() - t0
        rec["score_gpu_stage_thread_core_s"] = _core_s() - c0      # process-wide (other threads included)
        rec["score_nce"] = info["nce"]
        rec["score_cache"] = {k: info["cache_stats"][k] for k in ("ms2_miss", "rt_miss", "n_predicted_ms2",
                                                                 "predict_s")}
        rec["score_timing_s"] = info["timing_s"]
        rec["_post_future"] = fut
        return tbl, fut

    def collect(self, runs: dict):
        """Attach post-stage results to run records; returns a summary."""
        for name, rec in runs.items():
            fut = rec.pop("_post_future", None)
            if fut is not None:
                rec["post"] = fut.result()
        post = [r["post"] for r in runs.values() if "post" in r]
        g = [r.get("score_gpu_stage_wall_s", 0.0) for r in runs.values()]
        return {"setup": self.setup, "n_scored": len(post),
                "gpu_stage_wall_s_median": float(np.median(g)) if g else None,
                "gpu_stage_wall_s_sum": float(np.sum(g)),
                "post_wall_s_median": float(np.median([p["post_wall_s"] for p in post])) if post else None,
                "post_core_s_sum": float(sum(p["post_core_s"] for p in post)),
                "ids_1pct_sum": int(sum(p["ids_1pct"] for p in post)),
                "proteins_1pct_sum": int(sum(p["proteins_1pct"] for p in post)),
                "last_post_end": max((p["t_end"] for p in post), default=None)}

    def shutdown(self):
        if self.pool is not None:
            self.pool.shutdown(wait=True)


class _Done:
    def __init__(self, v):
        self.v = v

    def result(self):
        return self.v


class _FrozenScorer:
    """The frozen pre-optimization components (scripts/m3e/bench_prod.py path), for the before/after cohort."""

    def __init__(self, lib_root, model_path, instrument, index_dir):
        import sys
        from vulcan.predict.apd import APD
        from vulcan.predict.cache import PredCache
        from vulcan.rescore.features import FeatureBuilder
        from vulcan.rescore.static import StaticModel
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        sys.path.insert(0, os.path.join(root, "scripts", "m3"))
        import evaluate as EV
        import m3_data as MD
        self.EV, self.index_dir, self.lab = EV, index_dir, MD.labels(index_dir)
        t0 = time.perf_counter()
        self.pred = APD()
        self.t_predictor_load = time.perf_counter() - t0
        caches = {}
        st, vr = {"C": 57.021464}, {"M": 15.9949}

        def factory(nce):
            if nce not in caches:
                caches[nce] = PredCache(lib_root, self.pred, nce, instrument, st, vr)
            return caches[nce]
        self.factory = factory
        self.fb = FeatureBuilder(factory, st, vr)
        t1 = time.perf_counter()
        self.model = StaticModel.load(model_path)
        self.t_model_load = time.perf_counter() - t1

    def warm(self):
        t0 = time.perf_counter()
        c = self.factory(29.0)
        for k in ("rt", "ms2_z2", "ms2_z3"):
            c._tab(k)
        return time.perf_counter() - t0

    def score(self, l1, psm, sp):
        feat, info = self.fb.run(l1, psm, sp)
        t1 = time.perf_counter()
        s = self.model.apply(psm, feat)
        info["timing_s"]["gbm"] = time.perf_counter() - t1
        return s, feat, info

    def post(self, run, psm, s):
        from vulcan.stats._util import hash_tiebreak
        from vulcan.stats.tdc import double_competition
        t0, c0 = time.time(), _core_s()
        lab = self.lab
        r = double_competition(psm["spec_idx"].astype(np.int64), psm["pf_id"], lab["tdc_pair"][psm["pf_row"]], s,
                               lab["is_decoy"][psm["pf_row"]], seed=11, tiebreak=hash_tiebreak(psm["pf_id"], 11))
        t1 = time.time()
        srch = type("S", (), {})()
        srch.psm, srch.lab, srch.index, srch.unit = psm, lab, self.index_dir, psm["pf_id"]
        prot = self.EV.protein_level(srch, {"row": r.row, "score": r.score, "is_decoy": r.is_decoy, "q": r.q}, (0.01,))
        t2 = time.time()
        return {"run": run, "ids_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum()),
                "proteins_1pct": prot[0.01]["n_groups"], "compete_s": t1 - t0, "protein_s": t2 - t1,
                "post_wall_s": t2 - t0, "post_core_s": _core_s() - c0, "t_end": time.time()}
