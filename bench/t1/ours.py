"""T1 driver, stage `ours`: Vulcan production path on every replicate database (DESIGN §6.0 A1.1).

1. `python -m vulcan.cohort` (supervised decode workers -> staged Arrow stores -> resident GPU index) ingests the RAW
   files and runs the production search on prod_<first decoy seed>; the stores are kept for the other replicates.
2. The same engine (vulcan.gpu.engine.search, report_psms 5) searches every replicate index with the kept stores; the
   prod_<first> result is compared bit-for-bit with the cohort's output (same code path => must be identical).
3. L3 features (vulcan.rescore.features, AlphaPeptDeep through the content-hashed cache; NCE chosen per run and
   search, label-blind) and the FROZEN static_gbm_pred.pkl (sha256 checked) -> <rep>/<run>.scored.parquet.
Labels are never read here.
"""
from __future__ import annotations

import glob
import json
import os
import resource
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl")
MODEL_SHA256_PREFIX = "31bfe685"           # DESIGN §6.0 A1.1
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
REPORT = 5
PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt"]


def core_s():
    a = resource.getrusage(resource.RUSAGE_SELF)
    b = resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def run_cohort(raws, index_dir, out_dir, stage, workers, log) -> dict:
    stats = os.path.join(out_dir, "cohort_stats.json")
    if os.path.exists(stats):
        return json.load(open(stats))
    os.makedirs(out_dir, exist_ok=True)
    cmd = [sys.executable, "-m", "vulcan.cohort", "--index", index_dir, "--raws", *raws, "--out", out_dir,
           "--stage", stage, "--keep-stores", "--workers", str(workers), "--params", json.dumps({"report_psms": REPORT})]
    log("cohort: " + " ".join(cmd))
    with open(os.path.join(out_dir, "cohort.log"), "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT,
                       env={**os.environ, "PYTHONPATH": ROOT})
    return json.load(open(stats))


def search_replicates(stores, reps: dict, out_root: str, cohort_dir: str, first_prod: str, log) -> dict:
    """Engine search of every replicate index (report_psms 5). Returns per-replicate timing + the cohort identity check."""
    import pyarrow.parquet as pq
    import torch

    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(report_psms=REPORT)
    info = {}
    todo = [r for r in reps if not os.path.exists(os.path.join(out_root, r, "search_done.json"))]
    if not todo:
        return {r: json.load(open(os.path.join(out_root, r, "search_done.json"))) for r in reps}
    sps = {}
    t = time.perf_counter()
    for s in stores:
        sp = E.load_store(s)
        sps[sp.run] = sp
    log(f"loaded {len(sps)} stores in {time.perf_counter() - t:.1f}s")
    for r in todo:
        d = os.path.join(out_root, r)
        os.makedirs(d, exist_ok=True)
        t0, c0 = time.perf_counter(), core_s()
        l1 = L1Index(reps[r]["index"])
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        t_idx = time.perf_counter() - t0
        rec = {"index_build_s": t_idx, "runs": {}}
        for run, sp in sps.items():
            t1 = time.perf_counter()
            tbl, st, _ = E.search(sp, gi, l1, p)
            pq.write_table(tbl, os.path.join(d, f"{run}.psms.parquet"))
            rec["runs"][run] = {"search_s": time.perf_counter() - t1, "n_psms": tbl.num_rows,
                                "gpu_s": sum(st.gpu_ms.values()) / 1e3}
        if r == first_prod:
            same = {}
            for run in sps:
                a = pq.read_table(os.path.join(d, f"{run}.psms.parquet"))
                b = pq.read_table(os.path.join(cohort_dir, f"{run}.psms.parquet"))
                same[run] = bool(a.select(b.column_names).equals(b))
            rec["identical_to_cohort_output"] = same
            log(f"{r}: engine == cohort output: {same}")
        rec["wall_s"] = time.perf_counter() - t0
        rec["core_s"] = core_s() - c0
        del gi, l1
        torch.cuda.empty_cache()
        json.dump(rec, open(os.path.join(d, "search_done.json"), "w"), indent=1)
        info[r] = rec
        log(f"searched {r}: {rec['wall_s']:.1f}s")
    return {r: json.load(open(os.path.join(out_root, r, "search_done.json"))) for r in reps}


class Scorer:
    """Feature builder + frozen model, kept resident (one process for all replicates and the sage+pred comparator)."""

    def __init__(self, predcache: str, instrument: str, dataset: str):
        import hashlib

        from vulcan.predict.apd import APD, check_allowed
        from vulcan.predict.cache import PredCache
        from vulcan.rescore.features import FeatureBuilder
        from vulcan.rescore.static import StaticModel
        check_allowed(dataset)
        sha = hashlib.sha256(open(MODEL, "rb").read()).hexdigest()
        assert sha.startswith(MODEL_SHA256_PREFIX), f"frozen model hash mismatch: {sha}"
        self.model_sha256 = sha
        self.model = StaticModel.load(MODEL)
        self.pred = APD()
        self.caches = {}

        def factory(nce):
            if nce not in self.caches:
                self.caches[nce] = PredCache(predcache, self.pred, nce, instrument, STATIC, VAR)
            return self.caches[nce]
        self.fb = FeatureBuilder(factory, STATIC, VAR)

    def features(self, l1, psm, sp):
        return self.fb.run(l1, psm, sp)


def score_replicates(scorer: Scorer, stores, reps: dict, out_root: str, log) -> dict:
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch

    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    todo = [r for r in reps if not os.path.exists(os.path.join(out_root, r, "score_done.json"))]
    if todo:
        sps = {}
        for s in stores:
            sp = load_store(s)
            sps[sp.run] = sp
    for r in todo:
        d = os.path.join(out_root, r)
        t0, c0 = time.perf_counter(), core_s()
        l1 = L1Index(reps[r]["index"])
        rec = {"runs": {}}
        for run, sp in sps.items():
            t = pq.read_table(os.path.join(d, f"{run}.psms.parquet"))
            psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
            t1 = time.perf_counter()
            feats, finfo = scorer.features(l1, psm, sp)
            t2 = time.perf_counter()
            score = scorer.model.apply(psm, feats).astype(np.float64)
            pq.write_table(pa.table({"run": psm["run"], "spec_idx": psm["spec_idx"], "rank": psm["rank"],
                                     "pf_row": psm["pf_row"], "pf_id": psm["pf_id"], "charge": psm["charge"],
                                     "hyperscore": psm["hyperscore"], "score": score}),
                           os.path.join(d, f"{run}.scored.parquet"))
            rec["runs"][run] = {"features_s": t2 - t1, "apply_s": time.perf_counter() - t2, "nce": finfo["nce"],
                                "rt_mad_min": finfo["rt_mad_min"], "cache_stats": finfo["cache_stats"]}
        rec["wall_s"] = time.perf_counter() - t0
        rec["core_s"] = core_s() - c0
        rec["model_sha256"] = scorer.model_sha256
        json.dump(rec, open(os.path.join(d, "score_done.json"), "w"), indent=1, default=float)
        torch.cuda.empty_cache()
        log(f"scored {r}: {rec['wall_s']:.1f}s nce={[v['nce'] for v in rec['runs'].values()]}")
    return {r: json.load(open(os.path.join(out_root, r, "score_done.json"))) for r in reps}


def stores_of(stage: str, runs) -> list[str]:
    out = [os.path.join(stage, r + ".vspec.arrow") for r in runs]
    missing = [s for s in out if not os.path.exists(s)]
    assert not missing, f"missing kept stores: {missing}"
    return out


def run_names(cohort_dir: str) -> list[str]:
    return sorted(os.path.basename(f)[:-len(".psms.parquet")] for f in glob.glob(os.path.join(cohort_dir, "*.psms.parquet")))
