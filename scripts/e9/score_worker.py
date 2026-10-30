#!/usr/bin/env python3
"""E9 persistent scorer (pred env): L3 features + frozen default scorer for every finished stream batch.

    conda_envs/pred/bin/python scripts/e9/score_worker.py [--work /tmp/vulcan_e9] [--exit-when-idle 0]

Polls work/queue/bNNN.json (written by scripts/e9/stream.py).  For each run of the batch and each searched arm
(P = T+E_2001+decoys, C = T+Arabidopsis+decoys): loads the kept Arrow store, builds the label-blind L3 features
(vulcan.rescore.features, AlphaPeptDeep via the content-hashed cache; NCE chosen per run, label-blind), applies the
frozen static_gbm_pred model (D7) and writes work/out/<arm>/<run>.scored.parquet (spec_idx, rank, pf_row, pf_id,
charge, hyperscore, score, rt) plus <run>.feat.parquet; the PSM parquet is moved next to them.  Then the stores are
deleted and bNNN.done written with per-run timings.  Labels are never read here.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys
import time

t_start = time.perf_counter()
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

MODEL = os.path.join(ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl")
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/vulcan_e9")
    ap.add_argument("--index-p", default="/tmp/vulcan_e2/index/human_contam_E2001")
    ap.add_argument("--index-c", default="/tmp/vulcan_m3/index/human_contam_C")
    ap.add_argument("--predcache", default="/tmp/vulcan_e9/predcache")
    ap.add_argument("--exit-when-idle", type=float, default=0, help="exit after this many idle seconds (0 = never)")
    a = ap.parse_args()
    import torch

    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    from vulcan.rescore.static import StaticModel
    t_import = time.perf_counter() - t_start
    t0 = time.perf_counter()
    pred = APD()
    model = StaticModel.load(MODEL)
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PredCache(a.predcache, pred, nce, "QE", STATIC, VAR)
        return caches[nce]

    fb = FeatureBuilder(factory, STATIC, VAR, chunk=50_000)   # smaller chunks: the GPU is shared with the cohort
    l1 = {"P": L1Index(a.index_p), "C": L1Index(a.index_c)}
    t_setup = time.perf_counter() - t0
    q = os.path.join(a.work, "queue")
    with open(os.path.join(a.work, "timeline.jsonl"), "a") as fh:
        fh.write(json.dumps({"ev": "scorer_start", "t": time.time(), "import_s": t_import, "setup_s": t_setup,
                             "model": MODEL, "model_sha256": model.meta.get("sha256")}) + "\n")
    idle0 = time.time()
    while True:
        todo = sorted(f for f in glob.glob(os.path.join(q, "b*.json"))
                      if not os.path.exists(f[:-5] + ".done") and not os.path.exists(f[:-5] + ".claim"))
        claimed = None
        for f in todo:          # several workers may run: claim a batch atomically
            try:
                os.close(os.open(f[:-5] + ".claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                claimed = f
                break
            except FileExistsError:
                continue
        if claimed is None:
            if a.exit_when_idle and time.time() - idle0 > a.exit_when_idle:
                return
            time.sleep(2)
            continue
        todo = [claimed]
        job = json.load(open(todo[0]))
        bt0 = time.time()
        recs = []
        for store in job["stores"]:
            run = os.path.basename(store)[: -len(".vspec.arrow")]
            rec = {"run": run}
            if not os.path.exists(store):       # already scored by an interrupted worker (outputs kept) or lost
                rec["skipped"] = "store gone" + (" (scored)" if os.path.exists(
                    os.path.join(a.work, "out", "P", run + ".scored.parquet")) else " (NOT scored)")
                recs.append(rec)
                continue
            try:
                t1 = time.perf_counter()
                sp = load_store(store)
                rec["load_store_s"] = time.perf_counter() - t1
                rec["n_spectra"] = int(len(sp))
                for arm in ("P", "C"):
                    src = os.path.join(job["dir"], arm, run + ".psms.parquet")
                    if not os.path.exists(src):
                        rec[arm] = "missing"
                        continue
                    t2 = time.perf_counter()
                    t = pq.read_table(src)
                    psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
                    feats, info = fb.run(l1[arm], psm, sp)
                    t3 = time.perf_counter()
                    score = model.apply(psm, feats).astype(np.float32)
                    t_apply = time.perf_counter() - t3
                    od = os.path.join(a.work, "out", arm)
                    os.makedirs(od, exist_ok=True)
                    pq.write_table(pa.table({"spec_idx": psm["spec_idx"], "rank": psm["rank"],
                                             "pf_row": psm["pf_row"], "pf_id": psm["pf_id"],
                                             "charge": psm["charge"], "hyperscore": psm["hyperscore"],
                                             "score": score, "rt": psm["rt"]}),
                                   os.path.join(od, run + ".scored.parquet"))
                    pq.write_table(pa.table(feats), os.path.join(od, run + ".feat.parquet"))
                    shutil.move(src, os.path.join(od, run + ".psms.parquet"))
                    rec[arm] = {"n_psms": int(len(score)), "nce": info["nce"], "rt_mad_min": info["rt_mad_min"],
                                "features_s": info["timing_s"]["total"], "predict_s": info["timing_s"]["of_which_predict"],
                                "apply_s": t_apply, "total_s": time.perf_counter() - t2,
                                "cache_stats": info["cache_stats"]}
            except Exception as e:  # noqa: BLE001
                rec["error"] = repr(e)[:500]
            finally:
                for fp in (store, store + ".json"):
                    if os.path.exists(fp):
                        os.remove(fp)
            recs.append(rec)
            torch.cuda.empty_cache()
        out = {"batch": job["batch"], "wall_s": time.time() - bt0, "runs": recs, "worker_pid": os.getpid(),
               "gpu_max_mem_bytes": torch.cuda.max_memory_allocated()}
        json.dump(out, open(todo[0][:-5] + ".done", "w"), indent=1, default=float)
        with open(os.path.join(a.work, "timeline.jsonl"), "a") as fh:
            fh.write(json.dumps({"ev": "batch_scored", "t": time.time(), "batch": job["batch"], "n": len(recs),
                                 "wall_s": out["wall_s"]}) + "\n")
        idle0 = time.time()


if __name__ == "__main__":
    main()
