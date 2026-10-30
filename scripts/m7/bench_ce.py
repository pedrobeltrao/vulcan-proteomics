#!/usr/bin/env python3
"""Runtime of the frozen CE reranking stage (production top-5 PSM tables), persistent process, model resident.
Per run: spectrum top-63 extraction + tokenisation (host) and the GPU forward (CUDA events), D0 run 2 + the 10 D8
holdout runs (prod tables); first pass = warm-up (excluded), second pass timed under /tmp/vulcan_bench.lock (fcntl).
    python scripts/m7/bench_ce.py --out results/2026-10-28_M7_encoders/bench_ce.json
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402

import e6_eval as E6  # noqa: E402


def core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default=os.path.join(E6.RES, "models", "ce_s.pt"))
    a = ap.parse_args()
    import glob
    import pyarrow.parquet as pq
    import torch
    from vulcan.encoders.score import load_model, score_table
    from vulcan.index.build import L1Index
    t0 = time.perf_counter()
    model = load_model(a.model)
    t_model = time.perf_counter() - t0
    dirs = E6.psm_dirs()
    l1 = L1Index(dirs["D0_prod"][1])
    t_l1 = time.perf_counter() - t0 - t_model
    tables = [os.path.join(dirs["D0_prod"][0], "b1937_293T_proteinID_01B_QE3_122212.psms.parquet")]
    tables += sorted(glob.glob(os.path.join(dirs["D8_prod"][0], "*.psms.parquet")))
    data = []
    for f in tables:
        run = os.path.basename(f)[:-len(".psms.parquet")]
        t = pq.read_table(f, columns=["pf_row", "charge", "spec_idx"])
        data.append((run, {k: t.column(k).to_numpy() for k in ("pf_row", "charge", "spec_idx")}, E6.find_spectra(run)))
    for run, psm, sp in data:            # warm-up pass
        score_table(model, l1, psm, sp)
    rows = {}
    with open("/tmp/vulcan_bench.lock", "w") as lk:
        tw = time.perf_counter()
        fcntl.flock(lk, fcntl.LOCK_EX)
        waited = time.perf_counter() - tw
        load = os.getloadavg()
        for run, psm, sp in data:
            c0 = core_s()
            _, info = score_table(model, l1, psm, sp)
            info["core_s"] = core_s() - c0
            rows[run] = info
            print(run, {k: round(v, 3) if isinstance(v, float) else v for k, v in info.items()}, flush=True)
        fcntl.flock(lk, fcntl.LOCK_UN)
    g = np.array([r["gpu_s"] for r in rows.values()])
    w = np.array([r["wall_s"] for r in rows.values()])
    c = np.array([r["core_s"] for r in rows.values()])
    n = np.array([r["n_rows"] for r in rows.values()])
    ns = np.array([r["n_spectra"] for r in rows.values()])
    out = {"per_run": rows, "median_gpu_s": float(np.median(g)), "median_wall_s": float(np.median(w)),
           "median_core_s": float(np.median(c)), "candidates_per_gpu_s": float(n.sum() / g.sum()),
           "spectra_per_gpu_s": float(ns.sum() / g.sum()), "gpu_s_per_1e5_candidates": float(g.sum() / n.sum() * 1e5),
           "one_time": {"model_load_s": t_model, "l1_load_s": t_l1}, "lock_wait_s": waited, "loadavg": load,
           "model_sha256": model.sha256, "gpu_max_mem_gib": torch.cuda.max_memory_allocated() / 2 ** 30}
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "per_run"}, default=float))


if __name__ == "__main__":
    main()
