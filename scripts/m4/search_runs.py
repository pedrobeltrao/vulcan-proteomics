#!/usr/bin/env python3
"""M4 GPU searches with timing: one process, one resident index, several stores (closed / phospho / offset / open).

    python scripts/m4/search_runs.py --index DIR --var JSON --stores a.vspec.arrow ... --out OUTDIR
        [--mode closed|offset|open|offset+open] [--offsets-json FILE] [--report 1] [--repeat 2] [--params JSON]

Writes OUTDIR/<run>.psms.parquet and OUTDIR/stats.json (setup timings, per-run per-phase wall/GPU ms, core-s,
GPU memory, index size, hwinfo). Run it under `flock /tmp/vulcan_bench.lock` for timed numbers.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time

t_start = time.perf_counter()
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--var", required=True)
    ap.add_argument("--static", default='{"C": 57.021464}')
    ap.add_argument("--stores", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="closed")
    ap.add_argument("--offsets-json", default=None)
    ap.add_argument("--open-da", default="-150,500")
    ap.add_argument("--report", type=int, default=1)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--params", default="{}")
    a = ap.parse_args()
    import pyarrow.parquet as pq
    import torch

    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ptm import search as PS
    t_import = time.perf_counter() - t_start
    p = E.SearchParams(**json.loads(a.params))
    mp = PS.ModeParams()
    if "offset" in a.mode:
        mp.offsets = tuple(json.load(open(a.offsets_json)))
    if "open" in a.mode:
        lo, hi = map(float, a.open_da.split(","))
        mp.open_da = (lo, hi)
    t0 = time.perf_counter()
    torch.zeros(1, device="cuda")
    t_cuda = time.perf_counter() - t0
    t0 = time.perf_counter()
    l1 = L1Index(a.index)
    t_l1 = time.perf_counter() - t0
    t0 = time.perf_counter()
    gi = E.GpuIndex(l1, p, json.loads(a.static), json.loads(a.var))
    torch.cuda.synchronize()
    t_gi = time.perf_counter() - t0
    c_setup = core_s()
    os.makedirs(a.out, exist_ok=True)
    passes = []
    for rep in range(a.repeat):
        runs = []
        for s in a.stores:
            t0 = time.perf_counter()
            c0 = core_s()
            sp = E.load_store(s)
            t_load = time.perf_counter() - t0
            if a.mode == "closed":
                tbl, st, _ = E.search(sp, gi, l1, p, report=a.report)
            else:
                tbl, st = PS.search_modes(sp, gi, l1, p, mp, report=a.report)
            torch.cuda.synchronize()
            t1 = time.perf_counter()
            if rep == 0:
                pq.write_table(tbl, os.path.join(a.out, f"{sp.run}.psms.parquet"))
            runs.append({"run": sp.run, "load_store_s": t_load, "write_s": time.perf_counter() - t1,
                         "wall_s_total": time.perf_counter() - t0, "cpu_core_s": core_s() - c0,
                         "search_wall_s": st.wall_s, "gpu_ms": st.gpu_ms, "h2d_bytes": st.h2d_bytes,
                         "counts": st.counts, "spectra_per_s_search": st.counts["n_spectra"] / st.wall_s["total"]})
            print(json.dumps({"run": sp.run, "rep": rep, "wall": runs[-1]["wall_s_total"],
                              "search": st.wall_s["total"], "psms": tbl.num_rows}), flush=True)
        passes.append(runs)
    from vulcan.hw import fingerprint
    out = {"mode": a.mode, "params": p.to_dict(), "mode_params": {"offsets": list(mp.offsets), "open_da": mp.open_da,
           "open_tile": mp.open_tile}, "index": os.path.abspath(a.index),
           "index_meta": {k: v for k, v in l1.meta.items() if k != "fasta"}, "var": json.loads(a.var),
           "setup": {"python_import_s": t_import, "cuda_init_s": t_cuda, "l1_load_s": t_l1,
                     "gpu_index_build_s": t_gi, "gpu_index_bytes": gi.device_bytes, "n_postings": gi.n_postings,
                     "n_slices": len(gi.slices), "n_peptidoforms": len(l1), "cpu_core_s_setup_total": c_setup},
           "passes": passes, "process_wall_s": time.perf_counter() - t_start, "process_cpu_core_s": core_s(),
           "gpu_max_mem_bytes": torch.cuda.max_memory_allocated(), "hw": fingerprint()}
    with open(os.path.join(a.out, "stats.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({"process_wall_s": out["process_wall_s"], "setup": out["setup"]}, default=float))


if __name__ == "__main__":
    main()
