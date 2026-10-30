"""CLI: closed GPU search of Arrow stores against a built L1 index.

    python -m vulcan.gpu.run --index DIR --stores a.vspec.arrow [...] --out OUTDIR [--repeat 2] [--params JSON]

Writes OUTDIR/<run>.psms.parquet (label-free PSM features, top `report_psms` per spectrum) and OUTDIR/stats.json
(per phase wall / GPU ms, H2D bytes, CPU core-s, counts). --repeat > 1 re-runs the searches (warm steady state);
the PSM files come from the first pass, stats keep every pass.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import time

t_start = time.perf_counter()


def _core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--stores", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--params", default="{}", help="JSON overrides of SearchParams")
    ap.add_argument("--static", default='{"C": 57.021464}')
    ap.add_argument("--var", default='{"M": 15.9949}')
    a = ap.parse_args(argv)
    import torch
    import pyarrow.parquet as pq

    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    t_import = time.perf_counter() - t_start
    c_import = _core_s()
    p = E.SearchParams(**json.loads(a.params))
    t0 = time.perf_counter()
    torch.cuda.init()
    torch.zeros(1, device="cuda")
    t_cuda = time.perf_counter() - t0
    t0 = time.perf_counter()
    l1 = L1Index(a.index)
    t_l1 = time.perf_counter() - t0
    t0 = time.perf_counter()
    gi = E.GpuIndex(l1, p, json.loads(a.static), json.loads(a.var))
    t_gi = time.perf_counter() - t0
    c_setup = _core_s()
    os.makedirs(a.out, exist_ok=True)
    passes = []
    for rep in range(a.repeat):
        runs = []
        for s in a.stores:
            t0 = time.perf_counter()
            c0 = _core_s()
            sp = E.load_store(s)
            t_load = time.perf_counter() - t0
            tbl, st, _ = E.search(sp, gi, l1, p)
            t1 = time.perf_counter()
            if rep == 0:
                pq.write_table(tbl, os.path.join(a.out, f"{sp.run}.psms.parquet"))
            t_write = time.perf_counter() - t1
            runs.append({"run": sp.run, "load_store_s": t_load, "write_s": t_write,
                         "wall_s_total": time.perf_counter() - t0, "cpu_core_s": _core_s() - c0,
                         "search_wall_s": st.wall_s, "gpu_ms": st.gpu_ms, "gpu_s": sum(st.gpu_ms.values()) / 1e3,
                         "h2d_bytes": st.h2d_bytes, "counts": st.counts,
                         "spectra_per_s_search": st.counts["n_spectra"] / st.wall_s["total"]})
        passes.append(runs)
    out = {
        "params": p.to_dict(), "index": os.path.abspath(a.index), "index_meta": l1.meta,
        "setup": {"python_import_s": t_import, "cuda_init_s": t_cuda, "l1_load_s": t_l1, "gpu_index_build_s": t_gi,
                  "gpu_index_bytes": gi.device_bytes, "gpu_index_h2d_bytes": gi.h2d, "n_postings": gi.n_postings,
                  "n_bins": gi.n_bins, "n_slices": len(gi.slices), "cpu_core_s_import": c_import,
                  "cpu_core_s_setup_total": c_setup},
        "passes": passes,
        "process_wall_s": time.perf_counter() - t_start, "process_cpu_core_s": _core_s(),
        "gpu_max_mem_bytes": torch.cuda.max_memory_allocated(),
    }
    with open(os.path.join(a.out, "stats.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({"process_wall_s": out["process_wall_s"], "setup": out["setup"],
                      "pass_walls": [[r["wall_s_total"] for r in ps] for ps in passes]}, default=float))


if __name__ == "__main__":
    main()
