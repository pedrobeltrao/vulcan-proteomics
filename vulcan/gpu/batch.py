"""Pipelined multi-run search: the CPU reads/flattens store i+1 and writes the PSMs of run i-1 while the GPU
searches run i (DESIGN §2: work goes to whichever device finishes the pipeline soonest).

    python -m vulcan.gpu.batch --index DIR --stores a.vspec.arrow ... --out OUTDIR [--threads-io 2]
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import time
from concurrent.futures import ThreadPoolExecutor


def _core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def search_many(stores, gi, l1, p, out_dir, io_threads=2):
    import pyarrow.parquet as pq

    from . import engine as E
    os.makedirs(out_dir, exist_ok=True)
    t0, c0 = time.perf_counter(), _core_s()
    stats = []
    with ThreadPoolExecutor(io_threads) as ex:
        nxt = ex.submit(E.load_store, stores[0]) if stores else None
        writes = []
        for i in range(len(stores)):
            sp = nxt.result()
            nxt = ex.submit(E.load_store, stores[i + 1]) if i + 1 < len(stores) else None
            tbl, st, _ = E.search(sp, gi, l1, p)
            writes.append(ex.submit(pq.write_table, tbl, os.path.join(out_dir, f"{sp.run}.psms.parquet")))
            stats.append({"run": sp.run, "n_spectra": len(sp), "search_wall_s": st.wall_s["total"],
                          "gpu_ms": st.gpu_ms, "h2d_bytes": st.h2d_bytes})
        for w in writes:
            w.result()
    return {"wall_s": time.perf_counter() - t0, "cpu_core_s": _core_s() - c0, "runs": stats,
            "n_spectra": sum(s["n_spectra"] for s in stats)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--stores", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--io-threads", type=int, default=2)
    ap.add_argument("--params", default="{}")
    a = ap.parse_args(argv)
    t_start = time.perf_counter()
    import torch

    from . import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(**json.loads(a.params))
    torch.zeros(1, device="cuda")
    l1 = L1Index(a.index)
    gi = E.GpuIndex(l1, p, {"C": 57.021464}, {"M": 15.9949})
    t_setup = time.perf_counter() - t_start
    r = search_many(a.stores, gi, l1, p, a.out, a.io_threads)
    r.update(setup_s=t_setup, process_wall_s=time.perf_counter() - t_start, process_core_s=_core_s(),
             spectra_per_s_steady=r["n_spectra"] / r["wall_s"])
    with open(os.path.join(a.out, "batch_stats.json"), "w") as fh:
        json.dump(r, fh, indent=1)
    print(json.dumps({k: v for k, v in r.items() if k != "runs"}))


if __name__ == "__main__":
    main()
