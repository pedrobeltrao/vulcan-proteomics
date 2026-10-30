#!/usr/bin/env python3
"""M2b item 4: does H2D matter at cohort scale? Times the per-run host->device peak transfer of K0 as it is today
(pageable f32 m/z + f32 intensity + i64 lengths) vs pinned f32, and vs the DESIGN L0 compact wire format
(u32 log-ppm bin + f16 intensity, pinned; size only plus a pinned copy), next to the run's full GPU search time.

    python scripts/m2b/h2d_measure.py --stores a.vspec.arrow ... --out h2d.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stores", nargs="+", required=True)
    ap.add_argument("--index", default="/tmp/vulcan_e2/index/human_contam_prod")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import numpy as np
    import torch

    from vulcan.cohort import load_ms2
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    dev = torch.device("cuda")
    l1 = L1Index(a.index)
    p = E.SearchParams()
    gi = E.GpuIndex(l1, p, {"C": 57.021464}, {"M": 15.9949})

    def timed(fn):
        ts = []
        for _ in range(a.reps):
            torch.cuda.synchronize()
            t = time.perf_counter()
            fn()
            torch.cuda.synchronize()
            ts.append(time.perf_counter() - t)
        return float(np.median(ts))

    res = []
    for s in a.stores:
        sp = load_ms2(s)
        E.search(sp, gi, l1, p)          # warm
        t_search = timed(lambda: E.search(sp, gi, l1, p))
        lens = np.diff(sp.offsets)
        nbytes_f32 = sp.mz.nbytes + sp.intensity.nbytes + lens.nbytes
        t_page = timed(lambda: (torch.from_numpy(sp.mz).to(dev), torch.from_numpy(sp.intensity).to(dev),
                                torch.from_numpy(lens).to(dev)))
        pm, pi, pl = (torch.from_numpy(x).pin_memory() for x in (sp.mz, sp.intensity, lens))
        t_pin = timed(lambda: (pm.to(dev, non_blocking=True), pi.to(dev, non_blocking=True),
                               pl.to(dev, non_blocking=True)))
        # compact wire: u32 bin (log-ppm, 2^-16 sub-bin as in K1 keys) + f16 intensity (scaled per spectrum) + u16 len
        cb = torch.empty(len(sp.mz), dtype=torch.int32).pin_memory()
        ch = torch.empty(len(sp.mz), dtype=torch.float16).pin_memory()
        cl = torch.from_numpy(lens.astype(np.uint16).view(np.int16)).pin_memory()
        nbytes_compact = cb.numel() * 4 + ch.numel() * 2 + cl.numel() * 2
        t_cmp = timed(lambda: (cb.to(dev, non_blocking=True), ch.to(dev, non_blocking=True),
                               cl.to(dev, non_blocking=True)))
        t_pack = time.perf_counter()
        q = np.floor(np.log(np.maximum(sp.mz, 1e-3) / 50.0) * (1.0 / np.log1p(20e-6)) * 65536.0).astype(np.int64)
        _ = q.astype(np.uint32), sp.intensity.astype(np.float16)
        t_pack = time.perf_counter() - t_pack
        r = {"store": os.path.basename(s), "n_ms2": len(sp), "n_peaks_ms2": int(len(sp.mz)),
             "search_s": t_search, "h2d_bytes_f32": int(nbytes_f32), "h2d_pageable_s": t_page,
             "h2d_pinned_s": t_pin, "h2d_bytes_compact": int(nbytes_compact), "h2d_compact_pinned_s": t_cmp,
             "compact_pack_cpu_s": t_pack,
             "pageable_share_of_search": t_page / t_search, "pageable_GBps": nbytes_f32 / t_page / 1e9,
             "pinned_GBps": nbytes_f32 / t_pin / 1e9}
        print(json.dumps(r), flush=True)
        res.append(r)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
