#!/usr/bin/env python3
"""IT-CID GPU closed search (K0-K2, Da fragment tolerance) of Arrow stores against one family x kind index.

    python scripts/m3c/search.py --family human_tmt --kind prod --stores a.vspec.arrow ... --out DIR [--params JSON]

Writes DIR/<run>.psms.parquet (top-5 per spectrum, label-free) and DIR/stats.json (per-run wall / core-s).
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402


def core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--kind", default="prod")
    ap.add_argument("--stores", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", default="{}")
    a = ap.parse_args()
    import torch
    import pyarrow.parquet as pq
    from vulcan.gpu import engine as E
    par = dict(C.SEARCH, **json.loads(a.params))
    p = E.SearchParams(**{k: tuple(v) if isinstance(v, list) else v for k, v in par.items()})
    os.makedirs(a.out, exist_ok=True)
    t0 = time.perf_counter()
    torch.zeros(1, device="cuda")
    C.ensure_index(a.family, a.kind)
    l1 = C.load_l1(a.family, a.kind)
    st, var = C.mods(a.family)
    gi = E.GpuIndex(l1, p, st, var)
    setup = time.perf_counter() - t0
    runs = []
    for s in a.stores:
        t1, c1 = time.perf_counter(), core_s()
        sp = E.load_store(s)
        tbl, stt, _ = E.search(sp, gi, l1, p)
        pq.write_table(tbl, os.path.join(a.out, f"{sp.run}.psms.parquet"))
        runs.append({"run": sp.run, "wall_s": time.perf_counter() - t1, "core_s": core_s() - c1,
                     "search_wall_s": stt.wall_s, "counts": stt.counts})
        print(sp.run, round(runs[-1]["wall_s"], 2), stt.counts.get("n_psms"), flush=True)
    json.dump({"family": a.family, "kind": a.kind, "index": C.index_dir(a.family, a.kind), "params": p.to_dict(),
               "setup_s": setup, "runs": runs}, open(os.path.join(a.out, "stats.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
