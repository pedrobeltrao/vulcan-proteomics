#!/usr/bin/env python3
"""Idempotent M3 preparation: L1 indexes (dev entrapment arms P x3 seeds + C; training species prod) and training
stores.  python scripts/m3/prepare.py [--indexes] [--train-stores] [--jobs 4]"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402


def _idx(args):
    fam, kind = args
    return C.ensure_index(fam, kind)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indexes", action="store_true")
    ap.add_argument("--train-stores", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    if a.indexes:
        todo = [(f, k) for f in ("human_contam", "pfu_contam") for k in ("E2002", "E2003", "C")]
        todo += [("mouse_contam", "prod"), ("ecoli_contam", "prod")]
        for f, _ in todo:
            C.fasta_path(f)
        for f in ("human_contam", "pfu_contam"):
            C.foreign_list(f)
        with ProcessPoolExecutor(a.jobs) as ex:
            for d in ex.map(_idx, todo):
                print("index", d, flush=True)
    if a.train_stores:
        rows = [r for r in C.train_rows() if os.path.exists(C.train_raw(r))]
        with ProcessPoolExecutor(a.jobs) as ex:
            for d in ex.map(C.ensure_train_store, rows):
                print("store", d, flush=True)


if __name__ == "__main__":
    main()
