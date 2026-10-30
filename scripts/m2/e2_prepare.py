#!/usr/bin/env python3
"""Idempotently build everything E2 needs on local disk: stores, L1 indexes (prod, E2001), MGFs, Sage FASTAs."""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_common as C  # noqa: E402

for ds, d in C.DATASETS.items():
    for raw in d["raw"]:
        t = time.time(); C.ensure_store(raw); C.ensure_mgf(raw); print(ds, raw, "store+mgf", round(time.time() - t, 1), flush=True)
    for kind in ("prod", f"E{C.ENTRAP_SEED}"):
        t = time.time(); C.ensure_index(ds, kind); C.ensure_sage_fasta(ds, kind); print(ds, kind, "index+fasta", round(time.time() - t, 1), flush=True)
