#!/usr/bin/env python3
"""Ingest identity: stores written by the cohort path (persistent decode worker, RAW staged to local NVMe, rebuilt
VulcanRaw helper) must equal the one-shot `python -m vulcan.io.convert` stores used by E2 (every column, every row;
only run metadata such as created_utc / source timings may differ). Also checks the PSM tables are identical.

    python scripts/m2b/ingest_identity.py --out results/.../ingest_identity.json
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
RAWS = ["data/raw/PXD001077/Velos005137.raw", "data/raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw",
        "data/raw/PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw"]
E2_STORE = "/tmp/vulcan_e2/store"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=f"/tmp/vulcan_m2b_{os.getuid()}/identity")
    a = ap.parse_args()
    from vulcan.io.store import read_table, run_meta
    out_dir, stage = os.path.join(a.work, "out"), os.path.join(a.work, "stage")
    subprocess.run([sys.executable, "-m", "vulcan.cohort", "--index", "/tmp/vulcan_e2/index/human_contam_prod",
                    "--raws", *[os.path.join(ROOT, r) for r in RAWS], "--out", out_dir, "--stage", stage,
                    "--workers", "3", "--keep-stores"], check=True, cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT})
    res = []
    for r in RAWS:
        name = os.path.basename(r)[:-4]
        ref = os.path.join(E2_STORE, name + ".vspec.arrow")
        new = os.path.join(stage, name + ".vspec.arrow")
        ta = read_table(ref).replace_schema_metadata(None)
        tb = read_table(new).replace_schema_metadata(None)
        diff = [c for c in ta.column_names if not ta.column(c).equals(tb.column(c))]
        ma, mb = run_meta(ref), run_meta(new)
        meta_diff = sorted(k for k in set(ma) | set(mb) if ma.get(k) != mb.get(k))
        res.append({"run": name, "rows": ta.num_rows, "store_identical": ta.num_rows == tb.num_rows and not diff,
                    "diff_columns": diff, "run_meta_keys_differing": meta_diff})
        print(res[-1], flush=True)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
