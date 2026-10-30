#!/usr/bin/env python3
"""E2 speed benchmark: our K0-K2 GPU search vs Sage 0.14.7 at 32 and 4 threads, same spectra (store -> MGF for
Sage), same index (our L1 index -> peptide FASTA for Sage), same search parameters.

    python scripts/m2/e2_bench.py [--datasets D0 D1] [--kinds prod] [--reps 3]

Every timed command runs under scripts/baselines/timed.py --lock (flock /tmp/vulcan_bench.lock; GNU time -v), and
appends a record (with the vulcan.hw fingerprint) to results/2026-10-27_E2_gpu_search/timings.jsonl.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_common as C  # noqa: E402

from vulcan.gpu import engine as E  # noqa: E402
from vulcan.gpu import sage_bridge as SB  # noqa: E402
from vulcan.hw import n_threads  # noqa: E402

TIMED = [sys.executable, os.path.join(C.ROOT, "scripts", "baselines", "timed.py")]
LOG = os.path.join(C.RES, "timings.jsonl")


def timed(label, cmd, threads, stdout, meta):
    full = TIMED + ["--label", label, "--log", LOG, "--lock", "--threads", str(threads), "--stdout", stdout,
                    "--meta", json.dumps(meta), "--"] + cmd
    subprocess.run(full, check=True, cwd=C.ROOT, env={**os.environ, "PYTHONPATH": C.ROOT})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["D1", "D0"])
    ap.add_argument("--kinds", nargs="+", default=["prod"])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--repeat-inner", type=int, default=3)
    ap.add_argument("--tools", nargs="+", default=["ours", "sage32", "sage4"])
    a = ap.parse_args()
    os.makedirs(C.RES, exist_ok=True)
    os.makedirs(os.path.join(C.WORK, "bench"), exist_ok=True)
    p = E.SearchParams()
    for ds in a.datasets:
        raws = C.DATASETS[ds]["raw"]
        stores = [C.ensure_store(r) for r in raws]
        mgfs = [C.ensure_mgf(r) for r in raws]
        for kind in a.kinds:
            idx = C.ensure_index(ds, kind)
            fa = C.ensure_sage_fasta(ds, kind)
            for rep in range(a.reps):
                meta = {"dataset": ds, "index": kind, "rep": rep, "n_files": len(raws), "from": "store (ours) / MGF of the same store (Sage)",
                        "entrapment": kind != "prod"}
                if "ours" in a.tools:
                    out = os.path.join(C.WORK, "bench", f"ours_{ds}_{kind}_r{rep}")
                    timed(f"E2_ours_{ds}_{kind}", [sys.executable, "-m", "vulcan.gpu.run", "--index", idx, "--stores",
                                                    *stores, "--out", out, "--repeat", str(a.repeat_inner)],
                          n_threads(), out + ".log",
                          {**meta, "tool": "vulcan-gpu", "stats": os.path.join(out, "stats.json")})
                for t in (32, 4):
                    if f"sage{t}" not in a.tools:
                        continue
                    out = os.path.join(C.WORK, "bench", f"sage_t{t}_{ds}_{kind}_r{rep}")
                    os.makedirs(out, exist_ok=True)
                    cfg = SB.sage_config(fa, mgfs, out, p, C.STATIC, C.VAR, C.MAX_VAR)
                    SB.write_config(cfg, os.path.join(out, "config.json"))
                    timed(f"E2_sage_t{t}_{ds}_{kind}", [os.path.join(C.ROOT, "tools", "bin", "sage"),
                                                         os.path.join(out, "config.json"),
                                                         "--disable-telemetry-i-dont-want-to-improve-sage"],
                          t, os.path.join(out, "sage.log"), {**meta, "tool": "sage-0.14.7", "threads": t, "outdir": out})


if __name__ == "__main__":
    main()
