#!/usr/bin/env python3
"""PRIDE download throughput (E9 stream): k parallel curl streams, HTTPS vs FTP, fixed time window, NVMe sink.

Each stream downloads a distinct E9 manifest file for --window seconds (partial transfers are fine; bytes are counted
from what landed on disk), then the files are deleted.  Reports aggregate GB/s per (protocol, streams) and the implied
runs/hour at the manifest's mean run size.  Aspera: reported unavailable if no `ascp` binary is on PATH.

    python3 scripts/e9/dl_bench.py --streams 1 4 8 16 32 --protos https ftp --window 45 --out results/.../download.json
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def run(rows, proto, k, window, sink):
    os.makedirs(sink, exist_ok=True)
    procs, paths = [], []
    t0 = time.time()
    for i, r in enumerate(rows[:k]):
        url = r["url"] if proto == "https" else r["url"].replace("https://ftp.pride.ebi.ac.uk", "ftp://ftp.pride.ebi.ac.uk")
        p = os.path.join(sink, f"{proto}_{k}_{i}.part")
        paths.append(p)
        procs.append(subprocess.Popen(["curl", "-s", "-o", p, "--max-time", str(window), url]))
    for pr in procs:
        pr.wait()
    wall = time.time() - t0
    sizes = [os.path.getsize(p) if os.path.exists(p) else 0 for p in paths]
    for p in paths:
        if os.path.exists(p):
            os.remove(p)
    return {"proto": proto, "streams": k, "wall_s": wall, "bytes": sum(sizes),
            "gb_per_s": sum(sizes) / wall / 1e9, "per_stream_mb_s": [s / wall / 1e6 for s in sizes],
            "complete_files": sum(s == int(r["expected_size"]) for s, r in zip(sizes, rows[:k]))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=os.path.join(ROOT, "data", "manifests", "e9_d8.tsv"))
    ap.add_argument("--streams", type=int, nargs="+", default=[1, 4, 8, 16, 32])
    ap.add_argument("--protos", nargs="+", default=["https", "ftp"])
    ap.add_argument("--window", type=int, default=45)
    ap.add_argument("--offset", type=int, default=500, help="manifest row offset of the files used")
    ap.add_argument("--sink", default="/tmp/vulcan_e9/dlbench")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    # largest files first so a stream never finishes inside the window
    pool = sorted(rows[a.offset:], key=lambda r: -int(r["expected_size"]))
    mean_run = sum(int(r["expected_size"]) for r in rows) / len(rows)
    res = {"aspera_available": shutil.which("ascp") is not None, "mean_run_bytes": mean_run, "window_s": a.window,
           "host": os.uname().nodename, "t_start": time.strftime("%Y-%m-%dT%H:%M:%S"), "runs": []}
    shift = 0
    for proto in a.protos:
        for k in a.streams:
            sel = pool[shift:shift + k]
            shift = (shift + k) % max(1, len(pool) - max(a.streams))
            m = run(sel, proto, k, a.window, a.sink)
            m["implied_runs_per_h"] = m["gb_per_s"] * 1e9 * 3600 / mean_run
            res["runs"].append(m)
            print(json.dumps({x: m[x] for x in ("proto", "streams", "gb_per_s", "implied_runs_per_h")}), flush=True)
            json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
