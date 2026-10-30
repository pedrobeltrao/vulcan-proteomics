#!/usr/bin/env python3
"""Re-verify the direct Thermo reader on the current node and compare with E1 (17.8 core-s / 1e5 MS2, old node).

    python scripts/ingest/hw_ingest_check.py --raw data/raw/PXD001077/Velos005137.raw --out results/DIR/ingest.json \
        [--conc 1,8] [--stage /tmp/vulcan_hwchk]

Each case launches N concurrent `python -m vulcan.io.convert` processes (1 thread each, same input, distinct
outputs) under GNU time -v inside the benchmark lock; reports core-s per 1e5 MS2 per case, wall, and the hwinfo
fingerprint.  Staging dir is scratch (recreated).
"""
import argparse
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan import hw  # noqa: E402


def tv(path):
    t = open(path).read()
    g = lambda p: float(re.search(p, t).group(1))
    w = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", t).group(1)
    wall = sum(float(x) * 60 ** i for i, x in enumerate(reversed(w.split(":"))))
    return dict(user_s=g(r"User time \(seconds\): ([\d.]+)"), sys_s=g(r"System time \(seconds\): ([\d.]+)"), wall_s=wall,
                max_rss_mb=g(r"Maximum resident set size \(kbytes\): (\d+)") / 1024)


def run_case(raw, stage, n):
    env = dict(os.environ, **hw.thread_env(1), PYTHONPATH=ROOT, DOTNET_gcServer="0")
    procs = []
    t0 = time.perf_counter()
    for i in range(n):
        out = os.path.join(stage, f"o{i}.vspec.arrow")
        tf = os.path.join(stage, f"t{i}.txt")
        p = subprocess.Popen(["/usr/bin/time", "-v", "-o", tf, sys.executable, "-m", "vulcan.io.convert", raw, out],
                             env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        procs.append((p, tf, out))
    members = []
    for p, tf, out in procs:
        so, se = p.communicate()
        m = tv(tf)
        meta = json.loads(so.strip().splitlines()[-1])
        m.update(rc=p.returncode, n_ms2=meta["n_ms2"], n_spectra=meta["n_spectra"], timing_s=meta["timing_s"])
        m["core_s"] = m["user_s"] + m["sys_s"]
        m["core_s_per_1e5_ms2"] = m["core_s"] / m["n_ms2"] * 1e5
        members.append(m)
        os.remove(out)
    wall = time.perf_counter() - t0
    tot = sum(m["core_s"] for m in members)
    return dict(n_concurrent=n, batch_wall_s=wall, total_core_s=tot,
                core_s_per_1e5_ms2_mean=sum(m["core_s_per_1e5_ms2"] for m in members) / n,
                files_per_hour_per_node=n / wall * 3600, members=members)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--conc", default="1,8")
    ap.add_argument("--stage", default="/tmp/vulcan_hwchk")
    ap.add_argument("--reps", type=int, default=2, help="repeats of the 1-file case (first is page-cache-cold-ish)")
    a = ap.parse_args()
    shutil.rmtree(a.stage, ignore_errors=True)
    os.makedirs(a.stage)
    res = dict(raw=a.raw, e1_reference_core_s_per_1e5_ms2=17.8, e1_node="4-core EPYC 7742 (old node)", cases=[])
    lock = open(os.environ.get("VULCAN_BENCH_LOCK", "/tmp/vulcan_bench.lock"), "a")
    fcntl.flock(lock, fcntl.LOCK_EX)
    try:
        run_case(a.raw, a.stage, 1)  # warm-up (page cache, .NET/python import path)
        for n in [int(x) for x in a.conc.split(",")]:
            for r in range(a.reps if n == 1 else 1):
                c = run_case(a.raw, a.stage, n)
                c["rep"] = r
                res["cases"].append(c)
                print(f"n={n} rep={r}: core-s/1e5 MS2 {c['core_s_per_1e5_ms2_mean']:.2f}  wall {c['batch_wall_s']:.1f}s", flush=True)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN)
    res["hw"] = hw.fingerprint()
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
