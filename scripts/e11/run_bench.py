#!/usr/bin/env python3
"""E11 per-run benches, back to back (each run's repeats hold flock /tmp/vulcan_bench.lock, taken in bench_run.py):
  t8      8 threads, full chimeric v2_sep path (E7b protocol): MS1 on the CPU, then on the GPU
  t1_ms1  1 thread, MS1 stage only (the cohort decode-worker setting): CPU, then GPU
"""
import argparse
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_E11_gpu_ms1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib-root", default=os.path.join(ROOT, "data", "cache", "e7b_predlib"))
    ap.add_argument("--sets", nargs="+", default=["t8", "t1_ms1"])
    a = ap.parse_args()
    env = dict(os.environ, PYTHONPATH=ROOT)
    for s in a.sets:
        for b in ("cpu", "gpu"):
            name = f"bench_{b}_{s}"
            cmd = [sys.executable, os.path.join(ROOT, "scripts", "e11", "bench_run.py"), "--ms1-backend", b,
                   "--lib-root", a.lib_root, "--out", os.path.join(RES, name + ".json")]
            cmd += ["--threads", "8"] if s == "t8" else ["--threads", "1", "--ms1-only"]
            with open(os.path.join(RES, name + ".log"), "w") as fh:
                rc = subprocess.call(cmd, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT)
            print(name, "rc", rc, flush=True)


if __name__ == "__main__":
    main()
