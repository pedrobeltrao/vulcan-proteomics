#!/usr/bin/env python3
"""Run one command with full speed accounting (bench/README.md conventions).

    python scripts/baselines/timed.py --label NAME --log timings.jsonl [--lock] [--stdout F] -- cmd args...

* --lock      take the node-wide benchmark lock (`flock /tmp/vulcan_bench.lock`) for the
              duration of the command only (lock wait time is recorded separately, not counted).
* Threads: pinned through the environment (OMP/MKL/NUMBA/RAYON/OPENBLAS/VULCAN_THREADS) to --threads N,
  default vulcan.hw.n_threads() (affinity mask, or $VULCAN_THREADS). Use --threads 4 for M1-continuity runs.
* Every record carries the vulcan.hw fingerprint (key "hw") and the thread count (key "threads_env").
* Lock file: $VULCAN_BENCH_LOCK, default /tmp/vulcan_bench.lock (the lock file is not state; it may vanish).
* Measurement: GNU /usr/bin/time -v (wall, user, sys, max RSS of the largest process in the tree).
  core-s = user + sys.
Appends one JSON line per run to --log and exits with the command's exit code.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from vulcan import hw  # noqa: E402

LOCK = os.environ.get("VULCAN_BENCH_LOCK", "/tmp/vulcan_bench.lock")


def parse_time_v(text: str) -> dict:
    def grab(pat, conv=float):
        m = re.search(pat, text)
        return conv(m.group(1)) if m else None

    wall = None
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", text)
    if m:
        parts = [float(p) for p in m.group(1).split(":")]
        wall = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
    user = grab(r"User time \(seconds\): ([\d.]+)")
    sys_ = grab(r"System time \(seconds\): ([\d.]+)")
    return {
        "wall_s_timev": wall, "user_s": user, "sys_s": sys_,
        "core_s": (user or 0) + (sys_ or 0),
        "max_rss_mb": (grab(r"Maximum resident set size \(kbytes\): (\d+)", int) or 0) / 1024.0,
        "cpu_percent": grab(r"Percent of CPU this job got: (\d+)%", int),
        "fs_inputs": grab(r"File system inputs: (\d+)", int),
        "fs_outputs": grab(r"File system outputs: (\d+)", int),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--lock", action="store_true")
    ap.add_argument("--threads", type=int, default=None,
                    help="thread budget for the child (default: vulcan.hw.n_threads(): affinity or $VULCAN_THREADS)")
    ap.add_argument("--stdout", default=None, help="file for the command's stdout+stderr")
    ap.add_argument("--meta", default="{}", help="extra JSON merged into the record")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    nthr = a.threads or hw.n_threads()
    env = dict(os.environ, **hw.thread_env(nthr))
    tv = tempfile.NamedTemporaryFile("r", suffix=".timev", delete=False)
    inner = ["/usr/bin/time", "-v", "-o", tv.name, *cmd]
    out = open(a.stdout, "w") if a.stdout else None
    lock_fd = None
    t_wait0 = time.time()
    if a.lock:
        import fcntl
        lock_fd = open(LOCK, "a")
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
    la0 = os.getloadavg()
    t0 = time.time()
    start_iso = dt.datetime.now().isoformat(timespec="seconds")
    rc = subprocess.call(inner, env=env, stdout=out or sys.stdout, stderr=subprocess.STDOUT if out else sys.stderr)
    t1 = time.time()
    if lock_fd is not None:
        lock_fd.close()
    if out:
        out.close()
    rec = {"label": a.label, "cmd": cmd, "rc": rc, "locked": a.lock, "start": start_iso,
           "lock_wait_s": round(t0 - t_wait0, 2), "wall_s": round(t1 - t0, 3), "threads_env": nthr, "loadavg_before": la0, "loadavg_after": os.getloadavg(),
           "hw": hw.fingerprint()}
    rec.update(parse_time_v(open(tv.name).read()))
    rec.update(json.loads(a.meta))
    os.unlink(tv.name)
    with open(a.log, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps({k: rec[k] for k in ("label", "rc", "wall_s", "core_s", "max_rss_mb", "lock_wait_s")}),
          file=sys.stderr)
    sys.exit(rc)


if __name__ == "__main__":
    main()
