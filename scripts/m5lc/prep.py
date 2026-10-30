#!/usr/bin/env python3
"""M5 LiP/CF: RAW -> spectral store (copied to NVMe first; LESSONS) -> MS1 features + hypotheses.

    python scripts/m5lc/prep.py --set d7|d10cf [--jobs 8] [--ms1-threads 4]

Stores in $W/store, MS1 in $W/ms1, per-file timing JSON in $W/logs (W = $VULCAN_M5LC_WORK, default /tmp/vulcan_m5lc).
D7 stores already decoded by M3e (/tmp/vulcan_m3/store_null, same vulcan.io.convert) are copied instead of re-decoded.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

M3NULL = "/tmp/vulcan_m3/store_null"


def conv(raw, threads):
    base = os.path.basename(raw)[:-4]
    out = C.store(base)
    if os.path.exists(out) and os.path.exists(out + ".json"):
        return base, "exists", None
    old = os.path.join(M3NULL, base + ".vspec.arrow")
    if os.path.exists(old) and os.path.exists(old + ".json"):
        shutil.copyfile(old, out + ".partial")
        shutil.copyfile(old + ".json", out + ".json")
        os.replace(out + ".partial", out)
        return base, "copied_from_m3e", None
    stage = os.path.join(C.WORK, "rawstage")
    os.makedirs(stage, exist_ok=True)
    src = os.path.join(stage, os.path.basename(raw))
    t0 = time.perf_counter()
    shutil.copyfile(raw, src)
    t_copy = time.perf_counter() - t0
    t0 = time.perf_counter()
    r = subprocess.run(["/usr/bin/time", "-v", sys.executable, "-m", "vulcan.io.convert", src, os.path.dirname(out),
                        "--threads", str(threads)], cwd=C.ROOT, env={**os.environ, "PYTHONPATH": C.ROOT},
                       capture_output=True, text=True, timeout=7200)
    os.remove(src)
    rec = {"run": base, "copy_s": t_copy, "wall_s": time.perf_counter() - t0, "rc": r.returncode,
           "time_v": C.parse_time_v(r.stderr)}
    json.dump(rec, open(os.path.join(C.WORK, "logs", f"convert_{base}.json"), "w"), indent=1)
    return base, "ok" if r.returncode == 0 else "FAILED " + r.stderr[-400:], rec


def ms1(base, threads):
    outj = os.path.join(C.WORK, "ms1", base + ".ms1.json")
    if os.path.exists(outj):
        return base, "exists"
    r = subprocess.run(["/usr/bin/time", "-v", sys.executable, "-m", "vulcan.ms1", C.store(base),
                        os.path.join(C.WORK, "ms1"), "--threads", str(threads)], cwd=C.ROOT,
                       env={**os.environ, "PYTHONPATH": C.ROOT}, capture_output=True, text=True, timeout=7200)
    json.dump({"run": base, "rc": r.returncode, "time_v": C.parse_time_v(r.stderr)},
              open(os.path.join(C.WORK, "logs", f"ms1_{base}.time.json"), "w"), indent=1)
    return base, "ok" if r.returncode == 0 else "FAILED " + r.stderr[-400:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--ms1-threads", type=int, default=4)
    a = ap.parse_args()
    for d in ("store", "ms1", "logs"):
        os.makedirs(os.path.join(C.WORK, d), exist_ok=True)
    raws = C.d7_raws() if a.set == "d7" else C.cf_raws()
    with ThreadPoolExecutor(a.jobs) as ex:
        for base, st, _ in ex.map(lambda r: conv(r, a.threads), raws):
            print("convert", base, st, flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for base, st in ex.map(lambda r: ms1(os.path.basename(r)[:-4], a.ms1_threads), raws):
            print("ms1", base, st, flush=True)


if __name__ == "__main__":
    main()
