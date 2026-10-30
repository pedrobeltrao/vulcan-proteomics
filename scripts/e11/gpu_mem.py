#!/usr/bin/env python3
"""E11 GPU memory of the MS1 backend in one torch-free process (the cohort decode-worker setting): the idle CUDA
context after session creation, the peak during a run (nvidia-smi per-process samples every 20 ms) and what is left
after the run with VULCAN_MS1_GPU_TRIM=1.  Dev stores only."""
import argparse
import json
import os
import subprocess
import sys
import threading
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def used_mib(pid):
    o = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,used_memory", "--format=csv,noheader,nounits"],
                       capture_output=True, text=True).stdout
    for line in o.strip().splitlines():
        p, m = [x.strip() for x in line.split(",")]
        if int(p) == pid:
            return int(m)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stores", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.environ["VULCAN_MS1_GPU_TRIM"] = "1"
    from vulcan import hw
    hw.pin_thread_pools(1)
    from vulcan.ms1 import gpu
    from vulcan.ms1.features import FeatureParams, find_features, load_ms1
    pid = os.getpid()
    gpu.session()
    time.sleep(1)
    rec = {"context_idle_mib": used_mib(pid), "runs": {}}
    for st in a.stores:
        ms1 = load_ms1(st)
        peak = [0]
        stop = threading.Event()

        def samp():
            while not stop.is_set():
                peak[0] = max(peak[0], used_mib(pid))
                time.sleep(0.02)
        th = threading.Thread(target=samp)
        th.start()
        fr = find_features(ms1, FeatureParams(), backend="gpu", hill_peaks=False)
        stop.set()
        th.join()
        rec["runs"][os.path.basename(st)] = {"n_peaks": int(len(ms1.mz)), "peak_mib": peak[0],
                                             "after_trim_mib": used_mib(pid), "gpu_ms": fr.timings["gpu_ms"]}
        print(rec["runs"][os.path.basename(st)], flush=True)
    json.dump(rec, open(a.out, "w"), indent=1)
    print(rec["context_idle_mib"])


if __name__ == "__main__":
    main()
