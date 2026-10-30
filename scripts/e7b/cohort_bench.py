#!/usr/bin/env python3
"""E7b end-to-end cohort throughput in chimeric mode, RAW -> MS1 features/hypotheses -> chimeric search -> v2 scores
-> per-stratum TDC peptides -> picked protein groups (python -m vulcan.cohort --chimeric).  Design as 3i
(scripts/scoring_speed/cohort_bench.py): the 10 D8 holdout RAWs + the 2 D0 RAWs (local NVMe copies), human_contam_prod
index, the 3i predicted library (/tmp/vulcan_ss/predlib, GPU-resident), cycled to --repeat-to runs, 24 decode
workers, 3 post workers, one cohort process per config, each under `flock /tmp/vulcan_bench.lock`.

configs
  base12        the 12 base RAWs once, --chimeric --keep-stores, NOT timed: reference outputs + stores for
                scripts/e7b/equal_cohort.py (reference pipeline per run; every large-cohort run vs its base run)
  chim          --chimeric (fast MS1 + fast scoring/competition)            [after]
  chim_ms1ref   --chimeric with VULCAN_MS1_REFERENCE=1 (original MS1 kernels; scoring already fast)
  closed        3i production config (closed search + frozen v1 fast scoring), same session, for reference
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_E7b_fair_speed")
W = f"/tmp/vulcan_e7b_{os.getuid()}"
INDEX = "/tmp/vulcan_e2/index/human_contam_prod"
LIB = "/tmp/vulcan_ss/predlib"
RAWS = sorted(os.path.join("/tmp/vulcan_m3e/raw_d8", f) for f in os.listdir("/tmp/vulcan_m3e/raw_d8")
              if f.endswith(".raw")) + [os.path.join("/tmp/vulcan_ss/raw_d0", f) for f in
                                        sorted(os.listdir("/tmp/vulcan_ss/raw_d0"))]


def run(cfg, repeat_to, workers, post_workers, py):
    out = os.path.join(W, "out", cfg)
    stage = os.path.join(W, "stage", cfg)
    shutil.rmtree(out, ignore_errors=True)
    shutil.rmtree(stage, ignore_errors=True)
    os.makedirs(W, exist_ok=True)
    lst = os.path.join(W, f"rawlist_{cfg}.txt")
    open(lst, "w").write("\n".join(RAWS) + "\n")
    stats = os.path.join(RES, f"cohort_{cfg}.json")
    cmd = [py, "-m", "vulcan.cohort", "--index", INDEX, "--raw-list", lst, "--out", out, "--stage", stage,
           "--workers", str(workers), "--stats", stats, "--params", json.dumps({"report_psms": 5}),
           "--score-lib", LIB, "--score-post-workers", str(post_workers)]
    env = dict(os.environ, PYTHONPATH=ROOT)
    if cfg != "base12":
        cmd = ["flock", "/tmp/vulcan_bench.lock"] + cmd + ["--repeat-to", str(repeat_to)]
    if cfg in ("base12", "chim", "chim_ms1ref"):
        cmd += ["--chimeric"]
    if cfg == "base12":
        cmd += ["--keep-stores"]
    if cfg == "chim_ms1ref":
        env["VULCAN_MS1_REFERENCE"] = "1"
    t0 = time.time()
    la0 = os.getloadavg()
    with open(os.path.join(W, f"cohort_{cfg}.log"), "w") as fh:
        rc = subprocess.call(cmd, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT)
    rec = {"cfg": cfg, "rc": rc, "wall_s": time.time() - t0, "loadavg_start": la0, "loadavg_end": os.getloadavg(),
           "cmd": cmd, "out": out, "stats": stats}
    print(json.dumps(rec), flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", nargs="+", default=["base12", "chim", "chim_ms1ref", "closed"])
    ap.add_argument("--repeat-to", type=int, default=192)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--post-workers", type=int, default=3)
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()
    recs = [run(c, a.repeat_to, a.workers, a.post_workers, a.python) for c in a.configs]
    p = os.path.join(RES, "cohort_runs.jsonl")
    with open(p, "a") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()
