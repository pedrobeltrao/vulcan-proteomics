#!/usr/bin/env python3
"""E11 end-to-end chimeric cohort throughput, MS1 on the GPU vs the CPU (python -m vulcan.cohort --chimeric
[--ms1-backend gpu]).  Same design as E7b (scripts/e7b/cohort_bench.py): the 10 D8 holdout RAWs + the 2 D0 RAWs
(local NVMe copies), human_contam_prod index, the 3i predicted library (GPU-resident), cycled to --repeat-to runs,
24 decode workers, 3 post workers, one cohort process per config, each under `flock /tmp/vulcan_bench.lock`.

configs
  chim_gpu   --chimeric --ms1-backend gpu (E11)
  chim_cpu   --chimeric (E7b default CPU MS1), same session, for the before/after
RAWs are decoded in place (--no-stage-raw; they are already local NVMe copies) in both configs, unlike E7b.
Guards (LESSONS: /tmp filled to 100% in E7b): refuses to start below --min-free-gb free on /tmp, and a watcher
kills the cohort's process group if free space falls below --floor-gb (recorded as guard_killed).  Records load
average, the cohort's hw fingerprint (in the stats JSON), and GPU memory (nvidia-smi, every 2 s): node total, and per
process split into this cohort's process group vs every other process (the node is shared with other jobs).
"""
import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-28_E11_gpu_ms1")
W = f"/tmp/vulcan_e11_{os.getuid()}"
INDEX = "/tmp/vulcan_e2/index/human_contam_prod"
LIB = "/tmp/vulcan_ss/predlib"
RAWS = sorted(os.path.join("/tmp/vulcan_m3e/raw_d8", f) for f in os.listdir("/tmp/vulcan_m3e/raw_d8")
              if f.endswith(".raw")) + [os.path.join("/tmp/vulcan_ss/raw_d0", f) for f in
                                        sorted(os.listdir("/tmp/vulcan_ss/raw_d0"))]


def free_gb(path="/tmp"):
    return shutil.disk_usage(path).free / 1e9


def gpu_mem_mib():
    try:
        o = subprocess.run(["nvidia-smi", "--query-gpu=memory.used,memory.total", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=10).stdout.strip().split(",")
        return int(o[0]), int(o[1])
    except Exception:  # noqa: BLE001
        return None, None


def gpu_mem_by_group(pgid):
    """(MiB used by processes in our cohort's process group, MiB used by every other process) from nvidia-smi."""
    try:
        o = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,used_memory", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=10).stdout
    except Exception:  # noqa: BLE001
        return None, None
    ours = others = 0
    for line in o.strip().splitlines():
        try:
            pid, m = (int(x.strip()) for x in line.split(","))
        except ValueError:
            continue
        try:
            mine = os.getpgid(pid) == pgid
        except OSError:
            mine = False
        if mine:
            ours += m
        else:
            others += m
    return ours, others


def run(cfg, a):
    out = os.path.join(W, "out", cfg)
    stage = os.path.join(W, "stage", cfg)
    shutil.rmtree(out, ignore_errors=True)
    shutil.rmtree(stage, ignore_errors=True)
    os.makedirs(W, exist_ok=True)
    lst = os.path.join(W, f"rawlist_{cfg}.txt")
    open(lst, "w").write("\n".join(RAWS) + "\n")
    stats = os.path.join(RES, f"cohort_{cfg}.json")
    cmd = ["flock", a.lock, a.python, "-m", "vulcan.cohort", "--index", INDEX, "--raw-list", lst, "--out", out,
           "--stage", stage, "--workers", str(a.workers), "--stats", stats,
           "--params", json.dumps({"report_psms": 5}), "--score-lib", LIB, "--score-post-workers", str(a.post_workers),
           "--repeat-to", str(a.repeat_to), "--chimeric"]
    if cfg == "chim_gpu":
        cmd += ["--ms1-backend", "gpu", "--ms1-gpu-slots", str(a.gpu_slots)]
    if a.no_stage_raw:
        cmd += ["--no-stage-raw"]
    env = dict(os.environ, PYTHONPATH=ROOT)
    f0 = free_gb()
    if f0 < a.min_free_gb:
        rec = {"cfg": cfg, "rc": None, "refused": f"/tmp free {f0:.1f} GB < {a.min_free_gb} GB"}
        print(json.dumps(rec), flush=True)
        return rec
    t0 = time.time()
    la0 = os.getloadavg()
    samples, state = [], {"killed": None, "min_free_gb": f0}
    with open(os.path.join(W, f"cohort_{cfg}.log"), "w") as fh:
        pr = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)

        def watch():
            k = 0
            while pr.poll() is None:
                fg = free_gb()
                state["min_free_gb"] = min(state["min_free_gb"], fg)
                if fg < a.floor_gb and state["killed"] is None:
                    state["killed"] = f"/tmp free {fg:.1f} GB < floor {a.floor_gb} GB at {time.time() - t0:.0f}s"
                    os.killpg(pr.pid, signal.SIGTERM)
                    time.sleep(5)
                    if pr.poll() is None:
                        os.killpg(pr.pid, signal.SIGKILL)
                if k % 4 == 0:
                    u, tot = gpu_mem_mib()
                    ours, others = gpu_mem_by_group(pr.pid)
                    samples.append((round(time.time() - t0, 1), u, round(fg, 1), os.getloadavg()[0], ours, others))
                k += 1
                time.sleep(0.5)
        th = threading.Thread(target=watch, daemon=True)
        th.start()
        rc = pr.wait()
        th.join(5)
    _, tot = gpu_mem_mib()
    used = [s[1] for s in samples if s[1] is not None]
    rec = {"cfg": cfg, "rc": rc, "wall_s": time.time() - t0, "loadavg_start": la0, "loadavg_end": os.getloadavg(),
           "cmd": cmd, "out": out, "stats": stats, "tmp_free_gb_start": f0, "tmp_free_gb_min": state["min_free_gb"],
           "guard_floor_gb": a.floor_gb, "guard_killed": state["killed"], "gpu_mem_total_mib": tot,
           "gpu_mem_used_mib_max": max(used) if used else None,
           "gpu_mem_cohort_mib_max": max((x[4] for x in samples if x[4] is not None), default=None),
           "gpu_mem_others_mib_at_cohort_max": max(((x[4], x[5]) for x in samples if x[4] is not None),
                                                   default=(None, None))[1],
           "samples_t_gpumem_tmpfree_load1": samples}
    print(json.dumps({k: v for k, v in rec.items() if k != "samples_t_gpumem_tmpfree_load1"}), flush=True)
    shutil.rmtree(stage, ignore_errors=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", nargs="+", default=["chim_gpu", "chim_cpu"])
    ap.add_argument("--repeat-to", type=int, default=96)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--post-workers", type=int, default=3)
    ap.add_argument("--gpu-slots", type=int, default=2)
    ap.add_argument("--min-free-gb", type=float, default=50.0)
    ap.add_argument("--floor-gb", type=float, default=20.0)
    ap.add_argument("--stage-raw", dest="no_stage_raw", action="store_false",
                    help="copy each RAW to the stage dir before decoding (E7b). Default off in E11: the RAWs already sit "
                         "on local NVMe, and 24 concurrent copies (~55 GB) cannot be held above the /tmp floor")
    ap.add_argument("--lock", default=os.environ.get("VULCAN_BENCH_LOCK", "/tmp/vulcan_bench.lock"))
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()
    recs = [run(c, a) for c in a.configs]
    with open(os.path.join(RES, "cohort_runs.jsonl"), "a") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()
