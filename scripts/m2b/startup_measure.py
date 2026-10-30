#!/usr/bin/env python3
"""M2b item 2: process startup cost, Lustre env vs node-local staged env (scripts/m2b/stage_env.sh).

Each rep is a fresh process that imports numpy/pyarrow/torch, inits CUDA, loads the kernels .so, loads the L1 index
and builds the GPU index (everything a one-shot search pays before its first spectrum). Run under the bench lock.

    python scripts/m2b/startup_measure.py --out results/.../startup.json [--reps 3]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CHILD = r"""
import time, os, sys, json, resource
t0 = time.perf_counter()
import numpy, pyarrow
t1 = time.perf_counter()
import torch
t2 = time.perf_counter()
torch.zeros(1, device="cuda"); torch.cuda.synchronize()
t3 = time.perf_counter()
sys.path.insert(0, os.environ["ROOT"])
from vulcan.gpu import engine as E, _lib
from vulcan.index.build import L1Index
_lib.lib()
t4 = time.perf_counter()
l1 = L1Index(os.environ["INDEX"])
t5 = time.perf_counter()
gi = E.GpuIndex(l1, E.SearchParams(), {"C": 57.021464}, {"M": 15.9949}); torch.cuda.synchronize()
t6 = time.perf_counter()
r = resource.getrusage(resource.RUSAGE_SELF)
print("JSON" + json.dumps({"numpy_pyarrow_s": t1 - t0, "torch_import_s": t2 - t1, "cuda_init_s": t3 - t2,
      "vulcan_import_kernels_s": t4 - t3, "l1_load_s": t5 - t4, "gpu_index_build_s": t6 - t5,
      "in_process_total_s": t6 - t0, "core_s": r.ru_utime + r.ru_stime}))
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--index", default="/tmp/vulcan_e2/index/human_contam_prod")
    a = ap.parse_args()
    local = os.environ.get("VULCAN_LOCAL_ENV", f"/tmp/vulcan_env_{os.getuid()}")
    envs = {"lustre_env": os.path.join(ROOT, "conda_envs", "pe", "bin", "python"),
            "local_nvme_env": os.path.join(local, "pe", "bin", "python")}
    res = {}
    for name, py in envs.items():
        reps = []
        for _ in range(a.reps):
            t = time.perf_counter()
            p = subprocess.run([py, "-c", CHILD], capture_output=True, text=True,
                               env={**os.environ, "ROOT": ROOT, "INDEX": a.index, "VULCAN_THREADS": "4"})
            wall = time.perf_counter() - t
            line = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
            if not line:
                raise RuntimeError(p.stderr[-2000:])
            r = json.loads(line[0][4:])
            r["process_wall_s"] = wall
            r["interpreter_start_s"] = wall - r["in_process_total_s"]
            reps.append(r)
            print(name, json.dumps({k: round(v, 3) for k, v in r.items()}), flush=True)
        res[name] = {"reps": reps, "median": {k: sorted(x[k] for x in reps)[len(reps) // 2] for k in reps[0]}}
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
