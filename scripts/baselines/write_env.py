#!/usr/bin/env python3
"""Write results/2026-10-27_M1e_baselines/{env.txt, git_sha}."""
import json
import os
import platform
import sys
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-27_M1e_baselines")
sys.path.insert(0, ROOT)
from vulcan import hw as HW  # noqa: E402

WORK = os.environ.get("VULCAN_M1E_WORK", "/tmp/vulcan_m1e")


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120).stdout.strip()
    except Exception as e:  # noqa: BLE001
        return f"ERR {e}"


lines = [
    f"host: {platform.node()}  kernel: {platform.release()}",
    f"cpus_allowed: {sh('grep Cpus_allowed_list /proc/self/status')}",
    f"cpu: {sh('grep -m1 \"model name\" /proc/cpuinfo')}",
    f"mem: {sh('free -g | sed -n 2p')}",
    f"threads: {HW.n_threads()} (vulcan.hw.n_threads(): affinity mask or $VULCAN_THREADS; OMP/MKL/NUMBA/RAYON/OPENBLAS, Sage RAYON_NUM_THREADS, Comet num_threads; mokapot --max_workers 1)",
    f"hwinfo: {json.dumps(HW.fingerprint())}",
    "timed runs: flock /tmp/vulcan_bench.lock + GNU time -v (scripts/baselines/timed.py)",
    f"sage: {sh(ROOT + '/tools/bin/sage --version')} (release tag v0.14.7, commit 99407db6; crate string not bumped)",
    f"comet: {sh(ROOT + '/conda_envs/comet/bin/comet 2>&1 | grep -m1 version')}",
    f"ThermoRawFileParser: {sh(ROOT + '/conda_envs/trfp/bin/ThermoRawFileParser --version')}",
    f"mokapot env ({WORK}/envs/mokapot): {sh(WORK + '/envs/mokapot/bin/python -c \"import mokapot,numpy,pandas,sklearn;print(mokapot.__version__,numpy.__version__,pandas.__version__,sklearn.__version__)\"')} (mokapot numpy pandas sklearn)",
    f"FDRBench: 1.1.1 release zip, java: {sh(ROOT + '/conda_envs/java/bin/java -version 2>&1 | head -1')}",
    f"python (pe env): {sh(ROOT + '/conda_envs/pe/bin/python --version')}; numpy {sh(ROOT + '/conda_envs/pe/bin/python -c \"import numpy;print(numpy.__version__)\"')}",
]
with open(os.path.join(RES, "env.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
sha = sh(f"cd {ROOT} && git rev-parse HEAD") or "no-commit"
dirty = sh(f"cd {ROOT} && git status --porcelain | head -1")
with open(os.path.join(RES, "git_sha"), "w") as fh:
    fh.write(sha + (" (dirty working tree)" if dirty else "") + "\n")
print(open(os.path.join(RES, "env.txt")).read())
