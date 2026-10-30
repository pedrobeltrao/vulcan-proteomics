"""Hardware discovery and thread-pool pinning (decision D5).

Never size a pool from ``os.cpu_count()``: on cgroup/affinity-limited nodes it reports the whole
host (64 vs 32 usable here; 128 vs 4 on the old node).  Use :func:`n_threads`.

    n_threads()            affinity-mask size; env VULCAN_THREADS overrides (use 4 for M1 continuity)
    pin_thread_pools(n)    set OMP/MKL/OPENBLAS/NUMBA/RAYON/... env (+ torch/arrow if already imported)
    fingerprint()          JSON-able dict identifying the node (recorded with every benchmark result)
    thread_env(n)          the env-var dict alone (for subprocess launching)
"""
from __future__ import annotations

import os
import platform
import re
import shutil
import socket
import subprocess
import sys

THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS",
               "RAYON_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "VULCAN_THREADS")


def affinity_cpus() -> list[int]:
    try:
        return sorted(os.sched_getaffinity(0))
    except AttributeError:  # non-Linux
        return list(range(os.cpu_count() or 1))


def n_threads() -> int:
    """Thread budget: $VULCAN_THREADS if set (>=1), else the size of the affinity mask."""
    v = os.environ.get("VULCAN_THREADS", "").strip()
    if v:
        try:
            n = int(v)
            if n >= 1:
                return n
        except ValueError:
            pass
    return max(1, len(affinity_cpus()))


def thread_env(n: int | None = None) -> dict[str, str]:
    n = int(n) if n else n_threads()
    return {k: str(n) for k in THREAD_VARS}


def pin_thread_pools(n: int | None = None) -> int:
    """Pin every thread pool we know of to n (default n_threads()).  Call before importing numpy/numba
    for env-var pools to take effect; torch / pyarrow pools are also set live if already imported."""
    n = int(n) if n else n_threads()
    os.environ.update(thread_env(n))
    if "torch" in sys.modules:
        try:
            sys.modules["torch"].set_num_threads(n)
        except Exception:
            pass
    if "pyarrow" in sys.modules:
        try:
            sys.modules["pyarrow"].set_cpu_count(n)
        except Exception:
            pass
    if "numba" in sys.modules:
        try:
            import numba
            numba.set_num_threads(min(n, numba.config.NUMBA_NUM_THREADS))
        except Exception:
            pass
    return n


def _cpuinfo() -> tuple[str | None, list[str]]:
    model, flags = None, []
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if model is None and line.startswith("model name"):
                    model = line.split(":", 1)[1].strip()
                elif line.startswith("flags"):
                    flags = line.split(":", 1)[1].split()
                    break
    except OSError:
        model = platform.processor() or None
    return model, flags


def _physical_cores(cpus: list[int]) -> dict:
    """Distinct (package, core) pairs among affinity CPUs, from sysfs topology."""
    pairs, pkgs = set(), set()
    for c in cpus:
        base = f"/sys/devices/system/cpu/cpu{c}/topology/"
        try:
            pk = int(open(base + "physical_package_id").read())
            co = int(open(base + "core_id").read())
        except (OSError, ValueError):
            pairs.add(("?", c))
            continue
        pairs.add((pk, co))
        pkgs.add(pk)
    return {"physical_cores_in_affinity": len(pairs), "sockets_in_affinity": len(pkgs) or None}


def _meminfo_gb() -> float | None:
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal"):
                    return round(int(line.split()[1]) / 1024 ** 2, 1)
    except OSError:
        pass
    return None


def _cgroup_mem_gb() -> float | None:
    for p in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            s = open(p).read().strip()
            if s.isdigit() and int(s) < 1 << 50:
                return round(int(s) / 1024 ** 3, 1)
        except OSError:
            pass
    return None


def _gpu() -> dict:
    smi = shutil.which("nvidia-smi")
    out: dict = {"gpus": [], "driver": None, "cuda_driver": None, "nvcc": None, "torch_cuda": None}
    if smi:
        try:
            q = subprocess.run([smi, "--query-gpu=name,memory.total,driver_version,compute_cap",
                                "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=20)
            for line in q.stdout.strip().splitlines():
                p = [x.strip() for x in line.split(",")]
                if len(p) >= 3:
                    out["gpus"].append({"name": p[0], "mem_mib": int(float(p[1])), "compute_cap": p[3] if len(p) > 3 else None})
                    out["driver"] = p[2]
            m = re.search(r"CUDA Version: ([\d.]+)", subprocess.run([smi], capture_output=True, text=True, timeout=20).stdout)
            out["cuda_driver"] = m.group(1) if m else None
        except Exception:
            pass
    nvcc = shutil.which("nvcc")
    if nvcc:
        try:
            m = re.search(r"release ([\d.]+)", subprocess.run([nvcc, "--version"], capture_output=True, text=True, timeout=20).stdout)
            out["nvcc"] = m.group(1) if m else None
        except Exception:
            pass
    if "torch" in sys.modules:  # do not import torch just for a fingerprint (slow on cluster FS)
        t = sys.modules["torch"]
        try:
            out["torch_cuda"] = {"torch": t.__version__, "cuda": t.version.cuda, "available": t.cuda.is_available()}
        except Exception:
            pass
    return out


def fingerprint() -> dict:
    cpus = affinity_cpus()
    model, flags = _cpuinfo()
    fl = set(flags)
    avx512 = sorted(f for f in fl if f.startswith("avx512"))
    d = {
        "hostname": socket.gethostname(),
        "cpu_model": model,
        "os_cpu_count": os.cpu_count(),
        "affinity_count": len(cpus),
        "n_threads": n_threads(),
        "VULCAN_THREADS": os.environ.get("VULCAN_THREADS"),
        "avx2": "avx2" in fl,
        "avx512": avx512,
        "ram_gb": _meminfo_gb(),
        "cgroup_mem_limit_gb": _cgroup_mem_gb(),
        "kernel": platform.release(),
        "python": platform.python_version(),
    }
    d.update(_physical_cores(cpus))
    d.update(_gpu())
    gpus = d["gpus"]
    d["gpu_summary"] = (f"{gpus[0]['name']} {gpus[0]['mem_mib'] // 1024}GB drv {d['driver']}" if gpus else None)
    return d
