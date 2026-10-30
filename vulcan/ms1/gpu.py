"""E11 GPU MS1 backend (opt-in: `find_features(..., backend="gpu")`, `vulcan.cohort --ms1-backend gpu`).

ctypes binding of tools/vulcan_ms1_kernels/libvulcan_ms1.so (built from vulcan/ms1/ms1_kernels.cu by
scripts/e11/build_ms1_kernels.sh whenever the recorded source hash differs).  No torch: the library links cudart
statically, so CPU decode workers can use it without importing torch (LESSONS rule 4).

What runs where (results/2026-10-28_E11_gpu_ms1):
  GPU  hill linking, labelling, ordering, valley split, length filter, summaries; m/z argsort, per-hill norms,
       scan-block bins, first-isotope pre-screen, isotope chains + averagine prefix, nn >= 2 filter, charge ban,
       priority order; the per-peak feature map of the feature table.
  CPU  the sequential greedy acceptance (vulcan.ms1.isotopes._greedy), the per-feature table, hypotheses.
Every array equals the CPU path bitwise (tests/ms1/test_gpu_exact.py, scripts/e11/equal_gpu.py).  The CPU path stays
the reference and the default.
"""
from __future__ import annotations

import ctypes as C
import hashlib
import os
import subprocess
import threading
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "vulcan", "ms1", "ms1_kernels.cu")
OUT = os.path.join(ROOT, "tools", "vulcan_ms1_kernels")
SO = os.path.join(OUT, "libvulcan_ms1.so")
BUILD = os.path.join(ROOT, "scripts", "e11", "build_ms1_kernels.sh")
MAX_ISO_GPU = 16

P, I32, I64, F64, U64 = C.c_void_p, C.c_int, C.c_int64, C.c_double, C.c_uint64
_SIGS = {
    "ms1g_create": [P],
    "ms1g_release": [P],
    "ms1g_destroy": [P],
    "ms1g_set_pool_threshold": [U64],
    "ms1g_trim": [],
    "ms1g_mem_info": [P, P],
    "ms1g_hills": [P, I64, I64, P, P, I32, P, I32, P, P, I64, F64, I32, F64, P, P],
    "ms1g_hills_get": [P] + [P] * 12,
    "ms1g_isotopes": [P, P, I32, F64, I32, F64, I32, P, I32, I32, F64, F64, I32, I32, P, P],
    "ms1g_isotopes_get": [P] + [P] * 7,
    "ms1g_peak_map": [P, P, P, I64, P, P],
}

_lib = None
_lock = threading.Lock()
_session = None


def _src_hash() -> str:
    with open(SRC, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def lib():
    global _lib
    if _lib is not None:
        return _lib
    with _lock:
        if _lib is not None:
            return _lib
        stamp = os.path.join(OUT, "source.sha256")
        fresh = os.path.exists(SO) and os.path.exists(stamp) and open(stamp).read().strip() == _src_hash()
        if not fresh:
            subprocess.run(["bash", BUILD], check=True)
        L = C.CDLL(SO)
        for name, args in _SIGS.items():
            f = getattr(L, name)
            f.argtypes = args
            f.restype = C.c_int
        L.ms1g_gpu_ms.argtypes = [P, P]
        L.ms1g_gpu_ms.restype = None
        _lib = L
    return _lib


def available() -> bool:
    """True if the library loads and a CUDA device is usable (no torch import)."""
    try:
        return session() is not None
    except Exception:  # noqa: BLE001
        return False


def _chk(code: int, what: str):
    if code != 0:
        raise RuntimeError(f"vulcan.ms1.gpu {what}: error {code}")


def _p(a: np.ndarray | None):
    return None if a is None else a.ctypes.data


class Session:
    """One CUDA stream + the device buffers of the current run.  Not thread-safe: use one per thread (or the
    process-wide `session()` under `session_lock`)."""

    def __init__(self, pool_bytes: int | None = None):
        L = lib()
        h = C.c_void_p()
        _chk(L.ms1g_create(C.byref(h)), "create")
        self.h = h
        self.L = L
        if pool_bytes is None:
            pool_bytes = int(float(os.environ.get("VULCAN_MS1_GPU_POOL_MB", "4096")) * 2**20)
        _chk(L.ms1g_set_pool_threshold(int(pool_bytes)), "pool")
        self.trim_after_run = os.environ.get("VULCAN_MS1_GPU_TRIM", "0") == "1"

    def mem_info(self) -> tuple[int, int]:
        """(free, total) device bytes (cudaMemGetInfo)."""
        f, t = C.c_uint64(), C.c_uint64()
        _chk(self.L.ms1g_mem_info(C.byref(f), C.byref(t)), "mem_info")
        return f.value, t.value

    @staticmethod
    def bytes_needed(n_peaks: int) -> int:
        """Upper estimate of one run's device memory: measured peak 4.5 GiB at 79.9M peaks (E9_0000, 59 B/peak) and
        3.4 GiB at 50.3M (D0, 73 B/peak); 96 B/peak + 256 MiB covers both with margin (gpu_mem_worker.json)."""
        return 96 * int(n_peaks) + (256 << 20)

    def gpu_ms(self) -> dict:
        a = np.zeros(4, np.float32)
        self.L.ms1g_gpu_ms(self.h, a.ctypes.data)
        return {"hills": float(a[0]), "isotopes": float(a[1]), "peak_map": float(a[2])}

    def release(self):
        _chk(self.L.ms1g_release(self.h), "release")
        if self.trim_after_run:
            _chk(self.L.ms1g_trim(), "trim")

    def close(self):
        if self.h is not None:
            self.L.ms1g_destroy(self.h)
            self.h = None

    # ------------------------------------------------------------------ hills
    def hills(self, scan: np.ndarray, mz: np.ndarray, inten: np.ndarray, ms1_rt: np.ndarray, src: np.ndarray | None,
              *, htol_ppm: float, minlh: int, hvf: float, timings: dict | None = None, peaks: bool = True):
        """detect_hills(gap=0) on the GPU; same Hills arrays (bitwise).  peaks=False: the per-peak arrays (p_scan,
        p_mz, p_int, p_src) stay on the device (None on the host; fetch_peaks() before release) -- the feature table and
        hypotheses do not read them (production cohort)."""
        from .hills import Hills
        t = timings if timings is not None else {}
        t0 = time.perf_counter()
        n_rt = int(ms1_rt.shape[0])
        soff = np.searchsorted(scan, np.arange(n_rt + 1)).astype(np.int64)
        mz = np.ascontiguousarray(mz if mz.dtype in (np.float32, np.float64) else mz.astype(np.float64))
        inten = np.ascontiguousarray(inten if inten.dtype in (np.float32, np.float64) else inten.astype(np.float32))
        rt = np.ascontiguousarray(ms1_rt, dtype=np.float64)
        src = None if src is None else np.ascontiguousarray(src, dtype=np.int64)
        n = len(mz)
        nh, npk = C.c_int64(), C.c_int64()
        _chk(self.L.ms1g_hills(self.h, n, n_rt, _p(soff), _p(mz), int(mz.dtype == np.float64), _p(inten),
                               int(inten.dtype == np.float64), _p(src), _p(rt), n_rt, float(htol_ppm), int(minlh),
                               float(hvf), C.byref(nh), C.byref(npk)), "hills")
        nh, npk = nh.value, npk.value
        t1 = time.perf_counter()
        t["gpu_hills"] = t.get("gpu_hills", 0.0) + t1 - t0
        a = dict(off=np.empty(nh + 1, np.int64), p_scan=None, p_mz=None, p_int=None, p_src=None,
                 h_mz=np.empty(nh, np.float64),
                 h_start=np.empty(nh, np.int32), h_end=np.empty(nh, np.int32), h_apex_scan=np.empty(nh, np.int32),
                 h_apex_int=np.empty(nh, np.float32), h_area=np.empty(nh, np.float32), h_len=np.empty(nh, np.int32))
        _chk(self.L.ms1g_hills_get(self.h, *(_p(a[k]) for k in (
            "off", "p_scan", "p_mz", "p_int", "p_src", "h_mz", "h_start", "h_end", "h_apex_scan", "h_apex_int",
            "h_area", "h_len"))), "hills_get")
        if nh == 0:
            a["off"][0] = 0
        h = Hills(**a)
        if peaks:
            self.fetch_peaks(h)
        t["gpu_hills_d2h"] = t.get("gpu_hills_d2h", 0.0) + time.perf_counter() - t1
        return h

    def fetch_peaks(self, h) -> None:
        """Copy the device-resident per-peak hill arrays into h (p_scan, p_mz, p_int, p_src)."""
        npk = int(h.off[-1])
        h.p_scan, h.p_mz = np.empty(npk, np.int32), np.empty(npk, np.float64)
        h.p_int, h.p_src = np.empty(npk, np.float32), np.empty(npk, np.int64)
        _chk(self.L.ms1g_hills_get(self.h, None, _p(h.p_scan), _p(h.p_mz), _p(h.p_int), _p(h.p_src), *([None] * 7)),
             "hills_get")

    # ------------------------------------------------------------------ isotopes
    def isotopes(self, h, *, charges, itol_ppm, max_iso, min_trace_cos, min_overlap_for_cos, min_avg_cos,
                 min_explained, use_charge_ban, scan_bin, timings: dict | None = None) -> dict:
        """group_isotopes fast path on the GPU (hills of the last `hills` call) + the CPU greedy pass."""
        raw = self.isotopes_raw(charges=charges, itol_ppm=itol_ppm, max_iso=max_iso, min_trace_cos=min_trace_cos,
                                min_overlap_for_cos=min_overlap_for_cos, min_avg_cos=min_avg_cos,
                                min_explained=min_explained, use_charge_ban=use_charge_ban, scan_bin=scan_bin,
                                timings=timings)
        return finish_isotopes(raw, h, min_avg_cos=min_avg_cos, min_explained=min_explained, timings=timings)

    def isotopes_raw(self, *, charges, itol_ppm, max_iso, min_trace_cos, min_overlap_for_cos, min_avg_cos,
                     min_explained, use_charge_ban, scan_bin, timings: dict | None = None) -> dict:
        """Device part of group_isotopes: candidates (nn >= 2, after the charge ban) and their priority order, on
        the host.  The device may be released afterwards; finish_isotopes() runs the greedy pass."""
        from .averagine import get_table
        t = timings if timings is not None else {}
        t0 = time.perf_counter()
        table = np.ascontiguousarray(get_table(), dtype=np.float64)
        charges = np.ascontiguousarray(charges, dtype=np.int64)
        ncand, m = C.c_int64(), C.c_int64()
        _chk(self.L.ms1g_isotopes(self.h, _p(charges), len(charges), float(itol_ppm), int(max_iso),
                                  float(min_trace_cos), int(min_overlap_for_cos), _p(table), table.shape[0],
                                  table.shape[1], float(min_avg_cos), float(min_explained), int(bool(use_charge_ban)),
                                  int(scan_bin), C.byref(ncand), C.byref(m)), "isotopes")
        m = m.value
        r = {"ci": np.empty(m, np.int64), "cz": np.empty(m, np.int64), "hills": np.empty((m, max_iso), np.int64),
             "tcos": np.empty((m, max_iso), np.float32), "acos": np.empty(m, np.float32),
             "nn": np.empty(m, np.int32), "order": np.empty(m, np.int64), "n_candidates": int(ncand.value)}
        _chk(self.L.ms1g_isotopes_get(self.h, *(_p(r[k]) for k in ("ci", "cz", "hills", "tcos", "acos", "nn",
                                                                    "order"))), "isotopes_get")
        t["gpu_isotopes"] = t.get("gpu_isotopes", 0.0) + time.perf_counter() - t0
        return r

    def fetch_src(self, h) -> None:
        """Copy only p_src (the per-peak source index) into h: enough for the CPU peak map of the feature table."""
        h.p_src = np.empty(int(h.off[-1]), np.int64)
        _chk(self.L.ms1g_hills_get(self.h, None, None, None, None, _p(h.p_src), *([None] * 7)), "hills_get")

    def peak_map(self, hill_feat: np.ndarray, hill_iso: np.ndarray, n_all: int):
        """peak_feature / peak_iso (features._feature_table) from the device-resident hill peaks."""
        pf = np.empty(n_all, np.int32)
        pi = np.empty(n_all, np.int8)
        hf = np.ascontiguousarray(hill_feat, np.int32)
        hi = np.ascontiguousarray(hill_iso, np.int8)
        _chk(self.L.ms1g_peak_map(self.h, _p(hf), _p(hi), int(n_all), _p(pf), _p(pi)), "peak_map")
        return pf, pi


def finish_isotopes(raw: dict, h, *, min_avg_cos, min_explained, timings: dict | None = None) -> dict:
    """Host part of group_isotopes (no device needed): the sequential greedy acceptance, same numba kernel."""
    from .averagine import get_table
    from .isotopes import _greedy
    t = timings if timings is not None else {}
    t1 = time.perf_counter()
    table = np.ascontiguousarray(get_table(), dtype=np.float64)
    ci, cz, hills, tcos = raw["ci"], raw["cz"], raw["hills"], raw["tcos"]
    acc, acc_n, acc_cos = _greedy(raw["order"], hills, raw["nn"], cz, h.h_mz, h.h_apex_int, table, min_avg_cos,
                                  min_explained, h.n)
    sel = np.flatnonzero(acc)
    t["greedy"] = t.get("greedy", 0.0) + time.perf_counter() - t1
    return {
        "mono_hill": ci[sel].astype(np.int64),
        "charge": cz[sel].astype(np.int8),
        "n_isotopes": acc_n[sel],
        "iso_cos": acc_cos[sel],
        "hills": hills[sel],
        "trace_cos": tcos[sel],
        "n_candidates": raw["n_candidates"],
        "fast_path": True,
    }


session_lock = threading.Lock()
_gate = None


def set_gate(g) -> None:
    """Cross-process limit on concurrent GPU MS1 runs (vulcan.cohort decode workers: a multiprocessing Semaphore),
    so N workers never hold N runs' device buffers at once."""
    global _gate
    _gate = g


def gate():
    import contextlib
    return contextlib.nullcontext() if _gate is None else _gate


def session() -> Session:
    """Process-wide session (lazily created).  Callers hold `session_lock` for a whole run."""
    global _session
    if _session is None:
        _session = Session()
    return _session
