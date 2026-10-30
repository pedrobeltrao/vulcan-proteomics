"""Direct Thermo RAW -> store, in-process via pythonnet + the VulcanRaw C# helper.

The helper (vulcan/io/dotnet/VulcanRaw, built by scripts/ingest/build_vulcanraw.sh into
tools/vulcanraw/) wraps Thermo's RawFileReader (the DLLs shipped with ThermoRawFileParser in
conda_envs/trfp) and returns one columnar Batch per block of scans, so Python never loops
per spectrum. .NET runtime: the SDK's shared runtime in tools/dotnet-sdk if present, else a
synthetic DOTNET_ROOT built from TRFP's self-contained runtime (see STORE_SPEC.md).
"""
from __future__ import annotations

import ctypes
import json
import os
import time

import numpy as np

from . import schema as S
from .store import StoreWriter

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TRFP_BIN = os.environ.get("VULCAN_TRFP_BIN", os.path.join(ROOT, "conda_envs", "trfp", "bin"))
HELPER_DLL = os.environ.get("VULCAN_RAW_DLL", os.path.join(ROOT, "tools", "vulcanraw", "VulcanRaw.dll"))
SDK_ROOT = os.environ.get("VULCAN_DOTNET_ROOT", os.path.join(ROOT, "tools", "dotnet-sdk"))

_loaded = None


def _synthetic_dotnet_root() -> tuple[str, str]:
    """Make a DOTNET_ROOT layout from TRFP's self-contained runtime (no SDK needed)."""
    base = os.path.join(os.environ.get("TMPDIR", "/tmp"), f"vulcan_dotnet_{os.getuid()}")
    ver = "8.0.22"
    fx = os.path.join(base, "shared", "Microsoft.NETCore.App", ver)
    host = os.path.join(base, "host", "fxr", ver)
    if not os.path.exists(os.path.join(fx, "Microsoft.NETCore.App.deps.json")):
        os.makedirs(fx, exist_ok=True)
        os.makedirs(host, exist_ok=True)
        import shutil
        shutil.copy(os.path.join(TRFP_BIN, "libhostfxr.so"), host)
        for f in os.listdir(TRFP_BIN):
            if f.endswith((".dll", ".so")) and not os.path.exists(os.path.join(fx, f)):
                os.symlink(os.path.join(TRFP_BIN, f), os.path.join(fx, f))
        # hostfxr ignores a framework dir without <name>.deps.json; TRFP's lists the runtime pack
        shutil.copy(os.path.join(TRFP_BIN, "ThermoRawFileParser.deps.json"),
                    os.path.join(fx, "Microsoft.NETCore.App.deps.json"))
    return base, ver


def load_runtime():
    global _loaded
    if _loaded is not None:
        return _loaded
    # NB (M2b): RawFileReader takes a cross-process named mutex + shared memory per RAW ("Global\\CBaseIO_<path>")
    # under $HOME/.dotnet/corefx/Runtime/SharedMemory (NFS on this cluster; not moved by HOME or TMPDIR). Sporadic
    # decode hangs (futex wait, no CPU) were seen with 16-28 concurrent workers; vulcan.cohort supervises workers
    # (per-job timeout -> kill, respawn, retry) instead of relying on this path.
    os.environ.setdefault("DOTNET_SYSTEM_GLOBALIZATION_INVARIANT", "1")
    os.environ.setdefault("DOTNET_SYSTEM_GLOBALIZATION_PREDEFINED_CULTURES_ONLY", "0")
    os.environ.setdefault("DOTNET_gcServer", "0")
    os.environ.setdefault("DOTNET_GCHeapCount", "1")
    from clr_loader import get_coreclr
    from pythonnet import set_runtime
    shared = os.path.join(SDK_ROOT, "shared", "Microsoft.NETCore.App")
    if os.path.isdir(shared) and os.listdir(shared):
        root, ver = SDK_ROOT, sorted(os.listdir(shared))[-1]
    else:
        root, ver = _synthetic_dotnet_root()
    # per-process file, written atomically: concurrent decode workers sharing one path raced (a reader saw a
    # half-written file -> hostfxr InvalidConfigFile 0x80008093; M2b cohort at 24-28 workers)
    rc = os.path.join(os.environ.get("TMPDIR", "/tmp"), f"vulcan_raw_{os.getuid()}_{os.getpid()}.runtimeconfig.json")
    with open(rc + ".partial", "w") as fh:
        json.dump({"runtimeOptions": {"tfm": "net8.0", "framework": {"name": "Microsoft.NETCore.App", "version": ver},
                                      "configProperties": {"System.GC.Server": False, "System.GC.Concurrent": False,
                                                           "System.Globalization.Invariant": True,
                                                           "System.Globalization.PredefinedCulturesOnly": False}}}, fh)
    os.replace(rc + ".partial", rc)
    set_runtime(get_coreclr(runtime_config=rc, dotnet_root=root))
    import clr  # noqa: F401
    clr.AddReference(os.path.join(TRFP_BIN, "ThermoFisher.CommonCore.Data.dll"))
    clr.AddReference(os.path.join(TRFP_BIN, "ThermoFisher.CommonCore.RawFileReader.dll"))
    clr.AddReference(HELPER_DLL)
    from System.Runtime.InteropServices import GCHandle, GCHandleType
    from Vulcan.IO import RawExtractor
    try:
        os.remove(rc)
    except OSError:
        pass
    _loaded = (RawExtractor, GCHandle, GCHandleType)
    return _loaded


def _np(arr, dtype, n=None) -> np.ndarray:
    """Copy a .NET primitive array into numpy with one memmove (pinned GCHandle)."""
    _, GCHandle, GCHandleType = load_runtime()
    L = arr.Length if n is None else int(n)
    out = np.empty(L, dtype)
    if L:
        h = GCHandle.Alloc(arr, GCHandleType.Pinned)
        try:
            ctypes.memmove(out.ctypes.data, h.AddrOfPinnedObject().ToInt64(), L * out.itemsize)
        finally:
            h.Free()
    return out


_SCALARS = [  # (Batch field, store column, dtype)
    ("Scan", "scan", np.int32), ("MsLevel", "ms_level", np.uint8), ("RtMin", "rt_min", np.float64),
    ("InjTimeMs", "injection_time_ms", np.float32), ("PrecMz", "precursor_mz", np.float64),
    ("IsoTargetMz", "isolation_target_mz", np.float64), ("MonoMz", "mono_mz", np.float64),
    ("IsoLowerOffset", "isolation_lower_offset", np.float32), ("IsoUpperOffset", "isolation_upper_offset", np.float32),
    ("Charge", "precursor_charge", np.int16), ("MasterScan", "master_scan", np.int32),
    ("Activation", "activation", np.uint8), ("Nce", "nce", np.float32), ("SuppActivation", "supplemental_activation", np.uint8),
    ("SuppEnergy", "supplemental_energy", np.float32), ("Analyzer", "analyzer", np.uint8),
    ("Resolution", "resolution", np.float32), ("FaimsCv", "faims_cv", np.float32), ("Tic", "tic", np.float64),
    ("ScanLow", "scan_window_lower", np.float32), ("ScanHigh", "scan_window_upper", np.float32),
    ("Polarity", "polarity", np.int8), ("PeakSource", "peak_source", np.uint8),
]


def batch_to_cols(b) -> dict:
    n = b.N
    c = {col: _np(getattr(b, f), dt, n) for f, col, dt in _SCALARS}
    lvl = c["ms_level"]
    ms1 = lvl < 2
    c["precursor_charge_mask"] = c["precursor_charge"] <= 0
    c["master_scan_mask"] = c["master_scan"] < 0
    c["activation_mask"] = ms1
    c["supplemental_activation_mask"] = c["supplemental_activation"] == 0
    c["peak_offsets"] = _np(b.PeakOffsets, np.int64, n + 1)
    npk = int(b.NPeaks)
    c["mz"] = _np(b.Mz, np.float32, npk)
    c["intensity"] = _np(b.Intensity, np.float32, npk)
    fo = _np(b.FilterOffsets, np.int32, n + 1)
    if fo[-1] > 0:
        c["filter_offsets"] = fo
        c["filter_bytes"] = _np(b.FilterUtf8, np.uint8, int(fo[-1]))
    so = _np(b.SpsOffsets, np.int32, n + 1)
    if so[-1] > 0:
        c["sps_offsets"] = so
        c["sps"] = _np(b.Sps, np.float64)
    return c


def prefetch(path: str, chunk: int = 16 << 20) -> float:
    """Sequentially read a file into the page cache; returns seconds.

    RawFileReader memory-maps the RAW and touches it in a random-ish order; on Lustre with a
    cold cache every 4 KiB fault is an RPC (measured: 335 MB Velos file 53 s decode wall at
    ~5 core-s). One sequential pass (~1 s per 0.3-1 GB) removes that latency.
    """
    t = time.perf_counter()
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_SEQUENTIAL)
        buf = bytearray(chunk)
        mv = memoryview(buf)
        while os.readv(fd, [mv]) > 0:
            pass
    finally:
        os.close(fd)
    return time.perf_counter() - t


def convert(raw_path: str, out_path: str, run_name: str | None = None, run_id: int | None = None,
            batch_scans: int = S.BATCH_SPECTRA, compression: str | None = None,
            include_filter: bool = True, include_reference_peaks: bool = True,
            do_prefetch: bool = True, first_scan: int | None = None, last_scan: int | None = None,
            fast_trailer: bool = False, profile: bool = False, source_path: str | None = None) -> dict:
    c0 = time.process_time()
    t_pf = prefetch(raw_path) if do_prefetch else 0.0
    c_pf = time.process_time() - c0
    t0 = time.perf_counter(); c0 = time.process_time()
    RawExtractor, _, _ = load_runtime()
    t_rt = time.perf_counter() - t0; c_rt = time.process_time() - c0
    ex = RawExtractor(os.path.abspath(raw_path))
    ex.IncludeFilterString = include_filter
    ex.IncludeReferencePeaks = include_reference_peaks
    if hasattr(ex, "FastTrailer"):   # older helper builds lack these knobs
        ex.FastTrailer = bool(fast_trailer)
        ex.Profile = bool(profile)
    run_name = run_name or os.path.basename(raw_path).rsplit(".", 1)[0]
    meta = dict(run_name=run_name, vendor="Thermo", instrument_model=str(ex.InstrumentModel),
                instrument_serial=str(ex.InstrumentSerial), instrument_software=str(ex.SoftwareVersion),
                acquisition_date=str(ex.CreationDate), source_path=source_path or os.path.abspath(raw_path),
                source_size=os.path.getsize(raw_path), converter="vulcan.io.thermo_direct (RawFileReader via pythonnet)",
                converter_version="1", acquisition="unknown", peak_order="mz_ascending",
                include_reference_peaks=include_reference_peaks)
    if run_id is not None:
        meta["run_id"] = run_id
    w = StoreWriter(out_path, meta, compression=compression)
    t_read = t_write = c_read = c_write = 0.0
    lo = max(ex.FirstScan, first_scan or ex.FirstScan)
    hi = min(ex.LastScan, last_scan or ex.LastScan)
    for first in range(lo, hi + 1, batch_scans):
        a = time.perf_counter(); ac = time.process_time()
        b = ex.Read(first, min(first + batch_scans - 1, hi))
        cols = batch_to_cols(b)
        t_read += time.perf_counter() - a; c_read += time.process_time() - ac
        a = time.perf_counter(); ac = time.process_time()
        w.write(cols)
        t_write += time.perf_counter() - a; c_write += time.process_time() - ac
    prof = None
    if profile and hasattr(ex, "ProfTicks"):
        names = ["scan_stats", "scan_event_ms2", "trailer", "filter_obj", "filter_string", "peak_fetch", "peak_copy",
                 "other"]
        f = float(RawExtractor.ProfTicksPerSecond)
        prof = {n: int(ex.ProfTicks[i]) / f for i, n in enumerate(names)}
    ex.Dispose()
    a = time.perf_counter()
    extra = {"timing_s": {"prefetch": t_pf, "runtime_load": t_rt, "decode": t_read, "arrow_write": t_write},
             "cpu_s": {"prefetch": c_pf, "runtime_load": c_rt, "decode": c_read, "arrow_write": c_write},
             "fast_trailer": bool(fast_trailer)}
    if prof is not None:
        extra["dotnet_profile_s"] = prof
    final = w.close(extra)
    final["timing_s"]["close"] = time.perf_counter() - a
    return final
