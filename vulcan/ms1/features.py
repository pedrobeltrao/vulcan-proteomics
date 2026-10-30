"""MS1 feature finding on a Vulcan store: load MS1 centroids -> hills -> isotope features.

    res = find_features("/tmp/vulcan_store/e1/01A_direct.vspec.arrow", FeatureParams(), device="cuda")
    res.features            # pyarrow.Table, schema FEATURE_SCHEMA (vulcan/ms1/SCHEMA.md §2)
    res.peak_feature        # int32 per flat MS1 peak: feature id or -1 (used by hypotheses.py)
"""
from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from ..io import store
from .averagine import PROTON
from .hills import Hills, detect_hills
from .isotopes import group_isotopes

FEATURE_FINDER_VERSION = "vulcan.ms1 0.2"


@dataclass
class FeatureParams:
    htol_ppm: float = 8.0          # hill linking tolerance (biosaur2 -htol 8)
    itol_ppm: float = 8.0          # isotope tolerance (biosaur2 -itol 8)
    minlh: int = 2                 # min hill length in MS1 scans (biosaur2 -minlh 2)
    hvf: float = 1.3               # hill valley factor (biosaur2 -hvf 1.3)
    gap: int = 0                   # allowed missing MS1 scans inside a hill
    charges: tuple = (1, 2, 3, 4, 5, 6)
    max_iso: int = 10
    min_trace_cos: float = 0.6     # isotope trace vs mono trace cosine
    min_overlap_for_cos: int = 1   # overlaps shorter than this skip the trace-cosine test (1 = never skip;
                                   # 3 accepted 2-scan overlaps blindly: +50% noise features vs biosaur2, see NOTES)
    min_avg_cos: float = 0.8       # averagine cosine of isotope apex intensities
    min_explained: float = 0.5     # observed isotopes must hold >= this share of the averagine envelope (biosaur2)
    full_norm: bool = True         # trace cosine with full-hill norms (biosaur2) instead of overlap-only
    charge_ban: bool = True        # drop z if a multiple charge k*z of the same mono has >= as many isotopes
    min_intensity: float = 0.0
    mz_min: float = 0.0            # restrict MS1 peaks (biosaur2 default 350-1500; we keep all)
    mz_max: float = 1e9


FEATURE_SCHEMA = pa.schema([
    pa.field("feature_id", pa.uint32(), nullable=False),
    pa.field("run_id", pa.uint32(), nullable=False),
    pa.field("mono_mz", pa.float64(), nullable=False),
    pa.field("charge", pa.int8(), nullable=False),
    pa.field("neutral_mass", pa.float64(), nullable=False),
    pa.field("rt_apex", pa.float32(), nullable=False),
    pa.field("rt_start", pa.float32(), nullable=False),
    pa.field("rt_end", pa.float32(), nullable=False),
    pa.field("ms1_apex", pa.uint32(), nullable=False),     # MS1 ordinal of the apex
    pa.field("ms1_start", pa.uint32(), nullable=False),
    pa.field("ms1_end", pa.uint32(), nullable=False),
    pa.field("intensity_apex", pa.float32(), nullable=False),       # sum of isotope-hill apex intensities
    pa.field("intensity_mono_apex", pa.float32(), nullable=False),  # mono hill apex (biosaur2 intensityApex, -iuse 0)
    pa.field("intensity_area", pa.float32(), nullable=False),       # sum of isotope-hill areas (intensity x min)
    pa.field("n_isotopes", pa.uint8(), nullable=False),
    pa.field("iso_cos", pa.float32(), nullable=False),              # averagine cosine of isotope apex intensities
    pa.field("trace_cos", pa.float32(), nullable=False),            # mean isotope-vs-mono trace cosine
    pa.field("n_scans", pa.uint16(), nullable=False),               # mono hill length
])


@dataclass
class MS1Data:
    run_id: int
    ms1_index: np.ndarray   # store row index per MS1 ordinal
    ms1_scan: np.ndarray    # native scan number per MS1 ordinal
    ms1_rt: np.ndarray      # minutes
    mz: np.ndarray          # float32 flat peaks, mz ascending within a scan
    intensity: np.ndarray
    off: np.ndarray         # int64 (n_ms1+1)


@dataclass
class FeatureResult:
    features: pa.Table
    hills: Hills
    ms1: MS1Data
    peak_feature: np.ndarray   # int32 per flat MS1 peak
    peak_iso: np.ndarray       # int8 isotope index per flat MS1 peak (-1 if none)
    timings: dict = field(default_factory=dict)
    params: dict = field(default_factory=dict)


def load_ms1(path: str) -> MS1Data:
    meta = store.run_meta(path)
    t = store.read_table(path, columns=["index", "ms_level", "rt_min", "scan", "mz", "intensity"])
    t = t.filter(pc.equal(t.column("ms_level"), 1))
    pk = store.peaks_numpy(t)
    idx = t.column("index").to_numpy().astype(np.int64)
    scan = t.column("scan").to_numpy().astype(np.int64)
    rt = t.column("rt_min").to_numpy().astype(np.float64)
    mz, inten, off = pk.mz, pk.intensity, pk.offsets
    lens = np.diff(off)
    ordn = np.repeat(np.arange(len(lens), dtype=np.int64), lens)
    if meta.get("peak_order") != "mz_ascending" or (len(mz) > 1 and np.any(
            (np.diff(mz) < 0) & (ordn[1:] == ordn[:-1]))):
        o = np.lexsort((mz, ordn))
        mz, inten = mz[o], inten[o]
    return MS1Data(int(meta["run_id"]), idx, scan, rt, mz, inten, off)


MS1_BACKENDS = ("cpu", "gpu")


def gpu_supported(p: FeatureParams, device: str = "cpu", exact_reference: bool = False) -> bool:
    """The E11 GPU backend covers the default production configuration (gap 0, full-norm trace cosines)."""
    from .gpu import MAX_ISO_GPU
    return (not exact_reference) and device == "cpu" and p.gap == 0 and p.full_norm and 1 <= p.max_iso <= MAX_ISO_GPU


def find_features(path_or_ms1, params: FeatureParams | None = None, device: str = "cpu",
                  exact_reference: bool = False, backend: str = "cpu", hill_peaks: bool = True) -> FeatureResult:
    """exact_reference: original isotope kernels (vulcan.ms1.isotopes; the default fast path is bitwise equal).
    backend: "cpu" (default, the reference) or "gpu" (E11, vulcan.ms1.gpu: bitwise-equal tables; falls back to the
    CPU for configurations it does not cover, see gpu_supported; timings["backend"] records what ran).
    hill_peaks=False (GPU backend only): leave res.hills' per-peak arrays (p_scan, p_mz, p_int, p_src) as None instead of
    copying them to the host; nothing downstream (feature table, hypotheses, cohort) reads them."""
    if backend not in MS1_BACKENDS:
        raise ValueError(f"ms1 backend must be one of {MS1_BACKENDS}, got {backend!r}")
    p = params or FeatureParams()
    if backend == "gpu" and gpu_supported(p, device, exact_reference):
        return _find_features_gpu(path_or_ms1, p, hill_peaks)
    tm = {}
    t0 = time.perf_counter()
    ms1 = load_ms1(path_or_ms1) if isinstance(path_or_ms1, str) else path_or_ms1
    tm["load"] = time.perf_counter() - t0
    lens = np.diff(ms1.off)
    scan = np.repeat(np.arange(len(lens), dtype=np.int32), lens)
    sel = (ms1.mz >= p.mz_min) & (ms1.mz <= p.mz_max)
    src_all = np.flatnonzero(sel) if not sel.all() else None
    mz_s = ms1.mz[sel] if src_all is not None else ms1.mz
    in_s = ms1.intensity[sel] if src_all is not None else ms1.intensity
    sc_s = scan[sel] if src_all is not None else scan
    h = detect_hills(sc_s, mz_s, in_s, ms1.ms1_rt, htol_ppm=p.htol_ppm, minlh=p.minlh, hvf=p.hvf, gap=p.gap,
                     min_intensity=p.min_intensity, device=device, timings=tm, reference=exact_reference)
    if src_all is not None:
        h.p_src = src_all[h.p_src]
    t1 = time.perf_counter()
    g = group_isotopes(h, charges=p.charges, itol_ppm=p.itol_ppm, max_iso=p.max_iso,
                       min_trace_cos=p.min_trace_cos, min_overlap_for_cos=p.min_overlap_for_cos,
                       min_avg_cos=p.min_avg_cos, min_explained=p.min_explained, full_norm=p.full_norm,
                       use_charge_ban=p.charge_ban, exact_reference=exact_reference)
    tm["isotopes"] = time.perf_counter() - t1
    t2 = time.perf_counter()
    feats, peak_feature, peak_iso = _feature_table(ms1, h, g)
    tm["table"] = time.perf_counter() - t2
    tm["total"] = time.perf_counter() - t0
    return _result(p, ms1, h, g, feats, peak_feature, peak_iso, tm)


def _result(p, ms1, h, g, feats, peak_feature, peak_iso, tm) -> FeatureResult:
    params_d = asdict(p)
    params_d["charges"] = list(p.charges)
    params_d["version"] = FEATURE_FINDER_VERSION
    params_d["n_hills"] = int(h.n)
    params_d["n_candidates"] = g["n_candidates"]
    feats = feats.replace_schema_metadata({"vulcan_ms1": json.dumps(params_d)})
    return FeatureResult(feats, h, ms1, peak_feature, peak_iso, tm, params_d)


def _find_features_gpu(path_or_ms1, p: FeatureParams, hill_peaks: bool = True) -> FeatureResult:
    """E11 GPU backend: same steps and outputs as the CPU path (bitwise); see vulcan/ms1/gpu.py."""
    from . import gpu
    from . import isotopes as IS
    tm = {}
    t0 = time.perf_counter()
    ms1 = load_ms1(path_or_ms1) if isinstance(path_or_ms1, str) else path_or_ms1
    tm["load"] = time.perf_counter() - t0
    lens = np.diff(ms1.off)
    scan = np.repeat(np.arange(len(lens), dtype=np.int32), lens)
    sel = (ms1.mz >= p.mz_min) & (ms1.mz <= p.mz_max)
    src = np.flatnonzero(sel) if not sel.all() else None
    mz_s = ms1.mz[sel] if src is not None else ms1.mz
    in_s = ms1.intensity[sel] if src is not None else ms1.intensity
    sc_s = scan[sel] if src is not None else scan
    if p.min_intensity > 0:
        keep = in_s >= p.min_intensity
        src = (np.arange(len(mz_s), dtype=np.int64) if src is None else src)[keep]
        mz_s, in_s, sc_s = mz_s[keep], in_s[keep], sc_s[keep]
    tg = time.perf_counter()
    fallback = None
    early = False
    with gpu.session_lock, gpu.gate():
        tm["gpu_gate_wait"] = time.perf_counter() - tg
        S = gpu.session()
        free_b, _ = S.mem_info()
        if free_b < S.bytes_needed(len(mz_s)):      # not enough device memory now: the (bitwise-equal) CPU path
            fallback = f"gpu_mem: free {free_b >> 20} MiB < {S.bytes_needed(len(mz_s)) >> 20} MiB needed"
        else:
            try:
                h = S.hills(sc_s, mz_s, in_s, ms1.ms1_rt, src, htol_ppm=p.htol_ppm, minlh=p.minlh, hvf=p.hvf,
                            timings=tm, peaks=hill_peaks)
                t1 = time.perf_counter()
                contiguous = bool(np.all(h.h_len.astype(np.int64) == h.h_end.astype(np.int64)
                                         - h.h_start.astype(np.int64) + 1))
                kw = dict(charges=p.charges, itol_ppm=p.itol_ppm, max_iso=p.max_iso, min_trace_cos=p.min_trace_cos,
                          min_overlap_for_cos=p.min_overlap_for_cos, min_avg_cos=p.min_avg_cos,
                          min_explained=p.min_explained)
                early = gpu._gate is not None     # shared GPU (cohort workers): free the device before host work
                if contiguous:      # always for gap 0; otherwise the CPU reference kernels (as group_isotopes does)
                    raw = S.isotopes_raw(use_charge_ban=p.charge_ban, scan_bin=IS.SCAN_BIN, timings=tm, **kw)
                    g = None
                else:
                    if h.p_int is None:
                        S.fetch_peaks(h)
                    g = IS.group_isotopes(h, full_norm=p.full_norm, use_charge_ban=p.charge_ban, **kw)
                if early:
                    if h.p_src is None:
                        S.fetch_src(h)
                else:
                    if g is None:
                        g = gpu.finish_isotopes(raw, h, min_avg_cos=p.min_avg_cos, min_explained=p.min_explained,
                                                timings=tm)
                    tm["isotopes"] = time.perf_counter() - t1
                    t2 = time.perf_counter()
                    feats, peak_feature, peak_iso = _feature_table(ms1, h, g, peak_map=S.peak_map)
                    tm["table"] = time.perf_counter() - t2
                tm["gpu_ms"] = S.gpu_ms()
            except RuntimeError as e:   # e.g. cudaErrorMemoryAllocation (2) under contention: same result on the CPU
                if "error 2" not in str(e):
                    raise
                fallback = f"gpu_error: {e}"
            finally:
                S.release()
    tm["gpu_held"] = time.perf_counter() - tg - tm["gpu_gate_wait"]
    if fallback is None and early:          # device released: greedy + feature table (CPU peak map) on the host
        if g is None:
            g = gpu.finish_isotopes(raw, h, min_avg_cos=p.min_avg_cos, min_explained=p.min_explained, timings=tm)
        tm["isotopes"] = time.perf_counter() - t1
        t2 = time.perf_counter()
        feats, peak_feature, peak_iso = _feature_table(ms1, h, g)
        tm["table"] = time.perf_counter() - t2
        if not hill_peaks:
            h.p_src = None
    if fallback is not None:
        r = find_features(ms1, p)
        r.timings["backend"] = "cpu_fallback"
        r.timings["fallback_reason"] = fallback
        r.timings["gpu_gate_wait"] = tm["gpu_gate_wait"]
        return r
    tm["total"] = time.perf_counter() - t0
    tm["backend"] = "gpu"
    return _result(p, ms1, h, g, feats, peak_feature, peak_iso, tm)


def _feature_table(ms1: MS1Data, h: Hills, g: dict, peak_map=None):
    """peak_map: optional f(hill_feat, hill_iso, n_peaks) -> (peak_feature, peak_iso) (the GPU backend)."""
    n_f = len(g["mono_hill"])
    nn = g["n_isotopes"].astype(np.int64)
    hills = g["hills"].copy()
    kk = np.arange(hills.shape[1])[None, :]
    valid = kk < nn[:, None]
    hills[~valid] = -1
    hv = np.where(valid, hills, 0)
    mono = g["mono_hill"]
    z = g["charge"].astype(np.int64)
    apex_int = np.where(valid, h.h_apex_int[hv], 0.0)
    area = np.where(valid, h.h_area[hv], 0.0)
    # apex: apex of the most intense isotope hill
    top = hv[np.arange(n_f), np.argmax(apex_int, axis=1)] if n_f else np.empty(0, np.int64)
    apex_scan = h.h_apex_scan[top]
    start = h.h_start[mono]
    end = h.h_end[mono]
    tc = np.where(valid[:, 1:], g["trace_cos"][:, 1:], 0.0).sum(1) / np.maximum(nn - 1, 1)
    mono_mz = h.h_mz[mono]
    # order features by (apex rt, mono mz) for a deterministic, data-derived feature_id
    o = np.lexsort((mono_mz, apex_scan))
    cols = {
        "feature_id": np.arange(n_f, dtype=np.uint32),
        "run_id": np.full(n_f, ms1.run_id, np.uint32),
        "mono_mz": mono_mz[o],
        "charge": z[o].astype(np.int8),
        "neutral_mass": ((mono_mz - PROTON) * z)[o],
        "rt_apex": ms1.ms1_rt[apex_scan][o].astype(np.float32),
        "rt_start": ms1.ms1_rt[start][o].astype(np.float32),
        "rt_end": ms1.ms1_rt[end][o].astype(np.float32),
        "ms1_apex": apex_scan[o].astype(np.uint32),
        "ms1_start": start[o].astype(np.uint32),
        "ms1_end": end[o].astype(np.uint32),
        "intensity_apex": apex_int.sum(1)[o].astype(np.float32),
        "intensity_mono_apex": h.h_apex_int[mono][o].astype(np.float32),
        "intensity_area": area.sum(1)[o].astype(np.float32),
        "n_isotopes": nn[o].astype(np.uint8),
        "iso_cos": g["iso_cos"][o].astype(np.float32),
        "trace_cos": tc[o].astype(np.float32),
        "n_scans": np.minimum(h.h_len[mono][o], 65535).astype(np.uint16),
    }
    table = pa.table(cols, schema=FEATURE_SCHEMA)
    # per-peak feature / isotope index
    hill_feat = np.full(h.n, -1, np.int32)
    hill_iso = np.full(h.n, -1, np.int8)
    inv = np.empty(n_f, np.int64)
    inv[o] = np.arange(n_f)
    fi, ki = np.nonzero(valid)
    hill_feat[hills[fi, ki]] = inv[fi]
    hill_iso[hills[fi, ki]] = ki
    if peak_map is not None:
        peak_feature, peak_iso = peak_map(hill_feat, hill_iso, len(ms1.mz))
        return table, peak_feature, peak_iso
    lens = np.diff(h.off)
    p_hill = np.repeat(np.arange(h.n, dtype=np.int64), lens)
    peak_feature = np.full(len(ms1.mz), -1, np.int32)
    peak_iso = np.full(len(ms1.mz), -1, np.int8)
    peak_feature[h.p_src] = hill_feat[p_hill]
    peak_iso[h.p_src] = hill_iso[p_hill]
    return table, peak_feature, peak_iso


def write_features(table: pa.Table, path: str) -> None:
    with pa.OSFile(path + ".partial", "wb") as sink, pa.ipc.new_file(sink, table.schema) as w:
        w.write_table(table)
    import os
    os.replace(path + ".partial", path)
