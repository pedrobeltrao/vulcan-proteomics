"""Torch-free host-side pieces of the closed search (split out of vulcan.gpu.engine in M2b so that CPU decode
workers can build search-ready MS2 arrays and precursor hypotheses without importing torch). engine re-exports all of
these names, so `from vulcan.gpu.engine import SearchParams, Spectra, make_hypotheses` keeps working."""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

PROTON = np.float32(1.0072764)
NEUTRON = np.float32(1.00335)


@dataclass
class SearchParams:
    precursor_ppm: tuple = (-10.0, 10.0)
    fragment_ppm: tuple = (-20.0, 20.0)
    fragment_unit: str = "ppm"              # 3c: "ppm" (OT-HCD; default, bit-identical to E2) | "da" (IT-CID):
                                            # fragment_ppm then holds Da bounds (Sage Tolerance::Da), e.g. (-0.5, 0.5)
    isotope_errors: tuple = (0, 2)
    precursor_charge: tuple = (2, 4)        # tried when the instrument gives none
    max_peaks: int = 150
    min_peaks: int = 15
    fragment_min_mz: float = 150.0
    fragment_max_mz: float = 2000.0
    min_ion_index: int = 2
    max_fragment_charge: int | None = 1
    min_matched_peaks: int = 4
    topk: int = 500                         # E2 sweep: recall T/D 99.9/99.8% vs exhaustive on D0 (K=50: 96.8/93.3%)
    k1_threshold: int = 1
    pair_closure: bool = True
    report_psms: int = 1
    peak_mode: str = "global"               # "global" (Sage) | "window"
    window_th: float = 100.0
    window_k: int = 15
    # K1 layout
    chunk_size: int = 8192
    bin_ppm: float = 20.0
    tile: int = 1024
    k2_threads: int = 128
    max_slice_entries: int = 2 ** 31 - 1

    def to_dict(self):
        return asdict(self)


@dataclass
class Spectra:
    """MS2 spectra of one run (from the Arrow store), host side."""
    run: str
    scan: np.ndarray
    index: np.ndarray
    rt: np.ndarray
    precursor_mz: np.ndarray       # f64
    charge: np.ndarray             # i8, 0 = unknown
    mz: np.ndarray                 # f32 flat
    intensity: np.ndarray          # f32 flat
    offsets: np.ndarray            # i64 n+1

    def __len__(self):
        return len(self.scan)


def make_hypotheses(sp: Spectra, valid: np.ndarray, idx_mass32: np.ndarray, p: SearchParams):
    """Precursor hypotheses (spectrum, charge, isotope) and their candidate row windows (Sage f32 arithmetic)."""
    s_ids, zs, isos = [], [], []
    spec = np.nonzero(valid & np.isfinite(sp.precursor_mz))[0]
    zi = sp.charge[spec].astype(np.int64)
    known = zi > 0
    zr = np.arange(p.precursor_charge[0], p.precursor_charge[1] + 1)
    for iso in range(p.isotope_errors[0], p.isotope_errors[1] + 1):
        s_ids.append(spec[known]); zs.append(zi[known]); isos.append(np.full(known.sum(), iso))
        for z in zr:
            s_ids.append(spec[~known]); zs.append(np.full((~known).sum(), z)); isos.append(np.full((~known).sum(), iso))
    s = np.concatenate(s_ids)
    z = np.concatenate(zs)
    iso = np.concatenate(isos)
    o = np.lexsort((iso, z, s))
    s, z, iso = s[o], z[o], iso[o]
    mz = (sp.precursor_mz[s].astype(np.float32) - PROTON).astype(np.float32)
    pm = (mz * z.astype(np.float32)).astype(np.float32)
    M = (pm - (iso.astype(np.float32) * NEUTRON).astype(np.float32)).astype(np.float32)
    plo, phi = np.float32(p.precursor_ppm[0]), np.float32(p.precursor_ppm[1])
    lo = (M + (M * plo) / np.float32(1e6)).astype(np.float32)
    hi = (M + (M * phi) / np.float32(1e6)).astype(np.float32)
    rlo = np.searchsorted(idx_mass32, lo, side="left")
    rhi = np.searchsorted(idx_mass32, hi, side="right")
    mfc = p.max_fragment_charge
    maxfz = np.maximum(np.minimum(z, mfc + 1) if mfc else z, 2).astype(np.int8)
    return {"spec": s.astype(np.int32), "z": z.astype(np.int8), "iso": iso.astype(np.int8), "pm": pm, "M": M,
            "lo": rlo.astype(np.int64), "hi": rhi.astype(np.int64), "maxfz": maxfz}


def pre_hypotheses(sp: Spectra, gi, p: SearchParams):
    """Host-side precursor hypotheses, computable before K0 (e.g. in a loader thread; vulcan.cohort). Valid for the
    'global' peak view, where K0 keeps min(#peaks in [fragment_min_mz, fragment_max_mz], max_peaks) peaks.
    Returns (valid, hyp); search(pre_hyp=...) checks `valid` against K0's and recomputes on any mismatch."""
    assert p.peak_mode == "global"
    mz = sp.mz
    ok = (mz >= p.fragment_min_mz) & (mz <= p.fragment_max_mz) & ~np.isnan(mz)
    cs = np.zeros(len(ok) + 1, np.int64)
    np.cumsum(ok, out=cs[1:])
    n_in = cs[sp.offsets[1:]] - cs[sp.offsets[:-1]]
    valid = np.minimum(n_in, p.max_peaks) >= p.min_peaks
    return valid, make_hypotheses(sp, valid, gi.mass32_np, p)
