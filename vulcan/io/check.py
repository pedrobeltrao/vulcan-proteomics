"""Correctness checks for stores (used by scripts/ingest and tests/io).

compare_stores(a, b): field-by-field agreement of two stores of the same run (e.g. direct vs
TRFP-mzML path). check_against_pyteomics(store, mzml, n): independent re-read of sampled
spectra with pyteomics (lxml-based reference reader).
"""
from __future__ import annotations

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from . import store as St

_SCALAR_TOL = {  # field -> abs tolerance (None = exact)
    "ms_level": None, "scan": None, "precursor_charge": None, "master_scan": None, "activation": None,
    "rt_min": 1e-6, "injection_time_ms": 1e-3, "isolation_target_mz": 1e-6, "isolation_lower_offset": 1e-4,
    "isolation_upper_offset": 1e-4, "nce": 1e-4, "scan_window_lower": 1e-3, "scan_window_upper": 1e-3,
    "analyzer": None, "polarity": None, "faims_cv": 1e-3,
}


def _col(t, name):
    c = t.column(name).combine_chunks()
    out = c.to_numpy(zero_copy_only=False).astype(np.float64)  # nulls -> NaN
    if c.null_count and not pa.types.is_floating(c.type):
        out[np.asarray(c.is_null())] = -999.0
    return out


def compare_stores(path_a: str, path_b: str, fields=None) -> dict:
    ta, tb = St.read_table(path_a), St.read_table(path_b)
    res = {"n_a": ta.num_rows, "n_b": tb.num_rows, "fields": {}}
    if ta.num_rows != tb.num_rows:
        res["ok"] = False
        return res
    ok = True
    for f in (fields or _SCALAR_TOL):
        a, b = _col(ta, f), _col(tb, f)
        both_nan = np.isnan(a) & np.isnan(b)
        tol = _SCALAR_TOL.get(f)
        diff = np.abs(a - b)
        bad = ~both_nan & ~(diff <= (tol or 0))
        res["fields"][f] = int(bad.sum())
        ok &= not bad.any()
    # filter strings
    fa, fb = ta.column("filter_string"), tb.column("filter_string")
    if fa.null_count < len(fa) and fb.null_count < len(fb):
        res["fields"]["filter_string"] = int(pc.sum(pc.not_equal(fa, fb)).as_py() or 0)
        ok &= res["fields"]["filter_string"] == 0
    # precursor m/z: report agreement of b's selected ion with a's isolation target / mono
    ms2 = _col(ta, "ms_level") >= 2
    pa_, pb_ = _col(ta, "precursor_mz")[ms2], _col(tb, "precursor_mz")[ms2]
    iso_a = _col(ta, "isolation_target_mz")[ms2]
    res["precursor_mz_equal_frac"] = float(np.mean(np.abs(pa_ - pb_) < 1e-5))
    res["b_precursor_eq_a_isolation_target_frac"] = float(np.mean(np.abs(iso_a - pb_) < 1e-5))
    res["precursor_mz_max_abs_diff"] = float(np.nanmax(np.abs(pa_ - pb_))) if len(pa_) else 0.0
    # peaks
    A, B = St.peaks_numpy(ta), St.peaks_numpy(tb)
    res["peak_count_mismatch_spectra"] = int(np.sum(np.diff(A.offsets) != np.diff(B.offsets)))
    if len(A.mz) == len(B.mz):
        res["mz_max_abs_diff"] = float(np.max(np.abs(A.mz - B.mz))) if len(A.mz) else 0.0
        rel = np.abs(A.intensity - B.intensity) / np.maximum(np.abs(B.intensity), 1e-30)
        res["intensity_max_rel_diff"] = float(rel.max()) if len(rel) else 0.0
        ok &= res["mz_max_abs_diff"] == 0.0 and res["intensity_max_rel_diff"] <= 1e-6
    else:
        ok = False
    ok &= res["peak_count_mismatch_spectra"] == 0
    res["ok"] = bool(ok)
    return res


def check_against_pyteomics(store_path: str, mzml_path: str, n: int = 50, seed: int = 11) -> dict:
    """Re-read n sampled spectra (half MS1/half MS2) with pyteomics and compare."""
    from pyteomics import mzml as pmz

    t = St.read_table(store_path)
    P = St.peaks_numpy(t)
    lvl = _col(t, "ms_level")
    scans = _col(t, "scan").astype(int)
    rng = np.random.default_rng(seed)
    idx1 = np.nonzero(lvl == 1)[0]
    idx2 = np.nonzero(lvl >= 2)[0]
    pick = np.concatenate([rng.choice(idx1, min(n // 2, len(idx1)), replace=False),
                           rng.choice(idx2, min(n - n // 2, len(idx2)), replace=False)])
    prec = _col(t, "precursor_mz")
    iso = _col(t, "isolation_target_mz")
    charge = _col(t, "precursor_charge")
    rt = _col(t, "rt_min")
    nce = _col(t, "nce")
    bad = []
    with pmz.PreIndexedMzML(mzml_path) as r:
        for i in sorted(pick.tolist()):
            sid = f"controllerType=0 controllerNumber=1 scan={scans[i]}"
            s = r.get_by_id(sid)
            mz = np.asarray(s["m/z array"], np.float32)
            it = np.asarray(s["intensity array"], np.float32)
            a, b = P.offsets[i], P.offsets[i + 1]
            errs = []
            if s["ms level"] != lvl[i]:
                errs.append("ms_level")
            if len(mz) != b - a or not np.array_equal(mz, P.mz[a:b]) or not np.allclose(it, P.intensity[a:b], rtol=1e-6):
                errs.append("peaks")
            if abs(float(s["scanList"]["scan"][0]["scan start time"]) - rt[i]) > 1e-5:
                errs.append("rt")
            if lvl[i] >= 2:
                p = s["precursorList"]["precursor"][0]
                ion = p["selectedIonList"]["selectedIon"][0]
                sel = float(ion["selected ion m/z"])
                if abs(sel - prec[i]) > 1e-4 and abs(sel - iso[i]) > 1e-4:
                    errs.append("precursor_mz")
                z = int(ion.get("charge state", 0) or 0)
                if (z if z > 0 else -999) != charge[i]:
                    errs.append("charge")
                ce = p.get("activation", {}).get("collision energy")
                if ce is not None and abs(float(ce) - nce[i]) > 1e-3:
                    errs.append("nce")
            if errs:
                bad.append({"scan": int(scans[i]), "errors": errs})
    return {"n_checked": int(len(pick)), "n_bad": len(bad), "bad": bad[:10]}
