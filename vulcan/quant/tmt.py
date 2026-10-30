"""Isobaric (TMT) reporter quantification (DESIGN §4 L5): SPS-MS3 linkage, reporter extraction with a per-plex
(per-run) empirical tolerance, isolation purity from MS1 features, isotopic-impurity correction.

* Linkage: an MS3 row is linked to the MS2 whose native scan equals the MS3 `master_scan` (Thermo trailer
  "Master Scan Number"), and the link is verified by the filter string (the MS3 filter's first precursor m/z
  must equal the MS2 isolation target, 1e-3 Th).  For MS2-reporter methods the MS2 itself is the reporter scan.
* Extraction: per channel the most intense centroid within tol of the (offset-corrected) reporter m/z.  The tol
  is per plex: a first pass at +-20 ppm (well below the 6.3 mDa N/C spacing at m/z 127-131 only when combined
  with the closest-peak rule; it is narrowed after) estimates the run's reporter ppm offset and robust SD; the
  final window is offset +- max(3 SD, 3 ppm), capped at 15 ppm.
* Purity: the MS1 feature that explains the MS2 precursor (precursor-hypothesis table, the feature group with the
  PSM's charge and mass) gives `isolation_fraction` = feature intensity / total centroid intensity in the MS2
  isolation window at the master MS1 scan (vulcan.ms1).  Reported, and used as an optional filter.
* Impurity correction: lot-specific matrices are not deposited for PXD016766, so the OpenMS TMT11 default matrix
  (pyopenms TMTElevenPlexQuantitationMethod, a Thermo lot sheet) is used; per spectrum NNLS solve.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.optimize import nnls

TMT11 = {  # channel: (m/z, [-2, -1, +1, +2] % impurities, affected channel indices for -2,-1,+1,+2)
    "126": (126.127726, [0.0, 0.0, 8.6, 0.3], [-1, -1, 2, 4]),
    "127N": (127.124761, [0.0, 0.1, 7.8, 0.1], [-1, -1, 3, 5]),
    "127C": (127.131081, [0.0, 0.8, 6.9, 0.1], [-1, 0, 4, 6]),
    "128N": (128.128116, [0.0, 7.4, 7.4, 0.0], [-1, 1, 5, 7]),
    "128C": (128.134436, [0.0, 1.5, 6.2, 0.2], [0, 2, 6, 8]),
    "129N": (129.131471, [0.0, 1.5, 5.7, 0.1], [1, 3, 7, 9]),
    "129C": (129.137790, [0.0, 2.6, 4.8, 0.0], [2, 4, 8, 10]),
    "130N": (130.134825, [0.0, 2.2, 4.6, 0.0], [3, 5, 9, -1]),
    "130C": (130.141145, [0.0, 2.8, 4.5, 0.1], [4, 6, 10, -1]),
    "131N": (131.138180, [0.1, 2.9, 3.8, 0.0], [5, 7, -1, -1]),
    "131C": (131.144500, [0.0, 3.9, 2.8, 0.0], [6, 8, -1, -1]),
}
CHANNELS = list(TMT11)
REPORTER_MZ = np.array([v[0] for v in TMT11.values()])


def impurity_matrix(spec=TMT11) -> np.ndarray:
    """A[obs, src]: fraction of channel src's signal observed at channel obs (columns sum to <= 1)."""
    n = len(spec)
    A = np.zeros((n, n))
    for j, (_, imp, aff) in enumerate(spec.values()):
        A[j, j] = 1.0 - sum(imp) / 100.0
        for pct, k in zip(imp, aff):
            if k >= 0:
                A[k, j] += pct / 100.0
    return A


def correct(X: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Per-row non-negative least squares X_obs = A @ X_true."""
    out = np.zeros_like(X)
    for i in range(len(X)):
        x = X[i]
        if not np.isfinite(x).all() or x.sum() <= 0:
            out[i] = np.where(np.isfinite(x), x, 0)
            continue
        out[i] = nnls(A, x)[0]
    return out


_MS3_PARENT = re.compile(r"ms3\s+([0-9.]+)@")
_MS2_PARENT = re.compile(r"ms2\s+([0-9.]+)@")


def link_ms3(store_rows: pd.DataFrame) -> pd.DataFrame:
    """store_rows: columns index, scan, ms_level, master_scan, filter_string.  Returns MS2 rows with the store index
    of their reporter scan (`rep_index`; MS3 for SPS-MS3, the MS2 itself if the run has no MS3) and a
    `link_ok` flag (filter-string parent check)."""
    ms2 = store_rows[store_rows["ms_level"] == 2][["index", "scan", "filter_string"]].copy()
    ms3 = store_rows[store_rows["ms_level"] == 3][["index", "master_scan", "filter_string"]].copy()
    if not len(ms3):
        ms2["rep_index"] = ms2["index"]
        ms2["link_ok"] = True
        return ms2
    ms3 = ms3.dropna(subset=["master_scan"])
    ms3["master_scan"] = ms3["master_scan"].astype(np.int64)
    ms3["p3"] = ms3["filter_string"].str.extract(_MS3_PARENT)[0].astype(float)
    ms3 = ms3.drop_duplicates("master_scan")
    ms2["p2"] = ms2["filter_string"].str.extract(_MS2_PARENT)[0].astype(float)
    m = ms2.merge(ms3[["index", "master_scan", "p3"]].rename(columns={"index": "rep_index"}),
                  left_on="scan", right_on="master_scan", how="left")
    m["link_ok"] = m["rep_index"].notna() & (np.abs(m["p2"] - m["p3"]) < 1e-3)
    return m.drop(columns=["master_scan"])


def extract(mz_flat, int_flat, offsets, rows, tol_ppm, offset_ppm=0.0):
    """Reporter intensities (n_rows x 11) and ppm errors for the spectra `rows` (store row indices into the peak
    arrays).  Closest-in-intensity rule: the most intense centroid within the window."""
    n = len(rows)
    I = np.zeros((n, len(REPORTER_MZ)))
    E = np.full((n, len(REPORTER_MZ)), np.nan)
    target = REPORTER_MZ * (1 + offset_ppm * 1e-6)
    lo = target * (1 - tol_ppm * 1e-6)
    hi = target * (1 + tol_ppm * 1e-6)
    for i, r in enumerate(rows):
        a, b = offsets[r], offsets[r + 1]
        mz = mz_flat[a:b]
        it = int_flat[a:b]
        s = np.searchsorted(mz, lo)
        e = np.searchsorted(mz, hi, side="right")
        for c in range(len(target)):
            if e[c] > s[c]:
                k = s[c] + int(np.argmax(it[s[c]:e[c]]))
                I[i, c] = it[k]
                E[i, c] = (mz[k] - REPORTER_MZ[c]) / REPORTER_MZ[c] * 1e6
    return I, E


def calibrate_tolerance(E: np.ndarray, I: np.ndarray, k: float = 3.0, tmin: float = 3.0, tmax: float = 15.0):
    """Per-plex reporter mass offset and tolerance from a first-pass extraction (intensity-weighted)."""
    e = E[np.isfinite(E) & (I > 0)]
    if len(e) < 100:
        return 0.0, 10.0, {"n": int(len(e))}
    off = float(np.median(e))
    sd = float(1.4826 * np.median(np.abs(e - off)))
    tol = float(np.clip(k * sd, tmin, tmax))
    return off, tol, {"n": int(len(e)), "offset_ppm": off, "sd_ppm": sd, "tol_ppm": tol}
