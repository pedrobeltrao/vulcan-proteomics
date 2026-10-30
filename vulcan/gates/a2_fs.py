"""Amendment A2.4 fragment-shifted nulls (FS) and A2.2(e) single-precursor hypotheses.

FS (fixed here; A2.11 item 5):
  * donors = ALL stratum-A spectra (non-human, certified, identity-screened), after admission;
  * every MS2 peak m/z is shifted by delta = +2/7 Th, computed as float32(mz) + float32(2/7) (the store's dtype);
    intensities, peak order and peak count unchanged;
  * precursor = the CERTIFIED precursor (mono m/z, z) of the donor; RT = the donor's RT (the same A1.3 quantile map
    to the host is then applied in the mixture builder, as for every species null);
  * FS spectra get new scan numbers (donor scan + FS_SCAN_OFFSET) so they never collide with the donor in a mixture.
(e) single precursor hypothesis: an A or FS spectrum carries exactly one precursor hypothesis, (certified mono m/z,
  certified z, isotope error 0).  Every other spectrum (host, B, A1 offsets) keeps the production hypotheses.  The
  identity/family screens use the wider W (production isotope set), so the screened set covers the searched set.
"""
from __future__ import annotations

import numpy as np

FS_DELTA_TH = 2.0 / 7.0
FS_SCAN_OFFSET = 10_000_000


def fragment_shift(sp, idx, mz_mono, z, delta: float = FS_DELTA_TH):
    """FS copies of spectra ``idx`` of ``sp`` (vulcan.gpu.host.Spectra).  mz_mono / z: certified precursor per idx."""
    from vulcan.gpu.host import Spectra
    idx = np.asarray(idx, np.int64)
    starts, ends = sp.offsets[idx], sp.offsets[idx + 1]
    lens = ends - starts
    take = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)]) if idx.size else np.zeros(0, np.int64)
    mz = (sp.mz[take].astype(np.float32) + np.float32(delta)).astype(np.float32)
    off = np.zeros(idx.size + 1, np.int64)
    np.cumsum(lens, out=off[1:])
    return Spectra(run=sp.run + "_FS", scan=(sp.scan[idx].astype(np.int64) + FS_SCAN_OFFSET).astype(sp.scan.dtype),
                   index=sp.index[idx].copy(), rt=sp.rt[idx].copy(),
                   precursor_mz=np.asarray(mz_mono, np.float64).copy(), charge=np.asarray(z).astype(sp.charge.dtype),
                   mz=mz, intensity=sp.intensity[take].astype(np.float32).copy(), offsets=off)


def certified_precursor(sp, idx, mz_mono, z):
    """Copy of spectra ``idx`` with the precursor replaced by the certified (mono m/z, z) -- A spectra in mixtures."""
    from vulcan.gpu.host import Spectra
    idx = np.asarray(idx, np.int64)
    starts, ends = sp.offsets[idx], sp.offsets[idx + 1]
    take = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)]) if idx.size else np.zeros(0, np.int64)
    off = np.zeros(idx.size + 1, np.int64)
    np.cumsum(ends - starts, out=off[1:])
    return Spectra(run=sp.run, scan=sp.scan[idx].copy(), index=sp.index[idx].copy(), rt=sp.rt[idx].copy(),
                   precursor_mz=np.asarray(mz_mono, np.float64).copy(), charge=np.asarray(z).astype(sp.charge.dtype),
                   mz=sp.mz[take].copy(), intensity=sp.intensity[take].copy(), offsets=off)


def a2_hypotheses(sp, valid, idx_mass32, p, single_mask):
    """Production hypotheses for spectra with single_mask False; exactly (z, iso 0) at the given precursor for spectra
    with single_mask True (their charge must be known).  Same dict layout / sort order as host.make_hypotheses."""
    from dataclasses import replace

    from vulcan.gpu.host import make_hypotheses
    single_mask = np.asarray(single_mask, bool)
    if (sp.charge[single_mask] <= 0).any():
        raise ValueError("single-hypothesis spectra need a known (certified) charge")
    h1 = make_hypotheses(sp, valid & ~single_mask, idx_mass32, p)
    h2 = make_hypotheses(sp, valid & single_mask, idx_mass32, replace(p, isotope_errors=(0, 0)))
    out = {k: np.concatenate([h1[k], h2[k]]) for k in h1}
    o = np.lexsort((out["iso"], out["z"], out["spec"]))
    return {k: v[o] for k, v in out.items()}
