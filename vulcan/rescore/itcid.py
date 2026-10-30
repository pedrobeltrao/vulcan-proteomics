"""IT-CID instrument class (tree 3c): L3 feature assembly for low-resolution ion-trap CID MS2.  A separate class so
the OT-HCD path (vulcan.rescore.features / static, and the scoring-speed work on it) is untouched.

Differences from the OT-HCD FeatureBuilder (everything else, incl. label blindness, is inherited):
  * fragment matching tolerance in Da (default +-0.5 Da, as the K1/K2 search) instead of ppm, scanning up to
    `scan` centroids per window (IT windows are 1 Th wide);
  * predictor = Prosit 2020 CID / TMT-CID via Koina (vulcan.predict.koina), fixed CE 35 (no NCE grid: the CID
    model has no CE input, and the TMT model is queried with the CID CE a priori);
  * peptides Prosit cannot predict (length > 30, non-standard residues) get pred_missing = 1, fragment features 0
    and RT features set to the run median of the predicted PSMs (not 0, which would look like a perfect RT match).
"""
from __future__ import annotations

import time

import numpy as np
import torch

from .features import PRED_FEATURES, FeatureBuilder, SpecPeaks, calibration_mask, fragment_features, \
    fragment_mz, rt_calibrate

ITCID_PRED_FEATURES = PRED_FEATURES + ["pred_missing"]
CE = 35.0


def observed_da(peaks: SpecPeaks, spec, frag_mz, tol_da=0.5, scan=12):
    """Max observed intensity within +-tol_da Th of each fragment m/z (0 = unmatched)."""
    N, F, T = frag_mz.shape
    base = (spec.to(torch.float64) * 10000.0)[:, None, None]
    k = (base + frag_mz).reshape(-1)
    lo = torch.searchsorted(peaks.key, k - tol_da, side="left")
    best = torch.zeros_like(k, dtype=torch.float32)
    nk = peaks.key.numel()
    for j in range(scan):
        i = torch.clamp(lo + j, max=nk - 1)
        ok = (lo + j < nk) & (peaks.key[i] <= k + tol_da) & (peaks.key[i] >= k - tol_da)
        best = torch.where(ok, torch.maximum(best, peaks.it[i]), best)
    return best.reshape(N, F, T)


class ItcidFeatureBuilder(FeatureBuilder):
    def __init__(self, cache_factory, static_mods, var_mods, dev="cuda", tol_da=0.5, chunk=200_000):
        super().__init__(cache_factory, static_mods, var_mods, dev=dev, tol_ppm=0.0, chunk=chunk)
        self.tol_da = tol_da

    def _frag_block(self, l1, rows, charge, spec, peaks, cache):
        dev = self.dev
        Lmax = int(l1.length[rows].max())
        fmz, valid = fragment_mz(l1, rows, dev, self.static, self.var, Lmax)
        nt = float(self.var.get("^", 0.0))
        if nt:   # peptide-N-term mod (TMT, mask bit 63): b ions carry it, y ions do not
            has = torch.from_numpy(((l1.mask[rows] >> np.uint64(63)) & np.uint64(1)).astype(bool)).to(dev)
            add = torch.where(has, nt, 0.0)[:, None] * valid
            fmz[..., 0] += add
            fmz[..., 1] += add / 2
        obs = observed_da(peaks, torch.from_numpy(spec).to(dev), fmz, self.tol_da)
        off, flat = cache.lookup_ms2(l1, rows, charge)
        N = len(rows)
        lens = np.diff(off) // 4
        pred = torch.zeros(N, Lmax - 1, 4, dtype=torch.float32, device=dev)
        rr = np.repeat(np.arange(N), lens)
        ff = np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
        vals = torch.from_numpy(flat.astype(np.float32).reshape(-1, 4)).to(dev)
        pred[torch.from_numpy(rr).to(dev), torch.from_numpy(ff).to(dev)] = vals
        z = torch.from_numpy(np.asarray(charge, np.int64)).to(dev)
        typ_ok = torch.ones(N, 1, 4, dtype=torch.bool, device=dev)
        typ_ok[:, 0, 1] = z >= 3
        typ_ok[:, 0, 3] = z >= 3
        s = torch.from_numpy(spec).to(dev)
        in_rng = (fmz >= peaks.mzmin[s][:, None, None] - 0.5) & (fmz <= peaks.mzmax[s][:, None, None] + 0.5)
        vs = (valid[..., None] & typ_ok & in_rng).reshape(N, -1).float()
        f = fragment_features(pred.reshape(N, -1), obs.reshape(N, -1), vs, peaks.tic[s])
        f["pred_missing"] = (pred.reshape(N, -1).amax(1) <= 0).float()
        return f

    def run(self, l1, psm: dict, sp, nce=None):
        t0 = time.perf_counter()
        peaks = SpecPeaks(sp, self.dev)
        cache = self.cache_factory(CE)
        n = len(psm["pf_row"])
        out = {k: np.zeros(n, np.float32) for k in ITCID_PRED_FEATURES}
        t2 = time.perf_counter()
        st0 = dict(cache.stats)
        for b in range(0, n, self.chunk):
            sl = slice(b, min(n, b + self.chunk))
            f = self._frag_block(l1, psm["pf_row"][sl], psm["charge"][sl], psm["spec_idx"][sl], peaks, cache)
            for k, v in f.items():
                out[k][sl] = v.float().cpu().numpy()
        torch.cuda.synchronize()
        t_frag = time.perf_counter() - t2
        t3 = time.perf_counter()
        rtp = cache.lookup_rt(l1, psm["pf_row"]).astype(np.float64)
        have = np.isfinite(rtp)
        cal = calibration_mask(psm, 0.10) & have
        fit, mad = rt_calibrate(psm["rt"], np.where(have, rtp, 0.0), cal)
        res = np.abs(psm["rt"] - fit)
        med = float(np.median(res[have])) if have.any() else 0.0
        res = np.where(have, res, med)
        out["rt_abs_resid"] = res.astype(np.float32)
        out["rt_resid_z"] = (res / mad).astype(np.float32)
        out["rt_pred"] = np.where(have, rtp, np.nan).astype(np.float32)
        t_rt = time.perf_counter() - t3
        st = {k: cache.stats[k] - st0.get(k, 0) for k in cache.stats}
        info = {"nce": CE, "rt_calib_n": int(cal.sum()), "rt_mad_min": float(mad), "cache_key": cache.key,
                "cache_stats": st, "frac_pred_missing": float(out["pred_missing"].mean()),
                "timing_s": {"fragment_features": t_frag, "rt": t_rt, "total": time.perf_counter() - t0,
                             "of_which_predict": st["predict_s"]}}
        return out, info
