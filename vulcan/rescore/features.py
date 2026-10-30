"""L3 predicted-feature assembly on the GPU (DESIGN §4 L3).  LABEL-BLIND: inputs are the PSM table (index rows,
charges, spectrum indices), the index side (residues, mod masks) and the spectra; nothing from labels.arrow.

Per PSM:
  fragment features  predicted b/y (z1, z2 if precursor z >= 3) intensities vs the observed spectrum (all centroids,
                     max intensity within +-tol ppm): spectral angle (raw and sqrt), Pearson r, entropy similarity,
                     fraction of predicted intensity observed, top-6 predicted observed, observed/TIC.
  RT features        |RT - f(pred)| after a per-run label-blind calibration f (monotone binned-median fit on the
                     top-decile hyperscore rank-1 PSMs of BOTH labels, two outlier-rejection passes), and its
                     robust z.
  NCE                chosen per run from a grid by median spectral angle of the same label-blind calibration set.
"""
from __future__ import annotations

import time

import numpy as np
import torch

PROTON = 1.007276466812
H2O = 18.010564684
MONO = {"G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414, "T": 101.047679,
        "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578,
        "K": 128.094963, "E": 129.042593, "M": 131.040485, "H": 137.058912, "F": 147.068414, "R": 156.101111,
        "Y": 163.06332, "W": 186.079313, "U": 150.95363, "O": 237.147727}
NCE_GRID = (20.0, 23.0, 26.0, 29.0, 32.0, 35.0)

PRED_FEATURES = ["pred_sa", "pred_sa_sqrt", "pred_pearson", "pred_entropy", "pred_frac_obs", "pred_top6",
                 "pred_obs_tic", "rt_abs_resid", "rt_resid_z"]


def _aa_table(static_mods, dev):
    t = torch.zeros(128, dtype=torch.float64, device=dev)
    for a, m in MONO.items():
        t[ord(a)] = m + static_mods.get(a, 0.0)
    return t


class SpecPeaks:
    """All centroids of one run on the device, sorted by (spectrum, m/z), as float64 composite keys."""

    def __init__(self, sp, dev):
        n = len(sp)
        lens = np.diff(sp.offsets)
        seg = np.repeat(np.arange(n, dtype=np.int64), lens)
        mz = sp.mz.astype(np.float64)
        it = sp.intensity.astype(np.float32)
        ok = np.isfinite(mz) & (it > 0)
        seg, mz, it = seg[ok], mz[ok], it[ok]
        key = seg * 10000.0 + mz
        o = np.argsort(key, kind="stable")
        self.key = torch.from_numpy(key[o]).to(dev)
        self.it = torch.from_numpy(it[o]).to(dev)
        tic = np.bincount(seg, weights=it, minlength=n)
        mzmin = np.full(n, np.inf)
        mzmax = np.full(n, -np.inf)
        np.minimum.at(mzmin, seg, mz)
        np.maximum.at(mzmax, seg, mz)
        self.tic = torch.from_numpy(tic).to(dev)
        self.mzmin = torch.from_numpy(mzmin).to(dev)
        self.mzmax = torch.from_numpy(mzmax).to(dev)


def fragment_mz(l1, rows, dev, static_mods, var_mods, Lmax):
    """[N, Lmax-1, 4] fragment m/z (b1+, b2+, y1+, y2+) in alphabase row layout; 0 outside the peptide."""
    aa = _aa_table(static_mods, dev)
    var = torch.zeros(128, dtype=torch.float64, device=dev)
    for a, m in var_mods.items():
        var[ord(a)] = m
    N = len(rows)
    L = torch.from_numpy(l1.length[rows].astype(np.int64)).to(dev)
    off = torch.from_numpy(l1.seq_off[rows].astype(np.int64)).to(dev)
    res_all = torch.from_numpy(l1.residues).to(dev)
    pos = torch.arange(Lmax, device=dev)
    inpep = pos[None, :] < L[:, None]
    idx = torch.where(inpep, off[:, None] + pos[None, :], torch.zeros_like(off[:, None]))
    r = res_all[idx].long()
    mask = torch.from_numpy(l1.mask[rows].view(np.int64)).to(dev)
    bit = ((mask[:, None] >> pos[None, :]) & 1).bool()
    m = torch.where(inpep, aa[r] + torch.where(bit, var[r], 0.0), 0.0)
    cum = torch.cumsum(m, 1)                                  # cum[:, f] = sum residues 0..f
    tot = cum.gather(1, (L - 1)[:, None])
    b = cum[:, :Lmax - 1]
    y = tot - b + H2O
    out = torch.stack([b + PROTON, (b + 2 * PROTON) / 2, y + PROTON, (y + 2 * PROTON) / 2], dim=2)
    valid = (pos[None, :Lmax - 1] < (L - 1)[:, None])
    return out * valid[..., None], valid


def observed(peaks: SpecPeaks, spec, frag_mz, tol_ppm=20.0):
    """Max observed intensity within tol of each fragment m/z (0 = unmatched)."""
    N, F, T = frag_mz.shape
    base = (spec.to(torch.float64) * 10000.0)[:, None, None]
    k = (base + frag_mz).reshape(-1)
    tol = (frag_mz * tol_ppm * 1e-6).reshape(-1)
    lo = torch.searchsorted(peaks.key, k - tol, side="left")
    best = torch.zeros_like(k, dtype=torch.float32)
    nk = peaks.key.numel()
    for j in range(4):
        i = torch.clamp(lo + j, max=nk - 1)
        ok = (lo + j < nk) & (peaks.key[i] <= k + tol) & (peaks.key[i] >= k - tol)
        best = torch.where(ok, torch.maximum(best, peaks.it[i]), best)
    return best.reshape(N, F, T)


def _sa(p, o):
    num = (p * o).sum(1)
    den = torch.sqrt((p * p).sum(1) * (o * o).sum(1))
    cos = torch.where(den > 0, num / den.clamp_min(1e-30), torch.zeros_like(num)).clamp(-1, 1)
    return 1 - 2 * torch.arccos(cos) / np.pi


def _entropy(x):
    x = x / x.sum(1, keepdim=True).clamp_min(1e-30)
    return -(torch.where(x > 0, x * torch.log(x.clamp_min(1e-30)), torch.zeros_like(x))).sum(1)


def fragment_features(pred, obs, valid_slot, tic):
    """pred, obs: [N, S] (S = flattened fragment slots) with invalid slots zeroed."""
    p = pred * valid_slot
    o = obs * valid_slot
    f = {}
    f["pred_sa"] = _sa(p, o)
    f["pred_sa_sqrt"] = _sa(torch.sqrt(p), torch.sqrt(o))
    n = valid_slot.sum(1).clamp_min(1)
    pm = (p.sum(1) / n)[:, None]
    om = (o.sum(1) / n)[:, None]
    cov = ((p - pm) * (o - om) * valid_slot).sum(1)
    sd = torch.sqrt(((p - pm) ** 2 * valid_slot).sum(1) * ((o - om) ** 2 * valid_slot).sum(1))
    f["pred_pearson"] = torch.where(sd > 0, cov / sd.clamp_min(1e-30), torch.zeros_like(cov))
    # entropy similarity (Li et al. 2021) on the aligned slot vectors
    pn = p / p.sum(1, keepdim=True).clamp_min(1e-30)
    on = o / o.sum(1, keepdim=True).clamp_min(1e-30)
    hab = _entropy(pn + on)
    es = 1 - (2 * hab - _entropy(pn) - _entropy(on)) / np.log(4)
    f["pred_entropy"] = torch.where(o.sum(1) > 0, es, torch.zeros_like(es))
    f["pred_frac_obs"] = (p * (o > 0)).sum(1) / p.sum(1).clamp_min(1e-30)
    k = min(6, p.shape[1])
    top = torch.topk(p, k, dim=1).indices
    hit = (o.gather(1, top) > 0) & (p.gather(1, top) > 0)
    f["pred_top6"] = hit.float().sum(1) / (p.gather(1, top) > 0).float().sum(1).clamp_min(1)
    f["pred_obs_tic"] = o.sum(1) / tic.clamp_min(1e-30)
    return f


def rt_calibrate(rt_obs, rt_pred, calib_mask, n_bins=None, passes=2):
    """Label-blind monotone binned-median fit obs ~ f(pred) on calibration PSMs; returns (f(pred) for all, mad)."""
    x = rt_pred[calib_mask].astype(np.float64)
    y = rt_obs[calib_mask].astype(np.float64)
    keep = np.ones(len(x), bool)
    fx = fy = None
    mad = 1.0
    for _ in range(passes + 1):
        xs, ys = x[keep], y[keep]
        o = np.argsort(xs)
        xs, ys = xs[o], ys[o]
        nb = n_bins or int(np.clip(len(xs) // 100, 8, 60))
        edges = np.linspace(0, len(xs), nb + 1).astype(int)
        fx = np.array([np.median(xs[a:b]) for a, b in zip(edges[:-1], edges[1:]) if b > a])
        fy = np.array([np.median(ys[a:b]) for a, b in zip(edges[:-1], edges[1:]) if b > a])
        fy = np.maximum.accumulate(fy)                          # monotone
        res = y - _interp(x, fx, fy)
        mad = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep]))) + 1e-6
        keep = np.abs(res) <= 4 * mad
    return _interp(rt_pred.astype(np.float64), fx, fy), mad


def _interp(x, fx, fy):
    y = np.interp(x, fx, fy)
    if len(fx) >= 2:
        lo = (fy[1] - fy[0]) / max(fx[1] - fx[0], 1e-9)
        hi = (fy[-1] - fy[-2]) / max(fx[-1] - fx[-2], 1e-9)
        y = np.where(x < fx[0], fy[0] + (x - fx[0]) * lo, y)
        y = np.where(x > fx[-1], fy[-1] + (x - fx[-1]) * hi, y)
    return y


def calibration_mask(psm, top_frac=0.10, max_n=None):
    """Label-blind calibration PSMs: rank-1 PSMs in the top `top_frac` of hyperscore (both labels)."""
    r1 = psm["rank"] == 1
    hs = np.where(r1, psm["hyperscore"], -np.inf)
    thr = np.quantile(psm["hyperscore"][r1], 1 - top_frac)
    m = r1 & (hs >= thr)
    if max_n is not None and m.sum() > max_n:
        idx = np.flatnonzero(m)
        idx = idx[np.argsort(-psm["hyperscore"][idx], kind="stable")[:max_n]]
        m = np.zeros_like(m)
        m[idx] = True
    return m


class FeatureBuilder:
    def __init__(self, cache_factory, static_mods, var_mods, dev="cuda", tol_ppm=20.0, chunk=200_000):
        """cache_factory(nce) -> vulcan.predict.cache.PredCache."""
        self.cache_factory = cache_factory
        self.static, self.var = static_mods, var_mods
        self.dev = torch.device(dev)
        self.tol = tol_ppm
        self.chunk = chunk
        self.timing = {}

    def _frag_block(self, l1, rows, charge, spec, peaks, cache):
        dev = self.dev
        Lmax = int(l1.length[rows].max())
        fmz, valid = fragment_mz(l1, rows, dev, self.static, self.var, Lmax)
        obs = observed(peaks, torch.from_numpy(spec).to(dev), fmz, self.tol)
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
        return fragment_features(pred.reshape(N, -1), obs.reshape(N, -1), vs, peaks.tic[s])

    def choose_nce(self, l1, psm, sp, peaks, n_cal=2000):
        cal = np.flatnonzero(calibration_mask(psm, 0.10, n_cal))
        rows, z, spec = psm["pf_row"][cal], psm["charge"][cal], psm["spec_idx"][cal]
        med = {}
        for nce in NCE_GRID:
            f = self._frag_block(l1, rows, z, spec, peaks, self.cache_factory(nce))
            med[nce] = float(f["pred_sa"].median())
        return max(med, key=med.get), med

    def run(self, l1, psm: dict, sp, nce=None):
        """psm: dict of numpy columns (one run).  Returns (features dict of numpy arrays, info)."""
        t0 = time.perf_counter()
        peaks = SpecPeaks(sp, self.dev)
        t_peaks = time.perf_counter() - t0
        t1 = time.perf_counter()
        nce_grid = None
        if nce is None:
            nce, nce_grid = self.choose_nce(l1, psm, sp, peaks)
        t_nce = time.perf_counter() - t1
        cache = self.cache_factory(nce)
        n = len(psm["pf_row"])
        out = {k: np.zeros(n, np.float32) for k in PRED_FEATURES}
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
        rtp = cache.lookup_rt(l1, psm["pf_row"])
        cal = calibration_mask(psm, 0.10)
        fit, mad = rt_calibrate(psm["rt"], rtp, cal)
        res = psm["rt"] - fit
        out["rt_abs_resid"] = np.abs(res).astype(np.float32)
        out["rt_resid_z"] = (np.abs(res) / mad).astype(np.float32)
        out["rt_pred"] = rtp.astype(np.float32)
        t_rt = time.perf_counter() - t3
        st = {k: cache.stats[k] - st0.get(k, 0) for k in cache.stats}
        info = {"nce": nce, "nce_grid_median_sa": nce_grid, "rt_calib_n": int(cal.sum()), "rt_mad_min": float(mad),
                "cache_key": cache.key, "cache_stats": st,
                "timing_s": {"peaks": t_peaks, "nce_select": t_nce, "fragment_features": t_frag, "rt": t_rt,
                             "total": time.perf_counter() - t0, "of_which_predict": st["predict_s"]}}
        return out, info
