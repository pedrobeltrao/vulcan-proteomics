"""Fast production scoring for the FROZEN static_gbm_pred (DESIGN §6.0 A1.1): same features, same model, same
competition — only faster.  results/2026-10-27_scoring_speed has the equality record.

What is the same (bit-identical by construction, checked by tests/rescore/test_fast_equal.py and the D0/D1/D8 record):
* L3 features: the same `vulcan.rescore.features` functions (fragment_mz, observed, fragment_features,
  calibration_mask, rt_calibrate) on the same tensors with the same shapes (chunk 200k rows; per-NCE calibration
  blocks of the same size), the same label-blind per-run NCE choice over the same grid.  Changes:
  - `GpuPeaks`: the SpecPeaks arrays built on the GPU (float64 key, stable sort; TIC as the sequential float64 sum
    np.bincount computes, mz min/max exact) instead of a host argsort;
  - the NCE-independent half of a fragment block (fragment m/z, observed intensities, slot mask) is computed once
    per run instead of once per grid NCE (deterministic kernels on identical inputs);
  - predictions come from `vulcan.predict.library.LibCache` (memory-mapped consolidated library; identical values
    for hits) and are scattered into the padded layout on the GPU.
* GBM: `GpuHGB` evaluates the sklearn HistGradientBoosting trees on the GPU in float64 with sklearn's exact
  decision rule (NaN -> missing_go_to_left, else x <= threshold) and sums leaf values sequentially in tree order
  from the baseline (sklearn's `_raw_predict` order), so decision_function is reproduced bit for bit.
"""
from __future__ import annotations

import time

import numpy as np
import torch

from vulcan.rescore import features as F
from vulcan.rescore.static import feature_matrix

try:
    import numba as _nb

    @_nb.njit(cache=True, nogil=True)
    def _seg_stats_nb(mz, it, off, tic, mzmin, mzmax):        # pragma: no cover - jitted
        for s in range(off.shape[0] - 1):
            acc = 0.0
            lo = np.inf
            hi = -np.inf
            for k in range(off[s], off[s + 1]):
                m = np.float64(mz[k])
                w = it[k]
                if np.isfinite(m) and w > 0:
                    acc += np.float64(w)
                    if m < lo:
                        lo = m
                    if m > hi:
                        hi = m
            tic[s] = acc
            mzmin[s] = lo
            mzmax[s] = hi
except Exception:  # noqa: BLE001
    _seg_stats_nb = None


class GpuPeaks:
    """Same arrays as vulcan.rescore.features.SpecPeaks (key, it, tic, mzmin, mzmax), built on the device."""

    def __init__(self, sp, dev):
        n = len(sp)
        off = np.asarray(sp.offsets, np.int64)
        mz32 = np.ascontiguousarray(sp.mz)
        it32 = np.ascontiguousarray(sp.intensity, dtype=np.float32)
        lens = torch.from_numpy(np.diff(off)).to(dev)
        mz = torch.from_numpy(mz32).to(dev, non_blocking=False).double()
        it = torch.from_numpy(it32).to(dev)
        seg = torch.repeat_interleave(torch.arange(n, device=dev, dtype=torch.int64), lens)
        ok = torch.isfinite(mz) & (it > 0)
        seg, mz, it = seg[ok], mz[ok], it[ok]
        key = seg.double() * 10000.0 + mz
        key_s, o = torch.sort(key, stable=True)
        self.key = key_s
        self.it = it[o]
        tic = np.zeros(n, np.float64)
        mzmin = np.full(n, np.inf)
        mzmax = np.full(n, -np.inf)
        if _seg_stats_nb is not None and mz32.dtype == np.float32:
            _seg_stats_nb(mz32, it32, off, tic, mzmin, mzmax)
        else:                                            # reference path (exactly SpecPeaks)
            segc = np.repeat(np.arange(n, dtype=np.int64), np.diff(off))
            m64 = mz32.astype(np.float64)
            okc = np.isfinite(m64) & (it32 > 0)
            tic = np.bincount(segc[okc], weights=it32[okc], minlength=n)
            np.minimum.at(mzmin, segc[okc], m64[okc])
            np.maximum.at(mzmax, segc[okc], m64[okc])
        self.tic = torch.from_numpy(tic).to(dev)
        self.mzmin = torch.from_numpy(mzmin).to(dev)
        self.mzmax = torch.from_numpy(mzmax).to(dev)


class _ObsBlock:
    """NCE-independent part of FeatureBuilder._frag_block for one set of PSM rows."""

    def __init__(self, fb, l1, rows, charge, spec, peaks):
        dev = fb.dev
        self.rows, self.charge = rows, charge
        self.N = len(rows)
        self.Lmax = int(l1.length[rows].max())
        fmz, valid = F.fragment_mz(l1, rows, dev, fb.static, fb.var, self.Lmax)
        self.obs = F.observed(peaks, torch.from_numpy(spec).to(dev), fmz, fb.tol)
        z = torch.from_numpy(np.asarray(charge, np.int64)).to(dev)
        typ_ok = torch.ones(self.N, 1, 4, dtype=torch.bool, device=dev)
        typ_ok[:, 0, 1] = z >= 3
        typ_ok[:, 0, 3] = z >= 3
        s = torch.from_numpy(spec).to(dev)
        in_rng = (fmz >= peaks.mzmin[s][:, None, None] - 0.5) & (fmz <= peaks.mzmax[s][:, None, None] + 0.5)
        self.vs = (valid[..., None] & typ_ok & in_rng).reshape(self.N, -1).float()
        self.tic = peaks.tic[s]
        # padded-layout scatter indices of the flat prediction rows (row r, fragment f)
        lens = (l1.length[rows].astype(np.int64) - 1)
        lt = torch.from_numpy(lens).to(dev)
        self.rr = torch.repeat_interleave(torch.arange(self.N, device=dev), lt)
        start = torch.cumsum(lt, 0) - lt
        self.ff = torch.arange(int(lens.sum()), device=dev) - torch.repeat_interleave(start, lt)

    def features(self, fb, l1, cache):
        dev = fb.dev
        if getattr(cache, "device", None) is not None:          # library mirrored on the GPU
            off, flat_t = cache.lookup_ms2_dev(l1, self.rows, self.charge)
            vals = flat_t.float().reshape(-1, 4)
        else:
            off, flat = cache.lookup_ms2(l1, self.rows, self.charge)
            vals = torch.from_numpy(flat).to(dev).float().reshape(-1, 4)
        pred = torch.zeros(self.N, self.Lmax - 1, 4, dtype=torch.float32, device=dev)
        pred[self.rr, self.ff] = vals
        return F.fragment_features(pred.reshape(self.N, -1), self.obs.reshape(self.N, -1), self.vs, self.tic)


class FastFeatureBuilder(F.FeatureBuilder):
    """Drop-in for vulcan.rescore.features.FeatureBuilder (same outputs; see module doc)."""

    def choose_nce(self, l1, psm, sp, peaks, n_cal=2000):
        cal = np.flatnonzero(F.calibration_mask(psm, 0.10, n_cal))
        rows, z, spec = psm["pf_row"][cal], psm["charge"][cal], psm["spec_idx"][cal]
        blk = _ObsBlock(self, l1, rows, z, spec, peaks)
        med = {}
        for nce in F.NCE_GRID:
            f = blk.features(self, l1, self.cache_factory(nce))
            med[nce] = float(f["pred_sa"].median())
        return max(med, key=med.get), med

    def run(self, l1, psm: dict, sp, nce=None, primary=None):
        """primary (E7b, scorer v2 on chimeric PSM tables): optional bool mask of the primary (group 0) PSMs; the
        per-run label-blind calibrations (NCE choice, RT fit) then use those rows only, exactly as
        vulcan.rescore.features_v2.FeatureBuilderV2.run_v2 computes the v1 feature columns.  None: v1 semantics."""
        t0 = time.perf_counter()
        peaks = GpuPeaks(sp, self.dev)
        t_peaks = time.perf_counter() - t0
        t1 = time.perf_counter()
        nce_grid = None
        if primary is not None:
            pidx = np.flatnonzero(primary)
            psm_cal = {k: psm[k][pidx] for k in ("rank", "hyperscore", "pf_row", "charge", "spec_idx")}
        else:
            psm_cal = psm
        if nce is None:
            nce, nce_grid = self.choose_nce(l1, psm_cal, sp, peaks)
        t_nce = time.perf_counter() - t1
        cache = self.cache_factory(nce)
        n = len(psm["pf_row"])
        t2 = time.perf_counter()
        st0 = dict(cache.stats)
        parts = {}
        for b in range(0, n, self.chunk):
            sl = slice(b, min(n, b + self.chunk))
            blk = _ObsBlock(self, l1, psm["pf_row"][sl], psm["charge"][sl], psm["spec_idx"][sl], peaks)
            f = blk.features(self, l1, cache)
            for k, v in f.items():
                parts.setdefault(k, []).append(v.float())
        out = {k: np.zeros(n, np.float32) for k in F.PRED_FEATURES}      # same keys / order as FeatureBuilder.run
        if n:
            keys = list(parts)
            allf = torch.stack([torch.cat(parts[k]) for k in keys]).cpu().numpy()     # one D2H copy
            for i, k in enumerate(keys):
                out[k] = allf[i]
        t_frag = time.perf_counter() - t2
        t3 = time.perf_counter()
        if getattr(cache, "device", None) is not None:
            rtp = cache.lookup_rt_dev(l1, psm["pf_row"]).cpu().numpy()
        else:
            rtp = cache.lookup_rt(l1, psm["pf_row"])
        if primary is not None:
            cal = np.zeros(n, bool)
            cal[pidx[F.calibration_mask(psm_cal, 0.10)]] = True
        else:
            cal = F.calibration_mask(psm, 0.10)
        fit, mad = F.rt_calibrate(psm["rt"], rtp, cal)
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


class GpuHGB:
    """sklearn HistGradientBoosting{Classifier,Regressor} decision_function on the GPU, bit-identical (float64)."""

    def __init__(self, est, dev="cuda"):
        dev = torch.device(dev)
        preds = est._predictors
        if any(len(it) != 1 for it in preds):
            raise ValueError("only one tree per iteration (binary / regression) is supported")
        trees = [it[0].nodes for it in preds]
        if any(t["is_categorical"].any() for t in trees):
            raise ValueError("categorical splits are not supported")
        T, M = len(trees), max(len(t) for t in trees)
        fidx = np.zeros((T, M), np.int64)
        thr = np.zeros((T, M), np.float64)
        mgl = np.zeros((T, M), bool)
        left = np.zeros((T, M), np.int64)
        right = np.zeros((T, M), np.int64)
        leaf = np.ones((T, M), bool)
        val = np.zeros((T, M), np.float64)
        for i, t in enumerate(trees):
            k = len(t)
            fidx[i, :k], thr[i, :k], mgl[i, :k] = t["feature_idx"], t["num_threshold"], t["missing_go_to_left"] == 1
            left[i, :k], right[i, :k], leaf[i, :k], val[i, :k] = t["left"], t["right"], t["is_leaf"] == 1, t["value"]
        self.depth = int(max(t["depth"].max() for t in trees))
        g = lambda a: torch.from_numpy(a).to(dev)  # noqa: E731
        self.fidx, self.thr, self.mgl, self.left, self.right, self.leaf, self.val = map(
            g, (fidx, thr, mgl, left, right, leaf, val))
        self.T, self.M = T, M
        self.baseline = float(np.asarray(est._baseline_prediction).ravel()[0])
        self.n_features = int(est.n_features_in_)
        self.dev = dev

    def decision_function(self, X, chunk=100_000) -> np.ndarray:
        X = np.ascontiguousarray(X, dtype=np.float64)
        return np.concatenate([self._df(X[b:b + chunk]) for b in range(0, max(len(X), 1), chunk)]) \
            if len(X) else np.zeros(0)

    @torch.inference_mode()
    def _df(self, X) -> np.ndarray:
        X = torch.from_numpy(X).to(self.dev)
        N = X.shape[0]
        if X.shape[1] != self.n_features:
            raise ValueError("feature count mismatch")
        tix = torch.arange(self.T, device=self.dev)[None, :]
        node = torch.zeros(N, self.T, dtype=torch.int64, device=self.dev)
        for _ in range(self.depth + 1):
            f = self.fidx[tix, node]
            x = torch.gather(X, 1, f)
            go_left = torch.where(torch.isnan(x), self.mgl[tix, node], x <= self.thr[tix, node])
            nxt = torch.where(go_left, self.left[tix, node], self.right[tix, node])
            node = torch.where(self.leaf[tix, node], node, nxt)
        v = self.val[tix, node]                           # [N, T] float64 leaf values
        raw = torch.zeros(N, dtype=torch.float64, device=self.dev)
        raw += self.baseline                              # sklearn: zeros + baseline, then += tree by tree
        for t in range(self.T):
            raw += v[:, t]
        return raw.cpu().numpy()


class FastStatic:
    """Frozen StaticModel wrapper: same feature_matrix, model evaluated by GpuHGB (gbm) or on the host (linear)."""

    def __init__(self, model, dev="cuda"):
        self.m = model
        self.names = model.names
        self.gpu = GpuHGB(model.model, dev) if model.kind == "gbm" else None

    def apply(self, psm, feat):
        X = feature_matrix(psm, feat, self.names)
        return self.gpu.decision_function(X) if self.gpu is not None else self.m.score(X)
