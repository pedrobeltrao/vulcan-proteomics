"""Scorer v2 features (E7/v2; pre-registered in results/2026-10-28_E7_chimeric_v2/config.yaml).  LABEL-BLIND: the PSM
table (incl. the chimeric columns of vulcan.chimeric.search), the index side, the spectra and the MS1 feature table;
nothing from labels.arrow.

All v1 features (vulcan.rescore.features.FeatureBuilder: predicted fragment/RT features) are computed by the v1 code;
the per-run label-blind calibrations (NCE choice, RT fit, hyperscore run-z) use the PRIMARY (closed, group 0) PSMs
only, so for a closed PSM table the v1 features are exactly the v1 values.  Added:
  pred_sa_b, pred_sa_y, pred_pearson_b, pred_pearson_y   predicted vs observed per ion series (b: b+/b2+, y: y+/y2+)
  cleave_cov                                             fraction of backbone cleavages with any b/y observed
  longest_by_frac                                        longest run of consecutive covered cleavages / (len - 1)
  prec_ppm_resid, prec_ppm_resid_z                       signed precursor ppm (isotope-corrected) minus a per-run
                                                         label-blind calibration curve (binned median vs m/z on the
                                                         same calibration PSMs as the RT fit), and |.|/robust SD
  ms1_has_feature, ms1_iso_cos, ms1_trace_cos, ms1_n_iso, ms1_log_int_runz, ms1_isolation_fraction, ms1_rt_offset_rel
                                                         MS1 feature quality of the PSM's precursor group
  chim_rank, n_groups, n_features_in_window              chimeric context (label-free)
  z_inst_agree                                           PSM charge == instrument charge (0 if unknown)
"""
from __future__ import annotations

import numpy as np
import torch

from .features import FeatureBuilder, SpecPeaks, _sa, calibration_mask, fragment_features, fragment_mz, observed, \
    rt_calibrate, PRED_FEATURES

NEUTRON = 1.00335
V2_EXTRA = ["pred_sa_b", "pred_sa_y", "pred_pearson_b", "pred_pearson_y", "cleave_cov", "longest_by_frac",
            "prec_ppm_resid", "prec_ppm_resid_z", "ms1_has_feature", "ms1_iso_cos", "ms1_trace_cos", "ms1_n_iso",
            "ms1_log_int_runz", "ms1_isolation_fraction", "ms1_rt_offset_rel", "chim_rank", "n_groups",
            "n_features_in_window", "z_inst_agree"]


def _pearson(p, o, v):
    n = v.sum(1).clamp_min(1)
    pm = ((p * v).sum(1) / n)[:, None]
    om = ((o * v).sum(1) / n)[:, None]
    cov = ((p - pm) * (o - om) * v).sum(1)
    sd = torch.sqrt(((p - pm) ** 2 * v).sum(1) * ((o - om) ** 2 * v).sum(1))
    return torch.where(sd > 0, cov / sd.clamp_min(1e-30), torch.zeros_like(cov))


def _longest_run(x: torch.Tensor) -> torch.Tensor:
    """x: [N, F] bool -> longest run of True per row."""
    N, F = x.shape
    run = torch.zeros(N, dtype=torch.int32, device=x.device)
    best = torch.zeros_like(run)
    for f in range(F):
        run = torch.where(x[:, f], run + 1, torch.zeros_like(run))
        best = torch.maximum(best, run)
    return best


class FeatureBuilderV2(FeatureBuilder):
    def _frag_block(self, l1, rows, charge, spec, peaks, cache, extra=False):
        if not extra:
            return super()._frag_block(l1, rows, charge, spec, peaks, cache)
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
        vmask = valid[..., None] & typ_ok & in_rng                      # [N, F, 4]
        f = fragment_features(pred.reshape(N, -1), obs.reshape(N, -1), vmask.reshape(N, -1).float(), peaks.tic[s])
        for name, sl in (("b", slice(0, 2)), ("y", slice(2, 4))):
            vs = vmask[:, :, sl].reshape(N, -1).float()
            p = pred[:, :, sl].reshape(N, -1) * vs
            o = obs[:, :, sl].reshape(N, -1) * vs
            f[f"pred_sa_{name}"] = _sa(p, o)
            f[f"pred_pearson_{name}"] = _pearson(p, o, vs)
        cov = ((obs > 0) & vmask).any(2) & valid                        # [N, F] cleavage covered
        nclv = valid.sum(1).clamp_min(1).float()
        f["cleave_cov"] = cov.sum(1).float() / nclv
        f["longest_by_frac"] = _longest_run(cov).float() / nclv
        return f

    def run_v2(self, l1, psm: dict, sp, ms1_feat: dict | None, nce=None):
        """psm: dict of numpy columns of a (chimeric) PSM table of one run.  Returns (features dict, info)."""
        import time
        t0 = time.perf_counter()
        peaks = SpecPeaks(sp, self.dev)
        prim = psm["group"] == 0 if "group" in psm else np.ones(len(psm["rank"]), bool)
        pidx = np.flatnonzero(prim)
        psm_p = {k: v[pidx] for k, v in psm.items()}
        nce_grid = None
        if nce is None:
            nce, nce_grid = self.choose_nce(l1, psm_p, sp, peaks)
        cache = self.cache_factory(nce)
        n = len(psm["pf_row"])
        names = PRED_FEATURES + ["pred_sa_b", "pred_sa_y", "pred_pearson_b", "pred_pearson_y", "cleave_cov",
                                 "longest_by_frac"]
        out = {k: np.zeros(n, np.float32) for k in names}
        for b in range(0, n, self.chunk):
            sl = slice(b, min(n, b + self.chunk))
            f = self._frag_block(l1, psm["pf_row"][sl], psm["charge"][sl], psm["spec_idx"][sl], peaks, cache, extra=True)
            for k, v in f.items():
                out[k][sl] = v.float().cpu().numpy()
        # RT (v1 rule, calibration on primary PSMs only)
        rtp = cache.lookup_rt(l1, psm["pf_row"])
        cal_p = calibration_mask(psm_p, 0.10)
        cal = np.zeros(n, bool)
        cal[pidx[cal_p]] = True
        fit, mad = rt_calibrate(psm["rt"], rtp, cal)
        res = psm["rt"] - fit
        out["rt_abs_resid"] = np.abs(res).astype(np.float32)
        out["rt_resid_z"] = (np.abs(res) / mad).astype(np.float32)
        # precursor ppm vs a label-blind calibration curve (binned median vs m/z on the calibration PSMs)
        iso = psm["isotope_error"].astype(np.float64)
        em = psm["expmass"].astype(np.float64) - iso * NEUTRON
        cm = psm["calcmass"].astype(np.float64)
        ppm = (em - cm) / cm * 1e6
        mz = (psm["expmass"].astype(np.float64) + psm["charge"] * 1.007276466812) / psm["charge"]
        x, y = mz[cal], ppm[cal]
        o = np.argsort(x)
        x, y = x[o], y[o]
        nb = int(np.clip(len(x) // 200, 4, 20))
        edges = np.linspace(0, len(x), nb + 1).astype(int)
        fx = np.array([np.median(x[a:b]) for a, b in zip(edges[:-1], edges[1:]) if b > a])
        fy = np.array([np.median(y[a:b]) for a, b in zip(edges[:-1], edges[1:]) if b > a])
        resid = ppm - np.interp(mz, fx, fy)
        rc = resid[cal]
        sd = 1.4826 * np.median(np.abs(rc - np.median(rc))) + 1e-6
        out["prec_ppm_resid"] = resid.astype(np.float32)
        out["prec_ppm_resid_z"] = (np.abs(resid) / sd).astype(np.float32)
        # MS1 feature quality
        fid = psm["ms1_feature_id"].astype(np.int64)
        has = fid >= 0
        out["ms1_has_feature"] = has.astype(np.float32)
        out["ms1_iso_cos"] = np.where(has, psm["ms1_iso_cos"], 0).astype(np.float32)
        out["ms1_isolation_fraction"] = np.where(has, psm["ms1_isolation_fraction"], 0).astype(np.float32)
        tc = np.zeros(n, np.float32)
        ni = np.zeros(n, np.float32)
        rel = np.zeros(n, np.float32)
        li = np.zeros(n, np.float32)
        if ms1_feat is not None and has.any():
            f = fid[has]
            tc[has] = ms1_feat["trace_cos"][f]
            ni[has] = ms1_feat["n_isotopes"][f]
            width = np.maximum(ms1_feat["rt_end"][f] - ms1_feat["rt_start"][f], 1e-3)
            rel[has] = (psm["rt"][has] - ms1_feat["rt_apex"][f]) / width
            lint = np.log10(np.maximum(ms1_feat["intensity_apex"][f], 1.0))
            allint = np.log10(np.maximum(ms1_feat["intensity_apex"], 1.0))
            med = np.median(allint)
            iqr = max(np.quantile(allint, 0.75) - np.quantile(allint, 0.25), 1e-3)
            li[has] = (lint - med) / iqr
        out["ms1_trace_cos"], out["ms1_n_iso"], out["ms1_rt_offset_rel"], out["ms1_log_int_runz"] = tc, ni, rel, li
        out["chim_rank"] = np.minimum(psm["chim_rank"], 3).astype(np.float32)
        out["n_groups"] = psm["n_groups"].astype(np.float32)
        out["n_features_in_window"] = psm["n_features_in_window"].astype(np.float32)
        out["z_inst_agree"] = ((psm["inst_charge"] > 0) & (psm["inst_charge"] == psm["charge"])).astype(np.float32)
        info = {"nce": nce, "nce_grid_median_sa": nce_grid, "rt_mad_min": float(mad), "prec_ppm_sd": float(sd),
                "cache_stats": dict(cache.stats), "wall_s": time.perf_counter() - t0}
        return out, info


def read_ms1_features(path: str) -> dict:
    import pyarrow as pa
    import pyarrow.ipc as ipc
    with pa.OSFile(path, "rb") as f:
        t = ipc.open_file(f).read_all()
    return {c: t.column(c).to_numpy() for c in ("trace_cos", "n_isotopes", "rt_apex", "rt_start", "rt_end",
                                                  "intensity_apex")}


def model_psm(psm: dict) -> dict:
    """PSM dict for StaticModel.feature_matrix: hs_runz must be computed on primary rank-1 PSMs only (v1 semantics);
    feature_matrix uses rank == 1 rows, so secondary rows are moved to rank + 1000 (rank is not a model feature)."""
    if "group" not in psm:
        return psm
    q = dict(psm)
    q["rank"] = np.where(psm["group"] == 0, psm["rank"], psm["rank"] + 1000)
    return q
