"""Precursor-hypothesis table (DESIGN §4 L2): one row per (spectrum, charge, isotope error, precursor source).

Schema and semantics: vulcan/ms1/SCHEMA.md §3 (HYPOTHESIS_SCHEMA below is the source of truth; tests pin it).

For every MS2 the isolation window [target - lower, target + upper] is sliced from its master MS1 scan (the scan
the instrument picked the precursor from; the preceding MS1 if the store has none). Each in-window centroid that
belongs to an MS1 feature (via hill -> feature) adds its intensity to that feature. Features are ranked label-free by
in-window intensity (rank 1 = dominant precursor); `isolation_fraction` = feature in-window / total in-window.

Precursor groups (spectrum, source, charge, feature) come from three sources, in priority order:
  1. ms1_feature  every co-isolated feature with charge in `hyp_charges`, rank <= max_rank, fraction >= min_fraction
                  (a charge-1 feature only if the instrument also says 1). `instrument_agrees` marks the one feature
                  whose charge equals the instrument charge and whose mono equals the instrument m/z up to an isotope
                  offset (`instrument_isotope_offset`; the feature mono wins: it corrects instrument mono picks).
  2. instrument   instrument charge known and no agreeing feature: the instrument m/z and charge.
  3. fallback     no instrument charge and no feature group: z in `fallback_charges` at the instrument m/z.
Each group is expanded over `isotope_errors` (default -1..+3). Isotope error e means the observed m/z is the M+e
peak, so the searched neutral mass is (mz_cal - PROTON) * z - e * NEUTRON (Sage / vulcan.gpu convention).

Nothing here reads decoy, pair or entrapment labels (DESIGN §3 label blindness): inputs are the store and the
feature table only. Calibration (vulcan.calib) enters as a ppm offset applied to the m/z; the raw m/z is kept.
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass

import numba as nb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from ..io import store
from .averagine import NEUTRON, PROTON
from .features import FeatureResult

HYPOTHESIS_VERSION = "vulcan.ms1.hypotheses 1.0"
SOURCE = {"ms1_feature": 0, "instrument": 1, "fallback": 2}      # also the priority order
SOURCE_NAME = {v: k for k, v in SOURCE.items()}
ISOTOPE_ERRORS = (-1, 0, 1, 2, 3)

# isotope-error priors per source before a pilot pass supplies empirical ones (vulcan.calib.pilot)
DEFAULT_PRIORS = {
    "ms1_feature": {-1: 0.02, 0: 0.90, 1: 0.06, 2: 0.02, 3: 0.0},
    "instrument": {-1: 0.02, 0: 0.80, 1: 0.12, 2: 0.05, 3: 0.01},
    "fallback": {-1: 0.05, 0: 0.60, 1: 0.20, 2: 0.10, 3: 0.05},
}

HYPOTHESIS_SCHEMA = pa.schema([
    pa.field("spectrum_id", pa.uint64(), nullable=False),
    pa.field("hypothesis_id", pa.uint16(), nullable=False),       # 0..n-1 within the spectrum, priority order
    pa.field("precursor_group", pa.uint16(), nullable=False),     # (source, feature, charge) group within the spectrum
    pa.field("source", pa.uint8(), nullable=False),               # SOURCE code
    pa.field("charge", pa.int8(), nullable=False),
    pa.field("isotope_error", pa.int8(), nullable=False),
    pa.field("precursor_mz", pa.float64(), nullable=False),       # observed mono m/z (uncalibrated)
    pa.field("precursor_mz_cal", pa.float64(), nullable=False),   # after the run's precursor calibration (== raw if none)
    pa.field("neutral_mass", pa.float64(), nullable=False),       # (mz_cal - PROTON) * z - isotope_error * NEUTRON
    pa.field("isotope_error_prior", pa.float32(), nullable=False),
    pa.field("feature_id", pa.uint32()),                          # MS1 feature (FEATURE_SCHEMA.feature_id), null if none
    pa.field("feature_intensity", pa.float32()),                  # feature intensity_apex (sum of isotope apexes)
    pa.field("feature_window_intensity", pa.float32()),           # feature intensity inside the window at the master scan
    pa.field("isolation_fraction", pa.float32()),                 # feature_window_intensity / window total
    pa.field("chimeric_rank", pa.uint16()),                       # label-free rank of the feature in the window (1 = top)
    pa.field("feature_iso_cos", pa.float32()),
    pa.field("feature_rt_apex", pa.float32()),
    pa.field("instrument_agrees", pa.bool_(), nullable=False),
    pa.field("instrument_isotope_offset", pa.int8()),             # (instrument mono - feature mono) in isotope steps
    pa.field("n_features_in_window", pa.uint16(), nullable=False),  # all features (any charge) with peaks in the window
])
SCHEMA_FIELDS = [f.name for f in HYPOTHESIS_SCHEMA]


@dataclass
class HypothesisParams:
    hyp_charges: tuple = (2, 3, 4, 5, 6)   # feature charges that become hypotheses
    max_rank: int = 5                      # keep co-isolated features ranked <= max_rank ...
    min_fraction: float = 0.01             # ... and with isolation_fraction >= min_fraction
    agree_ppm: float = 20.0                # instrument vs feature mono agreement (after the isotope offset)
    max_agree_offset: int = 3
    fallback_charges: tuple = (2, 3, 4)
    isotope_errors: tuple = ISOTOPE_ERRORS
    default_half_window: float = 1.0       # Th, if the store has no isolation offsets
    master_fallback: bool = True           # if master_scan is missing, use the preceding MS1


# ---------------------------------------------------------------------------------------------- window slicing
@nb.njit(cache=True)
def _window_bounds(off, mz, s_ord, lo, hi, i0, i1):
    for m in range(s_ord.shape[0]):
        s = s_ord[m]
        if s < 0:
            i0[m] = 0
            i1[m] = 0
            continue
        a, b = off[s], off[s + 1]
        i0[m] = a + np.searchsorted(mz[a:b], lo[m])
        i1[m] = a + np.searchsorted(mz[a:b], hi[m], side="right")


@nb.njit(cache=True)
def _window_fill(i0, i1, inten, pfeat, out_start, out_spec, out_feat, out_int, total, nfeat):
    """Per MS2 m: total in-window intensity and one (m, feature, intensity) row per distinct feature."""
    for m in range(i0.shape[0]):
        k0 = out_start[m]
        k = k0
        tot = 0.0
        for i in range(i0[m], i1[m]):
            tot += inten[i]
            f = pfeat[i]
            if f < 0:
                continue
            found = False
            for r in range(k0, k):          # windows hold tens of peaks: linear merge is fine
                if out_feat[r] == f:
                    out_int[r] += inten[i]
                    found = True
                    break
            if not found:
                out_spec[k] = m
                out_feat[k] = f
                out_int[k] = inten[i]
                k += 1
        total[m] = tot
        nfeat[m] = k - k0


def window_features(ms1, peak_feature: np.ndarray, s_ord: np.ndarray, lo: np.ndarray, hi: np.ndarray):
    """(spec, feature, in-window intensity) rows, total window intensity and #features per MS2. Two-pass, exact
    allocation (the number of in-window peaks bounds the number of distinct features)."""
    n2 = len(s_ord)
    mz = np.ascontiguousarray(ms1.mz, dtype=np.float64)
    i0 = np.empty(n2, np.int64)
    i1 = np.empty(n2, np.int64)
    _window_bounds(ms1.off, mz, s_ord.astype(np.int64), lo.astype(np.float64), hi.astype(np.float64), i0, i1)
    cap = np.maximum(i1 - i0, 0)
    start = np.zeros(n2 + 1, np.int64)
    np.cumsum(cap, out=start[1:])
    o_spec = np.empty(start[-1], np.int64)
    o_feat = np.full(start[-1], -1, np.int64)
    o_int = np.zeros(start[-1], np.float64)
    total = np.zeros(n2)
    nfeat = np.zeros(n2, np.int64)
    _window_fill(i0, i1, ms1.intensity.astype(np.float64), peak_feature.astype(np.int64), start, o_spec, o_feat,
                 o_int, total, nfeat)
    keep = o_feat >= 0
    return o_spec[keep], o_feat[keep], o_int[keep], total, nfeat


# ---------------------------------------------------------------------------------------------- MS2 info
def ms2_info(path_or_table, ms1, p: HypothesisParams | None = None) -> dict:
    """Precursor fields of every MS2 in acquisition order, plus its master MS1 ordinal."""
    p = p or HypothesisParams()
    cols = ["spectrum_id", "index", "ms_level", "rt_min", "precursor_mz", "isolation_target_mz", "mono_mz",
            "isolation_lower_offset", "isolation_upper_offset", "precursor_charge", "master_scan"]
    t = store.read_table(path_or_table, columns=cols) if isinstance(path_or_table, str) else path_or_table.select(cols)
    t = t.filter(pc.equal(t.column("ms_level"), 2))
    g = lambda c: np.asarray(store.column_numpy(t, c), dtype=np.float64)
    idx = g("index").astype(np.int64)
    master = g("master_scan")
    master = np.where(np.isnan(master), -1, master).astype(np.int64)
    pos = np.searchsorted(ms1.ms1_scan, master)
    posc = np.clip(pos, 0, max(len(ms1.ms1_scan) - 1, 0))
    ok = (master >= 0) & (pos < len(ms1.ms1_scan)) & (ms1.ms1_scan[posc] == master) if len(ms1.ms1_scan) else \
        np.zeros(len(master), bool)
    s_ord = np.where(ok, pos, -1)
    prev = np.searchsorted(ms1.ms1_index, idx) - 1      # last MS1 before this MS2 in acquisition order
    if p.master_fallback:
        s_ord = np.where(s_ord < 0, prev, s_ord)
    target = g("isolation_target_mz")
    pmz = g("precursor_mz")
    target = np.where(np.isnan(target), pmz, target)
    pmz = np.where(np.isnan(pmz), target, pmz)
    lo_off = g("isolation_lower_offset")
    hi_off = g("isolation_upper_offset")
    lo_off = np.where(np.isnan(lo_off) | (lo_off <= 0), p.default_half_window, lo_off)
    hi_off = np.where(np.isnan(hi_off) | (hi_off <= 0), p.default_half_window, hi_off)
    z = g("precursor_charge")
    z = np.where(np.isnan(z), 0, z).astype(np.int64)
    return dict(spectrum_id=np.asarray(t.column("spectrum_id").to_numpy(), dtype=np.uint64), index=idx,
                rt=g("rt_min"), s_ord=s_ord, master_found=ok, target=target, lo=target - lo_off, hi=target + hi_off,
                pmz=pmz, mono=g("mono_mz"), z=z)


def _apply_ppm(mz: np.ndarray, rt: np.ndarray, calib) -> np.ndarray:
    if calib is None:
        return mz.copy()
    return calib.correct_precursor(mz, rt)


def _priors_array(priors: dict, isotope_errors) -> np.ndarray:
    out = np.zeros((3, len(isotope_errors)), np.float32)
    for name, code in SOURCE.items():
        d = priors[name]
        out[code] = [float(d.get(e, d.get(str(e), 0.0))) for e in isotope_errors]
    return out


# ---------------------------------------------------------------------------------------------- builder
def build_hypotheses(path_or_table, fr: FeatureResult, params: HypothesisParams | None = None,
                     priors: dict | None = None, calib=None):
    """Return (hypothesis table, per-MS2 stats dict). `calib`: vulcan.calib.calibration.Calibration or None."""
    p = params or HypothesisParams()
    pri = {k: dict(v) for k, v in DEFAULT_PRIORS.items()}
    for k, v in (priors or {}).items():
        pri[k] = {int(e): float(x) for e, x in v.items()}
    ms1 = fr.ms1
    info = ms2_info(path_or_table, ms1, p)
    n2 = len(info["s_ord"])
    w_spec, w_feat, w_int, total, nfeat = window_features(ms1, fr.peak_feature, info["s_ord"], info["lo"], info["hi"])
    F = fr.features
    f_z = F.column("charge").to_numpy().astype(np.int64)
    f_mz = F.column("mono_mz").to_numpy()
    f_int = F.column("intensity_apex").to_numpy()
    f_cos = F.column("iso_cos").to_numpy()
    f_rt = F.column("rt_apex").to_numpy()

    # label-free rank by in-window intensity (all features, any charge); tie-break on feature id (data-derived)
    o = np.lexsort((w_feat, -w_int, w_spec))
    w_spec, w_feat, w_int = w_spec[o], w_feat[o], w_int[o]
    k = len(w_spec)
    first = np.r_[0, np.flatnonzero(np.diff(w_spec)) + 1] if k else np.empty(0, np.int64)
    rank = (np.arange(k) - np.repeat(first, np.diff(np.r_[first, k]))).astype(np.int64) + 1 if k else \
        np.empty(0, np.int64)
    frac = w_int / np.maximum(total[w_spec], 1e-30) if k else np.empty(0)
    zi_w = info["z"][w_spec]
    keep = np.isin(f_z[w_feat], p.hyp_charges) | ((f_z[w_feat] == 1) & (zi_w == 1))
    keep &= (rank <= p.max_rank) & (frac >= p.min_fraction)
    fs, ff, fi, frc, frk = w_spec[keep], w_feat[keep], w_int[keep], frac[keep], rank[keep]

    # instrument agreement: same charge, mono equal up to an isotope offset
    zi, mi, dz = info["z"][fs], info["pmz"][fs], f_z[ff]
    koff = np.rint((mi - f_mz[ff]) * dz / NEUTRON).astype(np.int64)
    resid = np.abs(mi - (f_mz[ff] + koff * NEUTRON / np.maximum(dz, 1))) / mi * 1e6
    cand = (zi > 0) & (zi == dz) & (np.abs(koff) <= p.max_agree_offset) & (resid <= p.agree_ppm)
    agree = np.zeros(len(fs), bool)
    if cand.any():
        a = np.flatnonzero(cand)
        a = a[np.lexsort((frk[a], np.abs(koff[a]), fs[a]))]
        fst = np.r_[True, fs[a][1:] != fs[a][:-1]]
        agree[a[fst]] = True
    has_agree = np.zeros(n2, bool)
    has_agree[fs[agree]] = True
    has_feat = np.zeros(n2, bool)
    has_feat[fs] = True
    inst = np.flatnonzero((info["z"] > 0) & ~has_agree)
    fb = np.flatnonzero((info["z"] <= 0) & ~has_feat)
    fbz = np.asarray(p.fallback_charges, np.int64)
    fb_s, fb_z = np.repeat(fb, len(fbz)), np.tile(fbz, len(fb))

    # precursor groups
    n_f, n_i, n_b = len(fs), len(inst), len(fb_s)
    nan = lambda n: np.full(n, np.nan)
    G = dict(
        spec=np.concatenate([fs, inst, fb_s]),
        src=np.concatenate([np.zeros(n_f), np.ones(n_i), np.full(n_b, 2)]).astype(np.int64),
        z=np.concatenate([dz, info["z"][inst], fb_z]),
        mz=np.concatenate([f_mz[ff], info["pmz"][inst], info["pmz"][fb_s]]),
        fid=np.concatenate([ff, np.full(n_i + n_b, -1)]),
        fint=np.concatenate([f_int[ff], nan(n_i + n_b)]),
        fwin=np.concatenate([fi, nan(n_i + n_b)]),
        frac=np.concatenate([frc, nan(n_i + n_b)]),
        rank=np.concatenate([frk, np.full(n_i + n_b, -1)]),
        cos=np.concatenate([f_cos[ff], nan(n_i + n_b)]),
        frt=np.concatenate([f_rt[ff], nan(n_i + n_b)]),
        agree=np.concatenate([agree, np.zeros(n_i + n_b, bool)]),
        ioff=np.concatenate([np.where(agree, koff, -128), np.full(n_i + n_b, -128)]),
    )
    sub = np.where(G["src"] == 0, G["rank"], G["z"])
    og = np.lexsort((sub, G["src"], G["spec"]))
    G = {k_: v[og] for k_, v in G.items()}
    ng = len(og)
    gfirst = np.r_[0, np.flatnonzero(np.diff(G["spec"])) + 1] if ng else np.empty(0, np.int64)
    gid = (np.arange(ng) - np.repeat(gfirst, np.diff(np.r_[gfirst, ng]))) if ng else np.empty(0, np.int64)
    rt_g = info["rt"][G["spec"]]
    mz_cal = _apply_ppm(G["mz"], rt_g, calib)

    # expand isotope errors (group-major, so rows stay sorted by spectrum, group, isotope error)
    ie = np.asarray(p.isotope_errors, np.int64)
    ne = len(ie)
    R = np.repeat(np.arange(ng), ne)
    e = np.tile(ie, ng)
    n_rows = len(R)
    rspec = G["spec"][R]
    rfirst = np.r_[0, np.flatnonzero(np.diff(rspec)) + 1] if n_rows else np.empty(0, np.int64)
    hid = (np.arange(n_rows) - np.repeat(rfirst, np.diff(np.r_[rfirst, n_rows]))) if n_rows else np.empty(0, np.int64)
    if n_rows and hid.max() > 65535:
        raise ValueError("more than 65535 hypotheses for one spectrum")
    prior = _priors_array(pri, p.isotope_errors)[G["src"][R], np.tile(np.arange(ne), ng)]
    z_r = G["z"][R]
    mzc_r = mz_cal[R]
    mask_nan = lambda x: pa.array(x.astype(np.float32), pa.float32(), mask=np.isnan(x))
    fid_r = G["fid"][R]
    rk_r = G["rank"][R]
    io_r = G["ioff"][R]
    tbl = pa.table({
        "spectrum_id": pa.array(info["spectrum_id"][rspec], pa.uint64()),
        "hypothesis_id": pa.array(hid.astype(np.uint16), pa.uint16()),
        "precursor_group": pa.array(gid[R].astype(np.uint16), pa.uint16()),
        "source": pa.array(G["src"][R].astype(np.uint8), pa.uint8()),
        "charge": pa.array(z_r.astype(np.int8), pa.int8()),
        "isotope_error": pa.array(e.astype(np.int8), pa.int8()),
        "precursor_mz": pa.array(G["mz"][R], pa.float64()),
        "precursor_mz_cal": pa.array(mzc_r, pa.float64()),
        "neutral_mass": pa.array((mzc_r - PROTON) * z_r - e * NEUTRON, pa.float64()),
        "isotope_error_prior": pa.array(prior.astype(np.float32), pa.float32()),
        "feature_id": pa.array(np.clip(fid_r, 0, None).astype(np.uint32), pa.uint32(), mask=fid_r < 0),
        "feature_intensity": mask_nan(G["fint"][R]),
        "feature_window_intensity": mask_nan(G["fwin"][R]),
        "isolation_fraction": mask_nan(G["frac"][R]),
        "chimeric_rank": pa.array(np.clip(rk_r, 0, 65535).astype(np.uint16), pa.uint16(), mask=rk_r < 0),
        "feature_iso_cos": mask_nan(G["cos"][R]),
        "feature_rt_apex": mask_nan(G["frt"][R]),
        "instrument_agrees": pa.array(G["agree"][R], pa.bool_()),
        "instrument_isotope_offset": pa.array(np.clip(io_r, -128, 127).astype(np.int8), pa.int8(), mask=io_r == -128),
        "n_features_in_window": pa.array(np.minimum(nfeat[rspec], 65535).astype(np.uint16), pa.uint16()),
    }, schema=HYPOTHESIS_SCHEMA)
    meta = {"version": HYPOTHESIS_VERSION, "run_id": int(ms1.run_id), "source_codes": SOURCE_NAME,
            "isotope_errors": list(p.isotope_errors), "params": asdict(p),
            "priors": {k_: {str(e_): v_ for e_, v_ in d.items()} for k_, d in pri.items()},
            "calibration": None if calib is None else calib.summary(), "feature_params": fr.params}
    tbl = tbl.replace_schema_metadata({"vulcan_hypotheses": json.dumps(meta, default=str)})
    n_groups = np.bincount(G["spec"], minlength=n2) if ng else np.zeros(n2, np.int64)
    stats = dict(spectrum_id=info["spectrum_id"], rt=info["rt"], total=total, n_feat_in_window=nfeat,
                 master_found=info["master_found"], inst_z=info["z"], inst_mz=info["pmz"], has_agree=has_agree,
                 has_feat=has_feat, n_groups=n_groups, n_rows=n_groups * ne,
                 n_feature_groups=np.bincount(fs, minlength=n2) if n_f else np.zeros(n2, np.int64))
    return tbl, stats


def hypothesis_stats(tbl: pa.Table, st: dict) -> dict:
    """Summary numbers for metrics.json."""
    n2 = len(st["spectrum_id"])
    src = tbl.column("source").to_numpy()
    first_iso = tbl.column("isotope_error").to_numpy() == tbl.column("isotope_error").to_numpy()[0] if tbl.num_rows \
        else np.zeros(0, bool)
    g_src = src[first_iso]
    z = st["inst_z"]
    ng = st["n_groups"]
    ag = tbl.filter(pc.and_(tbl.column("instrument_agrees"), pc.equal(tbl.column("isotope_error"), 0)))
    ioff = ag.column("instrument_isotope_offset").to_numpy(zero_copy_only=False)
    top = tbl.filter(pc.and_(pc.equal(tbl.column("chimeric_rank"), 1), pc.equal(tbl.column("isotope_error"), 0)))
    top_z = top.column("charge").to_numpy()
    top_sid = top.column("spectrum_id").to_numpy()
    zmap = dict(zip(st["spectrum_id"].tolist(), z.tolist()))
    zi_top = np.array([zmap[s] for s in top_sid.tolist()], dtype=np.int64) if len(top_sid) else np.zeros(0, np.int64)
    both = zi_top > 0
    return {
        "n_ms2": int(n2),
        "n_rows": int(tbl.num_rows),
        "n_precursor_groups": int(ng.sum()),
        "rows_per_spectrum_mean": float(tbl.num_rows / max(n2, 1)),
        "groups_per_spectrum_mean": float(ng.mean()) if n2 else 0.0,
        "groups_per_spectrum_hist": {str(i): int(c) for i, c in enumerate(np.bincount(np.minimum(ng, 10)))},
        "groups_by_source": {SOURCE_NAME[c]: int((g_src == c).sum()) for c in SOURCE_NAME},
        "frac_ms2_master_found": float(st["master_found"].mean()) if n2 else 0.0,
        "frac_ms2_with_instrument_charge": float((z > 0).mean()) if n2 else 0.0,
        "frac_ms2_with_feature_group": float(st["has_feat"].mean()) if n2 else 0.0,
        "frac_ms2_instrument_agrees_with_a_feature": float(st["has_agree"].mean()) if n2 else 0.0,
        "frac_ms2_fallback_only": float(((z <= 0) & ~st["has_feat"]).mean()) if n2 else 0.0,
        "frac_ms2_ge2_feature_groups(chimeric)": float((st["n_feature_groups"] >= 2).mean()) if n2 else 0.0,
        "features_in_window_mean": float(st["n_feat_in_window"].mean()) if n2 else 0.0,
        "instrument_vs_feature_isotope_offset_hist": {str(int(k_)): int(v) for k_, v in
                                                      zip(*np.unique(ioff, return_counts=True))},
        "top_feature_charge_vs_instrument": {
            "n_both": int(both.sum()),
            "agree_frac": float((top_z[both] == zi_top[both]).mean()) if both.any() else None,
        },
    }


def write_table(table: pa.Table, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with pa.OSFile(path + ".partial", "wb") as sink, pa.ipc.new_file(sink, table.schema) as w:
        w.write_table(table)
    os.replace(path + ".partial", path)


def read_table(path: str) -> pa.Table:
    with pa.OSFile(path, "rb") as src:
        return pa.ipc.open_file(src).read_all()
