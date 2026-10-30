"""Label-free quantification from vulcan.ms1 features (DESIGN §4 L5).

Pipeline (per experiment = list of runs):
1. `link_psms`: each accepted PSM -> an MS1 feature of its run.  First through the precursor-hypothesis table
   (the spectrum's co-isolated `ms1_feature` groups with the PSM's charge whose feature mass matches the
   peptidoform mass within the run's empirical ppm window), then, for PSMs without a hypothesis link (instrument
   or fallback groups), a targeted lookup in the feature table (same charge, mass window, PSM RT inside the
   feature's elution bounds +-0.05 min).  The run's ppm offset/SD are estimated from the hypothesis links.
2. `precursor_table`: per run, (pf_id, charge) -> the feature linked by most PSMs (ties: larger area).
   Intensity = feature `intensity_area` (sum of isotope-hill areas).  Provenance 'ms2'.
3. `Alignment`: cross-run RT alignment to a consensus scale (all runs fitted to the median of the aligned apex
   RTs; 2 rounds), monotone piecewise-linear maps on quantile knots, with local robust residual SDs.  O(n_runs).
4. `match_between_runs` (opt-in): for every donor precursor (quantified in >= 1 run) missing in run r, the best
   feature in r inside the predicted (m/z, RT) window is scored with a FIXED a-priori score (mass, RT,
   isotope/trace shape, intensity consistency; nothing fitted on labels).  Transfer FDR is estimated PIP-ECHO
   style with three decoy sources:
     * RT-shifted decoy peaks: same m/z, window centre moved by a random 4-8 windows (>= 0.5 min);
     * m/z-shifted decoy peaks: neutral mass + 11 x 1.0005 Da at the predicted RT (IonQuant-style);
     * decoy-peptide donors: decoy peptidoforms that reached the ID score threshold are transferred exactly like
       targets; their accepted transfers estimate the propagated identification error.
   A target transfer is accepted at threshold s only if it beats both of its decoy peaks;
   FDR(s) = (max(D_rt(s), D_mz(s)) + 1) / T(s) + D_pep(s) / T(s), made monotone -> transfer q-value.
   Features already explained by an MS2-identified precursor (target or decoy) of the acceptor run are excluded
   from every window (target and decoy alike): MS2 evidence outranks inference, and on the D9 human-only control
   runs this conflict accounted for 58% of the false yeast/E. coli transfers (dev diagnostic, M6).
   q-values are computed per acceptor run (run-specific transfer FDR; pooled is an option).
   Transfers carry provenance 'mbr', the transfer q-value and the donor-run count, and are never IDs.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numba as nb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.ipc as ipc

PROTON = 1.007276466812
MZ_DECOY_SHIFT_DA = 11 * 1.0005


# ----------------------------------------------------------------------------------------------- loading
@dataclass
class RunMS1:
    run: str
    run_id: int
    feat: pd.DataFrame          # FEATURE_SCHEMA columns, row i == feature_id i
    hyp: pd.DataFrame           # ms1_feature hypothesis groups only (spectrum_id, charge, feature_id, fwi)
    ms2_store_index: np.ndarray # spec_idx -> store row index
    rt_range: tuple = (0.0, 0.0)
    ppm_offset: float = 0.0
    ppm_sd: float = 2.0
    stats: dict = field(default_factory=dict)


def _read_ipc(path, columns=None):
    with pa.memory_map(path) as src:
        t = ipc.open_file(src).read_all()
    if columns:
        t = t.select(columns)
    return t


def load_run_ms1(ms1_dir: str, store_path: str, run: str) -> RunMS1:
    ft = _read_ipc(os.path.join(ms1_dir, f"{run}.features.arrow")).to_pandas()
    assert (ft["feature_id"].to_numpy() == np.arange(len(ft))).all(), "feature table must be ordered by feature_id"
    ht = _read_ipc(os.path.join(ms1_dir, f"{run}.hypotheses.arrow"),
                   ["spectrum_id", "source", "charge", "feature_id", "feature_window_intensity"])
    h = ht.to_pandas()
    h = h[h["source"] == 0].drop_duplicates(["spectrum_id", "charge", "feature_id"]).copy()
    h["feature_id"] = h["feature_id"].astype(np.int64)
    h["charge"] = h["charge"].astype(np.int8)
    h["spectrum_id"] = h["spectrum_id"].astype(np.uint64)
    st = _read_ipc(store_path, ["index", "ms_level", "rt_min", "run_id"])
    lvl = st.column("ms_level").to_numpy()
    idx = st.column("index").to_numpy()
    rt = st.column("rt_min").to_numpy()
    run_id = int(st.column("run_id")[0].as_py())
    return RunMS1(run, run_id, ft, h[["spectrum_id", "charge", "feature_id", "feature_window_intensity"]],
                  idx[lvl == 2].astype(np.int64), (float(rt.min()), float(rt.max())))


# ----------------------------------------------------------------------------------------------- linking
def _robust(x):
    x = np.asarray(x, float)
    if len(x) == 0:
        return 0.0, np.nan
    m = float(np.median(x))
    return m, float(1.4826 * np.median(np.abs(x - m)))


def link_psms(psm: pd.DataFrame, r: RunMS1, tol_k: float = 4.0, tol_min: float = 3.0, tol_max: float = 10.0,
              rt_pad: float = 0.05) -> pd.DataFrame:
    """psm: one run's accepted PSMs with columns spec_idx, charge, calcmass, rt.  Returns psm with feature_id
    (-1 if none), link_ppm and link ('hyp' | 'targeted' | 'none'); sets r.ppm_offset / r.ppm_sd."""
    p = psm.copy()
    p["spectrum_id"] = (np.uint64(r.run_id) << np.uint64(32)) | r.ms2_store_index[p["spec_idx"].to_numpy()].astype(np.uint64)
    p["_i"] = np.arange(len(p))
    p["charge"] = p["charge"].astype(np.int8)
    m = p[["_i", "spectrum_id", "charge", "calcmass"]].merge(r.hyp, on=["spectrum_id", "charge"], how="inner")
    fm = r.feat["neutral_mass"].to_numpy()[m["feature_id"].to_numpy()]
    m["ppm"] = (fm - m["calcmass"].to_numpy()) / m["calcmass"].to_numpy() * 1e6
    off, sd = _robust(m["ppm"][np.abs(m["ppm"]) < 20])
    sd = sd if np.isfinite(sd) and sd > 0 else 1.0
    tol = float(np.clip(tol_k * sd, tol_min, tol_max))
    r.ppm_offset, r.ppm_sd = off, sd
    m = m[np.abs(m["ppm"] - off) <= tol].sort_values(["_i", "feature_window_intensity"], ascending=[True, False])
    m = m.drop_duplicates("_i")
    fid = np.full(len(p), -1, np.int64)
    ppm = np.full(len(p), np.nan)
    fid[m["_i"].to_numpy()] = m["feature_id"].to_numpy()
    ppm[m["_i"].to_numpy()] = m["ppm"].to_numpy()
    how = np.where(fid >= 0, "hyp", "none").astype(object)
    # targeted fallback
    todo = np.flatnonzero(fid < 0)
    if len(todo):
        F = r.feat
        key = F["charge"].to_numpy().astype(np.float64) * 1e5 + F["neutral_mass"].to_numpy()
        o = np.argsort(key, kind="stable")
        ks = key[o]
        cm = p["calcmass"].to_numpy()[todo] * (1 + off * 1e-6)
        z = p["charge"].to_numpy()[todo].astype(np.float64)
        rt = p["rt"].to_numpy()[todo]
        lo = np.searchsorted(ks, z * 1e5 + cm * (1 - tol * 1e-6), "left")
        hi = np.searchsorted(ks, z * 1e5 + cm * (1 + tol * 1e-6), "right")
        rs, re_, area, nm = (F["rt_start"].to_numpy(), F["rt_end"].to_numpy(), F["intensity_area"].to_numpy(),
                             F["neutral_mass"].to_numpy())
        for j in range(len(todo)):
            if hi[j] <= lo[j]:
                continue
            c = o[lo[j]:hi[j]]
            ok = (rs[c] - rt_pad <= rt[j]) & (rt[j] <= re_[c] + rt_pad)
            if ok.any():
                c = c[ok]
                b = c[np.argmax(area[c])]
                i = todo[j]
                fid[i] = b
                ppm[i] = (nm[b] - p["calcmass"].to_numpy()[i]) / p["calcmass"].to_numpy()[i] * 1e6
                how[i] = "targeted"
    p["feature_id"], p["link_ppm"], p["link"] = fid, ppm, how
    r.stats.update({"ppm_offset": off, "ppm_sd": sd, "ppm_tol": tol, "n_psms": int(len(p)),
                    "linked_hyp": int((how == "hyp").sum()), "linked_targeted": int((how == "targeted").sum()),
                    "unlinked": int((how == "none").sum())})
    return p.drop(columns=["_i"])


def precursor_table(linked: pd.DataFrame, r: RunMS1) -> pd.DataFrame:
    """(pf_id, charge) -> dominant linked feature; rows without a feature are kept with NaN intensity."""
    g = linked.groupby(["pf_id", "charge"], sort=False)
    base = g.agg(n_psms=("spec_idx", "size"), rt_psm=("rt", "median"), calcmass=("calcmass", "first"),
                 best_score=("score", "max")).reset_index()
    L = linked[linked["feature_id"] >= 0]
    if len(L):
        cnt = L.groupby(["pf_id", "charge", "feature_id"]).size().rename("k").reset_index()
        cnt["area"] = r.feat["intensity_area"].to_numpy()[cnt["feature_id"].to_numpy()]
        cnt = cnt.sort_values(["pf_id", "charge", "k", "area"], ascending=[True, True, False, False])
        best = cnt.drop_duplicates(["pf_id", "charge"])[["pf_id", "charge", "feature_id"]]
        base = base.merge(best, on=["pf_id", "charge"], how="left")
    else:
        base["feature_id"] = np.nan
    fid = base["feature_id"].fillna(-1).to_numpy().astype(np.int64)
    ok = fid >= 0
    F = r.feat
    for col, src in (("intensity", "intensity_area"), ("rt_apex", "rt_apex"), ("iso_cos", "iso_cos"),
                     ("trace_cos", "trace_cos"), ("feat_mass", "neutral_mass")):
        v = np.full(len(base), np.nan)
        v[ok] = F[src].to_numpy()[fid[ok]]
        base[col] = v
    base["feature_id"] = fid
    base["run"] = r.run
    base["provenance"] = "ms2"
    base["transfer_q"] = np.nan
    return base


# ----------------------------------------------------------------------------------------------- alignment
def _fit_monotone(x, y, n_knots=30, min_per=15):
    o = np.argsort(x)
    x, y = x[o], y[o]
    k = int(max(2, min(n_knots, len(x) // min_per)))
    edges = np.quantile(x, np.linspace(0, 1, k + 1))
    kx, ky, ks = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (x >= a) & (x <= b)
        if s.sum() < 3:
            continue
        kx.append(np.median(x[s]))
        ky.append(np.median(y[s] - x[s]))
    kx = np.array(kx)
    ky = np.array(ky) + kx
    ky = np.maximum.accumulate(ky)
    ky = ky + np.arange(len(ky)) * 1e-9          # strictly increasing for the inverse
    res = y - _apply(kx, ky, x)
    for a, b in zip(edges[:-1], edges[1:]):
        s = (x >= a) & (x <= b)
        ks.append(1.4826 * np.median(np.abs(res[s] - np.median(res[s]))) if s.sum() >= 3 else np.nan)
    ks = np.array(ks[:len(kx)], float)
    glob = 1.4826 * np.median(np.abs(res - np.median(res)))
    ks = np.where(np.isfinite(ks), ks, glob)
    return kx, ky, ks, float(glob)


def _apply(kx, ky, x):
    x = np.asarray(x, float)
    if len(kx) == 1:
        return x + (ky[0] - kx[0])
    y = np.interp(x, kx, ky)
    lo, hi = x < kx[0], x > kx[-1]
    s0 = (ky[1] - ky[0]) / max(kx[1] - kx[0], 1e-9)
    s1 = (ky[-1] - ky[-2]) / max(kx[-1] - kx[-2], 1e-9)
    y[lo] = ky[0] + (x[lo] - kx[0]) * s0
    y[hi] = ky[-1] + (x[hi] - kx[-1]) * s1
    return y


class Alignment:
    """Star alignment of every run to a consensus RT scale (2 rounds).  maps[run] = (kx, ky, local_sd, global_sd)."""

    def __init__(self, tables: dict, rounds: int = 2):
        rows = []
        for run, t in tables.items():
            q = t[np.isfinite(t["rt_apex"])]
            rows.append(pd.DataFrame({"run": run, "key": q["pf_id"].astype(np.uint64).astype(str) + "/" +
                                      q["charge"].astype(str), "rt": q["rt_apex"].to_numpy()}))
        L = pd.concat(rows)
        M = L.pivot_table(index="key", columns="run", values="rt", aggfunc="median")
        self.runs = list(tables)
        ref = M.notna().sum().idxmax()
        cons = M[ref].copy()
        self.maps = {}
        for _ in range(rounds):
            aligned = {}
            for run in self.runs:
                s = M[run].notna() & cons.notna()
                if s.sum() < 50:
                    raise ValueError(f"too few shared precursors to align {run}: {int(s.sum())}")
                self.maps[run] = _fit_monotone(M.loc[s, run].to_numpy(), cons[s].to_numpy())
                v = M[run].copy()
                ok = v.notna()
                v[ok] = _apply(self.maps[run][0], self.maps[run][1], v[ok].to_numpy())
                aligned[run] = v
            cons = pd.DataFrame(aligned).median(axis=1)
        self.consensus = cons
        self.n_anchor = {run: int((M[run].notna() & cons.notna()).sum()) for run in self.runs}

    def to_consensus(self, run, rt):
        kx, ky, _, _ = self.maps[run]
        return _apply(kx, ky, rt)

    def to_run(self, run, rt_cons):
        kx, ky, _, _ = self.maps[run]
        return _apply(ky, kx, rt_cons)

    def local_sd(self, run, rt_cons):
        kx, ky, ks, g = self.maps[run]
        return np.interp(rt_cons, ky, ks)

    def summary(self):
        return {run: {"n_anchors": self.n_anchor[run], "resid_sd_min": self.maps[run][3],
                      "n_knots": int(len(self.maps[run][0]))} for run in self.runs}


# ----------------------------------------------------------------------------------------------- MBR
@nb.njit(parallel=True, cache=True)
def _best_in_window(fkey, frt, fmass, fiso, ftr, flogi, fok, o, qz, qmass, qrt, qw, qsd_rt, qtol_ppm, qsd_ppm, qmu, sd_i,
                    w_iso, w_tr):
    n = len(qz)
    best = np.full(n, -np.inf)
    arg = np.full(n, -1, np.int64)
    for i in nb.prange(n):
        lo_k = qz[i] * 1e5 + qmass[i] * (1 - qtol_ppm[i] * 1e-6)
        hi_k = qz[i] * 1e5 + qmass[i] * (1 + qtol_ppm[i] * 1e-6)
        a = np.searchsorted(fkey, lo_k)
        b = np.searchsorted(fkey, hi_k, side="right")
        for j in range(a, b):
            f = o[j]
            if not fok[f]:
                continue
            d = frt[f] - qrt[i]
            if abs(d) > qw[i]:
                continue
            dp = (fmass[f] - qmass[i]) / qmass[i] * 1e6
            s = (-0.5 * (dp / qsd_ppm[i]) ** 2 - 0.5 * (d / qsd_rt[i]) ** 2 + w_iso * fiso[f] + w_tr * ftr[f]
                 - 0.5 * ((flogi[f] - qmu[i]) / sd_i) ** 2)
            if s > best[i]:
                best[i] = s
                arg[i] = f
    return best, arg


@dataclass
class MBRParams:
    rt_k: float = 4.0            # window half-width = max(rt_k * local SD, rt_min_w)
    rt_min_w: float = 0.05       # minutes
    ppm_k: float = 3.0           # mass half-width = max(ppm_k * run SD, ppm_min)
    ppm_min: float = 2.0
    w_iso: float = 2.0
    w_tr: float = 1.0
    rt_shift_lo: float = 4.0     # RT decoy shift in windows (and at least rt_shift_min minutes)
    rt_shift_hi: float = 8.0
    rt_shift_min: float = 0.5
    seed: int = 11
    exclude_claimed: bool = True  # same rule in target and decoy windows (exchangeability)
    per_run_q: bool = True        # transfer q-values per acceptor run (False: pooled over the experiment)


def match_between_runs(tables: dict, ms1: dict, align: Alignment, decoy_donor: np.ndarray,
                       params: MBRParams | None = None, acceptors=None, donor_keys=None) -> tuple[pd.DataFrame, dict]:
    """tables: run -> precursor_table (incl. decoy-peptide donors); decoy_donor: pf_ids that are decoy donors.
    Returns (transfers, info).  transfers: one row per (acceptor run, donor precursor) candidate with
    St/Sd_rt/Sd_mz, the chosen feature, intensity, transfer_q / transfer_q_peak (targets only; NaN for decoy-peptide
    rows).  donor_keys: optional set of (pf_id, charge) allowed as donors."""
    P = params or MBRParams()
    rng = np.random.default_rng(P.seed)
    allp = pd.concat([t.assign(run=run) for run, t in tables.items()])
    q = allp[allp["feature_id"] >= 0].copy()
    q["rt_cons"] = np.nan
    for run in tables:
        s = q["run"] == run
        q.loc[s, "rt_cons"] = align.to_consensus(run, q.loc[s, "rt_apex"].to_numpy())
    q["log2i"] = np.log2(q["intensity"].clip(lower=1))
    # run-level median centring of log intensities for the intensity-consistency term
    med = q.groupby("run")["log2i"].median()
    q["log2i_c"] = q["log2i"] - q["run"].map(med)
    don = q.groupby(["pf_id", "charge"]).agg(rt_cons=("rt_cons", "median"), mu=("log2i_c", "mean"),
                                              calcmass=("calcmass", "first"), n_donor_runs=("run", "nunique"))
    don = don.reset_index()
    resid = q["log2i_c"] - q.set_index(["pf_id", "charge"]).index.map(don.set_index(["pf_id", "charge"])["mu"]).to_numpy()
    multi = q.groupby(["pf_id", "charge"])["run"].transform("nunique").to_numpy() > 1
    sd_i = float(max(1.0, 1.4826 * np.median(np.abs(resid[multi])))) if multi.any() else 1.5
    don["is_decoy_donor"] = np.isin(don["pf_id"].to_numpy(), decoy_donor)
    if donor_keys is not None:      # e.g. donors restricted to a stricter identification level
        don = don[[k in donor_keys for k in zip(don["pf_id"], don["charge"])]].reset_index(drop=True)
    present = set(zip(allp["run"], allp["pf_id"], allp["charge"]))   # identified in run (with or w/o feature)
    out = []
    info = {"sd_log2_intensity": sd_i, "runs": {}}
    for run in (acceptors or list(tables)):
        r = ms1[run]
        miss = np.array([(run, a, b) not in present for a, b in zip(don["pf_id"], don["charge"])])
        D = don[miss].reset_index(drop=True)
        if not len(D):
            continue
        F = r.feat
        fkey_raw = F["charge"].to_numpy().astype(np.float64) * 1e5 + F["neutral_mass"].to_numpy()
        o = np.argsort(fkey_raw, kind="stable")
        fkey = fkey_raw[o]
        frt = F["rt_apex"].to_numpy().astype(np.float64)
        fmass = F["neutral_mass"].to_numpy()
        fiso = F["iso_cos"].to_numpy().astype(np.float64)
        ftr = F["trace_cos"].to_numpy().astype(np.float64)
        flogi = np.log2(np.clip(F["intensity_area"].to_numpy().astype(np.float64), 1, None)) - med[run]
        t_hat = align.to_run(run, D["rt_cons"].to_numpy())
        sd_rt = np.maximum(align.local_sd(run, D["rt_cons"].to_numpy()), 1e-3)
        W = np.maximum(P.rt_k * sd_rt, P.rt_min_w)
        sd_rt_eff = W / P.rt_k
        mass = D["calcmass"].to_numpy() * (1 + r.ppm_offset * 1e-6)
        tol = np.full(len(D), max(P.ppm_k * r.ppm_sd, P.ppm_min))
        sdp = np.full(len(D), max(r.ppm_sd, 0.5))
        z = D["charge"].to_numpy().astype(np.float64)
        mu = D["mu"].to_numpy()
        fok = np.ones(len(F), np.bool_)
        if P.exclude_claimed:     # a feature explained by an MS2-identified precursor of this run is not transferable
            cl = tables[run]["feature_id"].to_numpy()
            fok[cl[cl >= 0]] = False
        args = (fkey, frt, fmass, fiso, ftr, flogi, fok, o)
        St, At = _best_in_window(*args, z, mass, t_hat, W, sd_rt_eff, tol, sdp, mu, sd_i, P.w_iso, P.w_tr)
        # RT-shifted decoy window
        mag = np.maximum(rng.uniform(P.rt_shift_lo, P.rt_shift_hi, len(D)) * W, P.rt_shift_min)
        sign = np.where(rng.random(len(D)) < 0.5, -1.0, 1.0)
        lo_rt, hi_rt = r.rt_range
        t_d = t_hat + sign * mag
        bad = (t_d < lo_rt + W) | (t_d > hi_rt - W)
        t_d = np.where(bad, t_hat - sign * mag, t_d)
        bad2 = (t_d < lo_rt + W) | (t_d > hi_rt - W)
        Sr, Ar = _best_in_window(*args, z, mass, t_d, W, sd_rt_eff, tol, sdp, mu, sd_i, P.w_iso, P.w_tr)
        Sr[bad2] = -np.inf
        Sm, Am = _best_in_window(*args, z, mass + MZ_DECOY_SHIFT_DA, t_hat, W, sd_rt_eff, tol, sdp, mu, sd_i,
                                 P.w_iso, P.w_tr)
        T = pd.DataFrame({"run": run, "pf_id": D["pf_id"], "charge": D["charge"], "calcmass": D["calcmass"],
                          "n_donor_runs": D["n_donor_runs"], "is_decoy_donor": D["is_decoy_donor"],
                          "rt_pred": t_hat, "rt_window": W, "St": St, "Sd_rt": Sr, "Sd_mz": Sm, "feature_id": At})
        ok = At >= 0
        for col, src in (("intensity", "intensity_area"), ("rt_apex", "rt_apex"), ("feat_mass", "neutral_mass")):
            v = np.full(len(T), np.nan)
            v[ok] = F[src].to_numpy()[At[ok]]
            T[col] = v
        out.append(T)
        info["runs"][run] = {"n_candidates": int(len(T)), "n_with_target_peak": int(ok.sum()),
                             "median_rt_window_min": float(np.median(W)), "ppm_tol": float(tol[0])}
    T = pd.concat(out, ignore_index=True)
    T["transfer_q"], T["transfer_q_peak"] = np.nan, np.nan
    groups = [np.flatnonzero(T["run"].to_numpy() == r) for r in T["run"].unique()] if P.per_run_q else [np.arange(len(T))]
    for g in groups:     # run-specific transfer FDR: the null fraction differs between acceptor runs
        qa, qb = transfer_qvalues(T.iloc[g], P.seed)
        T.iloc[g, T.columns.get_loc("transfer_q")] = qa
        T.iloc[g, T.columns.get_loc("transfer_q_peak")] = qb
    T["provenance"] = "mbr"
    return T, info


def transfer_qvalues(T: pd.DataFrame, seed: int = 11) -> np.ndarray:
    """PIP-ECHO-style transfer q-values for target-donor rows (NaN for decoy-donor rows and for candidates that
    do not beat both decoy peaks).  Returns (q_total, q_peak): q_total includes the decoy-peptide (propagated
    identification error) term, q_peak only the RT/m-z decoy-peak (peak-matching error) term."""
    rng = np.random.default_rng(seed + 1)
    jit = rng.random(len(T)) * 1e-9                      # random tie-break, label independent
    St = T["St"].to_numpy() + jit
    Sr = T["Sd_rt"].to_numpy()
    Sm = T["Sd_mz"].to_numpy()
    tgt = ~T["is_decoy_donor"].to_numpy()
    win = tgt & np.isfinite(St) & (St > Sr) & (St > Sm)
    dr = tgt & np.isfinite(Sr) & (Sr >= St)                 # RT-decoy wins its competition
    dm = tgt & np.isfinite(Sm) & (Sm >= St)
    dp = ~tgt & np.isfinite(St) & (St > Sr) & (St > Sm)     # decoy-peptide transfer that would be accepted
    s_all = np.concatenate([St[win], Sr[dr], Sm[dm], St[dp]])
    lab = np.concatenate([np.zeros(win.sum()), np.ones(dr.sum()), 2 * np.ones(dm.sum()), 3 * np.ones(dp.sum())])
    o = np.argsort(-s_all, kind="stable")
    lab = lab[o]
    cT = np.cumsum(lab == 0)
    cR = np.cumsum(lab == 1)
    cM = np.cumsum(lab == 2)
    cP = np.cumsum(lab == 3)
    out = []
    for fdr in ((np.maximum(cR, cM) + 1) / np.maximum(cT, 1) + cP / np.maximum(cT, 1),
                (np.maximum(cR, cM) + 1) / np.maximum(cT, 1)):
        qv = np.minimum(np.minimum.accumulate(fdr[::-1])[::-1], 1.0)
        q_sorted = np.full(len(s_all), np.nan)
        q_sorted[o] = qv
        q = np.full(len(T), np.nan)
        q[np.flatnonzero(win)] = q_sorted[:win.sum()]
        out.append(q)
    return out[0], out[1]


# ----------------------------------------------------------------------------------------------- roll-up
def topn(M: pd.DataFrame, prot: pd.Series, n: int = 3) -> pd.DataFrame:
    """M: precursor x run linear intensities (NaN missing); prot: precursor -> protein key.  Per run, mean of the
    n most intense precursors of the protein."""
    L = M.copy()
    L["prot"] = prot.reindex(L.index).to_numpy()
    out = {}
    for run in M.columns:
        s = L[[run, "prot"]].dropna()
        s = s.sort_values(run, ascending=False).groupby("prot").head(n)
        out[run] = s.groupby("prot")[run].mean()
    return pd.DataFrame(out)


def maxlfq(M: pd.DataFrame, prot: pd.Series, min_ratios: int = 1) -> pd.DataFrame:
    """MaxLFQ (Cox 2014): per protein, pairwise run log-ratios = median over shared precursors (>= min_ratios),
    least-squares run profile per connected component, scaled so that the component's summed profile equals the
    summed precursor intensity.  Returns protein x run linear intensities."""
    L = np.log2(M)
    runs = list(M.columns)
    R = len(runs)
    res = {}
    for p, idx in prot.reindex(M.index).groupby(prot.reindex(M.index)).groups.items():
        X = L.loc[idx].to_numpy()
        obs = np.isfinite(X)
        A, b = [], []
        rat = np.full((R, R), np.nan)
        for j in range(R):
            for k in range(j + 1, R):
                s = obs[:, j] & obs[:, k]
                if s.sum() >= min_ratios:
                    rat[j, k] = np.median(X[s, j] - X[s, k])
        prof = np.full(R, np.nan)
        # connected components over runs with ratios
        adj = np.isfinite(rat) | np.isfinite(rat.T)
        seen = np.zeros(R, bool)
        for s0 in range(R):
            if seen[s0] or not obs[:, s0].any():
                continue
            comp, stack = [], [s0]
            seen[s0] = True
            while stack:
                u = stack.pop()
                comp.append(u)
                for v in np.flatnonzero(adj[u] & ~seen):
                    seen[v] = True
                    stack.append(v)
            comp = sorted(comp)
            if len(comp) == 1:
                prof[comp[0]] = np.log2(np.nansum(2 ** X[:, comp[0]]))
                continue
            ci = {c: i for i, c in enumerate(comp)}
            rows, rhs = [], []
            for j in comp:
                for k in comp:
                    if j < k and np.isfinite(rat[j, k]):
                        r = np.zeros(len(comp))
                        r[ci[j]], r[ci[k]] = 1, -1
                        rows.append(r)
                        rhs.append(rat[j, k])
            Am = np.vstack(rows + [np.ones(len(comp))])
            bm = np.array(rhs + [0.0])
            sol = np.linalg.lstsq(Am, bm, rcond=None)[0]
            tot = np.log2(np.nansum(2 ** X[:, comp]))
            shift = tot - np.log2(np.sum(2 ** sol))
            prof[comp] = sol + shift
        res[p] = 2 ** prof
    return pd.DataFrame(res, index=runs).T


def directlfq_rollup(M: pd.DataFrame, prot: pd.Series, workdir: str) -> pd.DataFrame:
    """directLFQ (Ammar 2023) on the precursor matrix (its own normalisation + protein estimation)."""
    import directlfq.lfq_manager as lfqm
    os.makedirs(workdir, exist_ok=True)
    df = M.copy()
    df.insert(0, "ion", [str(i) for i in df.index])
    df.insert(0, "protein", prot.reindex(M.index).to_numpy())
    df = df.dropna(subset=["protein"])
    path = os.path.join(workdir, "precursors.aq_reformat.tsv")
    df.fillna(0).to_csv(path, sep="\t", index=False)
    lfqm.run_lfq(path, num_cores=1)
    out = pd.read_csv(path.replace(".aq_reformat.tsv", ".aq_reformat.tsv.protein_intensities.tsv"), sep="\t")
    return out.set_index("protein")[list(M.columns)].replace(0, np.nan)


# ----------------------------------------------------------------------------------------------- XIC re-quant
NEUTRON = 1.0033548378


@dataclass
class MS1Peaks:
    rt: np.ndarray        # per MS1 scan
    off: np.ndarray       # n_scans + 1
    mz: np.ndarray        # flat, m/z ascending within a scan
    it: np.ndarray


def load_ms1_peaks(store_path: str) -> MS1Peaks:
    from vulcan.io import store as ST
    meta = ST.run_meta(store_path)
    assert meta.get("peak_order") == "mz_ascending", "XIC re-quant needs m/z-sorted MS1 peaks"
    t = ST.read_table(store_path, ["ms_level", "rt_min", "mz", "intensity"])
    lvl = t.column("ms_level").to_numpy()
    t = t.filter(pa.array(lvl == 1))
    pk = ST.peaks_numpy(t)
    return MS1Peaks(np.asarray(t.column("rt_min").to_numpy(), np.float64), pk.offsets.astype(np.int64),
                    pk.mz.astype(np.float64), pk.intensity.astype(np.float64))


@nb.njit(parallel=True, cache=True)
def _xic_iso(rt, off, mz, it, qmz, qz, qlo, qhi, tol_ppm, n_iso):
    """Per query: trapezoid area of each isotope's XIC (M+0..M+n_iso-1; per-scan most intense centroid within
    +-tol_ppm) over [qlo, qhi], and the number of MS1 scans with a mono signal."""
    n = len(qmz)
    out = np.zeros((n, n_iso))
    nsc = np.zeros(n, np.int64)
    for i in nb.prange(n):
        a = np.searchsorted(rt, qlo[i])
        b = np.searchsorted(rt, qhi[i], side="right")
        prev = np.zeros(n_iso)
        cur = np.zeros(n_iso)
        prev_t = 0.0
        for k in range(a, b):
            m0, m1 = off[k], off[k + 1]
            for j in range(n_iso):
                c = qmz[i] + j * NEUTRON / qz[i]
                lo = c * (1 - tol_ppm * 1e-6)
                hi = c * (1 + tol_ppm * 1e-6)
                p = m0 + np.searchsorted(mz[m0:m1], lo)
                v = 0.0
                while p < m1 and mz[p] <= hi:
                    if it[p] > v:
                        v = it[p]
                    p += 1
                cur[j] = v
            if cur[0] > 0:
                nsc[i] += 1
            if k > a:
                for j in range(n_iso):
                    out[i, j] += 0.5 * (cur[j] + prev[j]) * (rt[k] - prev_t)
            for j in range(n_iso):
                prev[j] = cur[j]
            prev_t = rt[k]
    return out, nsc


def combine_isotopes(areas: np.ndarray, neutral_mass: np.ndarray, iso: str = "avgproj") -> np.ndarray:
    """'sum': plain sum of the isotope areas (M6).  'avgproj' (M6d default): least-squares projection of the isotope
    areas onto the precursor's averagine envelope, rescaled to the envelope's total, i.e.
    sum_j a_j e_j / sum_j e_j^2 * sum_j e_j.  Weak isotopes get weight in proportion to their expected share, so
    intensity-dependent censoring of M+2 at the centroiding threshold (low-abundance side) no longer expands
    ratios (M6d diagnosis on D9: E. coli epsilon 0.162 -> 0.141, with the same CV)."""
    if iso == "sum":
        return areas.sum(1)
    from vulcan.ms1.averagine import get_table, STEP
    tab = get_table()
    k = areas.shape[1]
    x = np.clip(np.asarray(neutral_mass, float), 0, (len(tab) - 1) * STEP) / STEP
    i0 = np.minimum(np.floor(x).astype(np.int64), len(tab) - 2)
    f = (x - i0)[:, None]
    env = tab[i0, :k] * (1 - f) + tab[i0 + 1, :k] * f
    return (areas * env).sum(1) / (env * env).sum(1) * env.sum(1)


def _median_widths(tables, ms1):
    wid = []
    for run, t in tables.items():
        ok = t["feature_id"].to_numpy() >= 0
        F = ms1[run].feat
        w = np.full(len(t), np.nan)
        w[ok] = (F["rt_end"].to_numpy() - F["rt_start"].to_numpy())[t["feature_id"].to_numpy()[ok]]
        wid.append(pd.DataFrame({"pf_id": t["pf_id"].to_numpy(), "charge": t["charge"].to_numpy(), "w": w}))
    return pd.concat(wid).groupby(["pf_id", "charge"])["w"].median()


@dataclass
class XICParams:
    """Fixed-window XIC re-integration.  M6 = XICParams(iso='sum', width_mult=1.0, siblings=False).
    M6d defaults (chosen on the D9 selection runs A_REP1, A_REP2, B_REP1; confirmed on the held-out runs and D11):
    averagine projection of M+0..M+2, window 1.5 x the precursor's median mono-hill width, charge-sibling completion."""
    n_iso: int = 3
    iso: str = "avgproj"
    width_mult: float = 1.5
    min_half_width: float = 0.03
    tol_k: float = 3.0
    tol_min: float = 3.0
    siblings: bool = True
    sibling_min_scans: int = 3


def xic_requant(tables: dict, ms1: dict, stores: dict, align: Alignment, n_iso: int = 3, tol_k: float = 3.0,
                min_half_width: float = 0.03, transfers: pd.DataFrame | None = None,
                params: XICParams | None = None) -> dict:
    """Fixed-window XIC re-integration: for every (run, precursor) with a located MS1 feature (MS2-linked, or an
    MBR transfer), integrate the mono..+n_iso-1 isotope XICs over a window of the SAME width in every run
    (width_mult x the precursor's median mono-hill width across runs), centred on that run's feature apex, and
    combine the isotopes (params.iso).  Feature hills end where the signal reaches the noise floor, so hill areas
    shrink faster than linearly at low abundance (ratio expansion on D9: E. coli log2FC -2.28 vs -2.0 expected); a
    fixed window keeps the response linear.
    Returns {run: array of XIC areas aligned with tables[run] rows} and, if transfers is given, its 'xic' column.
    Without params the M6 behaviour (plain isotope sum, 1x width) is kept for the explicit n_iso/tol_k arguments."""
    P = params or XICParams(n_iso=n_iso, iso="sum", width_mult=1.0, min_half_width=min_half_width, tol_k=tol_k,
                            siblings=False)
    W = _median_widths(tables, ms1)
    out = {}
    for run, t in tables.items():
        Pk = load_ms1_peaks(stores[run])
        r = ms1[run]
        tol = max(P.tol_k * r.ppm_sd, P.tol_min)
        jobs = [("tab", t)]
        if transfers is not None:
            jobs.append(("tr", transfers[transfers["run"] == run]))
        for kind, x in jobs:
            fid = x["feature_id"].to_numpy().astype(np.int64)
            ok = fid >= 0
            z = x["charge"].to_numpy().astype(np.float64)
            nm = r.feat["neutral_mass"].to_numpy()[np.where(ok, fid, 0)]
            mz = nm / z + PROTON
            apex = x["rt_apex"].to_numpy().astype(np.float64)
            hw = W.reindex(pd.MultiIndex.from_arrays([x["pf_id"].to_numpy(), x["charge"].to_numpy()])).to_numpy() / 2
            hw = np.where(np.isfinite(hw), np.maximum(hw * P.width_mult, P.min_half_width), P.min_half_width)
            ar, _ = _xic_iso(Pk.rt, Pk.off, Pk.mz, Pk.it, mz, z, apex - hw, apex + hw, tol, P.n_iso)
            area = combine_isotopes(ar, nm, P.iso)
            area = np.where(ok & np.isfinite(apex), area, np.nan)
            if kind == "tab":
                out[run] = area
            else:
                transfers.loc[x.index, "xic"] = area
    return out


def charge_siblings(tables: dict, ms1: dict, stores: dict, params: XICParams | None = None,
                    exclude_pf=None) -> pd.DataFrame:
    """Charge-sibling completion (within run; M6d).  For every run and peptidoform quantified there (MS2-linked
    feature), integrate the peptidoform's other charge states that were identified in ANY run of the experiment but
    not in this run, in the same fixed window (anchor = the run's most intense quantified charge of the peptidoform:
    charge states co-elute) at the run-calibrated theoretical m/z.  The identity is this run's own MS2 ID of the
    peptidoform, so this is not a cross-run transfer (provenance 'ms2_sibling').  Values need a mono signal in
    >= sibling_min_scans MS1 scans.  Without it, charge-state sets differ between runs and summing charges to a
    peptidoform biases ratios (M6d: 7.6% of matched D9 peptidoforms, epsilon 0.30 vs 0.13)."""
    P = params or XICParams()
    W = _median_widths(tables, ms1)
    allc = pd.concat([t[["pf_id", "charge", "calcmass"]] for t in tables.values()]).drop_duplicates(["pf_id", "charge"])
    rows = []
    for run, t in tables.items():
        q = t[t["feature_id"] >= 0]
        if exclude_pf is not None:
            q = q[~q["pf_id"].isin(exclude_pf)]
        have = set(zip(t["pf_id"], t["charge"]))
        anc = q.sort_values("intensity", ascending=False).drop_duplicates("pf_id")
        aw = W.reindex(pd.MultiIndex.from_arrays([anc["pf_id"].to_numpy(), anc["charge"].to_numpy()])).to_numpy() / 2
        anc = anc.assign(hw=aw)[["pf_id", "rt_apex", "hw"]]
        c = allc.merge(anc, on="pf_id")
        c = c[[k not in have for k in zip(c["pf_id"], c["charge"])]].reset_index(drop=True)
        if not len(c):
            continue
        r = ms1[run]
        Pk = load_ms1_peaks(stores[run])
        tol = max(P.tol_k * r.ppm_sd, P.tol_min)
        z = c["charge"].to_numpy().astype(np.float64)
        nm = c["calcmass"].to_numpy() * (1 + r.ppm_offset * 1e-6)
        own = W.reindex(pd.MultiIndex.from_arrays([c["pf_id"].to_numpy(), c["charge"].to_numpy()])).to_numpy() / 2
        hw = np.where(np.isfinite(own), own, c["hw"].to_numpy())
        hw = np.where(np.isfinite(hw), np.maximum(hw * P.width_mult, P.min_half_width), P.min_half_width)
        apex = c["rt_apex"].to_numpy().astype(np.float64)
        ar, nsc = _xic_iso(Pk.rt, Pk.off, Pk.mz, Pk.it, nm / z + PROTON, z, apex - hw, apex + hw, tol, P.n_iso)
        v = combine_isotopes(ar, nm, P.iso)
        keep = (nsc >= P.sibling_min_scans) & np.isfinite(v) & (v > 0)
        c = c[keep]
        rows.append(pd.DataFrame({"run": run, "pf_id": c["pf_id"].to_numpy(), "charge": c["charge"].to_numpy(),
                                  "calcmass": c["calcmass"].to_numpy(), "rt_apex": c["rt_apex"].to_numpy(),
                                  "intensity": v[keep], "provenance": "ms2_sibling"}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["run", "pf_id", "charge", "calcmass", "rt_apex", "intensity", "provenance"])
