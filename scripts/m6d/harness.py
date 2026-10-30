"""M6d diagnostic harness: rebuild the M6 LFQ precursor tables (link + alignment) for any run set and re-integrate
MS1 XICs with variable isotope sets, windows, centring and peak rules, so that quant choices can be compared on the
same features.  Dev only (D9 / D11 / D12)."""
from __future__ import annotations

import os
import pickle
import sys

import numba as nb
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m6"))
from vulcan.quant import lfq as Q  # noqa: E402

NEUTRON = Q.NEUTRON
PROTON = Q.PROTON
ISO_OFFS = np.array([-1, 0, 1, 2, 3, 4], np.int64)     # column order of the per-isotope area matrix


def prepare(ids_dir, runs, ms1_dir, store_fn, alpha=0.01, cache=None):
    """Accepted PSMs -> linked precursor tables + alignment (exactly the M6 lfq_d9 path)."""
    if cache and os.path.exists(cache):
        return pickle.load(open(cache, "rb"))
    ids = pd.read_parquet(os.path.join(ids_dir, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(ids_dir, "peptides.parquet")).set_index("pf_id")
    acc = ids[ids[f"pass_{alpha}"] & ids["run"].isin(runs)].copy()
    ms1, tabs, linked = {}, {}, {}
    for run in runs:
        r = Q.load_run_ms1(ms1_dir, store_fn(run), run)
        ms1[run] = r
        L = Q.link_psms(acc[acc["run"] == run], r)
        linked[run] = L
        tabs[run] = Q.precursor_table(L, r)
    dec = pep.index[pep["is_decoy"]].to_numpy()
    al = Q.Alignment({r: t[~t["pf_id"].isin(dec)] for r, t in tabs.items()})
    obj = dict(ms1=ms1, tabs=tabs, linked=linked, align=al, pep=pep, dec=dec, runs=list(runs), ids=ids)
    if cache:
        pickle.dump(obj, open(cache, "wb"), protocol=5)
    return obj


@nb.njit(parallel=True, cache=True)
def xic_iso(rt, off, mz, it, qmz, qz, qlo, qhi, tol_ppm, offs, mode_sum):
    """Per query and isotope offset (offs, in neutron units from the mono m/z): trapezoid area of the per-scan max
    (mode_sum=False) or summed (True) centroid intensity within +-tol_ppm.  Also returns, per query, the number of
    scans with a mono signal and the RT of the scan with the largest summed (M..M+2) signal."""
    n = len(qmz)
    K = len(offs)
    out = np.zeros((n, K))
    nsc = np.zeros(n, np.int64)
    apx = np.full(n, np.nan)
    for i in nb.prange(n):
        a = np.searchsorted(rt, qlo[i])
        b = np.searchsorted(rt, qhi[i], side="right")
        prev = np.zeros(K)
        cur = np.zeros(K)
        prev_t = 0.0
        best = -1.0
        for k in range(a, b):
            m0, m1 = off[k], off[k + 1]
            for j in range(K):
                c = qmz[i] + offs[j] * NEUTRON / qz[i]
                lo = c * (1 - tol_ppm * 1e-6)
                hi = c * (1 + tol_ppm * 1e-6)
                p = m0 + np.searchsorted(mz[m0:m1], lo)
                v = 0.0
                while p < m1 and mz[p] <= hi:
                    if mode_sum:
                        v += it[p]
                    elif it[p] > v:
                        v = it[p]
                    p += 1
                cur[j] = v
            if k > a:
                for j in range(K):
                    out[i, j] += 0.5 * (cur[j] + prev[j]) * (rt[k] - prev_t)
            s3 = 0.0
            for j in range(K):
                if offs[j] >= 0 and offs[j] <= 2:
                    s3 += cur[j]
                if offs[j] == 0 and cur[j] > 0:
                    nsc[i] += 1
            if s3 > best:
                best = s3
                apx[i] = rt[k]
            for j in range(K):
                prev[j] = cur[j]
            prev_t = rt[k]
    return out, nsc, apx


_PEAKS = {}


def peaks(store_path):
    if store_path not in _PEAKS:
        _PEAKS[store_path] = Q.load_ms1_peaks(store_path)
    return _PEAKS[store_path]


def median_widths(obj):
    wid = []
    for run, t in obj["tabs"].items():
        ok = t["feature_id"].to_numpy() >= 0
        F = obj["ms1"][run].feat
        w = np.full(len(t), np.nan)
        w[ok] = (F["rt_end"].to_numpy() - F["rt_start"].to_numpy())[t["feature_id"].to_numpy()[ok]]
        wid.append(pd.DataFrame({"pf_id": t["pf_id"].to_numpy(), "charge": t["charge"].to_numpy(), "w": w}))
    return pd.concat(wid).groupby(["pf_id", "charge"])["w"].median()


def requant(obj, store_fn, rows: dict, centre="apex", width_mult=1.0, min_half_width=0.03, tol_k=3.0, tol_min=3.0,
            mode_sum=False, fixed_half_width=None, mz_source="feature"):
    """rows: run -> DataFrame with pf_id, charge, feature_id, rt_apex, calcmass (and rt_cons_pred if centre='pred').
    Returns run -> (iso_area[n, 6], n_scans, apex_rt)."""
    W = median_widths(obj)
    out = {}
    for run, x in rows.items():
        P = peaks(store_fn(run))
        r = obj["ms1"][run]
        tol = max(tol_k * r.ppm_sd, tol_min)
        z = x["charge"].to_numpy().astype(np.float64)
        fid = x["feature_id"].to_numpy().astype(np.int64)
        if mz_source == "feature":
            mz = r.feat["neutral_mass"].to_numpy()[np.where(fid >= 0, fid, 0)] / z + PROTON
        else:   # theoretical peptide m/z, run-calibrated
            mz = x["calcmass"].to_numpy() * (1 + r.ppm_offset * 1e-6) / z + PROTON
        if centre == "apex":
            c = x["rt_apex"].to_numpy().astype(np.float64)
        else:
            c = x["rt_pred"].to_numpy().astype(np.float64)
        if fixed_half_width is not None:
            hw = np.full(len(x), fixed_half_width)
        else:
            hw = W.reindex(pd.MultiIndex.from_arrays([x["pf_id"].to_numpy(), x["charge"].to_numpy()])).to_numpy() / 2
            if "hw" in x.columns:    # explicit half-width (e.g. sibling-charge rows: the identified charge's width)
                hw = np.where(np.isfinite(hw), hw, x["hw"].to_numpy())
            hw = np.where(np.isfinite(hw), np.maximum(hw * width_mult, min_half_width), min_half_width)
        a, n, ap = xic_iso(P.rt, P.off, P.mz, P.it, mz, z, c - hw, c + hw, tol, ISO_OFFS, mode_sum)
        bad = ~np.isfinite(c) | (fid < 0 if mz_source == "feature" else np.zeros(len(x), bool))
        a[bad] = np.nan
        out[run] = (a, n, ap)
    return out


def sibling_rows(obj, runs=None):
    """Charge-sibling completion (within run): for every run and peptidoform quantified there at charge z0 (MS2-linked
    feature), add the peptidoform's other charge states that were identified in ANY run of the experiment but not in
    this one.  Window = the z0 feature's apex in this run (charge states co-elute), m/z = calibrated theoretical m/z.
    Identity comes from this run's own MS2 ID of the peptidoform, so this is not a cross-run transfer."""
    runs = runs or obj["runs"]
    W = median_widths(obj)
    allc = pd.concat([obj["tabs"][r][["pf_id", "charge"]] for r in runs]).drop_duplicates()
    cm = pd.concat([obj["tabs"][r][["pf_id", "calcmass"]] for r in runs]).drop_duplicates("pf_id").set_index("pf_id")["calcmass"]
    out = {}
    for r in runs:
        t = obj["tabs"][r]
        t = t[t["feature_id"] >= 0]
        have = set(zip(obj["tabs"][r]["pf_id"], obj["tabs"][r]["charge"]))
        # anchor = the most intense quantified charge of the peptidoform in this run
        anc = t.sort_values("intensity", ascending=False).drop_duplicates("pf_id")
        anc_w = W.reindex(pd.MultiIndex.from_arrays([anc["pf_id"].to_numpy(), anc["charge"].to_numpy()])).to_numpy()
        anc = anc.assign(hw_anchor=anc_w / 2)[["pf_id", "rt_apex", "hw_anchor"]]
        c = allc.merge(anc, on="pf_id")
        c = c[[k not in have for k in zip(c["pf_id"], c["charge"])]].copy()
        c["calcmass"] = cm.reindex(c["pf_id"].to_numpy()).to_numpy()
        c["feature_id"] = -1
        c["hw"] = c["hw_anchor"]
        out[r] = c.reset_index(drop=True)
    return out
