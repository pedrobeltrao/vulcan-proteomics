#!/usr/bin/env python3
"""Feature-level comparison of vulcan.ms1 vs biosaur2 (and optional ID support from Sage PSMs).

    python scripts/m2/ms1_validate.py --ours F.arrow --biosaur2 B.features.tsv [--psms results.sage.tsv] --out m.json

Matching (any-match, not one-to-one): same charge, |neutral mass diff| <= ppm, RT intervals overlap (slack 0.05 min)
and |apex RT diff| <= rt_apex_tol. Both sides are restricted to the common m/z window [mzmin, mzmax] (biosaur2 default
350-1500). recall = fraction of biosaur2 features matched by ours; precision = fraction of ours matched by biosaur2.
ID support: target PSMs at q <= 1% (evaluation-only use of labels) whose (charge, calc mass, RT) falls on a feature.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow as pa

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
NEUTRON = 1.0033548378


def load_ours(path):
    with pa.OSFile(path, "rb") as f:
        t = pa.ipc.open_file(f).read_all()
    d = t.select(["mono_mz", "charge", "neutral_mass", "rt_apex", "rt_start", "rt_end", "intensity_mono_apex",
                  "intensity_apex", "n_isotopes", "n_scans"]).to_pandas()
    return d.rename(columns={"neutral_mass": "mass", "mono_mz": "mz", "intensity_mono_apex": "int_mono"})


def load_bs2(path):
    d = pd.read_csv(path, sep="\t", usecols=["massCalib", "rtApex", "intensityApex", "charge", "nIsotopes", "nScans",
                                             "mz", "rtStart", "rtEnd"])
    return d.rename(columns={"massCalib": "mass", "rtApex": "rt_apex", "rtStart": "rt_start", "rtEnd": "rt_end",
                             "intensityApex": "int_mono", "nIsotopes": "n_isotopes", "nScans": "n_scans"})


def pairs(A, B, ppm, rt_apex_tol, shift=0.0, slack=0.05):
    """All (i in A, j in B) with same charge, |mA + shift - mB| <= ppm, RT overlap and apex within tolerance."""
    oi, oj = [], []
    for z in np.intersect1d(A.charge.unique(), B.charge.unique()):
        a = A[A.charge == z]
        b = B[B.charge == z].sort_values("mass")
        bm = b.mass.to_numpy()
        am = a.mass.to_numpy() + shift
        lo = np.searchsorted(bm, am * (1 - ppm * 1e-6))
        hi = np.searchsorted(bm, am * (1 + ppm * 1e-6), side="right")
        n = hi - lo
        ia = np.repeat(np.arange(len(a)), n)
        jb = np.concatenate([np.arange(x, y) for x, y in zip(lo, hi)]) if len(a) else np.empty(0, int)
        ra0, ra1, rap = (a[c].to_numpy()[ia] for c in ("rt_start", "rt_end", "rt_apex"))
        rb0, rb1, rbp = (b[c].to_numpy()[jb] for c in ("rt_start", "rt_end", "rt_apex"))
        ok = (ra0 <= rb1 + slack) & (rb0 <= ra1 + slack) & (np.abs(rap - rbp) <= rt_apex_tol)
        oi.append(a.index.to_numpy()[ia[ok]])
        oj.append(b.index.to_numpy()[jb[ok]])
    return (np.concatenate(oi) if oi else np.empty(0, int)), (np.concatenate(oj) if oj else np.empty(0, int))


def compare(A, B, ppm=10.0, rt_apex_tol=0.5):
    ia, jb = pairs(A, B, ppm, rt_apex_tol)
    mA = np.zeros(len(A), bool)
    mB = np.zeros(len(B), bool)
    mA[A.index.get_indexer(np.unique(ia))] = True
    mB[B.index.get_indexer(np.unique(jb))] = True
    # one best pair per B feature (smallest |dppm|) for agreement stats
    dppm = (A.mass.to_numpy()[A.index.get_indexer(ia)] - B.mass.to_numpy()[B.index.get_indexer(jb)]) / \
        B.mass.to_numpy()[B.index.get_indexer(jb)] * 1e6
    df = pd.DataFrame({"i": ia, "j": jb, "d": dppm, "ad": np.abs(dppm)}).sort_values(["j", "ad"]).drop_duplicates("j")
    ai, bj = A.index.get_indexer(df.i), B.index.get_indexer(df.j)
    drt = A.rt_apex.to_numpy()[ai] - B.rt_apex.to_numpy()[bj]
    li = np.log10(np.maximum(A.int_mono.to_numpy()[ai], 1))
    lb = np.log10(np.maximum(B.int_mono.to_numpy()[bj], 1))
    # off-by-one-isotope mono disagreements among B features not matched exactly
    um = B[~mB]
    off = {}
    for s in (-1, 1):
        _, j2 = pairs(A, um, ppm, rt_apex_tol, shift=s * NEUTRON)
        off[f"{s:+d}"] = int(len(np.unique(j2)))
    # recall by biosaur2 intensity quartile / isotopes
    q = pd.qcut(B.int_mono.rank(method="first"), 4, labels=["Q1(low)", "Q2", "Q3", "Q4(high)"])
    rec_q = {str(k): float(mB[(q == k).to_numpy()].mean()) for k in q.cat.categories}
    q_a = pd.qcut(A.int_mono.rank(method="first"), 4, labels=["Q1(low)", "Q2", "Q3", "Q4(high)"])
    prec_q = {str(k): float(mA[(q_a == k).to_numpy()].mean()) for k in q_a.cat.categories}
    return {
        "n_ours": int(len(A)), "n_biosaur2": int(len(B)),
        "recall": float(mB.mean()), "precision": float(mA.mean()),
        "recall_by_biosaur2_intensity_quartile": rec_q, "precision_by_our_intensity_quartile": prec_q,
        "recall_biosaur2_nIsotopes>=3": float(mB[(B.n_isotopes >= 3).to_numpy()].mean()),
        "precision_ours_nIsotopes>=3": float(mA[(A.n_isotopes >= 3).to_numpy()].mean()),
        "unmatched_biosaur2_matched_with_mono_shift": off,
        "mass_diff_ppm": {"median": float(np.median(df.d)), "mad_sd": float(1.4826 * np.median(np.abs(df.d - np.median(df.d)))),
                          "p95_abs": float(np.quantile(np.abs(df.d), 0.95))},
        "rt_apex_diff_min": {"median": float(np.median(drt)), "mad_sd": float(1.4826 * np.median(np.abs(drt - np.median(drt)))),
                             "p95_abs": float(np.quantile(np.abs(drt), 0.95))},
        "log10_mono_apex_intensity_pearson": float(np.corrcoef(li, lb)[0, 1]) if len(li) > 2 else None,
        "params": {"ppm": ppm, "rt_apex_tol_min": rt_apex_tol, "matching": "any-match, same charge, RT overlap"},
    }


def id_support(feats: pd.DataFrame, psm: pd.DataFrame, ppm=10.0, rt_slack=0.1) -> dict:
    """Fraction of confident target PSMs whose (charge, calc mass, rt) lies on a feature (mono exact or any isotope)."""
    out = {}
    hit0 = np.zeros(len(psm), bool)
    hit_iso = np.zeros(len(psm), bool)
    for z in psm.charge.unique():
        f = feats[feats.charge == z].sort_values("mass")
        fm = f.mass.to_numpy()
        p = psm[psm.charge == z]
        for shift, arr in ((0.0, hit0), (None, hit_iso)):
            shifts = [0.0] if shift is not None else [-NEUTRON, NEUTRON, 2 * NEUTRON]
            for s in shifts:
                m = p.calcmass.to_numpy() + s
                lo = np.searchsorted(fm, m * (1 - ppm * 1e-6))
                hi = np.searchsorted(fm, m * (1 + ppm * 1e-6), side="right")
                n = hi - lo
                ip = np.repeat(np.arange(len(p)), n)
                jf = np.concatenate([np.arange(x, y) for x, y in zip(lo, hi)]) if len(p) else np.empty(0, int)
                rt = p.rt.to_numpy()[ip]
                ok = (f.rt_start.to_numpy()[jf] - rt_slack <= rt) & (rt <= f.rt_end.to_numpy()[jf] + rt_slack)
                arr[p.index.to_numpy()[np.unique(ip[ok])]] = True
    out["frac_psm_on_feature_mono"] = float(hit0.mean())
    out["frac_psm_on_feature_isotope_shifted_only"] = float((hit_iso & ~hit0).mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--biosaur2", required=True)
    ap.add_argument("--psms", default=None)
    ap.add_argument("--mzmin", type=float, default=350)
    ap.add_argument("--mzmax", type=float, default=1500)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    A = load_ours(a.ours)
    B = load_bs2(a.biosaur2)
    A = A[(A.mz >= a.mzmin) & (A.mz <= a.mzmax)].reset_index(drop=True)
    B = B[(B.mz >= a.mzmin) & (B.mz <= a.mzmax)].reset_index(drop=True)
    res = {"strict_10ppm": compare(A, B)}
    res["ours_nIso>=3_vs_biosaur2_nIso>=3"] = {k: v for k, v in compare(
        A[A.n_isotopes >= 3].reset_index(drop=True), B[B.n_isotopes >= 3].reset_index(drop=True)).items()
        if k in ("n_ours", "n_biosaur2", "recall", "precision")}
    if a.psms:
        p = pd.read_csv(a.psms, sep="\t", usecols=["label", "spectrum_q", "charge", "calcmass", "rt", "peptide"])
        p = p[(p.label == 1) & (p.spectrum_q <= 0.01)].reset_index(drop=True)
        mz = p.calcmass / p.charge + 1.007276
        p = p[(mz >= a.mzmin) & (mz <= a.mzmax)].reset_index(drop=True)
        res["id_support"] = {"n_psms": int(len(p)), "ours": id_support(A, p), "biosaur2": id_support(B, p)}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
