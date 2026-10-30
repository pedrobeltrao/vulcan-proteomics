"""DESIGN §6.8 quant metrics for two-condition multi-species benchmarks (ProteoBench HYE convention).

* epsilon (ProteoBench): per precursor, mean log2 intensity in A minus mean in B, minus the species' expected
  log2 fold change; summarised as median / mean |epsilon| per species and overall, with the number of
  precursors (a min-observation filter per condition is applied).  Contaminant and multi-species precursors are
  removed, as ProteoBench does.
* CV: per analyte and condition, SD / mean of the linear intensities (values present >= min_obs); the median is
  reported together with the number quantified and the matrix completeness.
* DE pAUC: protein-level Welch t-test A vs B; positives = changing species (YEAST, ECOLI), negatives = HUMAN.
  Standardised partial AUC over FPR in [0, 0.05] (area / 0.05) and TPR at empirical FDP <= 1% / 5% (FDP =
  human calls / all calls in the p-value ranking).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

SPECIES = ("HUMAN", "YEAST", "ECOLI")


def median_normalize(M: pd.DataFrame, ref_cols=None) -> tuple[pd.DataFrame, dict]:
    """Label-free run normalisation: shift each run so the median log2 ratio to the row-mean over precursors
    observed in all runs is 0 (no species information is used)."""
    L = np.log2(M)
    full = L.dropna()
    if len(full) < 50:
        full = L
    rowm = full.mean(axis=1)
    sh = {c: float(np.nanmedian(full[c] - rowm)) for c in M.columns}
    return M / pd.Series({c: 2 ** s for c, s in sh.items()}), sh


def epsilon(M: pd.DataFrame, species: pd.Series, A: list, B: list, expected: dict, min_obs: int = 1) -> dict:
    L = np.log2(M)
    nA, nB = L[A].notna().sum(axis=1), L[B].notna().sum(axis=1)
    keep = (nA >= min_obs) & (nB >= min_obs) & species.reindex(M.index).isin(expected).to_numpy()
    L = L[keep]
    sp = species.reindex(L.index)
    fc = L[A].mean(axis=1) - L[B].mean(axis=1)
    eps = fc - sp.map(expected)
    out = {"n": int(len(eps)), "median_abs_eps": float(np.median(np.abs(eps))) if len(eps) else np.nan,
           "mean_abs_eps": float(np.mean(np.abs(eps))) if len(eps) else np.nan, "per_species": {}}
    for s in expected:
        e = eps[sp == s]
        out["per_species"][s] = {"n": int(len(e)), "median_abs_eps": float(np.median(np.abs(e))) if len(e) else np.nan,
                                 "mean_abs_eps": float(np.mean(np.abs(e))) if len(e) else np.nan,
                                 "median_log2fc": float(np.median(fc[sp == s])) if len(e) else np.nan,
                                 "expected_log2fc": expected[s]}
    return out


def cv_stats(M: pd.DataFrame, A: list, B: list, min_obs: int = 3) -> dict:
    out = {}
    cvs = []
    for name, cols in (("A", A), ("B", B)):
        X = M[cols]
        ok = X.notna().sum(axis=1) >= min_obs
        X = X[ok]
        cv = X.std(axis=1, ddof=1) / X.mean(axis=1)
        out[name] = {"n": int(ok.sum()), "median_cv": float(np.median(cv)) if len(cv) else np.nan,
                     "frac_cv_lt_20pct": float((cv < 0.2).mean()) if len(cv) else np.nan}
        cvs.append(cv)
    allcv = pd.concat(cvs)
    out["median_cv"] = float(np.median(allcv)) if len(allcv) else np.nan
    return out


def completeness(M: pd.DataFrame, A: list, B: list) -> dict:
    X = M[A + B]
    X = X[X.notna().any(axis=1)]
    return {"n_any": int(len(X)), "n_complete": int(X.notna().all(axis=1).sum()),
            "n_ge1_per_condition": int(((X[A].notna().sum(axis=1) >= 1) & (X[B].notna().sum(axis=1) >= 1)).sum()),
            "n_3_per_condition": int(((X[A].notna().sum(axis=1) >= 3) & (X[B].notna().sum(axis=1) >= 3)).sum())
            if len(A) >= 3 else None,
            "completeness": float(X.notna().to_numpy().mean()) if len(X) else np.nan}


def de_pauc(P: pd.DataFrame, species: pd.Series, A: list, B: list, fpr_max: float = 0.05) -> dict:
    """P: protein x run linear intensities.  Proteins need >= 2 values in each condition."""
    L = np.log2(P)
    ok = (L[A].notna().sum(axis=1) >= 2) & (L[B].notna().sum(axis=1) >= 2)
    sp = species.reindex(L.index)
    ok &= sp.isin(SPECIES).to_numpy()
    L, sp = L[ok], sp[ok]
    t, p = stats.ttest_ind(L[A], L[B], axis=1, equal_var=False, nan_policy="omit")
    p = np.asarray(p, float)
    p = np.where(np.isfinite(p), p, 1.0)
    pos = (sp != "HUMAN").to_numpy()
    o = np.lexsort((np.arange(len(p)), p))
    tp = np.cumsum(pos[o])
    fp = np.cumsum(~pos[o])
    P_, N_ = max(pos.sum(), 1), max((~pos).sum(), 1)
    tpr, fpr = tp / P_, fp / N_
    fpr_c = np.concatenate([[0], fpr])
    tpr_c = np.concatenate([[0], tpr])
    m = fpr_c <= fpr_max
    xs = np.append(fpr_c[m], fpr_max)
    ys = np.append(tpr_c[m], tpr_c[m][-1])
    pauc = float(np.trapezoid(ys, xs) / fpr_max) if hasattr(np, "trapezoid") else float(np.trapz(ys, xs) / fpr_max)
    fdp = fp / np.maximum(tp + fp, 1)
    out = {"n_proteins": int(len(p)), "n_pos": int(pos.sum()), "n_neg": int((~pos).sum()),
           "pauc_fpr5pct_std": pauc, "auc": float(np.trapezoid(tpr_c, fpr_c)) if hasattr(np, "trapezoid") else None}
    for f in (0.01, 0.05):
        k = np.flatnonzero(fdp <= f)
        out[f"tpr_at_fdp{int(f * 100)}pct"] = float(tpr[k.max()]) if len(k) else 0.0
        out[f"tp_at_fdp{int(f * 100)}pct"] = int(tp[k.max()]) if len(k) else 0
    return out


def full_report(M: pd.DataFrame, prot_of: pd.Series, species_prec: pd.Series, protein_tables: dict,
                species_prot: pd.Series, A: list, B: list, expected: dict) -> dict:
    """M: precursor x run (normalised); protein_tables: name -> protein x run."""
    rep = {"precursor": {"completeness": completeness(M, A, B),
                         "epsilon_min1": epsilon(M, species_prec, A, B, expected, 1),
                         "epsilon_min3": epsilon(M, species_prec, A, B, expected, 3),
                         "cv_min3": cv_stats(M, A, B, 3)},
           "protein": {}}
    for name, P in protein_tables.items():
        rep["protein"][name] = {"completeness": completeness(P, A, B),
                                "epsilon_min2": epsilon(P, species_prot, A, B, expected, 2),
                                "cv_min3": cv_stats(P, A, B, 3),
                                "de": de_pauc(P, species_prot, A, B)}
    return rep
