"""G2 null test (DESIGN.md §6.4, v0.3: *differential* null test).

Pass condition (per instrument class): on the *same* null competitions, the learned
pipeline's null target-win rate may exceed the fixed-score pipeline's by at most
``margin_overall`` = 0.01 overall and ``margin_top`` = 0.02 within the top decile
of the fixed-score winning score.

Test implemented (paired, one-sided non-inferiority):

    d     = mean_i (B_i - A_i)          A_i, B_i = 1 if the target won competition i
                                         under the fixed (A) / learned (B) pipeline
    se    = sd(B - A) / sqrt(n)
    UCB   = d + z_conf * se             (z_conf = one-sided normal quantile)
    pass  = UCB <= margin

Default ``conf = 0.95`` (one-sided). This is the level implied by the v0.3 budget
(~110k null competitions per class, ~11k in the top decile at margin 0.02, 90%
power, worst-case independence of A and B); DESIGN does not state it explicitly
(see results/2026-10-26_stats_sims/design_issues.md).

Diagnostics (not gated): absolute target-win fraction with Wilson CI (the v0.2
equivalence band [0.49, 0.51]) overall and in the top decile, and a two-sample KS
test between the score distributions of null target winners and null decoy
winners.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats

__all__ = ["wilson_ci", "g2_differential", "g2_absolute", "ks_winners", "g2_power", "g2_sample_size", "leak_auc",
           "simulate_g2"]


def wilson_ci(k, n, conf: float = 0.99):
    z = stats.norm.ppf(0.5 + conf / 2)
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    p = k / np.maximum(n, 1)
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, c - h, c + h


def _top_mask(score, top_frac):
    s = np.asarray(score, float)
    thr = np.quantile(s, 1 - top_frac)
    return s >= thr


def _diff_test(a, b, margin, z):
    d = b.astype(np.float64) - a.astype(np.float64)
    n = d.size
    if n < 2:
        return {"n": int(n), "diff": float("nan"), "ucb": float("inf"), "pass": False, "margin": margin}
    m = float(d.mean())
    se = float(d.std(ddof=1) / math.sqrt(n))
    ucb = m + z * se
    return {"n": int(n), "diff": m, "se": se, "ucb": ucb, "margin": margin, "pass": bool(ucb <= margin),
            "p_fixed": float(a.mean()), "p_learned": float(b.mean())}


def g2_differential(win_fixed, win_learned, score_fixed, stratum=None, margin_overall: float = 0.01,
                    margin_top: float = 0.02, top_frac: float = 0.1, conf: float = 0.95) -> dict:
    """Differential G2 per stratum (instrument class) and overall verdict.

    Parameters
    ----------
    win_fixed, win_learned : bool per null competition (target won).
    score_fixed : fixed-pipeline winning score (defines the top decile).
    stratum : optional instrument-class label per competition.
    """
    a = np.asarray(win_fixed, bool)
    b = np.asarray(win_learned, bool)
    s = np.asarray(score_fixed, float)
    st = np.zeros(a.size, dtype=np.int64) if stratum is None else np.asarray(stratum)
    z = stats.norm.ppf(conf)
    out = {"conf_one_sided": conf, "strata": {}}
    ok = True
    for u in np.unique(st):
        m = st == u
        top = _top_mask(s[m], top_frac)
        r_all = _diff_test(a[m], b[m], margin_overall, z)
        r_top = _diff_test(a[m][top], b[m][top], margin_top, z)
        passed = r_all["pass"] and r_top["pass"]
        ok &= passed
        out["strata"][u.item() if hasattr(u, "item") else u] = {"overall": r_all, "top": r_top, "pass": passed}
    out["pass"] = bool(ok)
    out["fires"] = not out["pass"]
    return out


def g2_absolute(target_win, score, conf: float = 0.99, band=(0.49, 0.51), top_frac: float = 0.1) -> dict:
    """v0.2 absolute equivalence test, kept as a diagnostic only."""
    w = np.asarray(target_win, bool)
    top = _top_mask(score, top_frac)
    res = {}
    for name, m in (("overall", np.ones_like(w)), ("top", top)):
        k, n = int(w[m].sum()), int(m.sum())
        p, lo, hi = wilson_ci(k, n, conf)
        res[name] = {"n": n, "p": float(p), "ci": (float(lo), float(hi)),
                     "within_band": bool(lo >= band[0] and hi <= band[1])}
    res["pass"] = res["overall"]["within_band"] and res["top"]["within_band"]
    return res


def ks_winners(score, target_win) -> dict:
    """KS test: null target-winner vs null decoy-winner score distributions."""
    s = np.asarray(score, float)
    w = np.asarray(target_win, bool)
    if w.all() or (~w).all():
        return {"statistic": float("nan"), "pvalue": float("nan")}
    r = stats.ks_2samp(s[w], s[~w])
    return {"statistic": float(r.statistic), "pvalue": float(r.pvalue)}


def _pass_prob(n, margin, true_diff, var, z):
    se = math.sqrt(var / n)
    # P(dhat + z se <= margin), dhat ~ N(true_diff, se^2)
    return float(stats.norm.cdf((margin - z * se - true_diff) / se))


def _var_diff(p_fixed, p_learned, rho):
    """Var(B - A) for the coupling used in :func:`simulate_g2`."""
    p11 = rho * min(p_fixed, p_learned) + (1 - rho) * p_fixed * p_learned
    e = p_learned - p_fixed
    e2 = p_fixed + p_learned - 2 * p11
    return e2 - e * e


def g2_power(n_total: int, true_diff: float = 0.0, p_fixed: float = 0.5, rho: float = 0.0,
             margin_overall: float = 0.01, margin_top: float = 0.02, top_frac: float = 0.1,
             conf: float = 0.95) -> dict:
    """Normal-approximation probability that G2 *passes* (both clauses).

    The joint probability is reported as the product of the two clause
    probabilities, a lower bound because the clauses are positively dependent.
    ``rho`` = coupling between fixed and learned outcomes (0 = independent, worst
    case).
    """
    z = stats.norm.ppf(conf)
    var = _var_diff(p_fixed, p_fixed + true_diff, rho)
    pa = _pass_prob(n_total, margin_overall, true_diff, var, z)
    pt = _pass_prob(max(int(n_total * top_frac), 2), margin_top, true_diff, var, z)
    return {"p_pass_overall": pa, "p_pass_top": pt, "p_pass_joint_lb": pa * pt}


def g2_sample_size(power: float = 0.9, p_fixed: float = 0.5, rho: float = 0.0, margin_overall: float = 0.01,
                   margin_top: float = 0.02, top_frac: float = 0.1, conf: float = 0.95) -> int:
    """Smallest total number of null competitions per stratum such that a fair
    learned scorer (true difference 0) passes G2 with probability >= ``power``
    (joint, lower bound)."""
    lo, hi = 100, 100
    while g2_power(hi, 0.0, p_fixed, rho, margin_overall, margin_top, top_frac, conf)["p_pass_joint_lb"] < power:
        hi *= 2
        if hi > 10 ** 10:
            raise ValueError("unattainable")
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if g2_power(mid, 0.0, p_fixed, rho, margin_overall, margin_top, top_frac, conf)["p_pass_joint_lb"] >= power:
            hi = mid
        else:
            lo = mid
    return hi


def simulate_g2(n_total: int, n_rep: int, p_fixed: float = 0.5, p_learned: float = 0.5, rho: float = 0.0,
                top_frac: float = 0.1, seed=11, conf: float = 0.95, margin_overall: float = 0.01,
                margin_top: float = 0.02) -> dict:
    """Monte-Carlo pass rate of the differential G2 (exact multinomial sufficient
    statistics: the paired outcomes enter the test only through the counts of
    (A,B) in {00,01,10,11}, simulated separately for the top decile and the rest).

    Coupling: U ~ U(0,1), A = 1[U < p_fixed]; with probability rho the learned
    pipeline reuses U, else a fresh V: B = 1[U or V < p_learned].
    """
    rng = np.random.default_rng(seed)
    p11 = rho * min(p_fixed, p_learned) + (1 - rho) * p_fixed * p_learned
    p10 = p_fixed - p11
    p01 = p_learned - p11
    p00 = 1 - p11 - p10 - p01
    probs = np.array([p00, p01, p10, p11])
    n_top = int(round(n_total * top_frac))
    top = rng.multinomial(n_top, probs, size=n_rep)
    rest = rng.multinomial(n_total - n_top, probs, size=n_rep)
    z = stats.norm.ppf(conf)

    def ucb(c, n):
        n01, n10 = c[:, 1], c[:, 2]
        m = (n01 - n10) / n
        ex2 = (n01 + n10) / n
        var = (ex2 - m * m) * n / (n - 1)
        return m + z * np.sqrt(var / n)

    ok_all = ucb(top + rest, n_total) <= margin_overall
    ok_top = ucb(top, n_top) <= margin_top
    passed = ok_all & ok_top
    return {"n_total": n_total, "n_rep": n_rep, "p_fixed": p_fixed, "p_learned": p_learned, "rho": rho,
            "pass_rate": float(passed.mean()), "fire_rate": float(1 - passed.mean()),
            "pass_rate_overall_clause": float(ok_all.mean()), "pass_rate_top_clause": float(ok_top.mean())}


def leak_auc(score, target_win) -> dict:
    """G5 leak alarm (DESIGN v0.3 §6.4): AUC separating null target winners from
    null decoy winners by score, computed on G2 nulls; alarm if > 0.55."""
    from scipy.stats import rankdata

    s = np.asarray(score, float)
    w = np.asarray(target_win, bool)
    n1, n0 = int(w.sum()), int((~w).sum())
    if n1 == 0 or n0 == 0:
        return {"auc": float("nan"), "alarm": False}
    r = rankdata(s)
    auc = (r[w].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    return {"auc": float(auc), "alarm": bool(auc > 0.55)}
