"""Gate statistics for learned / predicted-feature scorers: DESIGN §6.4 (G1, G3, G4, G5) and Amendment A1 (§6.0:
tail-G2).  Pure numpy/scipy; inputs are per-null-spectrum arrays prepared by the caller.  Original PSM-level G2 and
G5's AUC are vulcan.stats.g2.{g2_differential, leak_auc} (unchanged).

Specifications fixed here before the A1 re-run (2026-10-27) -- where §6.4/A1 leave a detail open, the choice is the
one written below and is never revisited after results are seen:

tail-G2 (A1.2).  Unit = clean null spectrum (homolog-filtered, A1.3).  A spectrum is in the scorer's tail if its
  PSM-level winner (best of target or decoy under that scorer) scores >= the scorer's own 1% threshold in that run
  (lowest score of an accepted target at q <= 1%, production estimator).  n_tail = #tail spectra of the scorer
  under test.  Pass iff all of
    n_tail >= 200,
    one-sided 95% upper bound of the tail target fraction <= 0.60       (exact Clopper-Pearson),
    one-sided 95% upper bound of (frac_learned - frac_fixed) <= 0.10    (both on the same clean nulls; normal bound
        with the delta-method variance of the difference of two ratio estimators over the same spectra, so the
        overlap of the two tails is accounted for).
  n_tail < 200 -> "not-evaluable" (blocks FDR claims).
G1 (§6.4).  On null spectra: Spearman rho between the scorer's output and the sequence-only classifier logit of the
  winning candidate, restricted to the scorer's top score decile.  Pass iff the two-sided 99% CI (Fisher z,
  Bonett-Wright variance (1 + rho^2/2)/(n - 3)) lies inside (-0.05, 0.05).
G3 (§6.4).  Null spectra binned by the exhaustive best score (deciles of max(best target, best decoy)); per bin,
  recall of the exhaustive best target and of the exhaustive best decoy by the pruned (top-K) search.  Pass iff
  |recall_T - recall_D| <= 0.01 in every bin (point estimate; Newcombe CI reported).
G4 (§6.4).  Strata (instrument class, activation, charge 2/3/4+), partition fixed before estimation.  Decoys =
  rank-1 decoy-winning PSMs of dev/training runs.  Tolerance tau (fixed on dev data) = 95th percentile of the KS
  statistics between pairs of runs within the same stratum (runs with >= 200 decoys in the stratum).  Pass iff for
  every stratum KS(stratum decoys, decoys of all other strata) <= tau.
G5 (§6.4).  leak AUC <= 0.55 on null winners; reported with a 95% CI (Hanley-McNeil) and in the top score decile.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats

TAIL_N_MIN = 200
TAIL_MARGIN_FRAC = 0.60
TAIL_MARGIN_DIFF = 0.10
TAIL_CONF = 0.95


# ------------------------------------------------------------------------------------------------ tail-G2
def cp_upper(k: int, n: int, conf: float = TAIL_CONF) -> float:
    """One-sided exact (Clopper-Pearson) upper bound of a binomial proportion."""
    if n == 0:
        return 1.0
    if k >= n:
        return 1.0
    return float(stats.beta.ppf(conf, k + 1, n - k))


def tail_threshold(score, is_decoy, q, alpha: float = 0.01) -> float:
    """The method's own threshold: lowest score of an accepted target (q <= alpha); +inf if none."""
    acc = (np.asarray(q) <= alpha) & ~np.asarray(is_decoy, bool)
    return float(np.min(np.asarray(score, float)[acc])) if acc.any() else float("inf")


def tail_g2(tail_learned, target_learned, tail_fixed, target_fixed, n_min: int = TAIL_N_MIN,
            margin_frac: float = TAIL_MARGIN_FRAC, margin_diff: float = TAIL_MARGIN_DIFF,
            conf: float = TAIL_CONF) -> dict:
    """Tail-G2 on aligned clean null spectra (one row per spectrum).

    tail_* : bool, the spectrum's winner under that pipeline is at/above its own 1% threshold.
    target_* : bool, that winner is a target.
    """
    a = np.asarray(tail_learned, bool)
    b = a & np.asarray(target_learned, bool)
    c = np.asarray(tail_fixed, bool)
    d = c & np.asarray(target_fixed, bool)
    nL, kL, nF, kF = int(a.sum()), int(b.sum()), int(c.sum()), int(d.sum())
    N = a.size
    pL = kL / nL if nL else float("nan")
    pF = kF / nF if nF else float("nan")
    out = {"n_null_clean": int(N), "n_tail": nL, "n_tail_target": kL, "frac_target": pL,
           "frac_target_ucb": cp_upper(kL, nL, conf), "n_tail_fixed": nF, "n_tail_fixed_target": kF,
           "frac_target_fixed": pF, "n_tail_both": int((a & c).sum()), "conf_one_sided": conf,
           "margins": {"n_min": n_min, "frac": margin_frac, "diff": margin_diff}}
    if nL and nF:
        infl = (b - pL * a) / nL - (d - pF * c) / nF
        var = N / max(N - 1, 1) * float(np.sum(infl ** 2))
        diff = pL - pF
        out.update(diff=diff, diff_se=math.sqrt(var), diff_ucb=diff + stats.norm.ppf(conf) * math.sqrt(var))
    else:
        out.update(diff=float("nan"), diff_se=float("nan"), diff_ucb=float("inf"))
    if nL < n_min:
        out["verdict"] = "not-evaluable"
    else:
        ok = out["frac_target_ucb"] <= margin_frac and out["diff_ucb"] <= margin_diff
        out["verdict"] = "pass" if ok else "fail"
    out["frac_clause_pass"] = bool(out["frac_target_ucb"] <= margin_frac)
    out["diff_clause_pass"] = bool(out["diff_ucb"] <= margin_diff)
    return out


def simulate_tail_g2(n_rep: int, n_null: int, q_tail: float, p_fixed: float = 0.5, p_learned: float = 0.5,
                     overlap: float = 0.9, seed=11, **kw) -> dict:
    """Monte-Carlo verdict rates of :func:`tail_g2`.

    n_null clean nulls; each is in the fixed tail w.p. q_tail; the learned tail keeps a fixed-tail spectrum w.p.
    ``overlap`` and adds others so that its expected size is also q_tail * n_null.  Tail winners are targets w.p.
    p_fixed / p_learned; a spectrum in both tails keeps the same label (coupled) when p_learned == p_fixed, else the
    learned label is redrawn with prob. |p_learned - p_fixed| / (1 - min(p)) towards the biased side.
    """
    rng = np.random.default_rng(seed)
    v = {"pass": 0, "fail": 0, "not-evaluable": 0}
    nt = []
    for _ in range(n_rep):
        c = rng.random(n_null) < q_tail
        keep = c & (rng.random(n_null) < overlap)
        n_add = rng.binomial(n_null, q_tail * (1 - overlap))
        cand = np.flatnonzero(~c)
        add = rng.choice(cand, size=min(n_add, cand.size), replace=False)
        a = keep.copy()
        a[add] = True
        tf = rng.random(n_null) < p_fixed
        tl = tf.copy()
        if p_learned > p_fixed:
            flip = ~tf & (rng.random(n_null) < (p_learned - p_fixed) / (1 - p_fixed))
            tl[flip] = True
        elif p_learned < p_fixed:
            flip = tf & (rng.random(n_null) < (p_fixed - p_learned) / p_fixed)
            tl[flip] = False
        r = tail_g2(a, tl, c, tf, **kw)
        v[r["verdict"]] += 1
        nt.append(r["n_tail"])
    return {k: x / n_rep for k, x in v.items()} | {"mean_n_tail": float(np.mean(nt)), "n_rep": n_rep,
                                                   "p_fixed": p_fixed, "p_learned": p_learned}


# ------------------------------------------------------------------------------------------------ G1
def spearman_ci(x, y, conf: float = 0.99) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    if n < 10:
        return {"n": int(n), "rho": float("nan"), "ci": (float("nan"), float("nan"))}
    rho = float(stats.spearmanr(x, y).statistic)
    z = math.atanh(max(min(rho, 0.999999), -0.999999))
    se = math.sqrt((1 + rho * rho / 2) / (n - 3))
    h = stats.norm.ppf(0.5 + conf / 2) * se
    return {"n": int(n), "rho": rho, "ci": (math.tanh(z - h), math.tanh(z + h)), "conf": conf}


def g1_gate(score, seq_logit, top_frac: float = 0.1, bound: float = 0.05, conf: float = 0.99) -> dict:
    s = np.asarray(score, float)
    thr = np.quantile(s, 1 - top_frac)
    top = s >= thr
    r = spearman_ci(s[top], np.asarray(seq_logit, float)[top], conf)
    r["pass"] = bool(np.isfinite(r["rho"]) and r["ci"][0] > -bound and r["ci"][1] < bound)
    r["all_nulls"] = spearman_ci(s, seq_logit, conf)
    r["bound"] = bound
    return r


# ------------------------------------------------------------------------------------------------ G3
def newcombe_diff(k1, n1, k2, n2, conf=0.95):
    from vulcan.stats.g2 import wilson_ci
    p1, l1, u1 = wilson_ci(k1, n1, conf)
    p2, l2, u2 = wilson_ci(k2, n2, conf)
    d = p1 - p2
    lo = d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    hi = d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return float(d), (float(lo), float(hi))


def g3_gate(best_t_exh, best_d_exh, best_t_k, best_d_k, n_bins: int = 10, tol: float = 0.01) -> dict:
    """Arrays per null spectrum: exhaustive best target / decoy score (-inf if none) and the pruned search's."""
    bt, bd = np.asarray(best_t_exh, float), np.asarray(best_d_exh, float)
    kt, kd = np.asarray(best_t_k, float), np.asarray(best_d_k, float)
    best = np.maximum(bt, bd)
    fin = np.isfinite(best)
    edges = np.quantile(best[fin], np.linspace(0, 1, n_bins + 1))
    b = np.clip(np.searchsorted(edges, best, side="right") - 1, 0, n_bins - 1)
    bins, ok = [], True
    for i in range(n_bins):
        m = fin & (b == i)
        ht, hd = m & np.isfinite(bt), m & np.isfinite(bd)
        kT, nT = int((kt[ht] == bt[ht]).sum()), int(ht.sum())
        kD, nD = int((kd[hd] == bd[hd]).sum()), int(hd.sum())
        d, ci = newcombe_diff(kT, max(nT, 1), kD, max(nD, 1))
        rec = {"bin": i, "score_lo": float(edges[i]), "score_hi": float(edges[i + 1]), "n_target": nT, "n_decoy": nD,
               "recall_target": kT / max(nT, 1), "recall_decoy": kD / max(nD, 1), "diff": d, "diff_ci95": ci,
               "pass": bool(abs(d) <= tol)}
        ok &= rec["pass"]
        bins.append(rec)
    return {"bins": bins, "tol": tol, "pass": bool(ok), "n_null": int(fin.sum())}


# ------------------------------------------------------------------------------------------------ G4
def ks(a, b) -> float:
    return float(stats.ks_2samp(a, b).statistic)


def g4_gate(score, stratum, run, min_n: int = 200, pct: float = 95.0) -> dict:
    s = np.asarray(score, float)
    st = np.asarray(stratum)
    rn = np.asarray(run)
    within = []
    per = {}
    for u in np.unique(st):
        m = st == u
        runs = [r for r in np.unique(rn[m]) if (m & (rn == r)).sum() >= min_n]
        for i in range(len(runs)):
            for j in range(i + 1, len(runs)):
                within.append(ks(s[m & (rn == runs[i])], s[m & (rn == runs[j])]))
        per[str(u)] = {"n": int(m.sum()), "n_runs_used_for_tau": len(runs)}
    tau = float(np.percentile(within, pct)) if within else float("nan")
    ok = np.isfinite(tau)
    for u in np.unique(st):
        m = st == u
        k = ks(s[m], s[~m]) if m.sum() and (~m).sum() else float("nan")
        q = np.quantile(s[m], [0.5, 0.9, 0.99]).tolist() if m.sum() else []
        per[str(u)].update(ks_vs_rest=k, quantiles_50_90_99=q, pass_=bool(np.isfinite(k) and k <= tau))
        ok &= per[str(u)]["pass_"]
    return {"tau": tau, "tau_rule": f"{pct}th pct of within-stratum between-run KS (runs with >= {min_n} decoys)",
            "n_within_pairs": len(within), "within_ks_quantiles_10_50_90": (np.percentile(within, [10, 50, 90]).tolist()
                                                                          if within else []),
            "strata": per, "pass": bool(ok)}


# ------------------------------------------------------------------------------------------------ G5
def auc_ci(score, target_win, conf: float = 0.95) -> dict:
    from vulcan.stats.g2 import leak_auc
    w = np.asarray(target_win, bool)
    r = leak_auc(score, w)
    a = r["auc"]
    n1, n0 = int(w.sum()), int((~w).sum())
    if n1 and n0 and np.isfinite(a):
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        se = math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
        h = stats.norm.ppf(0.5 + conf / 2) * se
        r.update(se=se, ci=(a - h, a + h))
    r.update(n_target=n1, n_decoy=n0)
    return r


def g5_gate(score, target_win, top_frac: float = 0.1, alarm: float = 0.55) -> dict:
    s = np.asarray(score, float)
    w = np.asarray(target_win, bool)
    allr = auc_ci(s, w)
    top = s >= np.quantile(s, 1 - top_frac)
    topr = auc_ci(s[top], w[top])
    return {"all": allr, "top_decile": topr, "alarm_at": alarm, "pass": bool(allr["auc"] <= alarm)}
