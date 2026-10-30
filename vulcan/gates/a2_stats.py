"""Amendment A2 gate statistics (DESIGN §6.0b = reviews/amendment_A2_draft.md): tail-G2 with the project-cluster
bootstrap max rule and the hyperscore pool-validity control (A2.5), G1 with n_min and the hyperscore negative control
(A2.6), G3 per-bin TOST with n_min and bottom-half merging (A2.8).

The A1 functions in vulcan.gates.stats are NOT changed and are reused where A2 keeps A1's rule (CP bound, delta-method
difference bound, Spearman/Fisher-z CI, Newcombe interval).  A1 margins are imported, never redefined (P3).

Implementation choices fixed here before DN2 is listed (A2.11 item 5; stricter reading where the text is open; each
is also listed in scripts/a2/IMPLEMENTATION_NOTES.md and config/a2.yaml):

tail-G2 (A2.5)
  * unit / tail membership / thresholds exactly as A1 (vulcan.gates.stats docstring).
  * cluster = DN2 project (accession).  Bootstrap: B = 10,000 resamples of the J projects with replacement,
    numpy.random.default_rng(83); one (B, J) index matrix drawn once and used for BOTH the absolute and the
    difference statistic (paired: a project brings its learned AND hyperscore tails).
  * absolute statistic per resample = sum(k_learned)/sum(n_learned) over the resampled projects' tails; a resample
    with no learned-tail spectrum is assigned 1.0 (stricter than dropping it).  UCB_cluster = the 95th percentile
    (numpy 'linear' quantile, i.e. type 7).  UCB_abs = max(exact CP one-sided 95%, UCB_cluster).
  * difference per resample = p_learned - p_fixed on the resampled projects; a resample in which either tail is
    empty is assigned +1.0 (stricter).  UCB_diff = max(A1 delta-method bound, 95th percentile).
  * J < 2 projects: the cluster bound is undefined -> +inf (stricter: the gate cannot pass).
  * pool validity: two-sided 95% Wilson CI of the hyperscore tail target fraction must contain 0.5 (closed
    interval); an empty hyperscore tail is invalid (stricter).  Invalid -> "not-evaluable" (invalid pool) whatever
    the clauses say (A2.1: the control "can only turn a verdict into not evaluable"); the clause results are kept
    in the output for reporting.
  * verdict precedence: n_tail < 200 -> not-evaluable; else invalid pool -> not-evaluable (invalid pool);
    else pass iff UCB_abs <= 0.60 and UCB_diff <= 0.10, else fail.
G1 (A2.6)
  * A1 statistic unchanged (vulcan.gates.stats.g1_gate: top score decile = score >= 90% quantile; Spearman with the
    frozen classifier logit of the winning candidate; two-sided 99% Fisher-z CI with Bonett-Wright variance inside
    (-0.05, 0.05)).
  * n_top = the number of spectra in the learned scorer's top decile.  n_top < 7,200 -> not-evaluable.
  * negative control: the same statistic for hyperscore (its own top decile of its own winners on the same pool,
    logit of ITS winning candidate); its 99% CI must contain 0 (closed), else not-evaluable (invalid pool).  If the
    hyperscore top decile itself has < 10 spectra the control is undefined -> invalid (stricter).
G3 (A2.8)
  * bins = deciles of the exhaustive best score max(best target, best decoy) over pool spectra with a finite
    exhaustive best (same edges / searchsorted rule as A1).
  * counts per bin: n_T = spectra with a finite exhaustive best target, n_D = with a finite exhaustive best decoy.
  * top five deciles (bins 5..9): each needs n_T >= 2000 AND n_D >= 2000, else not-evaluable.
  * bottom half (bins 0..4), merge order: walk upward from bin 0, accumulating consecutive deciles into a group
    until the group has n_T >= 2000 AND n_D >= 2000, then close it and start the next group.  If the last group
    (the one containing bin 4) is still short, there is no higher bottom-half decile to merge with and merging may
    not cross into the top half -> not-evaluable (stricter: never merge downward, never test an under-sized group).
  * per group: Delta = recall_T - recall_D, 90% two-sided Newcombe hybrid-score interval (A1's newcombe_diff with
    conf 0.90); pass iff lo >= -0.01 and hi <= +0.01 in EVERY group.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats

from vulcan.gates import stats as GS
from vulcan.stats.g2 import wilson_ci

# A1 margins (unchanged, P3) -- imported, not redefined
TAIL_N_MIN = GS.TAIL_N_MIN
TAIL_MARGIN_FRAC = GS.TAIL_MARGIN_FRAC
TAIL_MARGIN_DIFF = GS.TAIL_MARGIN_DIFF
TAIL_CONF = GS.TAIL_CONF
G1_BOUND = 0.05
G1_CONF = 0.99
G3_TOL = 0.01
# A2 additions
BOOT_B = 10_000
BOOT_SEED = 83            # SEEDS[3]
VALIDITY_CONF = 0.95
G1_N_MIN = 7_200
G3_N_MIN = 2_000
G3_CONF = 0.90            # two one-sided 95% bounds (TOST)
NOT_EVALUABLE = "not-evaluable"
INVALID = "not-evaluable (invalid pool)"


# ------------------------------------------------------------------------------------------------ tail-G2
def _project_sums(cluster, *bools):
    """Per-project sums of boolean arrays.  Projects are sorted (np.unique) so the result is order-free."""
    cl = np.asarray(cluster)
    u, inv = np.unique(cl, return_inverse=True)
    return u, [np.bincount(inv, weights=np.asarray(b, bool).astype(float), minlength=u.size) for b in bools]


def cluster_bootstrap_tail(cluster, tail_learned, target_learned, tail_fixed, target_fixed, B: int = BOOT_B,
                           seed: int | np.random.Generator = BOOT_SEED, conf: float = TAIL_CONF) -> dict:
    """Project-cluster percentile bootstrap of the learned tail target fraction and of the paired difference
    (learned - fixed).  See the module docstring for the degenerate-resample rules."""
    a = np.asarray(tail_learned, bool)
    c = np.asarray(tail_fixed, bool)
    u, (nL, kL, nF, kF) = _project_sums(cluster, a, a & np.asarray(target_learned, bool), c,
                                        c & np.asarray(target_fixed, bool))
    J = u.size
    out = {"n_projects": int(J), "B": int(B), "conf": conf}
    if J < 2:
        out.update(ucb_abs=float("inf"), ucb_diff=float("inf"), note="J < 2: cluster bound undefined")
        return out
    g = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    idx = g.integers(0, J, size=(B, J))
    sL, sKL, sF, sKF = nL[idx].sum(1), kL[idx].sum(1), nF[idx].sum(1), kF[idx].sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pL = np.where(sL > 0, sKL / np.maximum(sL, 1), 1.0)
        pF = sKF / np.maximum(sF, 1)
        d = np.where((sL > 0) & (sF > 0), pL - pF, 1.0)
    out.update(ucb_abs=float(np.quantile(pL, conf)), ucb_diff=float(np.quantile(d, conf)),
               n_degenerate_abs=int((sL == 0).sum()), n_degenerate_diff=int(((sL == 0) | (sF == 0)).sum()))
    return out


def pool_validity(tail_fixed, target_fixed, conf: float = VALIDITY_CONF) -> dict:
    """A2.5 clause 4: hyperscore tail target fraction, two-sided Wilson CI must contain 0.5."""
    c = np.asarray(tail_fixed, bool)
    n = int(c.sum())
    k = int((c & np.asarray(target_fixed, bool)).sum())
    if n == 0:
        return {"n": 0, "k": 0, "ci": (float("nan"), float("nan")), "valid": False, "note": "empty hyperscore tail"}
    p, lo, hi = wilson_ci(k, n, conf)
    lo, hi = float(lo), float(hi)
    return {"n": n, "k": k, "frac": float(p), "ci": (lo, hi), "conf": conf, "valid": bool(lo <= 0.5 <= hi)}


def tail_g2_a2(cluster, tail_learned, target_learned, tail_fixed, target_fixed, B: int = BOOT_B,
               seed: int | np.random.Generator = BOOT_SEED, n_min: int = TAIL_N_MIN,
               margin_frac: float = TAIL_MARGIN_FRAC, margin_diff: float = TAIL_MARGIN_DIFF,
               conf: float = TAIL_CONF, validity_conf: float = VALIDITY_CONF) -> dict:
    """A2.5 tail-G2 on aligned stratum-A null spectra (one row per spectrum; ``cluster`` = project accession)."""
    a1 = GS.tail_g2(tail_learned, target_learned, tail_fixed, target_fixed, n_min=n_min, margin_frac=margin_frac,
                    margin_diff=margin_diff, conf=conf)
    cb = cluster_bootstrap_tail(cluster, tail_learned, target_learned, tail_fixed, target_fixed, B=B, seed=seed,
                                conf=conf)
    val = pool_validity(tail_fixed, target_fixed, validity_conf)
    ucb_abs = max(a1["frac_target_ucb"], cb["ucb_abs"])
    ucb_diff = max(a1["diff_ucb"], cb["ucb_diff"])
    clause_abs = bool(ucb_abs <= margin_frac)
    clause_diff = bool(ucb_diff <= margin_diff)
    if a1["n_tail"] < n_min:
        verdict = NOT_EVALUABLE
    elif not val["valid"]:
        verdict = INVALID
    else:
        verdict = "pass" if (clause_abs and clause_diff) else "fail"
    return {"verdict": verdict, "n_tail": a1["n_tail"], "ucb_abs": ucb_abs, "ucb_diff": ucb_diff,
            "ucb_abs_cp": a1["frac_target_ucb"], "ucb_abs_cluster": cb["ucb_abs"],
            "ucb_diff_delta": a1["diff_ucb"], "ucb_diff_cluster": cb["ucb_diff"],
            "clause_abs_pass": clause_abs, "clause_diff_pass": clause_diff, "pool_validity": val,
            "a1_statistics": a1, "cluster_bootstrap": cb,
            "margins": {"n_min": n_min, "frac": margin_frac, "diff": margin_diff, "conf_one_sided": conf}}


# ------------------------------------------------------------------------------------------------ G1
def g1_gate_a2(score_learned, logit_learned, score_fixed, logit_fixed, n_min: int = G1_N_MIN,
               top_frac: float = 0.1, bound: float = G1_BOUND, conf: float = G1_CONF) -> dict:
    """A2.6 G1 on stratum-A null winners.  *_learned: learned score / classifier logit of the learned winner;
    *_fixed: hyperscore / logit of the hyperscore winner, same spectra."""
    r = GS.g1_gate(score_learned, logit_learned, top_frac=top_frac, bound=bound, conf=conf)
    ctl = GS.g1_gate(score_fixed, logit_fixed, top_frac=top_frac, bound=bound, conf=conf)
    n_top = int(r["n"])
    lo, hi = ctl["ci"]
    ctl_valid = bool(ctl["n"] >= 10 and np.isfinite(lo) and lo <= 0.0 <= hi)
    if n_top < n_min:
        verdict = NOT_EVALUABLE
    elif not ctl_valid:
        verdict = INVALID
    else:
        verdict = "pass" if r["pass"] else "fail"
    return {"verdict": verdict, "n_top": n_top, "n_min": n_min, "rho": r["rho"], "ci": r["ci"],
            "ci_inside_bound": bool(r["pass"]), "control_hyperscore": {"n_top": int(ctl["n"]), "rho": ctl["rho"],
                                                                        "ci": ctl["ci"], "contains_0": ctl_valid},
            "bound": bound, "conf": conf, "learned_detail": r}


def fisher_pooled_rho(groups, conf: float = G1_CONF) -> dict:
    """Reported-only (A2.6 item 5): Fisher-z pooled Spearman rho over groups [(x, y), ...] (e.g. within target winners
    and within decoy winners), weights n - 3, Bonett-Wright variance per group."""
    zs, ws, vs = [], [], []
    for x, y in groups:
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        if x.size < 10:
            continue
        rho = float(stats.spearmanr(x, y).statistic)
        zs.append(math.atanh(max(min(rho, 0.999999), -0.999999)))
        ws.append(x.size - 3)
        vs.append((1 + rho * rho / 2) / (x.size - 3))
    if not zs:
        return {"rho": float("nan"), "ci": (float("nan"), float("nan"))}
    w = np.asarray(ws, float)
    z = float(np.sum(w * np.asarray(zs)) / w.sum())
    se = math.sqrt(float(np.sum(w * w * np.asarray(vs))) / w.sum() ** 2)
    h = stats.norm.ppf(0.5 + conf / 2) * se
    return {"rho": math.tanh(z), "ci": (math.tanh(z - h), math.tanh(z + h)), "n_groups": len(zs), "conf": conf}


# ------------------------------------------------------------------------------------------------ G3
def tost_bin(kT: int, nT: int, kD: int, nD: int, tol: float = G3_TOL, conf: float = G3_CONF) -> dict:
    d, (lo, hi) = GS.newcombe_diff(kT, max(nT, 1), kD, max(nD, 1), conf)
    return {"n_target": int(nT), "n_decoy": int(nD), "recall_target": kT / max(nT, 1), "recall_decoy": kD / max(nD, 1),
            "diff": d, "ci90": (lo, hi), "pass": bool(nT > 0 and nD > 0 and lo >= -tol and hi <= tol)}


def g3_groups(nT_bin, nD_bin, n_min: int = G3_N_MIN, n_bins: int = 10):
    """A2.8 bin grouping.  Returns (groups as lists of decile indices, None) or (None, reason)."""
    half = n_bins // 2
    for b in range(half, n_bins):
        if nT_bin[b] < n_min or nD_bin[b] < n_min:
            return None, f"top-half decile {b}: n_T {int(nT_bin[b])}, n_D {int(nD_bin[b])} < {n_min}"
    groups, cur, cT, cD = [], [], 0, 0
    for b in range(half):
        cur.append(b)
        cT += int(nT_bin[b])
        cD += int(nD_bin[b])
        if cT >= n_min and cD >= n_min:
            groups.append(cur)
            cur, cT, cD = [], 0, 0
    if cur:
        return None, (f"bottom-half deciles {cur} reach only n_T {cT}, n_D {cD} < {n_min} and may not merge into "
                      f"the top half")
    return groups + [[b] for b in range(half, n_bins)], None


def g3_gate_a2(best_t_exh, best_d_exh, best_t_k, best_d_k, n_bins: int = 10, n_min: int = G3_N_MIN,
               tol: float = G3_TOL, conf: float = G3_CONF) -> dict:
    """A2.8 G3 on one pool.  Arrays per pool spectrum: exhaustive best target / decoy score (-inf if none) and the
    pruned search's best target / decoy score."""
    bt, bd = np.asarray(best_t_exh, float), np.asarray(best_d_exh, float)
    kt, kd = np.asarray(best_t_k, float), np.asarray(best_d_k, float)
    best = np.maximum(bt, bd)
    fin = np.isfinite(best)
    if fin.sum() < n_bins:
        return {"verdict": NOT_EVALUABLE, "reason": "fewer pool spectra than bins", "n_null": int(fin.sum())}
    edges = np.quantile(best[fin], np.linspace(0, 1, n_bins + 1))
    b = np.clip(np.searchsorted(edges, best, side="right") - 1, 0, n_bins - 1)
    ht, hd = fin & np.isfinite(bt), fin & np.isfinite(bd)
    nT = np.array([int((ht & (b == i)).sum()) for i in range(n_bins)])
    nD = np.array([int((hd & (b == i)).sum()) for i in range(n_bins)])
    groups, why = g3_groups(nT, nD, n_min, n_bins)
    base = {"n_null": int(fin.sum()), "edges": edges.tolist(), "n_target_by_decile": nT.tolist(),
            "n_decoy_by_decile": nD.tolist(), "n_min": n_min, "tol": tol, "conf": conf}
    if groups is None:
        return base | {"verdict": NOT_EVALUABLE, "reason": why}
    res, ok = [], True
    for g in groups:
        m = np.isin(b, g)
        mt, md = m & ht, m & hd
        r = tost_bin(int((kt[mt] == bt[mt]).sum()), int(mt.sum()), int((kd[md] == bd[md]).sum()), int(md.sum()),
                     tol, conf)
        r["deciles"] = g
        ok &= r["pass"]
        res.append(r)
    return base | {"verdict": "pass" if ok else "fail", "groups": res}


# ------------------------------------------------------------------------------------------------ MC helpers
def simulate_tail_g2_a2(n_rep: int, n_tail_mean: float, J: int = 10, p_fixed: float = 0.5, p_learned: float = 0.5,
                        between_sd: float = 0.0, B: int = 2000, seed: int = 131, validity: bool = True) -> dict:
    """Verdict rates of :func:`tail_g2_a2` on synthetic projects (independent learned and fixed tails, the
    conservative case for the difference clause; Dirichlet(3) project sizes; Beta between-project heterogeneity of
    the learned fraction).  validity=False reports the rate ignoring the pool-validity control (for comparison with
    the A1-margin table S1 of reviews/checks/a2_power)."""
    g = np.random.default_rng(seed)
    v = {"pass": 0, "fail": 0, NOT_EVALUABLE: 0, INVALID: 0}
    for _ in range(n_rep):
        w = g.dirichlet(np.full(J, 3.0))
        nL = g.multinomial(g.poisson(n_tail_mean), w)
        nF = g.multinomial(g.poisson(n_tail_mean), w)
        if between_sd > 0:
            s = p_learned * (1 - p_learned) / between_sd ** 2 - 1
            pl = g.beta(p_learned * s, (1 - p_learned) * s, size=J)
        else:
            pl = np.full(J, p_learned)
        kL = g.binomial(nL, pl)
        kF = g.binomial(nF, p_fixed)
        cl, tl, gl, tf, gf = [], [], [], [], []
        for j in range(J):
            m = nL[j] + nF[j]
            cl.append(np.full(m, j))
            tl.append(np.r_[np.ones(nL[j], bool), np.zeros(nF[j], bool)])
            gl.append(np.r_[np.arange(nL[j]) < kL[j], np.zeros(nF[j], bool)])
            tf.append(np.r_[np.zeros(nL[j], bool), np.ones(nF[j], bool)])
            gf.append(np.r_[np.zeros(nL[j], bool), np.arange(nF[j]) < kF[j]])
        r = tail_g2_a2(np.concatenate(cl), np.concatenate(tl), np.concatenate(gl), np.concatenate(tf),
                       np.concatenate(gf), B=B, seed=g)
        verdict = r["verdict"]
        if not validity and verdict == INVALID:
            verdict = "pass" if (r["clause_abs_pass"] and r["clause_diff_pass"]) else "fail"
        v[verdict] += 1
    return {k: x / n_rep for k, x in v.items()} | {"n_rep": n_rep, "n_tail_mean": n_tail_mean, "J": J,
                                                   "p_fixed": p_fixed, "p_learned": p_learned}


def g1_pass_rate_mc(n_top: int, rho: float, n_rep: int, seed: int = 131, bound: float = G1_BOUND,
                    conf: float = G1_CONF) -> float:
    """MC pass rate of the A2.6 CI rule at a given n_top and true correlation (bivariate normal; Spearman of normal
    data has Pearson-to-Spearman attenuation (6/pi) asin(r/2), so the latent r is chosen to give Spearman rho)."""
    g = np.random.default_rng(seed)
    r_lat = 2 * math.sin(math.pi * rho / 6)
    ok = 0
    for _ in range(n_rep):
        x = g.standard_normal(n_top)
        y = r_lat * x + math.sqrt(1 - r_lat * r_lat) * g.standard_normal(n_top)
        c = GS.spearman_ci(x, y, conf)
        ok += bool(c["ci"][0] > -bound and c["ci"][1] < bound)
    return ok / n_rep
