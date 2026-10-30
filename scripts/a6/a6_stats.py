"""Amendment A6 (DESIGN §6.0e) gate statistics for v1 on HN1 (A6.6.1, A6.3.5, Appendix A).

Reuse first: every analytic component and the cluster bootstrap of tail-G2 are vulcan.gates.a2_stats /
vulcan.gates.stats / vulcan.stats.g2 called UNCHANGED; reviews/checks/a4_power.cap_factors is copied verbatim.
New code (A6.8.3): the outer envelope over project and cluster units for every bound-based row, the G1 jackknife-t
hulls, the TOST project/cluster bootstrap, the G2-original bootstraps, the cap subsample, the IUT precedence, the
stricter tail-G2 difference clause and the subset check.

Levels (A6.2 / A5.10.2, alpha/2 for the re-tested v1 rows):
  tail-G2   one-sided 0.975 (analytic CP / delta, project bootstrap, cluster bootstrap: UCB = max of the three)
  G1        two-sided 0.995 (Fisher-z Bonett-Wright, project jackknife-t, cluster jackknife-t: hull in z space)
  TOST (G3) two-sided 0.95  (Newcombe, project bootstrap, cluster bootstrap: hull)
  G2-original, G5: unchanged levels (G2-original one-sided 0.95 with margins 0.01 overall / 0.02 top decile;
            G5 alarm 0.55).  G2-original UCB = max(§6.4 UCB, project bootstrap, cluster bootstrap) (A6.6.1).
Validity controls unchanged: hyperscore tail Wilson two-sided 95% CI must contain 0.5; G1 hyperscore control
two-sided 99% CI must contain 0 (A6.6.1 text "99% CI ∋ 0"; computed with the same vulcan.gates.stats.g1_gate).

Implementation choices fixed here (A6.1.4: stricter reading where the text is open; recorded under A6.9):
  * Bootstrap: B = 10,000, numpy.random.default_rng(83).  ONE generator per row call: the cluster component is drawn
    first (inside a2_stats.tail_g2_a2 for tail-G2), then the project component from the same generator, unless a
    separate project seed is passed (used only by the acceptance test to reproduce scripts/c9/c9_power.tail_iut).
  * Degenerate resamples are assigned the failing value (as a2_stats): tail 1.0 / diff +1.0; TOST and G2-original
    difference +1.0 (upper) and -1.0 (lower) i.e. the resample is counted as outside the margin on both sides.
  * Fewer than 10 clusters (or fewer than 2 projects) -> every clustered row is not evaluable (A6.6.4).
  * Stricter tail difference clause: max(UCB_diff, UCB_abs - min(p_hat_hyperscore, 0.5)) <= 0.10 (A6.6.1).
  * G1 jackknife: the top decile is fixed on the full pool (vulcan.gates.stats.g1_gate rule); leave-one-unit-out
    Spearman on the remaining top spectra; z-space; se_jack = sqrt((J-1)/J * sum (z_-j - mean)^2); half-width
    t_{J-1}(0.5 + conf/2) * se_jack around the full-pool z (reviews/checks/a4_power.g1_a4_pass, extended to projects).
  * G3: grouping and per-group Newcombe exactly as a2_stats.g3_gate_a2 at conf 0.95 (A2.8 grouping, n_min 2,000 per
    label); each group's interval is hulled with the project and cluster percentile bootstrap intervals (two-sided
    0.95) of recall_T - recall_D (a4_power.tost_hull_pass rule, extended to projects).
  * G2-original: analytic = vulcan.stats.g2.g2_differential (single stratum, as scripts/a2/gates_a2.evaluate); the
    top decile is fixed on the full pool by the fixed winning score; unit bootstraps of the paired mean difference,
    95th percentile, overall and top decile; pass iff every UCB <= its margin.
  * Subset check (A6.6.1 / A5.10.3 item 7): per instrument class, cluster and project with n_tail >= 200, point
    estimates p_L <= 0.60 and max(p_L - p_F, p_L - min(p_F, 0.5)) <= 0.10 (the stricter clause; p_F undefined -> 0.5).
    An exceedance is a FAIL.
  * IUT precedence over pools (and over rows for the status): invalid > fail > not evaluable > pass.
No function here reads any file.  Inputs are per-spectrum arrays produced at step 5 by gates_hn1.py.
"""
from __future__ import annotations

import hashlib
import math

import numpy as np
from scipy import stats

from vulcan.gates import a2_stats as A2S
from vulcan.gates import stats as GS
from vulcan.stats import g2 as G2

# ------------------------------------------------------------------------------------------------ constants
TAIL_CONF = 0.975
G1_CONF = 0.995
G1_CONTROL_CONF = 0.99
TOST_CONF = 0.95
G2O_CONF = 0.95
G2O_MARGINS = (0.01, 0.02)
G5_ALARM = 0.55
FS_G5_MIN_WINNERS = 200
BOOT_B = 10_000
BOOT_SEED = 83
TAIL_N_MIN = 200
G1_N_MIN = 11_100
G3_N_MIN = 2_000
J_MIN = 10
CAP_SHARE = 0.10
CAP_HASH_SEED = 47
TAIL_MARGIN_FRAC = A2S.TAIL_MARGIN_FRAC     # 0.60 (A1, unchanged)
TAIL_MARGIN_DIFF = A2S.TAIL_MARGIN_DIFF     # 0.10
G1_BOUND = A2S.G1_BOUND                     # 0.05
G3_TOL = A2S.G3_TOL                         # 0.01

PASS, FAIL = "pass", "fail"
NOT_EVALUABLE = A2S.NOT_EVALUABLE           # "not-evaluable"
INVALID = A2S.INVALID                       # "not-evaluable (invalid pool)"
PRECEDENCE = (INVALID, FAIL, NOT_EVALUABLE, PASS)


def _gen(seed):
    return seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)


# ------------------------------------------------------------------------------------------------ IUT
def iut(*verdicts) -> str:
    """Intersection-union combination with precedence invalid > fail > not evaluable > pass (A6.3.5)."""
    vs = list(verdicts)
    for o in PRECEDENCE[:-1]:
        if o in vs:
            return o
    if not vs or any(v != PASS for v in vs):
        raise ValueError(f"unknown verdict in {vs}")
    return PASS


# ------------------------------------------------------------------------------------------------ cap (A6.3.5)
def cap_factors(c, share):
    """A4.3.9 cluster cap: per-cluster retention factors f_j in (0, 1] such that no cluster holds more than `share`
    of the capped total of c (here the label-collapsed tail proxy).  Applied as a seeded uniform subsample of that
    cluster's A spectra (and, with B's own shares, of its B spectra), so it is label-free (P1/P2)."""
    c = np.asarray(c, float)
    if c.sum() <= 0 or c.size * share < 1:
        return np.ones_like(c)
    x = c.copy()
    for _ in range(200):
        L = share * x.sum()
        if (x <= L * (1 + 1e-9)).all():
            break
        x = np.minimum(x, L)
    return np.where(c > 0, x / np.maximum(c, 1e-12), 1.0)


def cap_factors3(proxy, a_cand, a_all, share: float = CAP_SHARE):
    """Minimum of the three water-filling factor sets (tail proxy, A_cand, |A| per cluster).  None if J*share < 1
    (HN1 not evaluable, A6.3.5)."""
    proxy, a_cand, a_all = (np.asarray(x, float) for x in (proxy, a_cand, a_all))
    if proxy.size * share < 1:
        return None
    return np.minimum.reduce([cap_factors(proxy, share), cap_factors(a_cand, share), cap_factors(a_all, share)])


def cap_u(run: str, store_index: int, seed: int = CAP_HASH_SEED) -> float:
    """u(sha256(run | store index | 47)) in [0, 1): first 16 hex digits / 2**64."""
    h = hashlib.sha256(f"{run}|{int(store_index)}|{seed}".encode()).hexdigest()
    return int(h[:16], 16) / float(1 << 64)


def capped_membership(runs, store_idx, cluster, factors: dict, seed: int = CAP_HASH_SEED) -> np.ndarray:
    """bool per A spectrum: kept in the capped pool iff u(sha256(run|store index|47)) < f_{cluster}."""
    return np.array([cap_u(r, i, seed) < factors[c] for r, i, c in zip(runs, store_idx, cluster)], bool)


# ------------------------------------------------------------------------------------------------ unit helpers
def _units(x):
    u, inv = np.unique(np.asarray(x), return_inverse=True)
    return u, inv


def n_units(x) -> int:
    return int(np.unique(np.asarray(x)).size)


# ------------------------------------------------------------------------------------------------ tail-G2
def tail_g2_a6(project, cluster, tail_learned, target_learned, tail_fixed, target_fixed, B: int = BOOT_B,
               seed=BOOT_SEED, seed_project=None, conf: float = TAIL_CONF, n_min: int = TAIL_N_MIN,
               j_min: int = J_MIN) -> dict:
    """A6.6.1 tail-G2.  The analytic (CP / delta) bounds, the cluster bootstrap and the Wilson validity control are
    vulcan.gates.a2_stats.tail_g2_a2 called unchanged with cluster = the A6.3.3 component; the project bootstrap is
    a2_stats.cluster_bootstrap_tail with cluster = project.  UCB = max of the three."""
    g = _gen(seed)
    a2 = A2S.tail_g2_a2(cluster, tail_learned, target_learned, tail_fixed, target_fixed, B=B, seed=g, conf=conf,
                        n_min=n_min)
    gp = g if seed_project is None else _gen(seed_project)
    pb = A2S.cluster_bootstrap_tail(project, tail_learned, target_learned, tail_fixed, target_fixed, B=B, seed=gp,
                                    conf=conf)
    a1 = a2["a1_statistics"]
    ucb_abs = max(a2["ucb_abs"], pb["ucb_abs"])
    ucb_diff = max(a2["ucb_diff"], pb["ucb_diff"])
    pF = a1["frac_target_fixed"]
    pF_ref = 0.5 if not np.isfinite(pF) else min(pF, 0.5)
    diff_stat = max(ucb_diff, ucb_abs - pF_ref)
    clause_abs = bool(ucb_abs <= TAIL_MARGIN_FRAC)
    clause_diff = bool(diff_stat <= TAIL_MARGIN_DIFF)
    J, Jp = n_units(cluster), n_units(project)
    if a2["n_tail"] < n_min or J < j_min:
        verdict = NOT_EVALUABLE
    elif not a2["pool_validity"]["valid"]:
        verdict = INVALID
    else:
        verdict = PASS if (clause_abs and clause_diff) else FAIL
    return {"verdict": verdict, "n_tail": a2["n_tail"], "n_clusters": J, "n_projects": Jp,
            "ucb_abs": ucb_abs, "ucb_diff": ucb_diff, "diff_clause_stat": diff_stat,
            "ucb_abs_cp": a2["ucb_abs_cp"], "ucb_abs_cluster": a2["ucb_abs_cluster"], "ucb_abs_project": pb["ucb_abs"],
            "ucb_diff_delta": a2["ucb_diff_delta"], "ucb_diff_cluster": a2["ucb_diff_cluster"],
            "ucb_diff_project": pb["ucb_diff"], "p_hat_learned": a1["frac_target"], "p_hat_fixed": pF,
            "clause_abs_pass": clause_abs, "clause_diff_pass": clause_diff,
            "pool_validity": a2["pool_validity"], "a2_component": a2, "project_bootstrap": pb,
            "conf_one_sided": conf, "B": B}


# ------------------------------------------------------------------------------------------------ G1
def _spearman_z(x, y):
    r = float(stats.spearmanr(x, y).statistic)
    return r, math.atanh(max(min(r, 0.999999), -0.999999))


def jackknife_t_interval(x, y, unit, conf: float = G1_CONF) -> dict:
    """Leave-one-unit-out jackknife-t interval of the Spearman rho in Fisher-z space, centred at the full-sample z."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    u, inv = _units(unit)
    J = u.size
    r, z = _spearman_z(x, y)
    if J < 2:
        return {"J": int(J), "z": z, "lo_z": -math.inf, "hi_z": math.inf, "ci": (-1.0, 1.0)}
    zs = np.empty(J)
    for j in range(J):
        m = inv != j
        zs[j] = _spearman_z(x[m], y[m])[1] if m.sum() >= 10 else (math.inf if j % 2 else -math.inf)
    if not np.all(np.isfinite(zs)):
        return {"J": int(J), "z": z, "lo_z": -math.inf, "hi_z": math.inf, "ci": (-1.0, 1.0)}
    se = math.sqrt((J - 1) / J * float(((zs - zs.mean()) ** 2).sum()))
    h = float(stats.t.ppf(0.5 + conf / 2, J - 1)) * se
    return {"J": int(J), "z": z, "se_jack": se, "lo_z": z - h, "hi_z": z + h,
            "ci": (math.tanh(z - h), math.tanh(z + h))}


def g1_a6(project, cluster, score_learned, logit_learned, score_fixed, logit_fixed, n_min: int = G1_N_MIN,
          conf: float = G1_CONF, control_conf: float = G1_CONTROL_CONF, bound: float = G1_BOUND,
          j_min: int = J_MIN) -> dict:
    """A6.6.1 G1: Fisher-z (a2_stats.g1_gate_a2 at conf 0.995, unchanged) hulled with project and cluster
    jackknife-t; pass iff the hull is inside (-0.05, 0.05); n_top >= n_min; hyperscore control 99% CI contains 0."""
    a2 = A2S.g1_gate_a2(score_learned, logit_learned, score_fixed, logit_fixed, n_min=n_min, bound=bound, conf=conf)
    ctl = GS.g1_gate(score_fixed, logit_fixed, bound=bound, conf=control_conf)
    clo, chi = ctl["ci"]
    ctl_valid = bool(ctl["n"] >= 10 and np.isfinite(clo) and clo <= 0.0 <= chi)
    s = np.asarray(score_learned, float)
    top = s >= np.quantile(s, 0.9)                      # the vulcan.gates.stats.g1_gate rule
    x, y = s[top], np.asarray(logit_learned, float)[top]
    pr, cl = np.asarray(project)[top], np.asarray(cluster)[top]
    out = {"n_top": a2["n_top"], "rho": a2["rho"], "ci_fisher": a2["ci"], "a2_component": a2,
           "control_hyperscore": {"n_top": int(ctl["n"]), "rho": ctl["rho"], "ci": ctl["ci"], "conf": control_conf,
                                  "contains_0": ctl_valid}, "conf": conf, "bound": bound}
    J = n_units(cluster)
    if a2["n_top"] < n_min or J < j_min or x.size < 10:
        out.update(verdict=NOT_EVALUABLE, n_clusters=J)
        return out
    jp = jackknife_t_interval(x, y, pr, conf)
    jc = jackknife_t_interval(x, y, cl, conf)
    z = math.atanh(max(min(a2["rho"], 0.999999), -0.999999))
    zf_lo, zf_hi = math.atanh(a2["ci"][0]), math.atanh(a2["ci"][1])
    lo = min(zf_lo, jp["lo_z"], jc["lo_z"])
    hi = max(zf_hi, jp["hi_z"], jc["hi_z"])
    hull = (math.tanh(lo), math.tanh(hi))
    inside = bool(hull[0] > -bound and hull[1] < bound)
    verdict = INVALID if not ctl_valid else (PASS if inside else FAIL)
    out.update(verdict=verdict, z=z, ci_hull=hull, jack_project=jp, jack_cluster=jc, n_clusters=J,
               n_projects=n_units(project), hull_inside_bound=inside)
    return out


# ------------------------------------------------------------------------------------------------ bootstraps
def unit_bootstrap_diff(unit, k1, n1, k2, n2, B: int, gen, conf: float, two_sided: bool = True) -> dict:
    """Percentile bootstrap over units of (sum k1 / sum n1) - (sum k2 / sum n2); per-unit count arrays computed
    from per-spectrum indicator arrays.  Degenerate resample (an empty denominator) -> +1 for the upper / -1 for the
    lower quantile (counted as outside on both sides)."""
    u, inv = _units(unit)
    J = u.size
    sums = [np.bincount(inv, weights=np.asarray(a, float), minlength=J) for a in (k1, n1, k2, n2)]
    if J < 2:
        return {"J": int(J), "lo": -math.inf, "hi": math.inf}
    idx = gen.integers(0, J, size=(B, J))
    K1, N1, K2, N2 = (s[idx].sum(1) for s in sums)
    ok = (N1 > 0) & (N2 > 0)
    d = K1 / np.maximum(N1, 1) - K2 / np.maximum(N2, 1)
    dhi = np.where(ok, d, 1.0)
    dlo = np.where(ok, d, -1.0)
    if two_sided:
        a = (1 - conf) / 2
        return {"J": int(J), "lo": float(np.quantile(dlo, a)), "hi": float(np.quantile(dhi, 1 - a)),
                "n_degenerate": int((~ok).sum())}
    return {"J": int(J), "lo": -math.inf, "hi": float(np.quantile(dhi, conf)), "n_degenerate": int((~ok).sum())}


# ------------------------------------------------------------------------------------------------ G3
def g3_a6(project, cluster, best_t_exh, best_d_exh, best_t_k, best_d_k, B: int = BOOT_B, seed=BOOT_SEED,
          conf: float = TOST_CONF, n_min: int = G3_N_MIN, tol: float = G3_TOL, j_min: int = J_MIN) -> dict:
    """A6.6.1 G3: a2_stats.g3_gate_a2 at conf 0.95 (grouping + Newcombe, unchanged), each group's interval hulled
    with the project and cluster bootstrap (two-sided 0.95) percentile intervals."""
    a2 = A2S.g3_gate_a2(best_t_exh, best_d_exh, best_t_k, best_d_k, n_min=n_min, tol=tol, conf=conf)
    J = n_units(cluster)
    out = {"a2_component": a2, "n_clusters": J, "conf": conf, "tol": tol}
    if a2["verdict"] == NOT_EVALUABLE or J < j_min:
        out["verdict"] = NOT_EVALUABLE
        return out
    g = _gen(seed)
    bt, bd = np.asarray(best_t_exh, float), np.asarray(best_d_exh, float)
    kt, kd = np.asarray(best_t_k, float), np.asarray(best_d_k, float)
    best = np.maximum(bt, bd)
    fin = np.isfinite(best)
    edges = np.asarray(a2["edges"])
    b = np.clip(np.searchsorted(edges, best, side="right") - 1, 0, len(edges) - 2)
    ht, hd = fin & np.isfinite(bt), fin & np.isfinite(bd)
    groups, ok = [], True
    for r in a2["groups"]:
        m = np.isin(b, r["deciles"]) & fin
        mt, md = m & ht, m & hd
        kT = mt & (kt == bt)
        kD = md & (kd == bd)
        cb = unit_bootstrap_diff(np.asarray(cluster), kT, mt, kD, md, B, g, conf)
        pb = unit_bootstrap_diff(np.asarray(project), kT, mt, kD, md, B, g, conf)
        lo = min(r["ci90"][0], cb["lo"], pb["lo"])
        hi = max(r["ci90"][1], cb["hi"], pb["hi"])
        p = bool(r["n_target"] > 0 and r["n_decoy"] > 0 and lo >= -tol and hi <= tol)
        ok &= p
        groups.append(r | {"ci_newcombe": r["ci90"], "ci_cluster": (cb["lo"], cb["hi"]),
                           "ci_project": (pb["lo"], pb["hi"]), "ci_hull": (lo, hi), "pass_hull": p})
    out.update(verdict=PASS if ok else FAIL, groups=groups)
    return out


# ------------------------------------------------------------------------------------------------ G2-original, G5
def g2_original_a6(project, cluster, win_fixed, win_learned, score_fixed, B: int = BOOT_B, seed=BOOT_SEED,
                   conf: float = G2O_CONF, margins=G2O_MARGINS, top_frac: float = 0.1, j_min: int = J_MIN) -> dict:
    """A6.6.1 G2-original: vulcan.stats.g2.g2_differential (unchanged, single stratum) with each UCB replaced by
    max(analytic, project bootstrap, cluster bootstrap)."""
    a = np.asarray(win_fixed, bool)
    lw = np.asarray(win_learned, bool)
    s = np.asarray(score_fixed, float)
    an = G2.g2_differential(a, lw, s, margin_overall=margins[0], margin_top=margins[1], top_frac=top_frac, conf=conf)
    st = an["strata"][0]
    J = n_units(cluster)
    out = {"analytic": an, "n_clusters": J, "conf_one_sided": conf, "margins": list(margins)}
    if a.size == 0 or J < j_min:
        out["verdict"] = NOT_EVALUABLE
        return out
    g = _gen(seed)
    top = s >= np.quantile(s, 1 - top_frac)
    ones = np.ones(a.size, bool)
    res = {}
    for name, m, marg, an_ucb in (("overall", ones, margins[0], st["overall"]["ucb"]),
                                  ("top", top, margins[1], st["top"]["ucb"])):
        cb = unit_bootstrap_diff(np.asarray(cluster), lw & m, m, a & m, m, B, g, conf, two_sided=False)
        pb = unit_bootstrap_diff(np.asarray(project), lw & m, m, a & m, m, B, g, conf, two_sided=False)
        ucb = max(an_ucb, cb["hi"], pb["hi"])
        res[name] = {"ucb_analytic": an_ucb, "ucb_cluster": cb["hi"], "ucb_project": pb["hi"], "ucb": ucb,
                     "margin": marg, "pass": bool(ucb <= marg)}
    out.update(res, verdict=PASS if (res["overall"]["pass"] and res["top"]["pass"]) else FAIL)
    return out


def g5_a6(score, win_t, alarm: float = G5_ALARM, min_winners: int | None = None) -> dict:
    """G5 unchanged (vulcan.gates.stats.g5_gate); FS: evaluable iff >= 200 target and >= 200 decoy winners."""
    w = np.asarray(win_t, bool)
    nT, nD = int(w.sum()), int((~w).sum())
    if (min_winners is not None and (nT < min_winners or nD < min_winners)) or nT == 0 or nD == 0:
        return {"verdict": NOT_EVALUABLE, "nT": nT, "nD": nD}
    r = GS.g5_gate(np.asarray(score, float), w, alarm=alarm)
    return {"verdict": PASS if r["pass"] else FAIL, "nT": nT, "nD": nD, "detail": r}


# ------------------------------------------------------------------------------------------------ subset check
def subset_check(groups: dict, tail_learned, target_learned, tail_fixed, target_fixed,
                 n_min: int = TAIL_N_MIN) -> dict:
    """groups: {kind: label array} (e.g. instrument class, cluster, project).  For each subset with n_tail >= 200:
    p_L <= 0.60 and max(p_L - p_F, p_L - min(p_F, 0.5)) <= 0.10 (point estimates); an exceedance is a FAIL."""
    tl, gl = np.asarray(tail_learned, bool), np.asarray(target_learned, bool)
    tf, gf = np.asarray(tail_fixed, bool), np.asarray(target_fixed, bool)
    rows, fails = {}, []
    for kind, lab in groups.items():
        lab = np.asarray(lab)
        for u in np.unique(lab):
            m = lab == u
            nL = int((tl & m).sum())
            if nL < n_min:
                continue
            pL = float((tl & gl & m).sum()) / nL
            nF = int((tf & m).sum())
            pF = float((tf & gf & m).sum()) / nF if nF else float("nan")
            ref = 0.5 if not np.isfinite(pF) else min(pF, 0.5)
            d = max(pL - pF if np.isfinite(pF) else -math.inf, pL - ref)
            ok = bool(pL <= TAIL_MARGIN_FRAC and d <= TAIL_MARGIN_DIFF)
            key = f"{kind}:{u}"
            rows[key] = {"n_tail": nL, "p_learned": pL, "p_fixed": pF, "diff_stat": d, "pass": ok}
            if not ok:
                fails.append(key)
    return {"verdict": FAIL if fails else PASS, "fails": fails, "subsets": rows}


# ------------------------------------------------------------------------------------------------ pool + v1
def evaluate_pool(D: dict, m, B: int = BOOT_B, seed: int = BOOT_SEED, n_min: dict | None = None) -> dict:
    """All v1 rows on one pool (capped or uncapped).  D: per-null-spectrum arrays aligned over A and FS:
      'pool' ('A'/'FS'), 'project', 'cluster', 'instrument';
      'fixed' and 'v1': dicts with 'tail1' (bool), 'win_t' (bool), 'score', 'logit';
      'g3': dict with 'bt','bd','kt','kd' (per spectrum; -inf where none; meaningful on A).
    m: bool mask (pool membership).  Each row uses a fresh default_rng(seed) (one generator per row call)."""
    nm = {"tail": TAIL_N_MIN, "g1": G1_N_MIN, "g3": G3_N_MIN} | (n_min or {})
    pool = np.asarray(D["pool"])
    pr, cl = np.asarray(D["project"]), np.asarray(D["cluster"])
    F, L = D["fixed"], D["v1"]
    m = np.asarray(m, bool)
    mA, mFS = m & (pool == "A"), m & (pool == "FS")

    def tail(mm):
        return tail_g2_a6(pr[mm], cl[mm], L["tail1"][mm], L["win_t"][mm], F["tail1"][mm], F["win_t"][mm], B=B,
                          seed=seed, n_min=nm["tail"])

    def g1(mm):
        return g1_a6(pr[mm], cl[mm], L["score"][mm], L["logit"][mm], F["score"][mm], F["logit"][mm],
                     n_min=nm["g1"])

    r = {"tail_G2": tail(mA), "G1": g1(mA)}
    G = D["g3"]
    r["G3"] = g3_a6(pr[mA], cl[mA], G["bt"][mA], G["bd"][mA], G["kt"][mA], G["kd"][mA], B=B, seed=seed,
                    n_min=nm["g3"])
    r["G2_original"] = g2_original_a6(pr[mA], cl[mA], F["win_t"][mA], L["win_t"][mA], F["score"][mA], B=B, seed=seed)
    r["G5"] = g5_a6(L["score"][mA], L["win_t"][mA])
    r["FS"] = {"tail_G2": tail(mFS) if mFS.any() else {"verdict": NOT_EVALUABLE},
               "G1": g1(mFS) if mFS.any() else {"verdict": NOT_EVALUABLE},
               "G5": g5_a6(L["score"][mFS], L["win_t"][mFS], min_winners=FS_G5_MIN_WINNERS)}
    groups = {"instrument": np.asarray(D["instrument"])[mA], "cluster": cl[mA], "project": pr[mA]}
    r["subset_check"] = subset_check(groups, L["tail1"][mA], L["win_t"][mA], F["tail1"][mA], F["win_t"][mA],
                                     n_min=nm["tail"])
    return r


def fs_outcome(fs: dict) -> str:
    """FS rows (A2.4 rule 3): an evaluable fail blocks (-> fail); not evaluable / invalid does not block (-> pass)."""
    return FAIL if any(v.get("verdict") == FAIL for v in fs.values()) else PASS


def pool_row_verdicts(r: dict) -> dict:
    return {"tail_G2": r["tail_G2"]["verdict"], "G1": r["G1"]["verdict"], "G3": r["G3"]["verdict"],
            "G2_original": r["G2_original"]["verdict"], "G5": r["G5"]["verdict"], "FS": fs_outcome(r["FS"]),
            "subset_check": r["subset_check"]["verdict"]}


def v1_verdict(res_capped: dict, res_uncapped: dict) -> dict:
    """v1 verdict (A6.6.1): every row on both pools (IUT per row), then the overall status with the same precedence.
    Returns the per-row IUT verdicts and the a6_status enum value (PASS / FAILED / not evaluable / invalid pool)."""
    vc, vu = pool_row_verdicts(res_capped), pool_row_verdicts(res_uncapped)
    rows = {k: iut(vc[k], vu[k]) for k in vc}
    o = iut(*rows.values())
    status = {PASS: "PASS", FAIL: "FAILED", NOT_EVALUABLE: "not evaluable", INVALID: "invalid pool"}[o]
    return {"rows": rows, "rows_capped": vc, "rows_uncapped": vu, "overall": o, "status": status}
