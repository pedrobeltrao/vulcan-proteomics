"""A6.8.3 acceptance test of scripts/a6/a6_stats.py (SYNTHETIC inputs only; no data file is read).

(i)   a6_stats reproduces vulcan.gates.a2_stats / vulcan.gates.stats / vulcan.stats.g2 EXACTLY where it calls them
      (tail-G2 analytic + cluster bootstrap + validity; G1 Fisher-z; G3 grouping + Newcombe; G2-original analytic).
(ii)  a6_stats.cap_factors == reviews/checks/a4_power.cap_factors (bitwise) and the min-of-three cap == c9_power.capped3;
      on the c9_power.tail_iut synthetic states (clean, margin tau, concentrated leak, residual contamination ->
      invalid pool) with the same seeds, the a2 component of a6 tail-G2 reproduces c9_power's per-pool verdicts
      (c9_power uses tail_g2_a2 with cluster bootstrap only) and the IUT precedence reproduces c9_power's combined
      verdict; capped/uncapped n_tail identical.  c9_power has no project bootstrap and no stricter difference
      clause: for those, the test asserts a6 is never more lenient (a6 pass => c9 pass) and that the bounds of the
      shared components agree to 1e-12 against a direct a2_stats call on the same arrays and generator state.
(iii) one-sided coverage >= 0.95 of the NEW interval components at their nominal levels, 1,000 synthetic replicates
      each (bootstrap B = 1,000 in the coverage replicates for speed; production B = 10,000):
        G1 cluster / project jackknife-t (two-sided 0.995), TOST cluster / project bootstrap (two-sided 0.95),
        G2-original cluster / project bootstrap UCB (one-sided 0.95).
      Truth: true parameter 0; two scenarios per component: 'hom' (no between-unit heterogeneity) and 'het'
      (unit-level random effects with mean 0: G1 rho_j ~ N(0, 0.03), TOST gap_j ~ N(0, 0.01), G2 d_j ~ N(0, 0.01));
      J = 37 clusters (A6.10 J p50) / 74 projects, log-normal unit sizes (log-sd 1.0, as a4_power.cluster_weights).
      PASS criterion (A6.8.3, read strictly): the component's own empirical one-sided coverage >= 0.95 on each side,
      in each scenario.  'within_2se_of_0.95' and the envelopes actually used in the gate are reported, not used.
Results are written to results/2026-10-28_A6_hn1/a6_stats_acceptance.json by the last test.
    conda_envs/pe/bin/python -m pytest tests/a6/test_a6_stats.py -q
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):   # shared node: do not oversubscribe
    os.environ.setdefault(_v, "4")

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
for p in (ROOT, os.path.join(ROOT, "scripts", "a6"), os.path.join(ROOT, "scripts", "c9"),
          os.path.join(ROOT, "reviews", "checks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import a6_stats as S  # noqa: E402
from vulcan.gates import a2_stats as A2S  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402
from vulcan.stats import g2 as G2  # noqa: E402

OUT = os.path.join(ROOT, "results", "2026-10-28_A6_hn1", "a6_stats_acceptance.json")
RES: dict = {"items": {}, "t0": time.time()}
N_COV = int(os.environ.get("A6_COV_REPS", "1000"))
B_COV = 1000
J_COV = 37
COV_MIN = 0.95


def _rec(key, **kw):
    RES["items"].setdefault(key, {}).update(kw)


# ------------------------------------------------------------------------------------------------ synthetic states
def tail_state(kind, seed, J=30, n=1600):
    """Per-spectrum tail arrays with J clusters (2 projects each)."""
    g = np.random.default_rng(seed)
    w = g.lognormal(0, 1.0, J)
    w /= w.sum()
    nL = g.multinomial(n, w)
    nF = g.multinomial(n, w)
    pf = 0.5
    if kind == "clean":
        pl = np.full(J, 0.5)
    elif kind == "margin":
        pl = np.clip(g.normal(0.60, 0.10, J), 0, 1)
    elif kind == "leak":
        pl = np.full(J, 0.5)
        pl[g.permutation(J)[:max(1, round(0.25 * J))]] = 0.9
    elif kind == "invalid":
        pf, pl = 0.58, np.full(J, 0.62)
    cl, pr, tl, gl, tf, gf = [], [], [], [], [], []
    for j in range(J):
        kL, kF = g.binomial(nL[j], pl[j]), g.binomial(nF[j], pf)
        m = nL[j] + nF[j]
        cl.append(np.full(m, j))
        pr.append(np.where(np.arange(m) % 2 == 0, 2 * j, 2 * j + 1))
        tl.append(np.r_[np.ones(nL[j], bool), np.zeros(nF[j], bool)])
        gl.append(np.r_[np.arange(nL[j]) < kL, np.zeros(nF[j], bool)])
        tf.append(np.r_[np.zeros(nL[j], bool), np.ones(nF[j], bool)])
        gf.append(np.r_[np.zeros(nL[j], bool), np.arange(nF[j]) < kF])
    return tuple(np.concatenate(x) for x in (pr, cl, tl, gl, tf, gf))


STATES = ("clean", "margin", "leak", "invalid")


# ------------------------------------------------------------------------------------------------ (i)
@pytest.mark.parametrize("kind", STATES)
def test_i_tail_reproduces_a2_stats(kind):
    pr, cl, tl, gl, tf, gf = tail_state(kind, 7)
    r = S.tail_g2_a6(pr, cl, tl, gl, tf, gf, B=2000, seed=83)
    g = np.random.default_rng(83)
    a2 = A2S.tail_g2_a2(cl, tl, gl, tf, gf, B=2000, seed=g, conf=0.975)
    pb = A2S.cluster_bootstrap_tail(pr, tl, gl, tf, gf, B=2000, seed=g, conf=0.975)
    for k in ("ucb_abs_cp", "ucb_abs_cluster", "ucb_diff_delta", "ucb_diff_cluster", "n_tail"):
        assert r[k] == a2[k], k
    assert r["pool_validity"] == a2["pool_validity"]
    assert r["ucb_abs_project"] == pb["ucb_abs"] and r["ucb_diff_project"] == pb["ucb_diff"]
    assert r["ucb_abs"] == max(a2["ucb_abs"], pb["ucb_abs"])
    # never more lenient than A2 at the same level
    if r["verdict"] == S.PASS:
        assert a2["verdict"] == "pass"
    if a2["verdict"] == S.INVALID:
        assert r["verdict"] == S.INVALID
    _rec("i_exact_a2_stats", **{f"tail_{kind}": {"a6": r["verdict"], "a2_component": a2["verdict"],
                                                 "ucb_abs": r["ucb_abs"], "diff_stat": r["diff_clause_stat"]}})


def test_i_state_verdicts_sane():
    v = {k: S.tail_g2_a6(*tail_state(k, 7), B=2000)["verdict"] for k in STATES}
    assert v["clean"] == S.PASS
    assert v["margin"] == S.FAIL
    assert v["leak"] == S.FAIL
    assert v["invalid"] == S.INVALID
    _rec("i_exact_a2_stats", state_verdicts=v)


def _g1_data(seed, n=60_000, J=30, rho=0.0, tau=0.0):
    g = np.random.default_rng(seed)
    w = g.lognormal(0, 1, J)
    cl = g.choice(J, size=n, p=w / w.sum())
    pr = 2 * cl + (g.random(n) < 0.5)
    rj = np.clip(rho + g.normal(0, tau, J), -0.9, 0.9) if tau else np.full(J, rho)
    x = g.standard_normal(n)
    y = rj[cl] * x + np.sqrt(1 - rj[cl] ** 2) * g.standard_normal(n)
    xf = g.standard_normal(n)
    yf = g.standard_normal(n)
    return pr, cl, x, y, xf, yf


def test_i_g1_reproduces_a2_stats():
    pr, cl, x, y, xf, yf = _g1_data(3, n=120_000)
    r = S.g1_a6(pr, cl, x, y, xf, yf)
    a2 = A2S.g1_gate_a2(x, y, xf, yf, n_min=11_100, conf=0.995)
    assert r["rho"] == a2["rho"] and tuple(r["ci_fisher"]) == tuple(a2["ci"]) and r["n_top"] == a2["n_top"]
    ctl = GS.g1_gate(xf, yf, conf=0.99)
    assert r["control_hyperscore"]["ci"] == ctl["ci"]
    lo, hi = r["ci_hull"]
    assert lo <= a2["ci"][0] + 1e-15 and hi >= a2["ci"][1] - 1e-15
    assert r["verdict"] in (S.PASS, S.FAIL)            # clean but J = 30 log-normal clusters: the hull may be wide
    assert (r["verdict"] == S.PASS) == (lo > -0.05 and hi < 0.05)
    r2 = S.g1_a6(*_g1_data(4, n=120_000, rho=0.08)[:4], xf, yf)
    assert r2["verdict"] == S.FAIL
    r3 = S.g1_a6(pr[:50_000], cl[:50_000], x[:50_000], y[:50_000], xf[:50_000], yf[:50_000])
    assert r3["verdict"] == S.NOT_EVALUABLE            # n_top 5,000 < 11,100
    _rec("i_exact_a2_stats", g1={"clean": r["verdict"], "rho_0.08": r2["verdict"], "small": r3["verdict"],
                                 "hull": r["ci_hull"], "fisher": a2["ci"]})


def _g3_data(seed, n=60_000, J=30, gap=0.0):
    g = np.random.default_rng(seed)
    w = g.lognormal(0, 1, J)
    cl = g.choice(J, size=n, p=w / w.sum())
    pr = 2 * cl + (g.random(n) < 0.5)
    bt = g.normal(10, 3, n)
    bd = g.normal(10, 3, n)
    kt = np.where(g.random(n) < 0.97 - gap, bt, bt - 1)
    kd = np.where(g.random(n) < 0.97, bd, bd - 1)
    return pr, cl, bt, bd, kt, kd


def test_i_g3_reproduces_a2_stats():
    pr, cl, bt, bd, kt, kd = _g3_data(5, n=80_000)
    r = S.g3_a6(pr, cl, bt, bd, kt, kd, B=1000)
    a2 = A2S.g3_gate_a2(bt, bd, kt, kd, n_min=2000, conf=0.95)
    assert [g["ci90"] for g in r["groups"]] == [g["ci90"] for g in a2["groups"]]
    assert all(g["ci_hull"][0] <= g["ci90"][0] and g["ci_hull"][1] >= g["ci90"][1] for g in r["groups"])
    assert r["verdict"] in (S.PASS, S.FAIL)
    r2 = S.g3_a6(*_g3_data(6, n=80_000, gap=0.05), B=1000)
    assert r2["verdict"] == S.FAIL
    r3 = S.g3_a6(*_g3_data(7, n=10_000), B=200)
    assert r3["verdict"] == S.NOT_EVALUABLE
    _rec("i_exact_a2_stats", g3={"clean": r["verdict"], "gap_0.05": r2["verdict"], "small": r3["verdict"]})


def test_i_g2_original_reproduces_g2():
    g = np.random.default_rng(9)
    n, J = 300_000, 30                                   # the §6.4 budget is ~110k per class; top decile 30k
    cl = g.integers(0, J, n)
    pr = 2 * cl + (g.random(n) < 0.5)
    a, b, s = g.random(n) < 0.5, g.random(n) < 0.5, g.normal(size=n)
    r = S.g2_original_a6(pr, cl, a, b, s, B=1000)
    an = G2.g2_differential(a, b, s)
    assert r["analytic"]["strata"][0]["overall"]["ucb"] == an["strata"][0]["overall"]["ucb"]
    assert r["overall"]["ucb"] >= an["strata"][0]["overall"]["ucb"]
    assert r["verdict"] == S.PASS
    b2 = np.where(g.random(n) < 0.04, True, b)
    assert S.g2_original_a6(pr, cl, a, b2, s, B=1000)["verdict"] == S.FAIL
    _rec("i_exact_a2_stats", g2_original={"clean": r["verdict"], "overall_ucb": r["overall"]["ucb"]})


# ------------------------------------------------------------------------------------------------ (ii)
def test_ii_cap_factors_verbatim():
    import a4_power as AP
    import c9_power as C9
    g = np.random.default_rng(1)
    n_eq = 0
    for t in range(500):
        J = int(g.integers(1, 60))
        c = g.lognormal(0, 1.5, J) * (g.random(J) > 0.1)
        assert np.array_equal(S.cap_factors(c, 0.10), AP.cap_factors(c, 0.10))
        ca, cA = g.lognormal(5, 1, J), g.lognormal(7, 1, J)
        mine = S.cap_factors3(c, ca, cA)
        ref = C9.capped3(c, ca, cA)
        assert (mine is None) == (ref is None)
        if ref is not None:
            assert np.array_equal(mine, ref[3])
            n_eq += 1
    u = [S.cap_u("run_x", i) for i in range(10_000)]
    assert 0 <= min(u) and max(u) < 1 and abs(np.mean(u) - 0.5) < 0.02
    _rec("ii_c9_a4_reproduction", cap_factors={"n_cases": 500, "n_min_of_three_cases": n_eq, "bitwise_equal": True})


def _iut_c9(vu, vc):
    for o in (A2S.INVALID, "fail", A2S.NOT_EVALUABLE):
        if o in (vu, vc):
            return o
    return "pass"


def test_ii_iut_precedence():
    vs = (S.PASS, S.FAIL, S.NOT_EVALUABLE, S.INVALID)
    for a in vs:
        for b in vs:
            assert S.iut(a, b) == _iut_c9(a, b)
    _rec("ii_c9_a4_reproduction", iut_precedence={"n_combinations": 16, "equal": True})


def _snap(seed, J=36):
    g = np.random.default_rng(seed)
    cu = g.lognormal(0, 1.0, J)
    cu = np.round(cu / cu.sum() * 1400)
    cc = cu * S.cap_factors(cu, 0.10)
    return {"cp_uncapped": cu.tolist(), "cp_capped": cc.tolist()}


C9_SCEN = {"clean": ("clean", None), "margin_tau": ("tau", (0.60, 0.10)), "leak25": ("concentrated", (0.25, 0.90)),
           "contam_0.05": ("contam", 0.05)}


def _c9_arrays(snap, rng, kappa, scen):
    """Exactly c9_power.tail_iut's generation (same calls, same order) up to the gate calls."""
    cu = np.asarray(snap["cp_uncapped"], float)
    cc = np.asarray(snap["cp_capped"], float)
    J = cu.size
    kind, par = scen
    p_fixed = 0.5
    if kind == "clean":
        pl = np.full(J, 0.5)
    elif kind == "tau":
        mu, tau = par
        s = mu * (1 - mu) / tau ** 2 - 1
        pl = rng.beta(mu * s, (1 - mu) * s, size=J)
    elif kind == "concentrated":
        frac, phi = par
        pl = np.full(J, 0.5)
        pl[rng.permutation(J)[:max(1, int(round(frac * J)))]] = phi
    elif kind == "contam":
        p_fixed = 0.5 + par
        pl = np.full(J, 0.5 + 1.5 * par)
    nL = rng.poisson(kappa * cu)
    nF = rng.poisson(kappa * cu)
    kL = rng.binomial(nL, pl)
    kF = rng.binomial(nF, p_fixed)
    f = np.where(cu > 0, cc / np.maximum(cu, 1e-12), 1.0)
    nLc = rng.binomial(nL, f)
    kLc = rng.hypergeometric(np.maximum(kL, 0), np.maximum(nL - kL, 0), nLc)
    nFc = rng.binomial(nF, f)
    kFc = rng.hypergeometric(np.maximum(kF, 0), np.maximum(nF - kF, 0), nFc)

    def arrays(nl, kl, nf, kf):
        cl = np.concatenate([np.full(nl[j] + nf[j], j) for j in range(J)])
        tl = np.concatenate([np.r_[np.ones(nl[j], bool), np.zeros(nf[j], bool)] for j in range(J)])
        gl = np.concatenate([np.r_[np.arange(nl[j]) < kl[j], np.zeros(nf[j], bool)] for j in range(J)])
        tf = np.concatenate([np.r_[np.zeros(nl[j], bool), np.ones(nf[j], bool)] for j in range(J)])
        gf = np.concatenate([np.r_[np.zeros(nl[j], bool), np.arange(nf[j]) < kf[j]] for j in range(J)])
        return cl, tl, gl, tf, gf
    return arrays(nL, kL, nF, kF), arrays(nLc, kLc, nFc, kFc), int(nLc.sum())


def test_ii_c9_tail_iut_reproduction():
    import c9_power as C9
    B = 400
    rows, n_agree_final, n_more_lenient = {}, 0, 0
    for name, scen in C9_SCEN.items():
        for kappa in (1.0, 0.5):
            for seed in range(6):
                snap = _snap(1000 + seed)
                ref_v, ref_n = C9.tail_iut(snap, np.random.default_rng(seed), kappa, scen, B=B)
                rng = np.random.default_rng(seed)
                (u, c, n_c) = _c9_arrays(snap, rng, kappa, scen)
                # a6 with project == cluster, cluster component drawn from the shared rng (as c9), project bootstrap
                # from a separate generator so the shared stream is consumed exactly as in c9_power.
                ru = S.tail_g2_a6(u[0], u[0], *u[1:], B=B, seed=rng, seed_project=83)
                st = rng.bit_generator.state
                rc = S.tail_g2_a6(c[0], c[0], *c[1:], B=B, seed=rng, seed_project=83)
                # direct a2_stats call on the capped arrays from the same generator state: bounds to 1e-12
                g2 = np.random.default_rng()
                g2.bit_generator.state = st
                a2c = A2S.tail_g2_a2(c[0], *c[1:], B=B, seed=g2, conf=0.975)
                for k in ("ucb_abs_cp", "ucb_abs_cluster", "ucb_diff_delta", "ucb_diff_cluster"):
                    assert abs(rc[k] - a2c[k]) <= 1e-12 or rc[k] == a2c[k], k
                comp = S.iut(ru["a2_component"]["verdict"], rc["a2_component"]["verdict"])
                assert comp == ref_v, (name, kappa, seed, comp, ref_v)
                assert n_c == ref_n
                fin = S.iut(ru["verdict"], rc["verdict"])
                n_agree_final += fin == ref_v
                if fin == S.PASS and ref_v != "pass":
                    n_more_lenient += 1
                rows[f"{name}|k{kappa}|s{seed}"] = {"c9": ref_v, "a6_a2_component_iut": comp, "a6_final": fin}
    assert n_more_lenient == 0
    _rec("ii_c9_a4_reproduction", c9_tail_iut={"n_states": len(rows), "a2_component_iut_equal_c9": True,
                                               "n_capped_tail_equal": True, "bounds_tol": 1e-12,
                                               "a6_final_equal_c9": int(n_agree_final),
                                               "a6_final_more_lenient_than_c9": n_more_lenient,
                                               "B": B, "states": rows,
                                               "note": "c9_power.tail_iut has no project bootstrap and no stricter "
                                                       "difference clause; a6's final verdict may only be stricter"})


# ------------------------------------------------------------------------------------------------ (iii)
def _weights(g, J):
    w = g.lognormal(0, 1.0, J)
    return w / w.sum()


def _cov_summary(lo_ok, hi_ok, nominal_one_sided):
    lo_ok, hi_ok = np.asarray(lo_ok), np.asarray(hi_ok)
    return {"n_rep": int(lo_ok.size), "nominal_one_sided": nominal_one_sided,
            "coverage_lower": float(lo_ok.mean()), "coverage_upper": float(hi_ok.mean()),
            "mc_se": float(math.sqrt(0.95 * 0.05 / lo_ok.size)),
            "pass": bool(lo_ok.mean() >= COV_MIN and hi_ok.mean() >= COV_MIN),
            "within_2se_of_0.95": bool(min(lo_ok.mean(), hi_ok.mean()) >= COV_MIN - 2 * math.sqrt(0.95 * 0.05 /
                                                                                                    lo_ok.size))}


@pytest.mark.parametrize("het", ["het", "hom"])
@pytest.mark.parametrize("unit", ["cluster", "project"])
def test_iii_g1_jackknife_coverage(unit, het):
    """True Spearman 0 (unit-level rho_j ~ N(0, 0.03), symmetric); n_top 4,000 per replicate; J = 37 clusters
    (project unit: 74 projects with project-level heterogeneity)."""
    g = np.random.default_rng((2001 if unit == "cluster" else 2002) + (10 if het == "hom" else 0))
    tau = 0.03 if het == "het" else 0.0
    J = J_COV if unit == "cluster" else 2 * J_COV
    n = 4000
    lo_ok, hi_ok, env_lo, env_hi = [], [], [], []
    for _ in range(N_COV):
        u = g.choice(J, size=n, p=_weights(g, J))
        rj = g.normal(0, tau, J) if tau else np.zeros(J)
        rl = 2 * np.sin(np.pi * rj / 6)
        x = g.standard_normal(n)
        y = rl[u] * x + np.sqrt(1 - rl[u] ** 2) * g.standard_normal(n)
        r = S.jackknife_t_interval(x, y, u, S.G1_CONF)
        lo_ok.append(r["lo_z"] <= 0)
        hi_ok.append(r["hi_z"] >= 0)
        fz = GS.spearman_ci(x, y, S.G1_CONF)["ci"]
        env_lo.append(min(r["lo_z"], math.atanh(fz[0])) <= 0)
        env_hi.append(max(r["hi_z"], math.atanh(fz[1])) >= 0)
    s = _cov_summary(lo_ok, hi_ok, 0.9975)
    s["envelope_with_fisher"] = _cov_summary(env_lo, env_hi, 0.9975)
    s["tau"] = tau
    _rec("iii_coverage", **{f"G1_jackknife_{unit}_{het}": s})
    assert s["pass"], s


@pytest.mark.parametrize("het", ["het", "hom"])
@pytest.mark.parametrize("unit", ["cluster", "project"])
def test_iii_tost_bootstrap_coverage(unit, het):
    """recall_T - recall_D, true 0: per-unit gap_j ~ N(0, 0.01) on recall 0.95; n = 4,000 per label; two-sided 0.95."""
    g = np.random.default_rng((3001 if unit == "cluster" else 3002) + (10 if het == "hom" else 0))
    tau = 0.01 if het == "het" else 0.0
    J = J_COV if unit == "cluster" else 2 * J_COV
    n = 4000
    lo_ok, hi_ok, env_lo, env_hi = [], [], [], []
    for _ in range(N_COV):
        w = _weights(g, J)
        nT, nD = g.multinomial(n, w), g.multinomial(n, w)
        gap = g.normal(0, tau, J) if tau else np.zeros(J)
        kT, kD = g.binomial(nT, np.clip(0.95 + gap, 0, 1)), g.binomial(nD, 0.95)
        uT, uD = np.repeat(np.arange(J), nT), np.repeat(np.arange(J), nD)
        unit_arr = np.r_[uT, uD]
        hitT = np.r_[np.concatenate([np.arange(nT[j]) < kT[j] for j in range(J)]), np.zeros(n, bool)]
        isT = np.r_[np.ones(n, bool), np.zeros(n, bool)]
        hitD = np.r_[np.zeros(n, bool), np.concatenate([np.arange(nD[j]) < kD[j] for j in range(J)])]
        isD = ~isT
        r = S.unit_bootstrap_diff(unit_arr, hitT, isT, hitD, isD, B_COV, g, S.TOST_CONF)
        lo_ok.append(r["lo"] <= 0)
        hi_ok.append(r["hi"] >= 0)
        _, (nl, nh) = GS.newcombe_diff(int(kT.sum()), n, int(kD.sum()), n, S.TOST_CONF)
        env_lo.append(min(r["lo"], nl) <= 0)
        env_hi.append(max(r["hi"], nh) >= 0)
    s = _cov_summary(lo_ok, hi_ok, 0.975)
    s["envelope_with_newcombe"] = _cov_summary(env_lo, env_hi, 0.975)
    s["tau"] = tau
    _rec("iii_coverage", **{f"TOST_bootstrap_{unit}_{het}": s})
    assert s["pass"], s


@pytest.mark.parametrize("het", ["het", "hom"])
@pytest.mark.parametrize("unit", ["cluster", "project"])
def test_iii_g2_original_bootstrap_coverage(unit, het):
    """paired mean(learned - fixed) target-win, true 0: fixed ~ Bern(0.5), learned ~ Bern(0.5 + d_j),
    d_j ~ N(0, 0.01); n = 20,000; one-sided 0.95 UCB."""
    g = np.random.default_rng((4001 if unit == "cluster" else 4002) + (10 if het == "hom" else 0))
    tau = 0.01 if het == "het" else 0.0
    J = J_COV if unit == "cluster" else 2 * J_COV
    n = 20_000
    hi_ok, env_hi = [], []
    for _ in range(N_COV):
        u = g.choice(J, size=n, p=_weights(g, J))
        d = g.normal(0, tau, J) if tau else np.zeros(J)
        a = g.random(n) < 0.5
        b = g.random(n) < 0.5 + d[u]
        ones = np.ones(n, bool)
        r = S.unit_bootstrap_diff(u, b, ones, a, ones, B_COV, g, S.G2O_CONF, two_sided=False)
        hi_ok.append(r["hi"] >= 0)
        an = G2.g2_differential(a, b, g.normal(size=n))["strata"][0]["overall"]["ucb"]
        env_hi.append(max(r["hi"], an) >= 0)
    s = _cov_summary(np.ones(N_COV, bool), hi_ok, 0.95)
    s["coverage_lower"] = None
    s["envelope_with_analytic"] = _cov_summary(np.ones(N_COV, bool), env_hi, 0.95)
    s["envelope_with_analytic"]["coverage_lower"] = None
    s["tau"] = tau
    _rec("iii_coverage", **{f"G2_original_bootstrap_{unit}_{het}": s})
    assert s["pass"], s


# ------------------------------------------------------------------------------------------------ v1 combiner
def test_v1_verdict_combiner():
    base = {"tail_G2": {"verdict": "pass"}, "G1": {"verdict": "pass"}, "G3": {"verdict": "pass"},
            "G2_original": {"verdict": "pass"}, "G5": {"verdict": "pass"},
            "FS": {"tail_G2": {"verdict": S.NOT_EVALUABLE}, "G1": {"verdict": S.INVALID}, "G5": {"verdict": "pass"}},
            "subset_check": {"verdict": "pass"}}
    assert S.v1_verdict(base, base)["status"] == "PASS"
    bad = json.loads(json.dumps(base))
    bad["FS"]["G1"]["verdict"] = "fail"
    assert S.v1_verdict(base, bad)["status"] == "FAILED"
    sub = json.loads(json.dumps(base))
    sub["subset_check"]["verdict"] = "fail"
    assert S.v1_verdict(sub, base)["status"] == "FAILED"
    inv = json.loads(json.dumps(base))
    inv["tail_G2"]["verdict"] = S.INVALID
    assert S.v1_verdict(inv, bad)["status"] == "invalid pool"
    ne = json.loads(json.dumps(base))
    ne["G3"]["verdict"] = S.NOT_EVALUABLE
    assert S.v1_verdict(base, ne)["status"] == "not evaluable"
    _rec("v1_combiner", ok=True)


def test_evaluate_pool_smoke():
    """End-to-end on a small synthetic pool (n_min lowered for the smoke run only; production uses the floor)."""
    g = np.random.default_rng(11)
    n, J = 30_000, 20
    cl = g.integers(0, J, n)
    pr = 2 * cl + (g.random(n) < 0.5)
    pool = np.where(g.random(n) < 0.85, "A", "FS")

    def pipe():
        sc = g.normal(size=n)
        return {"tail1": sc > 1.0, "win_t": g.random(n) < 0.5, "score": sc, "logit": g.normal(size=n)}
    D = {"pool": pool, "project": pr, "cluster": cl, "instrument": np.where(cl % 2, "QE", "HF"),
         "fixed": pipe(), "v1": pipe(),
         "g3": {"bt": g.normal(10, 3, n), "bd": g.normal(10, 3, n)}}
    D["g3"]["kt"], D["g3"]["kd"] = D["g3"]["bt"].copy(), D["g3"]["bd"].copy()
    r = S.evaluate_pool(D, np.ones(n, bool), B=300, n_min={"g1": 2000, "g3": 500})
    v = S.pool_row_verdicts(r)
    assert set(v) == {"tail_G2", "G1", "G3", "G2_original", "G5", "FS", "subset_check"}
    _rec("v1_combiner", evaluate_pool_smoke=v)


def test_zz_write_results():
    RES["wall_s"] = time.time() - RES.pop("t0")
    cov = RES["items"].get("iii_coverage", {})
    RES["summary"] = {
        "i_exact_a2_stats": "i_exact_a2_stats" in RES["items"],
        "ii_c9_a4_reproduction": "c9_tail_iut" in RES["items"].get("ii_c9_a4_reproduction", {}),
        "iii_coverage_all_components_ge_0.95": bool(cov) and all(v["pass"] for v in cov.values()),
        "n_coverage_components": len(cov), "coverage_reps": N_COV, "coverage_B": B_COV, "J": J_COV,
        "failed_coverage_components": sorted(k for k, v in cov.items() if not v["pass"])}
    RES["summary"]["acceptance_pass"] = bool(RES["summary"]["i_exact_a2_stats"]
                                             and RES["summary"]["ii_c9_a4_reproduction"]
                                             and RES["summary"]["iii_coverage_all_components_ge_0.95"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(RES, open(OUT, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
