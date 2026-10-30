"""A2 gate statistics (vulcan.gates.a2_stats): tail-G2 max rule + pool validity (A2.5), G1 n_min + control (A2.6),
G3 TOST per bin + n_min + merging (A2.8).  Synthetic data only."""
import numpy as np
import pytest

from vulcan.gates import a2_stats as S
from vulcan.gates import stats as GS


def _tails(rng, J=10, n=600, pL=0.5, pF=0.5, n_null=200_000):
    proj = rng.integers(0, J, n_null)
    tl = np.zeros(n_null, bool)
    tl[rng.choice(n_null, n, replace=False)] = True
    tf = np.zeros(n_null, bool)
    tf[rng.choice(n_null, n, replace=False)] = True
    gl = rng.random(n_null) < pL
    gf = rng.random(n_null) < pF
    return proj, tl, gl, tf, gf


def test_a1_margins_imported_unchanged():
    assert (S.TAIL_N_MIN, S.TAIL_MARGIN_FRAC, S.TAIL_MARGIN_DIFF, S.TAIL_CONF) == (200, 0.60, 0.10, 0.95)
    assert (S.BOOT_B, S.BOOT_SEED, S.G1_N_MIN, S.G3_N_MIN, S.G3_CONF) == (10_000, 83, 7_200, 2_000, 0.90)


def test_tail_not_evaluable_below_200():
    rng = np.random.default_rng(1)
    r = S.tail_g2_a2(*_tails(rng, n=150), B=500)
    assert r["verdict"] == S.NOT_EVALUABLE


def test_tail_fair_passes_and_max_rule():
    rng = np.random.default_rng(2)
    r = S.tail_g2_a2(*_tails(rng, n=800), B=2000)
    assert r["verdict"] == "pass", r
    assert r["ucb_abs"] == max(r["ucb_abs_cp"], r["ucb_abs_cluster"])
    assert r["ucb_diff"] == max(r["ucb_diff_delta"], r["ucb_diff_cluster"])


def test_tail_biased_fails():
    rng = np.random.default_rng(3)
    r = S.tail_g2_a2(*_tails(rng, n=800, pL=0.72), B=2000)
    assert r["verdict"] == "fail" and not r["clause_abs_pass"] and not r["clause_diff_pass"]


@pytest.mark.parametrize("pF", [0.30, 0.70])
def test_tail_invalid_pool_either_direction(pF):
    rng = np.random.default_rng(4)
    r = S.tail_g2_a2(*_tails(rng, n=800, pF=pF), B=500)
    assert not r["pool_validity"]["valid"] and r["verdict"] == S.INVALID


def test_tail_empty_hyperscore_tail_invalid():
    rng = np.random.default_rng(5)
    proj, tl, gl, tf, gf = _tails(rng, n=800)
    r = S.tail_g2_a2(proj, tl, gl, np.zeros_like(tf), gf, B=200)
    assert r["verdict"] == S.INVALID


def test_tail_cluster_bound_deterministic_and_order_free():
    rng = np.random.default_rng(6)
    proj, tl, gl, tf, gf = _tails(rng, n=700)
    a = S.cluster_bootstrap_tail(proj, tl, gl, tf, gf, B=3000)
    b = S.cluster_bootstrap_tail(proj, tl, gl, tf, gf, B=3000)
    assert a == b
    perm = rng.permutation(proj.size)
    c = S.cluster_bootstrap_tail(proj[perm], tl[perm], gl[perm], tf[perm], gf[perm], B=3000)
    assert a["ucb_abs"] == c["ucb_abs"] and a["ucb_diff"] == c["ucb_diff"]


def test_tail_single_project_cannot_pass():
    rng = np.random.default_rng(7)
    proj, tl, gl, tf, gf = _tails(rng, n=800)
    r = S.tail_g2_a2(np.zeros_like(proj), tl, gl, tf, gf, B=200)
    assert r["ucb_abs"] == float("inf") and r["verdict"] == "fail"


def test_tail_heterogeneous_projects_cluster_term_binds():
    # one project carries a pro-target tail (the M3e S. aureus pattern): the cluster term must exceed CP
    rng = np.random.default_rng(8)
    J, per = 10, 80
    proj = np.repeat(np.arange(J), per)
    n = proj.size
    tl = np.ones(n, bool)
    gl = rng.random(n) < np.where(proj == 0, 0.95, 0.47)
    tf = np.ones(n, bool)
    gf = rng.random(n) < 0.5
    r = S.tail_g2_a2(proj, tl, gl, tf, gf, B=4000)
    assert r["ucb_abs_cluster"] > r["ucb_abs_cp"]


def test_g1_not_evaluable_below_nmin():
    rng = np.random.default_rng(9)
    n = 50_000                               # top decile 5,000 < 7,200
    r = S.g1_gate_a2(rng.random(n), rng.standard_normal(n), rng.random(n), rng.standard_normal(n))
    assert r["n_top"] < 7200 and r["verdict"] == S.NOT_EVALUABLE


def test_g1_pass_and_fail_and_invalid_control():
    rng = np.random.default_rng(10)
    n = 120_000
    s, lg = rng.standard_normal(n), rng.standard_normal(n)
    sf, lf = rng.standard_normal(n), rng.standard_normal(n)
    assert S.g1_gate_a2(s, lg, sf, lf)["verdict"] == "pass"
    lg2 = lg + 0.6 * s                        # strong correlation in the top decile
    assert S.g1_gate_a2(s, lg2, sf, lf)["verdict"] == "fail"
    lf2 = lf + 0.6 * sf                       # hyperscore itself correlated: pool invalid
    assert S.g1_gate_a2(s, lg, sf, lf2)["verdict"] == S.INVALID


def test_fisher_pooled_rho_zero():
    rng = np.random.default_rng(11)
    r = S.fisher_pooled_rho([(rng.random(5000), rng.random(5000)), (rng.random(3000), rng.random(3000))])
    assert r["ci"][0] < 0 < r["ci"][1]


# ------------------------------------------------------------------------------------------------ G3
def _g3_data(rng, n_per_bin_label, rT, rD, rD_top=None, n_bins=10):
    n = n_per_bin_label * n_bins
    best = rng.random(n) * 100
    is_t = rng.random(n) < 0.5
    bt = np.where(is_t | (rng.random(n) < 0.9), best - rng.random(n) * is_t, -np.inf)
    bd = np.where(~is_t | (rng.random(n) < 0.9), best - rng.random(n) * ~is_t, -np.inf)
    top = best >= np.quantile(best, 0.9)
    recD = np.where(top & (rD_top is not None), rD_top if rD_top is not None else rD, rD)
    kt = np.where(rng.random(n) < rT, bt, bt - 1)
    kd = np.where(rng.random(n) < recD, bd, bd - 1)
    return bt, bd, kt, kd


def test_g3_groups_logic():
    g, why = S.g3_groups([2500] * 10, [2500] * 10)
    assert why is None and g == [[i] for i in range(10)]
    g, why = S.g3_groups([1200, 1200, 900, 2100, 2500] + [3000] * 5, [1500, 800, 2000, 2100, 2500] + [3000] * 5)
    assert why is None and g[:3] == [[0, 1], [2, 3], [4]]
    g, why = S.g3_groups([2500, 2500, 2500, 2500, 1000] + [3000] * 5, [2500] * 10)
    assert g is None and "bottom-half" in why                          # leftover containing decile 4
    g, why = S.g3_groups([2500] * 9 + [1999], [2500] * 10)
    assert g is None and "top-half decile 9" in why


def test_g3_not_evaluable_small_bins():
    rng = np.random.default_rng(12)
    r = S.g3_gate_a2(*_g3_data(rng, 300, 0.999, 0.999))
    assert r["verdict"] == S.NOT_EVALUABLE


def test_g3_pass_symmetric_and_fail_deficit():
    rng = np.random.default_rng(13)
    r = S.g3_gate_a2(*_g3_data(rng, 6000, 0.999, 0.999))
    assert r["verdict"] == "pass", [x["ci90"] for x in r.get("groups", [])]
    r = S.g3_gate_a2(*_g3_data(rng, 6000, 0.999, 0.999, rD_top=0.97))
    assert r["verdict"] == "fail"


def test_tost_never_passes_when_a1_point_rule_fails():
    rng = np.random.default_rng(14)
    for _ in range(3000):
        n1, n2 = rng.integers(50, 5000, 2)
        k1 = rng.binomial(n1, rng.uniform(0.95, 1.0))
        k2 = rng.binomial(n2, rng.uniform(0.95, 1.0))
        t = S.tost_bin(int(k1), int(n1), int(k2), int(n2))
        if t["pass"]:
            assert abs(k1 / n1 - k2 / n2) <= 0.01


def test_tost_interval_is_newcombe_90():
    d, ci = GS.newcombe_diff(1990, 2000, 1985, 2000, 0.90)
    t = S.tost_bin(1990, 2000, 1985, 2000)
    assert t["ci90"] == ci and t["diff"] == d
