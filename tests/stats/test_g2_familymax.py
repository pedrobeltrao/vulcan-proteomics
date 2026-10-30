import numpy as np
import pytest

from vulcan.stats.familymax import ReplicateSpec, reconstruct_competition, replicate_specs, top_k_per_family
from vulcan.stats.g2 import (g2_absolute, g2_differential, g2_power, g2_sample_size, ks_winners, simulate_g2,
                             wilson_ci)


def test_wilson():
    p, lo, hi = wilson_ci(50, 100, 0.95)
    assert p == 0.5 and lo == pytest.approx(0.4038, abs=1e-4) and hi == pytest.approx(0.5962, abs=1e-4)


def test_sample_size_and_power():
    n = g2_sample_size()
    assert 90_000 < n < 140_000
    assert g2_power(n)["p_pass_joint_lb"] >= 0.9 > g2_power(n - 1000)["p_pass_joint_lb"]
    # correlated pipelines need fewer competitions
    assert g2_sample_size(rho=0.5) < n


def test_g2_differential_pass_and_fire():
    rng = np.random.default_rng(11)
    n = 120_000
    s = rng.normal(size=n)
    a = rng.random(n) < 0.5
    same = g2_differential(a, a, s)
    assert same["pass"]
    b = np.where(rng.random(n) < 0.5, a, rng.random(n) < 0.5)   # fair, partly coupled
    assert g2_differential(a, b, s)["pass"]
    c = np.where(rng.random(n) < 0.5, a, rng.random(n) < 0.55)  # learned excess ~ +0.025
    r = g2_differential(a, c, s)
    assert r["fires"] and r["strata"][0]["overall"]["diff"] > 0.01
    strat = np.where(np.arange(n) < n // 2, "OT-HCD", "IT-CID")
    r = g2_differential(a, a, s, stratum=strat)
    assert set(r["strata"]) == {"OT-HCD", "IT-CID"} and r["pass"]


def test_g2_absolute_and_ks():
    rng = np.random.default_rng(23)
    n = 200_000
    s = rng.normal(size=n)
    w = rng.random(n) < 0.5
    assert g2_absolute(w, s)["overall"]["within_band"]
    w2 = rng.random(n) < 0.55
    assert not g2_absolute(w2, s)["pass"]
    k = ks_winners(s, w)
    assert 0 <= k["statistic"] < 0.02 and k["pvalue"] > 1e-4
    k2 = ks_winners(s + 0.3 * w, w)
    assert k2["pvalue"] < 1e-10


def test_g2_mc_power_and_firing():
    """DESIGN v0.3 §6.4: >= 90% power at p = 0.5 and the gate fires at +0.02."""
    # the helper's n for 90% power sits exactly on the boundary; budget for 92%
    n = g2_sample_size(power=0.92)
    ok = simulate_g2(n, 10_000, 0.5, 0.5, rho=0.0, seed=11)
    assert ok["pass_rate"] >= 0.9
    assert abs(simulate_g2(g2_sample_size(), 10_000, seed=23)["pass_rate"] - 0.9) < 0.015
    fire = simulate_g2(n, 10_000, 0.5, 0.52, rho=0.0, seed=23)
    assert fire["fire_rate"] >= 0.99
    biased = simulate_g2(n, 10_000, 0.5, 0.55, rho=0.0, seed=47)
    assert biased["fire_rate"] == 1.0


def _random_candidates(rng, n_spec=150, n_hyp=2, n_per_fam=4):
    fams = ["T", "E_2001", "E_2002", "dec_1001(T)", "dec_1002(T)", "dec_1001(E_2001)", "dec_1002(E_2001)",
            "dec_1001(E_2002)", "dec_1002(E_2002)"]
    rows = []
    for sp in range(n_spec):
        for h in range(n_hyp):
            for f in fams:
                for _ in range(n_per_fam):
                    rows.append((sp, h, f))
    sp, hy, fa = map(np.array, zip(*rows))
    score = np.round(rng.normal(size=sp.size), 1)   # coarse -> many exact ties
    tiebreak = rng.permutation(sp.size)            # stands in for hash(sequence, seed)
    return sp, hy, fa, score, tiebreak


def test_familymax_reconstruction_equals_separate_competition():
    rng = np.random.default_rng(83)
    sp, hy, fa, score, tb = _random_candidates(rng)
    keep = top_k_per_family(sp, hy, fa, score, tb, k=2)
    assert keep.size == 150 * 2 * 9 * 2
    specs = replicate_specs([1001, 1002], [2001, 2002])
    assert len(specs) == 4
    for spec in specs.values():
        full = reconstruct_competition(sp, fa, score, tb, spec)          # "separate search"
        fm = reconstruct_competition(sp[keep], fa[keep], score[keep], tb[keep], spec)
        np.testing.assert_array_equal(keep[fm.row], full.row)
        np.testing.assert_array_equal(fm.runner_up, full.runner_up)
        np.testing.assert_array_equal(fm.is_decoy, full.is_decoy)
        np.testing.assert_array_equal(fm.is_entrapment, full.is_entrapment)
    # top-1 per family is NOT enough for the runner-up score
    k1 = top_k_per_family(sp, hy, fa, score, tb, k=1)
    spec = ReplicateSpec(("T",), (), ("dec_1001(T)",))
    full = reconstruct_competition(sp, fa, score, tb, spec)
    fm1 = reconstruct_competition(sp[k1], fa[k1], score[k1], tb[k1], spec)
    np.testing.assert_array_equal(k1[fm1.row], full.row)
    assert not np.array_equal(fm1.runner_up, full.runner_up)


def test_leak_auc():
    from vulcan.stats.g2 import leak_auc

    rng = np.random.default_rng(5)
    s = rng.normal(size=20000)
    w = rng.random(20000) < 0.5
    assert not leak_auc(s, w)["alarm"]
    assert leak_auc(s + 0.5 * w, w)["alarm"]
    assert leak_auc([1.0, 2.0], [False, True])["auc"] == 1.0
