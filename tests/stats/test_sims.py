"""Monte-Carlo estimator verification (DESIGN.md §6.7).

Default runs use reduced replicate counts (total well under a minute). The
``slow`` variants use >= 1e4 replicates (``pytest --runslow``); the same
simulators produce results/2026-10-26_stats_sims/metrics.json via
scripts/run_stats_sims.py.
"""
import numpy as np
import pytest

from vulcan.stats import sim
from vulcan.stats.g2 import g2_absolute, g2_differential, g2_sample_size, simulate_g2

GATED = ("0.01", "0.05")
SLOW_R = 10_000


def _check(res, alphas=GATED):
    for a in alphas:
        r = res[a]
        # mean FDP <= alpha up to Monte-Carlo error (99% one-sided)
        assert r["mean_fdp"] - 2.326 * r["se"] <= float(a), (a, r)


# --------------------------------------------------------------------------- TDC
@pytest.mark.parametrize("R", [2000, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_tdc_controls(R):
    _check(sim.sim_tdc(R))
    tight = sim.sim_tdc(R, m=4000, mu=3.0, true_always_win=True)
    _check(tight)
    assert tight["0.05"]["mean_fdp"] > 0.045  # the tight regime really is tight


@pytest.mark.parametrize("R", [500, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_double_competition_controls(R):
    res = sim.sim_double_competition(R)
    _check(res["peptide_double_competition"])
    assert res["peptide_double_competition"]["0.01"]["mean_discoveries"] > 100


@pytest.mark.parametrize("R", [1000, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_group_walk_controls(R):
    res = sim.sim_groupwalk(R, frac_true=(0.05, 0.3, 0.8), mus=(2.0, 3.0, 4.0))
    _check(res["group_walk"])
    _check(res["pooled_tdc"])
    # group structure helps here
    assert res["group_walk"]["0.01"]["mean_discoveries"] > res["pooled_tdc"]["0.01"]["mean_discoveries"]


@pytest.mark.parametrize("mode", ["single", "two"])
@pytest.mark.parametrize("R", [300, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_reset_controls(R, mode):
    res = sim.sim_reset(R, mode=mode)
    assert res["frac_training_not_started"] < 0.05   # the learned path is what is being tested
    _check(res)


@pytest.mark.parametrize("R", [1000, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_picked_protein_controls(R):
    res = sim.sim_picked_protein(R)
    _check(res["picked"])
    _check(res["classic"])


@pytest.mark.parametrize("R", [50, pytest.param(2000, marks=pytest.mark.slow)])
def test_picked_group_general_controls(R):
    res = sim.sim_picked_group_general(R)
    _check(res["picked_group"])


@pytest.mark.parametrize("R", [1000, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_multi_decoy_controls(R):
    _check(sim.sim_multi_decoy(R))


def test_entrapment_estimators_bracket_true_fdp():
    res = sim.sim_entrapment_estimators(500)
    for a in GATED:
        r = res[a]
        assert r["lb"] <= r["true"] <= r["paired"] <= r["combined"] + 1e-12


# ------------------------------------------------------------- biased null p=0.55
@pytest.mark.parametrize("R", [1000, pytest.param(SLOW_R, marks=pytest.mark.slow)])
def test_biased_null_tdc_inflates_and_g2_fires(R):
    fair = sim.sim_tdc(R, m=4000, mu=3.0, true_always_win=True)
    biased = sim.sim_tdc(R, m=4000, mu=3.0, true_always_win=True, p_null_target=0.55)
    for a in GATED:
        b = biased[a]
        assert b["mean_fdp"] - 2.576 * b["se"] > float(a)          # TDC no longer controls
        ratio = b["mean_fdp"] / fair[a]["mean_fdp"]
        assert abs(ratio - 0.55 / 0.45) < 0.1                       # ~ p/(1-p) inflation
    # differential G2 fires (fixed pipeline fair, learned pipeline p = 0.55)
    n = g2_sample_size()
    assert simulate_g2(n, 10_000, 0.5, 0.55, seed=11)["fire_rate"] == 1.0
    # and on explicit per-competition data, also via the absolute diagnostic
    rng = np.random.default_rng(131)
    s = rng.normal(size=n)
    a = rng.random(n) < 0.5
    b = rng.random(n) < 0.55
    assert g2_differential(a, b, s)["fires"]
    assert not g2_absolute(b, s)["pass"]
