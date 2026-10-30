import numpy as np
import pytest

from vulcan.stats.reset import RESET_CONSTANTS, LDALearner, LinearSVMLearner, reset, two_decoy_competition


def _manual_q(score, is_dec, factor):
    o = np.argsort(-score)
    D = np.cumsum(is_dec[o])
    T = np.cumsum(~is_dec[o])
    f = factor * (D + 1) / np.maximum(T, 1)
    q = np.minimum.accumulate(np.minimum(f, 1)[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = q
    return out


def test_constants():
    assert RESET_CONSTANTS["single"] == {"c_train": 3.0, "c_final": 2.0}
    assert RESET_CONSTANTS["two"] == {"c_train": 2.0, "c_final": 1.0}


@pytest.mark.parametrize("mode,cf", [("single", 2.0), ("two", 1.0)])
def test_final_estimator_explicit(mode, cf):
    rng = np.random.default_rng(11)
    n = 600
    s = rng.normal(size=n) + np.r_[np.full(200, 3.0), np.zeros(n - 200)]
    isd = np.r_[np.zeros(200, bool), rng.random(n - 200) < 0.5]
    r = reset(s[:, None], isd, s, mode=mode, seed=83, n_iter=0)
    keep = ~r.is_training_decoy
    assert np.all(np.isnan(r.q[~keep]))
    assert not np.any(r.is_training_decoy & ~isd)  # only decoys can be training decoys
    # FDR_final = c * (D_est + 1) / T over targets + estimating decoys
    np.testing.assert_allclose(r.q[keep], _manual_q(s[keep], isd[keep], cf))
    # roughly half the decoys go to training
    frac = r.is_training_decoy.sum() / isd.sum()
    assert 0.4 < frac < 0.6


def test_training_positive_set_explicit():
    rng = np.random.default_rng(23)
    n = 3000
    s = rng.normal(size=n) + np.r_[np.full(1200, 4.0), np.zeros(n - 1200)]
    isd = np.r_[np.zeros(1200, bool), rng.random(n - 1200) < 0.5]
    X = np.column_stack([s, rng.normal(size=n)])
    r = reset(X, isd, s, mode="single", seed=47, n_iter=1)
    tr = r.is_training_decoy
    qt = _manual_q(s, tr, 3.0)  # pseudo-targets vs training decoys, factor 3
    assert r.n_positives[0] == int(((~tr) & (qt <= 0.01)).sum())


def test_seeded_split_reproducible():
    rng = np.random.default_rng(1)
    s = rng.normal(size=200)
    isd = rng.random(200) < 0.5
    a = reset(s[:, None], isd, s, seed=131, n_iter=0)
    b = reset(s[:, None], isd, s, seed=131, n_iter=0)
    np.testing.assert_array_equal(a.is_training_decoy, b.is_training_decoy)


def test_learners_recover_direction():
    rng = np.random.default_rng(5)
    Xp = rng.normal(size=(300, 3)) + np.array([2.0, 0.0, -1.0])
    Xn = rng.normal(size=(300, 3))
    for L in (LDALearner(), LinearSVMLearner(grid=(1.0,), cv=2)):
        L.fit(Xp, Xn)
        sp, sn = L.decision_function(Xp), L.decision_function(Xn)
        assert sp.mean() > sn.mean() + 1.0


def test_reset_learning_improves_power():
    rng = np.random.default_rng(83)
    m = 2000
    null = np.arange(m) >= 800
    XT = rng.normal(size=(m, 3)); XT[~null] += [3.5, 2.0, 1.5]
    XD = rng.normal(size=(m, 3))
    dwin = XD[:, 0] > XT[:, 0]
    X = np.where(dwin[:, None], XD, XT)
    r = reset(X, dwin, X[:, 0], seed=11)
    r0 = reset(X, dwin, X[:, 0], seed=11, n_iter=0)
    assert r.accepted(0.01, dwin).sum() > r0.accepted(0.01, dwin).sum()
    # best-feature start (Percolator style) also works
    r2 = reset(X, dwin, None, seed=11)
    assert r2.accepted(0.01, dwin).sum() > 0


def test_two_decoy_competition():
    which, w = two_decoy_competition([1.0, 5.0, 2.0], [3.0, 1.0, 2.5], [2.0, 0.0, 4.0], seed=1)
    np.testing.assert_array_equal(which, [1, 0, 2])
    np.testing.assert_allclose(w, [3.0, 5.0, 4.0])


def test_bad_mode():
    with pytest.raises(ValueError):
        reset(np.zeros((3, 1)), np.zeros(3, bool), np.zeros(3), mode="triple")
