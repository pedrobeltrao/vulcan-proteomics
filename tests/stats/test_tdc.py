import numpy as np
import pytest

from vulcan.stats.tdc import (compete, count_discoveries, double_competition, multi_decoy_competition,
                              multi_decoy_qvalues, psm_competition, tdc_fdr_curve, tdc_qvalues)


def test_tdc_qvalues_hand():
    # scores desc, labels T T D T T D
    s = np.array([10, 9, 8, 7, 6, 5.0])
    d = np.array([0, 0, 1, 0, 0, 1], bool)
    fdr, T, D = tdc_fdr_curve(s, d)
    np.testing.assert_allclose(fdr, [1.0, 0.5, 1.0, 2 / 3, 0.5, 0.75])
    np.testing.assert_array_equal(T, [1, 2, 2, 3, 4, 4])
    np.testing.assert_array_equal(D, [0, 0, 1, 1, 1, 2])
    q = tdc_qvalues(s, d)
    np.testing.assert_allclose(q, [0.5, 0.5, 0.5, 0.5, 0.5, 0.75])
    # input order does not matter
    p = np.array([3, 0, 5, 1, 4, 2])
    np.testing.assert_allclose(tdc_qvalues(s[p], d[p]), q[p])
    assert count_discoveries(q, d, 0.5) == 4
    assert count_discoveries(q, d, 0.49) == 0


def test_tdc_factor_and_no_plus_one():
    s = np.array([10, 9, 8, 7, 6, 5.0])
    d = np.array([0, 0, 1, 0, 0, 1], bool)
    fdr2, _, _ = tdc_fdr_curve(s, d, factor=2.0)
    np.testing.assert_allclose(fdr2, [2.0, 1.0, 2.0, 4 / 3, 1.0, 1.5])
    q2 = tdc_qvalues(s, d, factor=2.0)
    np.testing.assert_allclose(q2, [1.0] * 6)  # capped at 1
    q0 = tdc_qvalues(s, d, plus_one=False)
    np.testing.assert_allclose(q0, [0, 0, 1 / 4, 1 / 4, 1 / 4, 0.5])


def test_tdc_ties_decoys_counted_first():
    # target and decoy tied at the top: the target's threshold includes the decoy
    s = np.array([3, 3, 2, 1.0])
    d = np.array([0, 1, 0, 0], bool)
    fdr, T, D = tdc_fdr_curve(s, d)
    np.testing.assert_array_equal(T, [1, 1, 2, 3])
    np.testing.assert_array_equal(D, [1, 1, 1, 1])
    np.testing.assert_allclose(fdr, [2.0, 2.0, 1.0, 2 / 3])
    q = tdc_qvalues(s, d)
    assert q[0] == q[1] == pytest.approx(2 / 3)
    # tied targets share a threshold and a q-value
    s = np.array([5, 5, 5, 1.0])
    d = np.array([0, 0, 0, 1], bool)
    q = tdc_qvalues(s, d)
    np.testing.assert_allclose(q[:3], 1 / 3)


def test_tdc_blocks_independent():
    rng = np.random.default_rng(11)
    s = rng.normal(size=300)
    d = rng.random(300) < 0.3
    blk = np.repeat([0, 1, 2], 100)
    qb = tdc_qvalues(s, d, block=blk)
    for b in range(3):
        m = blk == b
        np.testing.assert_allclose(qb[m], tdc_qvalues(s[m], d[m]))


def test_qvalue_equals_largest_set_rule():
    rng = np.random.default_rng(23)
    for _ in range(50):
        n = 200
        s = rng.normal(size=n) + np.r_[np.full(80, 2.0), np.zeros(n - 80)]
        d = rng.random(n) < 0.35
        q = tdc_qvalues(s, d)
        for a in (0.01, 0.05, 0.2):
            o = np.argsort(-s)
            D = np.cumsum(d[o])
            T = np.cumsum(~d[o])
            ok = np.flatnonzero((D + 1) / np.maximum(T, 1) <= a)
            n_ref = int(T[ok[-1]]) if ok.size else 0
            assert count_discoveries(q, d, a) == n_ref


def test_compete_winner_and_seeded_ties():
    g = np.array([0, 0, 1, 1, 1, 2])
    s = np.array([1.0, 2.0, 5.0, 5.0, 3.0, 0.0])
    w = compete(g, s, seed=11)
    assert w[0] == 1 and w[2] == 5 and w[1] in (2, 3)
    np.testing.assert_array_equal(w, compete(g, s, seed=11))
    picks = [compete(g, s, seed=k)[1] for k in range(400)]
    frac = np.mean(np.array(picks) == 2)
    assert 0.4 < frac < 0.6
    # explicit tiebreak key: smallest wins
    w = compete(g, s, tiebreak=np.array([0, 0, 9, 1, 0, 0]))
    assert w[1] == 3
    # NaN treated as -inf
    w = compete(np.array([0, 0]), np.array([np.nan, -5.0]), seed=1)
    assert w[0] == 1


def test_psm_competition_all_hypotheses():
    # spectrum 7 has two charge hypotheses each with target+decoy candidates
    spec = np.array([7, 7, 7, 7, 8, 8])
    score = np.array([1.0, 4.0, 3.0, 2.0, 0.5, 0.7])
    w = psm_competition(spec, score, seed=1)
    np.testing.assert_array_equal(w, [1, 5])


def test_double_competition_hand():
    # spectra -> candidates (unit = peptide, pair id shared by target & decoy)
    spec = np.array([0, 0, 1, 1, 2, 2, 3, 3])
    unit = np.array(["A", "rA", "A", "B", "rB", "C", "rA", "C"])
    pair = np.array(["A", "A", "A", "B", "B", "C", "A", "C"])
    dec = np.array([0, 1, 0, 0, 1, 0, 1, 0], bool)
    score = np.array([5.0, 1.0, 3.0, 2.0, 4.0, 1.0, 6.0, 2.0])
    r = double_competition(spec, unit, pair, score, dec, seed=11)
    # PSM winners: s0->A(5), s1->A(3), s2->rB(4), s3->rA(6)
    np.testing.assert_array_equal(np.sort(r.psm_winners), [0, 2, 4, 6])
    # best per peptide: A=5, rB=4, rA=6 ; pair A: rA(6) beats A(5); pair B: rB(4)
    got = dict(zip(r.pair.tolist(), zip(r.score.tolist(), r.is_decoy.tolist())))
    assert got == {"A": (6.0, True), "B": (4.0, True)}


def test_multi_decoy_d1_is_tdc():
    rng = np.random.default_rng(47)
    t = rng.normal(size=500)
    D = rng.normal(size=(500, 1))
    w, lab = multi_decoy_competition(t, D, 1, 1, seed=1)
    np.testing.assert_allclose(w, np.maximum(t, D[:, 0]))
    np.testing.assert_array_equal(lab, np.where(t > D[:, 0], 1, -1))
    q = multi_decoy_qvalues(w, lab, 1, 1)
    np.testing.assert_allclose(q, tdc_qvalues(w, lab < 0))


def test_multi_decoy_rank_labels():
    t = np.array([5.0, 3.0, 1.0, 2.5])
    D = np.array([[1.0, 2.0, 0.0], [4.0, 2.0, 1.0], [4.0, 3.0, 2.0], [3.0, 2.0, 4.0]])
    w, lab = multi_decoy_competition(t, D, 2, 2, seed=1)
    # ranks: 1, 2, 4, 3 -> target, target, decoy, decoy (d+1 = 4, k_t = k_d = 2)
    np.testing.assert_array_equal(lab, [1, 1, -1, -1])
    np.testing.assert_allclose(w, [5, 4, 4, 4])
    with pytest.raises(ValueError):
        multi_decoy_competition(t, D, 3, 2)


def test_hash_tiebreak_is_order_free():
    from vulcan.stats._util import hash_tiebreak

    rng = np.random.default_rng(131)
    n = 2000
    g = rng.integers(0, 300, n)
    s = np.round(rng.normal(size=n), 0)  # many ties
    seq = np.array([f"PEPT{i}IDEK" for i in range(n)])
    w = compete(g, s, tiebreak=hash_tiebreak(seq, 1001))
    p = rng.permutation(n)
    wp = compete(g[p], s[p], tiebreak=hash_tiebreak(seq[p], 1001))
    np.testing.assert_array_equal(np.sort(seq[w]), np.sort(seq[p][wp]))
    ti = hash_tiebreak(np.arange(5), 3)
    assert ti.dtype == np.uint64 and len(set(ti.tolist())) == 5
    assert not np.array_equal(hash_tiebreak(np.arange(5), 3), hash_tiebreak(np.arange(5), 4))
