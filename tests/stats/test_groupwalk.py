import numpy as np
import pytest

from vulcan.stats.groupwalk import (check_pair_invariant, group_walk, group_walk_batch, label_blind_pair_counts,
                                    merge_adjacent, merge_rare_to_other, pool_chimeric_rank)
from vulcan.stats.tdc import tdc_qvalues


def test_group_walk_hand_K1():
    # group 0 ascending: 1(D), 3(T); group 1: 2(T), 4(T)
    s = np.array([1.0, 3.0, 2.0, 4.0])
    d = np.array([1, 0, 0, 0], bool)
    g = np.array([0, 0, 1, 1])
    q, fr = group_walk(s, d, g, K=1, tie_break="first", return_frontier=True)
    np.testing.assert_allclose(q, [2 / 3, 1 / 3, 1 / 3, 1 / 3])
    np.testing.assert_array_equal(fr, [[0, 0], [1, 0], [1, 1], [2, 1], [2, 2]])


def test_group_walk_factor():
    s = np.array([1.0, 3.0, 2.0, 4.0])
    d = np.array([1, 0, 0, 0], bool)
    g = np.array([0, 0, 1, 1])
    q = group_walk(s, d, g, K=1, tie_break="first", factor=2.0)
    np.testing.assert_allclose(q, [1.0, 2 / 3, 2 / 3, 2 / 3])


def test_single_group_equals_tdc():
    rng = np.random.default_rng(11)
    s = rng.normal(size=400) + np.r_[np.full(150, 2.5), np.zeros(250)]
    d = rng.random(400) < 0.4
    q = group_walk(s, d, np.zeros(400, int), seed=1)
    np.testing.assert_allclose(q[~d], tdc_qvalues(s, d)[~d])


def test_batch_matches_scalar():
    rng = np.random.default_rng(23)
    R, sizes, K = 15, (60, 80, 50), 5
    G, L = len(sizes), max(sizes)
    labs = np.zeros((R, G, L), bool)
    nul = np.zeros((R, G, L), bool)
    alphas = (0.05, 0.1, 0.2)
    ref = {a: ([], []) for a in alphas}
    for r in range(R):
        S, Dl, Nl, Gl = [], [], [], []
        for g, m in enumerate(sizes):
            nt = m // 3
            t = rng.normal(size=m); t[:nt] += 2.5 + g * 0.5
            dd = rng.normal(size=m)
            w = np.maximum(t, dd); isd = dd > t
            isn = np.arange(m) >= nt
            o = np.argsort(w)
            labs[r, g, :m] = isd[o]; nul[r, g, :m] = isn[o]
            S.append(w); Dl.append(isd); Nl.append(isn & ~isd); Gl.append(np.full(m, g))
        S, Dl, Nl, Gl = map(np.concatenate, (S, Dl, Nl, Gl))
        q = group_walk(S, Dl, Gl, K=K, tie_break="first")
        for a in alphas:
            acc = (q <= a) & ~Dl
            ref[a][0].append(acc.sum()); ref[a][1].append((acc & Nl).sum())
    out = group_walk_batch(labs, nul, np.asarray(sizes), alphas, K=K)
    for a in alphas:
        np.testing.assert_array_equal(out[a][0], ref[a][0])
        np.testing.assert_array_equal(out[a][1], ref[a][1])


def test_merge_adjacent_hand():
    order = ["a", "b", "c", "d", "e"]
    counts = {"a": 500, "b": 50, "c": 120, "d": 300, "e": 90}
    m = merge_adjacent(order, counts, 200)
    assert m["a"] == ("a",)
    assert m["b"] == m["c"] == m["d"] == m["e"] == ("b", "c", "d", "e")
    # equal-size neighbours: merge with the preceding one
    m = merge_adjacent([1, 2, 3], {1: 250, 2: 10, 3: 250}, 200)
    assert m[1] == m[2] == (1, 2) and m[3] == (3,)
    # everything too small -> one pooled block (pooled TDC fallback)
    m = merge_adjacent([1, 2, 3], {1: 10, 2: 10, 3: 10}, 200)
    assert len(set(m.values())) == 1
    # deterministic
    assert merge_adjacent(order, counts, 200) == merge_adjacent(order, counts, 200)


def test_merge_rare_and_rank_pooling():
    m = merge_rare_to_other({"phospho": 5000, "acetyl": 150, "GG": 30, "ox": 900}, 200)
    assert m == {"phospho": "phospho", "acetyl": "other", "GG": "other", "ox": "ox"}
    np.testing.assert_array_equal(pool_chimeric_rank([1, 2, 3, 4, 7]), [1, 2, 3, 3, 3])


def test_label_blind_counts_and_pair_invariant():
    pair = np.array([1, 1, 2, 3, 3, 4])
    strat = np.array(["x", "x", "x", "y", "y", "y"])
    assert label_blind_pair_counts(pair, strat) == {"x": 2, "y": 2}
    with pytest.raises(ValueError):
        check_pair_invariant(np.array([1, 1]), np.array(["x", "y"]))
