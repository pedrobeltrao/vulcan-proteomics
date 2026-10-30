import numpy as np

from vulcan.stats.protein import (classic_protein_fdr, picked_group_fdr, picked_protein_fdr_simple, score_groups,
                                  subset_grouping)


def test_subset_grouping_hand():
    pep2prot = {"p1": ["A", "B"], "p2": ["A", "B"], "p3": ["A"], "p4": ["C"], "p5": ["C", "D"], "p6": ["D"]}
    groups = subset_grouping(pep2prot)
    assert sorted(map(sorted, groups)) == [["A", "B"], ["C"], ["D"]]
    scores, leading = score_groups(groups, pep2prot, {k: float(i) for i, k in enumerate(pep2prot)})
    gi = {tuple(sorted(g)): i for i, g in enumerate(groups)}
    assert scores[gi[("A", "B")]] == 2.0       # p1..p3 unique to the group, best = p3
    assert scores[gi[("C",)]] == 3.0           # p5 shared C/D -> discarded
    assert scores[gi[("D",)]] == 5.0
    assert leading[gi[("A", "B")]] == ("A",)   # A has the most unique peptides


def _db(extra=None):
    pep2prot = {"t1": ["A"], "t2": ["A"], "t3": ["B"], "d1": ["rA"], "d2": ["rB"]}
    pair = {"A": "A", "rA": "A", "B": "B", "rB": "B"}
    isd = {"A": False, "B": False, "rA": True, "rB": True}
    return pep2prot, pair, isd


def test_picked_competition_hand():
    p2p, pair, isd = _db()
    score = {"t1": 5.0, "t2": 4.0, "t3": 2.0, "d1": 6.0, "d2": 1.0}
    r = picked_group_fdr(p2p, score, pair, isd, rescue=False, seed=11)
    kept = dict(zip([g[0] for g in r.groups], r.score))
    assert kept == {"rA": 6.0, "B": 2.0}       # A loses to rA; rB loses to B
    np.testing.assert_allclose(r.q, [1.0, 1.0])
    c = classic_protein_fdr(p2p, score, pair, isd, seed=11)
    assert sorted(g[0] for g in c.groups) == ["A", "B", "rA", "rB"]


def test_rescued_subsets_hand():
    pep2prot = {"b1": ["B"], "x1": ["X1"], "x2": ["X2"], "x3": ["X3"], "c1": ["C"],
                "s1": ["A", "A-2"], "s2": ["A", "A-2"], "uA": ["A"], "uA2": ["A-2"],
                "rb1": ["rB"], "ruA": ["rA"], "rz1": ["rZ"]}
    score = {"b1": 10.0, "x1": 9.5, "x2": 9.4, "x3": 9.3, "c1": 6.5, "s1": 8.0, "s2": 7.0, "uA": 0.5, "uA2": 0.4,
             "rb1": 0.6, "ruA": 0.3, "rz1": 0.45}
    prots = ["B", "X1", "X2", "X3", "C", "A", "A-2", "Z"]
    pair = {**{p: p for p in prots}, **{"r" + p: p for p in prots}}
    isd = {**{p: False for p in prots}, **{"r" + p: True for p in prots}}
    r0 = picked_group_fdr(pep2prot, score, pair, isd, rescue=False, seed=1)
    assert ("A",) in r0.groups and ("A-2",) in r0.groups     # split by low-scoring unique peptides
    r = picked_group_fdr(pep2prot, score, pair, isd, rescue=True, rescue_alpha=0.25, seed=1)
    merged = [g for g in r.groups if set(g) == {"A", "A-2"}]
    assert len(merged) == 1
    i = r.groups.index(merged[0])
    assert r.score[i] == 8.0 and not r.is_decoy[i] and not r.is_obsolete[i]
    assert ("rA",) not in r.groups                            # picked out by A
    assert ("rB",) not in r.groups                            # picked out by B
    assert ("rZ",) in r.groups                                # unpaired decoy competes
    assert merged[0] in r.accepted(0.5)


def test_picked_simple_vectorised():
    s, isd, q, pi, bb = picked_protein_fdr_simple([5.0, 3.0], [4.0, 6.0], seed=1)
    assert sorted(zip(s.tolist(), isd.tolist())) == [(5.0, False), (6.0, True)]
    s, isd, q, pi, bb = picked_protein_fdr_simple([5.0, 3.0], [4.0, 6.0], picked=False)
    assert len(s) == 4
    s, isd, q, pi, bb = picked_protein_fdr_simple([5.0, -np.inf], [4.0, -np.inf], seed=1)
    assert len(s) == 1
