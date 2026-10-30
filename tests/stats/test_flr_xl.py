import numpy as np
import pytest

from vulcan.stats.flr import acceptor_counts, decoy_aa_flr, flr_curve, flr_strata
from vulcan.stats.xl import is_self_link, merge_level, xl_level_fdr, xl_qvalues

PEPS = ["PESTAK", "SAAYK", "LIGHTK"]


def test_acceptor_counts():
    assert acceptor_counts(PEPS, "A") == (5, 3)
    assert acceptor_counts(PEPS, "L") == (5, 2)   # I counted as L


def test_decoy_aa_flr_hand():
    r = decoy_aa_flr(["S", "A", "T"], [0, 1, 2], PEPS, "A")
    assert (r["n_pA"], r["n_pSTY"], r["N_STY"], r["N_A"]) == (1, 2, 5, 3)
    assert r["flr"] == pytest.approx(1 * (5 / 3) / 2)
    r = decoy_aa_flr(["S", "I", "T"], [0, 2, 2], PEPS, "L")
    assert (r["n_pA"], r["n_pSTY"], r["N_STY"], r["N_A"]) == (1, 2, 3, 2)
    assert r["flr"] == pytest.approx(0.75)
    assert decoy_aa_flr(["S"], [0], PEPS, "A")["flr"] == 0.0


def test_flr_curve_hand():
    f = flr_curve([0.9, 0.8, 0.7], ["S", "A", "T"], [0, 1, 2], PEPS, "A")
    np.testing.assert_allclose(f, [0.0, 5 / 6, 5 / 6])
    f = flr_curve([0.9, 0.8, 0.7], ["S", "A", "T"], [0, 1, 2], PEPS, "A", monotone=False)
    np.testing.assert_allclose(f, [0.0, 1.0, 5 / 6])


def test_flr_strata():
    st = flr_strata([1, 2, 1], [2, 0, 4], PEPS, [0, 1, 2], "A")
    np.testing.assert_array_equal(st["phospho_class"], ["mono", "multi", "mono"])
    # PESTAK pos2 S: T at 3 and A at 4 nearby; SAAYK pos0 S: A at 1; LIGHTK pos4 T: no acceptor within 2
    np.testing.assert_array_equal(st["close_site"], [True, True, False])


def test_xl_qvalues_hand():
    s = np.array([10, 9, 8, 7, 6, 5.0])
    c = np.array([0, 0, 1, 0, 2, 1])
    q1, q2 = xl_qvalues(s, c)
    np.testing.assert_allclose(q1, [0, 0, 0, 0, 0, 1 / 3])
    np.testing.assert_allclose(q2, [0, 0, 1 / 3, 1 / 3, 1 / 3, 2 / 3])
    q1p, _ = xl_qvalues(s, c, plus_one=True)
    np.testing.assert_allclose(q1p, [1 / 3] * 5 + [2 / 3])  # raw: 1, 1/2, 1, 2/3, 1/3, 2/3
    q1m, _ = xl_qvalues(s, c, min_tt=3)
    np.testing.assert_allclose(q1m, [0, 0, 0, 0, 0, 1 / 3])
    q1m, _ = xl_qvalues(s, c, min_tt=4)
    np.testing.assert_allclose(q1m, [1] * 6)


def test_xl_merge_and_self_hetero():
    base = {"rA": "A", "rB": "B"}
    pa = np.array(["A", "A", "A", "B", "rA", "rA"])
    pb = np.array(["A", "rA", "B", "A", "B", "rB"])
    da = np.array([0, 0, 0, 0, 1, 1], bool)
    db = np.array([0, 1, 0, 0, 0, 1], bool)
    s = np.array([9.0, 3.0, 5.0, 7.0, 4.0, 2.0])
    np.testing.assert_array_equal(is_self_link(pa, pb, base), [True, True, False, False, False, False])
    ent = merge_level("ppi", s, pa, pb, da, db)
    got = dict(zip(ent["key"].tolist(), zip(ent["score"].tolist(), ent["n_decoy"].tolist())))
    assert got == {"A~A": (9.0, 0), "A~rA": (3.0, 1), "A~B": (7.0, 0), "B~rA": (4.0, 1), "rA~rB": (2.0, 2)}
    out = xl_level_fdr("ppi", s, pa, pb, da, db, base_of=base)
    assert sorted(out["self"]["key"].tolist()) == ["A~A", "A~rA"]
    assert sorted(out["hetero"]["key"].tolist()) == ["A~B", "B~rA", "rA~rB"]
    h = out["hetero"]
    # hetero: A~B (7, TT), B~rA (4, TD), rA~rB (2, DD): TD-DD at 2 -> 0/1
    np.testing.assert_allclose(h["q_tddd"], [0, 0, 0])
    np.testing.assert_allclose(h["q_td"], [0, 1, 1])
    rp = merge_level("residue_pair", s, pa, pb, da, db, site_a=[1, 1, 2, 5, 3, 3], site_b=[7, 2, 5, 2, 1, 1])
    assert "A:2~B:5" in rp["key"].tolist() and (rp["key"] == "A:2~B:5").sum() == 1
