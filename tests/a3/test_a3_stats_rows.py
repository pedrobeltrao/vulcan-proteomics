"""Synthetic pools for the new A3 rows (V6, V8, V9c, V12) and the tier logic (vulcan.gates.a3_stats)."""
import numpy as np
import pytest

from vulcan.gates import a3_stats as S


def _v6_pool(n, drop_t=0.118, drop_d=0.118, frac_sec=0.5, seed=0):
    g = np.random.default_rng(seed)
    t = g.random(n) < 0.5
    rank = np.where(g.random(n) < frac_sec, 2, 1)
    p = np.where(t, drop_t, drop_d)
    kept = g.random(n) >= p
    score = g.normal(size=n)
    return kept, t, score, rank


def test_v6_symmetric_passes_and_classes():
    kept, t, s, r = _v6_pool(1_200_000)
    out = S.nf_symmetry(kept, t, s, r, n_min=22_600)
    for c in ("class1", "class2"):
        assert out[c]["bins"]["top_decile"]["n_target"] >= 22_600 or out[c]["verdict"] == S.NE
    assert out["class1"]["verdict"] == S.PASS and out["class2"]["verdict"] == S.PASS
    assert abs(out["class1"]["drop_rate"] - 0.118) < 0.005


def test_v6_asymmetric_fails_and_small_not_evaluable():
    kept, t, s, r = _v6_pool(1_200_000, drop_t=0.10, drop_d=0.13)
    out = S.nf_symmetry(kept, t, s, r, n_min=22_600)
    assert out["class1"]["verdict"] == S.FAIL and out["class2"]["verdict"] == S.FAIL
    kept, t, s, r = _v6_pool(100_000)
    out = S.nf_symmetry(kept, t, s, r, n_min=22_600)
    assert out["class1"]["verdict"] == S.NE and out["class2"]["verdict"] == S.NE


def test_v6_top_decile_asymmetry_only_is_caught():
    g = np.random.default_rng(3)
    n = 1_500_000
    t = g.random(n) < 0.5
    s = g.normal(size=n)
    top = s >= np.quantile(s, 0.9)
    p = np.where(top & t, 0.08, 0.118)          # targets dropped less in the top decile only
    kept = g.random(n) >= p
    out = S.nf_symmetry_class(kept, t, s, n_min=22_600)
    assert out["bins"]["top_decile"]["pass"] is False and out["verdict"] == S.FAIL
    assert out["bins"]["rest"]["pass"] is True


def test_v12_truncation_symmetry():
    g = np.random.default_rng(5)
    n = 400_000
    bt, bd = g.normal(size=n), g.normal(size=n)
    kt, kd = g.random(n) >= 0.0954, g.random(n) >= 0.0954
    r = S.truncation_symmetry(bt, bd, kt, kd, n_min=18_800)
    assert r["verdict"] == S.PASS, r
    r = S.truncation_symmetry(bt, bd, kt, g.random(n) >= 0.14, n_min=18_800)
    assert r["verdict"] == S.FAIL
    r = S.truncation_symmetry(bt[:50_000], bd[:50_000], kt[:50_000], kd[:50_000], n_min=18_800)
    assert r["verdict"] == S.NE


def test_v12_missing_label_is_not_counted():
    bt = np.array([1.0, -np.inf, 2.0] * 1000)
    bd = np.array([0.5, 0.7, -np.inf] * 1000)
    r = S.truncation_symmetry(bt, bd, np.ones(3000, bool), np.ones(3000, bool), n_min=10)
    assert r["verdict"] in (S.PASS, S.NE)
    assert sum(r["n_target_by_decile"]) == 2000 and sum(r["n_decoy_by_decile"]) == 2000


def test_g2_original_n_min_and_g5_min_winners():
    g = np.random.default_rng(1)
    n = 50_000
    wf = g.random(n) < 0.5
    wl = np.where(g.random(n) < 0.95, wf, g.random(n) < 0.5)
    sf = g.normal(size=n)
    assert S.g2_original(wf, wl, sf, n_min=107_570)["verdict"] == S.NE
    assert S.g2_original(wf, wl, sf)["verdict"] in (S.PASS, S.FAIL)
    assert S.g5(g.normal(size=300), np.r_[np.ones(150, bool), np.zeros(150, bool)], min_winners=200)["verdict"] == S.NE
    assert S.g5(g.normal(size=2000), g.random(2000) < 0.5, min_winners=200)["verdict"] == S.PASS
    w = g.random(4000) < 0.5
    assert S.g5(np.where(w, 1.0, 0.0) + g.normal(size=4000), w)["verdict"] == S.FAIL


def test_tier_outcome_rules():
    t = S.tier_outcome({"V1": S.PASS, "V2": S.PASS, "V10": S.NE}, fs_keys=("V10",))
    assert t["outcome"] == S.PASS                                   # FS not evaluable never blocks
    assert S.tier_outcome({"V1": S.PASS, "V10": S.FAIL}, fs_keys=("V10",))["outcome"] == S.FAIL
    assert S.tier_outcome({"V1": S.NE, "V2": S.PASS})["outcome"] == S.NE
    assert S.tier_outcome({"V1": S.NE, "V2": S.FAIL})["outcome"] == S.FAIL     # a fail is final, even with NE
    assert S.tier_outcome({"V3": S.NE_IMPL})["outcome"] == S.NE
    assert S.tier_outcome({"V1": S.NE_INVALID})["outcome"] == S.NE
    lic = S.licensed({"outcome": S.PASS}, {"outcome": S.NE})
    assert "stratum-1" in lic["licensed_list"]
    assert S.licensed({"outcome": S.PASS}, {"outcome": S.PASS})["licensed_list"].startswith("full")
    assert S.licensed({"outcome": S.NE}, {"outcome": S.PASS})["licensed_list"] == "none"


@pytest.mark.parametrize("p,lo,hi", [(0.118, 20_000, 26_000), (0.0954, 16_500, 21_500)])
def test_n_min_symmetric_matches_analytic(p, lo, hi):
    r = S.n_min_symmetric(p, n_rep=3000)
    assert lo <= r["n_min"] <= hi and r["pass_prob_at_n_min"] >= 0.90
