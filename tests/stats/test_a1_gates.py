"""Amendment A1 gate statistics (vulcan.gates): tail-G2 MC (fires on a biased simulation, passes a fair one),
homolog filter vs brute force, RT quantile map, G1/G3/G4/G5 helpers, CI conventions of evaluate_arm."""
import random

import numpy as np
import pytest

from vulcan.gates import stats as GS
from vulcan.gates.homolog import HomologFilter, hamming_le2_bruteforce, shares_kmer_bruteforce
from vulcan.gates.nulls import quantile_map


def test_cp_upper():
    assert GS.cp_upper(0, 0) == 1.0
    u = GS.cp_upper(100, 200)
    assert 0.555 < u < 0.565        # exact CP; normal approx 0.558
    assert GS.cp_upper(200, 200) == 1.0


def test_tail_g2_not_evaluable_below_200():
    rng = np.random.default_rng(1)
    n = 10000
    t = rng.random(n) < 0.01                      # ~100 tail spectra
    tg = rng.random(n) < 0.5
    r = GS.tail_g2(t, tg, t, tg)
    assert r["n_tail"] < 200 and r["verdict"] == "not-evaluable"


def test_tail_g2_mc_fair_passes():
    # fair scorer: tail target fraction 0.5 for both pipelines, n_tail ~ 400
    r = GS.simulate_tail_g2(400, 40000, 0.01, 0.5, 0.5, overlap=0.9, seed=11)
    assert r["pass"] >= 0.95, r


@pytest.mark.parametrize("p_learned", [0.70, 0.75])
def test_tail_g2_mc_fires_on_biased_simulation(p_learned):
    # biased learned scorer: tail target fraction 0.70+ (fixed pipeline fair at 0.5): the gate must fire
    r = GS.simulate_tail_g2(400, 40000, 0.01, 0.5, p_learned, overlap=0.9, seed=23)
    assert r["fail"] >= 0.95, r


def test_tail_g2_mc_diff_clause_fires():
    # both pipelines biased (fixed 0.45 -> learned 0.62): the frac clause alone (0.60) is borderline, the
    # differential clause (+0.17 > 0.10) must fire
    r = GS.simulate_tail_g2(400, 60000, 0.01, 0.45, 0.62, overlap=0.9, seed=47)
    assert r["fail"] >= 0.9, r


def test_tail_g2_diff_variance_accounts_for_overlap():
    rng = np.random.default_rng(5)
    n = 50000
    t = rng.random(n) < 0.01
    tg = rng.random(n) < 0.5
    r_same = GS.tail_g2(t, tg, t, tg)
    assert r_same["diff"] == 0 and r_same["diff_se"] == 0


def test_homolog_filter_matches_bruteforce():
    rng = random.Random(3)
    AA = "ACDEFGHIKLMNPQRSTVWY"
    prots = ["".join(rng.choice(AA) for _ in range(rng.randint(40, 200))) for _ in range(30)]
    hf = HomologFilter(prots)
    peps = []
    for _ in range(300):
        L = rng.randint(7, 24)
        if rng.random() < 0.6:
            s = rng.choice(prots)
            a = rng.randint(0, len(s) - L)
            w = list(s[a:a + L])
            for _ in range(rng.randint(0, 4)):
                w[rng.randrange(L)] = rng.choice(AA)
            peps.append("".join(w))
        else:
            peps.append("".join(rng.choice(AA) for _ in range(L)))
    got = hf.near(peps)
    ref = np.array([hamming_le2_bruteforce(p, prots) or shares_kmer_bruteforce(p, prots) for p in peps])
    assert np.array_equal(got, ref)
    # I/L collapse
    hf2 = HomologFilter(["MKAAAILLLKR"])
    assert hf2.near(["AAALLLL"])[0]


def test_quantile_map():
    src = np.array([5.0, 10, 20, 40, 80])
    host = np.linspace(30, 200, 1001)
    m = quantile_map(src, host)
    assert np.all(np.diff(m) > 0) and m.min() >= 30 and m.max() <= 200
    big = np.random.default_rng(1).uniform(0, 100, 20000)
    mm = quantile_map(big, host)
    assert abs(np.median(mm) - np.median(host)) < 1.0


def test_g1_gate():
    rng = np.random.default_rng(2)
    s = rng.normal(size=200000)
    assert GS.g1_gate(s, rng.normal(size=s.size))["pass"]
    assert not GS.g1_gate(s, s + rng.normal(size=s.size) * 3)["pass"]


def test_g3_gate():
    rng = np.random.default_rng(3)
    n = 50000
    bt, bd = rng.normal(size=n), rng.normal(size=n)
    kt, kd = bt.copy(), bd.copy()
    kd[rng.random(n) < 0.05] = -np.inf             # decoys lost 5% -> asymmetric
    assert GS.g3_gate(bt, bd, bt, bd)["pass"]
    assert not GS.g3_gate(bt, bd, kt, kd)["pass"]


def test_g4_gate():
    rng = np.random.default_rng(4)
    runs = np.repeat(np.arange(6), 3000)
    st = np.tile(np.repeat([2, 3, 4], 1000), 6)
    s = rng.normal(size=runs.size) + 0.05 * rng.normal(size=6)[runs]
    assert GS.g4_gate(s, st, runs)["pass"]
    s2 = s + (st == 4) * 1.0
    assert not GS.g4_gate(s2, st, runs)["pass"]


def test_g5_auc_ci():
    rng = np.random.default_rng(6)
    s = rng.normal(size=20000)
    w = rng.random(20000) < 0.5
    r = GS.g5_gate(s, w)
    assert r["pass"] and r["all"]["ci"][0] < 0.5 < r["all"]["ci"][1]


def test_evaluate_arm_ci_conventions():
    from vulcan.stats.entrapment import Replicate, evaluate_arm, paired_counts
    rng = np.random.default_rng(7)
    reps = []
    for s in (2001, 2002, 2003):
        m = 3000
        pair = np.r_[np.arange(m), np.arange(m)]
        ent = np.r_[np.zeros(m, bool), np.ones(m, bool)]
        sc = rng.normal(size=2 * m)
        acc = np.r_[rng.random(m) < 0.8, rng.random(m) < 0.01]
        reps.append(Replicate(paired_counts(pair, ent, sc, acc), s, 1001))
    a = evaluate_arm(reps, 0.01, n_boot=300, ci_convention="a")
    b = evaluate_arm(reps, 0.01, n_boot=300, ci_convention="b")
    assert a["ci95"] == a["ci95_a"] and b["ci95"] == b["ci95_b"]
    assert a["ci95_a"] == b["ci95_a"]
    assert (a["ci95_a"][1] - a["ci95_a"][0]) >= (a["ci95_b"][1] - a["ci95_b"][0])
    with pytest.raises(ValueError):
        evaluate_arm(reps, 0.01, ci_convention="c")
