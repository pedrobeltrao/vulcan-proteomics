import numpy as np
import pytest

from vulcan.stats.entrapment import (Replicate, combine_arms, combined, combined_counts, evaluate_arm, lower_bound,
                                     paired, paired_counts, protein_group_is_target, protein_paired_counts, verdict)


def _hand_table():
    # pair, is_entrapment, score, accepted
    rows = [
        ("P1", 0, 10, 1), ("P1", 1, 12, 1),   # both accepted, E > T
        ("P2", 0, 9, 1), ("P2", 1, 8, 1),     # both accepted, E < T
        ("P3", 0, 1, 0), ("P3", 1, 7, 1),     # E accepted, T not
        ("P4", 0, 6, 1),
        ("P5", 0, 11, 1), ("P5", 1, 5, 1),    # both, E < T
        ("P6", 0, 4, 1), ("P6", 1, 2, 0),
        ("P7", 0, 3, 1), ("P8", 0, 3.5, 1),
    ]
    p, e, s, a = zip(*rows)
    return np.array(p), np.array(e, bool), np.array(s, float), np.array(a, bool)


def test_paired_counts_hand():
    c = paired_counts(*_hand_table())
    n_t, n_e, n_eo, n_egt = c.totals()
    assert (n_t, n_e, n_eo, n_egt) == (7, 4, 1, 1)
    assert paired(n_t, n_e, n_eo, n_egt) == pytest.approx(7 / 11)
    assert combined(n_t, n_e, 1.0) == pytest.approx(8 / 11)
    assert combined(n_t, n_e, 2.0) == pytest.approx(6 / 11)
    assert lower_bound(n_t, n_e) == pytest.approx(4 / 11)
    assert lower_bound(0, 0) == 0.0


def _fdrbench_port(pair, ent, score, acc):
    """Literal port of FDRBench FDPCalc.java (peptide level) for one cutoff."""
    order = np.argsort(-score, kind="stable")
    order = order[acc[order]]  # the list up to the cutoff i
    n_ent = n_t = 0
    all_hits = set()
    for k in order:
        all_hits.add((pair[k], ent[k]))
        if ent[k]:
            n_ent += 1
        else:
            n_t += 1
    cur = set()
    n_p_t_s = n_p_s_t = 0
    for k in order[::-1]:  # from low to high confidence
        cur.add((pair[k], ent[k]))
        if ent[k]:
            if (pair[k], False) in cur:
                n_p_t_s += 1
            elif (pair[k], False) not in all_hits:
                n_p_s_t += 1
    return (n_ent + 2.0 * n_p_t_s + n_p_s_t) / (n_ent + n_t)


def test_paired_matches_fdrbench_port():
    rng = np.random.default_rng(11)
    for _ in range(30):
        m = 300
        ts = rng.normal(size=m) + np.r_[np.full(100, 2.0), np.zeros(m - 100)]
        es = rng.normal(size=m)
        pair = np.r_[np.arange(m), np.arange(m)]
        ent = np.r_[np.zeros(m, bool), np.ones(m, bool)]
        s = np.r_[ts, es]
        acc = s >= np.quantile(s, 0.8)
        c = paired_counts(pair, ent, s, acc)
        assert paired(*c.totals()) == pytest.approx(_fdrbench_port(pair, ent, s, acc))


def test_verdict_table():
    assert verdict(0.01, 0.008, 0.012, 0.003, 5000) == "supported"
    assert verdict(0.01, 0.008, 0.016, 0.003, 5000) == "supported-noisy"
    assert verdict(0.01, 0.03, 0.04, 0.012, 5000) == "violated"
    assert verdict(0.01, 0.012, 0.02, 0.005, 5000) == "inconclusive"
    assert verdict(0.01, 0.0, 0.0, 0.0, 1999) == "not-evaluable"
    assert verdict(0.01, 0.0, 0.0, 0.0, 2000) == "supported"


def test_evaluate_arm_and_bootstrap():
    rng = np.random.default_rng(23)
    reps = []
    for j in range(3):
        m = 20000
        ts = rng.normal(size=m) + np.r_[np.full(8000, 3.0), np.zeros(m - 8000)]
        es = rng.normal(size=m)
        pair = np.r_[np.arange(m), np.arange(m)]
        ent = np.r_[np.zeros(m, bool), np.ones(m, bool)]
        s = np.r_[ts, es]
        acc = s >= 3.3
        reps.append(Replicate(paired_counts(pair, ent, s, acc), entrap_seed=2001 + j, decoy_seed=1001))
    res = evaluate_arm(reps, 0.01, "paired", n_boot=300, seed=11)
    manual = np.mean([paired(*r.counts.totals()) for r in reps])
    assert res["mean_fdp"] == pytest.approx(manual)
    lo, hi = res["ci95"]
    assert lo < res["mean_fdp"] < hi
    assert res["gated"] is True
    assert res["verdict"] in ("supported", "supported-noisy", "inconclusive", "violated", "not-evaluable")
    res2 = evaluate_arm(reps, 0.001, "lower_bound", n_boot=100)
    assert res2["gated"] is False
    # combined arm with item-level bootstrap units
    c = combined_counts(np.arange(10), np.r_[np.zeros(8, bool), np.ones(2, bool)], np.ones(10, bool))
    assert c.totals()[:2] == (8, 2)


def test_combine_arms():
    assert combine_arms({"P": "supported", "C": "supported"}, "static_low_capacity") == "supported"
    assert combine_arms({"P": "supported", "C": "supported-noisy"}, "static_low_capacity") == "supported-noisy"
    assert combine_arms({"P": "supported", "C": "supported"}, "sequence_capacity") == "not-evaluable"
    assert combine_arms({"P": "supported", "C": "supported", "M": "inconclusive"},
                        "sequence_capacity") == "inconclusive"
    assert combine_arms({"P": "supported", "C": "violated"}, "fixed") == "violated"
    assert combine_arms({"C": "supported"}, "fixed") == "supported"


def test_protein_counting_rule_and_paired():
    is_e = {"A": False, "B": False, "E_A": True, "E_B": True}
    groups = [["A", "E_B"], ["E_A"], ["B"]]
    np.testing.assert_array_equal(protein_group_is_target(groups, is_e), [True, False, True])
    groups = [["A"], ["E_A"], ["E_B"], ["B"]]
    score = [10.0, 12.0, 5.0, 3.0]
    acc = [True, True, True, False]
    c = protein_paired_counts(groups, score, acc, is_e, {"E_A": "A", "E_B": "B"})
    assert c.totals() == (1, 2, 1, 1)
