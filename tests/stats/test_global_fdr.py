"""L6 global FDR (vulcan.stats.global_fdr): Monte-Carlo control vs N and the naive-merge inflation it fixes."""
import numpy as np

from vulcan.stats import global_fdr as G
from vulcan.stats._util import hash_tiebreak


def _sim(n_runs, rng, n_true=1500, n_null=30000, p_detect=0.6, n_null_spec=4000):
    """Pooled PSM winners of n_runs runs. Units: 0..n_true-1 present targets; null pairs n_true..; unit id
    2*k = target of pair k, 2*k+1 = its decoy. Null spectra hit a random null pair, target or decoy with p 1/2."""
    run, unit, pair, score, dec = [], [], [], [], []
    for r in range(n_runs):
        det = np.flatnonzero(rng.random(n_true) < p_detect)
        run.append(np.full(len(det), r)); pair.append(det); unit.append(2 * det)
        score.append(rng.normal(3.0, 1.0, len(det))); dec.append(np.zeros(len(det), bool))
        k = n_true + rng.integers(0, n_null, n_null_spec)
        d = rng.random(n_null_spec) < 0.5
        run.append(np.full(n_null_spec, r)); pair.append(k); unit.append(2 * k + d)
        score.append(rng.normal(0.0, 1.0, n_null_spec)); dec.append(d)
    c = lambda x: np.concatenate(x)  # noqa: E731
    return c(run), c(unit), c(pair), c(score), c(dec), n_true


def _fdp(units, n_true):
    t = units // 2
    return float((t >= n_true).sum() / max(len(units), 1))


def test_global_controls_naive_inflates():
    rng = np.random.default_rng(11)
    alpha = 0.01
    res = {}
    for n_runs in (1, 20):
        g_fdp, n_fdp = [], []
        for _ in range(30):
            run, unit, pair, score, dec, n_true = _sim(n_runs, rng)
            tb = hash_tiebreak(unit, 11)
            glob = G.global_competition(run, unit, pair, score, dec, tb)
            per = G.run_competitions(run, unit, pair, score, dec, tb)
            g_fdp.append(_fdp(glob.unit[glob.accepted(alpha)], n_true))
            u, s, nr = G.naive_union(per, alpha)
            n_fdp.append(_fdp(u, n_true))
        res[n_runs] = (np.mean(g_fdp), np.mean(n_fdp))
    # global: mean FDP <= alpha (MC slack); naive at N=20 clearly inflated
    assert res[1][0] <= alpha * 1.5
    assert res[20][0] <= alpha * 1.5
    assert res[20][1] > 3 * alpha
    assert res[20][1] > res[1][1]


def test_run_competitions_match_single_run_global():
    rng = np.random.default_rng(23)
    run, unit, pair, score, dec, _ = _sim(3, rng, n_true=300, n_null=3000, n_null_spec=800)
    tb = hash_tiebreak(unit, 11)
    per = G.run_competitions(run, unit, pair, score, dec, tb)
    for r in range(3):
        m = run == r
        g = G.global_competition(run[m], unit[m], pair[m], score[m], dec[m], tb[m])
        pm = per.run == r
        a = np.argsort(per.unit[pm]); b = np.argsort(g.unit)
        assert np.array_equal(per.unit[pm][a], g.unit[b])
        assert np.allclose(per.q[pm][a], g.q[b])


def test_permutation_invariance():
    rng = np.random.default_rng(47)
    run, unit, pair, score, dec, _ = _sim(4, rng, n_true=200, n_null=2000, n_null_spec=500)
    score = np.round(score, 1)              # force ties
    tb = hash_tiebreak(unit, 11)
    g1 = G.global_competition(run, unit, pair, score, dec, tb)
    p = rng.permutation(len(run))
    g2 = G.global_competition(run[p], unit[p], pair[p], score[p], dec[p], tb[p])
    a, b = np.argsort(g1.unit), np.argsort(g2.unit)
    assert np.array_equal(g1.unit[a], g2.unit[b]) and np.array_equal(g1.q[a], g2.q[b])


def test_comparability_identical_runs_small_ks():
    rng = np.random.default_rng(83)
    run = np.repeat(np.arange(5), 2000)
    s = rng.normal(size=len(run))
    out = G.score_comparability(run, s, np.ones(len(run), bool))
    assert out["summary"]["ks_max"] < 0.06
    s2 = s + (run == 0) * 1.0
    out2 = G.score_comparability(run, s2, np.ones(len(run), bool))
    assert max(r["ks"] for r in out2["rows"] if r["run"] == "0") > 0.3
