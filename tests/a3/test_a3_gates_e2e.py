"""End-to-end synthetic-label run of scripts/a3/gates_a3 (A3.1.10): synthetic chimeric mixtures (host, A, B, FS, offsets),
A-syn tables, labels, V12 and G3 inputs -> mixture_units (the production competition code: stage-1 winners, CONGA
neighbour filter, double competition, per-stratum q_sep, per-pipeline tau) -> concat_units -> evaluate_v2sep / evaluate_ce_s
-> overall.  The n_min of the power-sized rows are scaled down in a COPY of the config (plumbing test, not power)."""
import copy

import numpy as np
import pytest

import gates_a3 as G
from vulcan.gates import a2_stats as A2S
from vulcan.gates import a3_config
from vulcan.gates import a3_stats as S

CFG = a3_config.load()


def small_cfg():
    """Real config except: bootstrap B (speed), V6 n_min (its power needs ~450k winners per class; the clean run uses
    non-overlapping ion sets so the filter keeps everything and V6 is exact) and V8 n_min (g2_sample_size ~108k)."""
    c = copy.deepcopy(CFG)
    c["rows"]["V1"]["bootstrap_B"] = 500
    c["rows"]["V6"]["n_min_per_label_per_bin"] = 200
    c["rows"]["V8"]["n_min"] = 2000
    return c


class Ions:
    """overlap=True: 20 random bins of 3,000 per peptidoform (random collisions -> neighbour-filter drops);
    overlap=False: disjoint bins (the filter keeps every winner)."""

    def __init__(self, overlap):
        self.overlap = overlap
        self.c = {}

    def sets(self, rows):
        out = []
        for r in rows:
            r = int(r)
            if r not in self.c:
                self.c[r] = frozenset(np.random.default_rng(r).integers(0, 3000, 20).tolist() if self.overlap
                                      else range(20 * r, 20 * r + 20))
            out.append(self.c[r])
        return out


N_PAIRS = 300_000


def make_mixture(g, leak=0.0, n_host=8000, n_A=9000, n_B=5000, n_FS=600, n_OFF=400, K=5, projects=12):
    prov = np.r_[np.zeros(n_host, int), np.full(n_OFF, 1), np.full(n_A, 20), np.full(n_B, 21), np.full(n_FS, 22)]
    n = prov.size
    project = np.where(prov == 0, "host", np.array([f"PXD{k:06d}" for k in g.integers(0, projects, n)], dtype=object))
    ng = np.where(np.isin(prov, [20, 22]), 1, g.choice([1, 2, 3], size=n, p=[0.45, 0.35, 0.2]))
    rows = {k: [] for k in ("spec_idx", "group", "chim_rank", "pf_row", "hyperscore", "v2", "v2_noRT")}
    for s in range(n):
        for gr in range(int(ng[s])):
            pairs = g.integers(0, N_PAIRS, K)
            side = g.random(K) < 0.5
            pf = 2 * pairs + side
            hs = g.gumbel(0.3 * gr, 1.0, K)
            if prov[s] == 0 and g.random() < (0.7 if gr == 0 else 0.5):   # host: real identifications, all ranks
                pf[0] = 2 * pairs[0]
                hs[0] = g.normal(10.0, 1.5)
            v2 = 0.8 * hs + 0.6 * g.normal(size=K)
            if prov[s] != 0:
                v2 = v2 + leak * (~side)
            for k, v in (("spec_idx", np.full(K, s)), ("group", np.full(K, gr)), ("chim_rank", np.full(K, gr + 1)),
                         ("pf_row", pf), ("hyperscore", hs), ("v2", v2), ("v2_noRT", v2 + 0.01 * g.normal(size=K))):
                rows[k].append(v)
    chim = {k: np.concatenate(v) for k, v in rows.items()}
    chim["pf_id"] = chim["pf_row"].astype(np.int64)
    a = np.flatnonzero(prov == 20)
    ar = {k: [] for k in ("spec_idx", "group", "chim_rank", "pf_row", "hyperscore", "v2", "v2_noRT")}
    for s in a:
        for rk in (2, 3):
            pairs = g.integers(0, N_PAIRS, K)
            side = g.random(K) < 0.5
            hs = g.gumbel(0.3 * (rk - 1), 1.0, K)
            v2 = 0.8 * hs + 0.6 * g.normal(size=K) + leak * (~side)
            for k, v in (("spec_idx", np.full(K, s)), ("group", np.full(K, rk - 1)), ("chim_rank", np.full(K, rk)),
                         ("pf_row", 2 * pairs + side), ("hyperscore", hs), ("v2", v2),
                         ("v2_noRT", v2 + 0.01 * g.normal(size=K))):
                ar[k].append(v)
    asyn = {k: np.concatenate(v) for k, v in ar.items()}
    asyn["pf_id"] = asyn["pf_row"].astype(np.int64)
    return chim, asyn, prov, project


def run(leak, seed=7, n_mix=12, overlap=False):
    g = np.random.default_rng(seed)
    is_decoy = (np.arange(2 * N_PAIRS) % 2).astype(bool)
    pair = np.arange(2 * N_PAIRS) // 2
    ions = Ions(overlap)
    parts = []
    for _ in range(n_mix):
        chim, asyn, prov, project = make_mixture(g, leak)
        asyn["logit_fn"] = lambda w: g.normal(size=len(w))
        parts.append(G.mixture_units(chim, asyn, is_decoy, pair, prov, project, np.ones(prov.size, bool), ions,
                                     logit_fn=lambda w: g.normal(size=len(w))))
    D = G.concat_units(parts)
    n = 260_000
    bt, bd = g.normal(size=n), g.normal(size=n)
    v12 = {"best_t": bt, "best_d": bd, "kept_t": g.random(n) >= 0.0954, "kept_d": g.random(n) >= 0.0954}
    X = [g.normal(size=200_000) for _ in range(2)]
    g3a = A2S.g3_gate_a2(X[0], X[1], np.where(g.random(200_000) < 0.999, X[0], -np.inf),
                         np.where(g.random(200_000) < 0.999, X[1], -np.inf))
    g3 = {"A": g3a, "AB_all": g3a, "AB_secondary": g3a}
    return D, G.evaluate_v2sep(D, small_cfg(), g3, copy.deepcopy(g3a), v12), parts


@pytest.fixture(scope="module")
def clean():
    return run(0.0)


def test_e2e_structure(clean):
    D, r, parts = clean
    U = D["units"]
    assert set(np.unique(U["pool"])) == {"A", "B", "FS", "OFF"}
    assert (U["rank"][U["pool"] == "A"] == 1).all() and (U["rank"][U["pool"] == "FS"] == 1).all()
    assert (U["rank"][U["pool"] == "B"] >= 2).any()
    assert set(np.unique(D["asyn"]["rank"])) == {2, 3}
    th = parts[0]["thresholds"]["v2"]
    assert set(th["rank_to_merged"]) == {1, 2, 3} and np.isfinite(th["tau"][th["rank_to_merged"][1]])
    assert D["v6"]["v2"]["kept"].size > 0 and D["v6"]["v2"]["kept"].all()          # disjoint ions: nothing dropped


def test_e2e_clean_scorer_passes_both_tiers(clean):
    _, r, _ = clean
    assert r["V3"]["byte_identical_to_v1"] is True
    bad = {k: r[k]["verdict"] for k in CFG["tiers"]["tier1"] + CFG["tiers"]["tier2"] if r[k]["verdict"] != S.PASS}
    assert r["tier1"]["outcome"] == S.PASS, (r["tier1"], bad)
    assert r["tier2"]["outcome"] == S.PASS, (r["tier2"], bad)
    assert r["licensed"]["licensed_list"].startswith("full")
    assert r["V1"]["n_tail"] >= 200 and r["V13"]["n_tail"] >= 200


def test_e2e_v3_mismatch_is_not_evaluable_implementation(clean):
    D, _, _ = clean
    g = np.random.default_rng(1)
    X = [g.normal(size=50_000) for _ in range(2)]
    g3 = A2S.g3_gate_a2(X[0], X[1], X[0], X[1], n_min=100)
    other = dict(g3, n_null=g3["n_null"] + 1)
    r = G.evaluate_v2sep(D, small_cfg(), {"A": g3, "AB_all": g3, "AB_secondary": g3}, other, None)
    assert r["V3"]["verdict"] == S.NE_IMPL and r["V12"]["verdict"] == S.NE
    assert r["tier1"]["outcome"] == S.NE and r["licensed"]["licensed_list"] == "none"
    real = copy.deepcopy(CFG)                     # real n_min for V6 -> not evaluable at this size (blocks tier 1)
    real["rows"]["V1"]["bootstrap_B"] = 500
    r2 = G.evaluate_v2sep(D, real, {"A": g3, "AB_all": g3, "AB_secondary": g3}, g3, None)
    assert r2["V6.class1"]["verdict"] == S.NE and r2["tier1"]["outcome"] in (S.NE, S.FAIL)


def test_e2e_leaking_scorer_fails_tier1():
    D, r, _ = run(1.5, seed=11, n_mix=3, overlap=True)
    assert 0 < 1 - D["v6"]["v2"]["kept"].mean() < 0.5                       # random ion collisions -> filter drops
    assert r["tier1"]["outcome"] == S.FAIL
    assert "V1" in r["tier1"]["failed_rows"] and r["licensed"]["licensed_list"] == "none"


def test_e2e_ce_s_and_overall(clean):
    g = np.random.default_rng(3)
    n = 260_000
    pool = np.array(["A"] * 200_000 + ["B"] * 59_900 + ["FS"] * 100, dtype=object)
    proj = np.array([f"PXD{k:06d}" for k in g.integers(0, 12, n)], dtype=object)
    wt = g.random(n) < 0.5
    hs = g.gumbel(size=n)
    P = {"hyperscore": {"win_t": wt, "score": hs, "tail1": hs >= np.quantile(hs, 0.99), "logit": g.normal(size=n)}}
    ce = 0.8 * hs + 0.6 * g.normal(size=n)
    P["ce_s"] = {"win_t": np.where(g.random(n) < 0.95, wt, g.random(n) < 0.5), "score": ce,
                 "tail1": ce >= np.quantile(ce, 0.99), "logit": g.normal(size=n)}
    bt, bd = g.normal(size=n), g.normal(size=n)
    v12 = {"best_t": bt, "best_d": bd, "kept_t": g.random(n) >= 0.0954, "kept_d": g.random(n) >= 0.0954}
    g3 = {"A": {"verdict": "pass"}, "AB": {"verdict": "pass"}}
    r = G.evaluate_ce_s({"pool": pool, "project": proj}, P, small_cfg(), g3, v12)
    assert r["verdict"]["outcome"] == S.PASS, r["verdict"]
    assert r["licenses"].startswith("nothing")
    ov = G.overall(clean[1], r)
    assert ov["v2_sep"]["tier1"]["outcome"] == S.PASS and ov["ce_s"]["outcome"] == S.PASS


def test_v12_inputs_best_and_kept():
    spec_a = np.array([10, 11, 12])
    cs = np.array([10, 10, 10, 11, 11, 12])
    sc = np.array([1.0, 3.0, 2.0, 5.0, 4.0, 0.5])
    dec = np.array([False, False, True, True, False, True])
    k5 = np.array([True, False, True, True, True, False])
    r = G.v12_inputs(spec_a, cs, sc, dec, k5)
    assert r["best_t"].tolist() == [3.0, 4.0, -np.inf] and r["kept_t"].tolist() == [False, True, False]
    assert r["best_d"].tolist() == [2.0, 5.0, 0.5] and r["kept_d"].tolist() == [True, True, False]


def test_v12_inputs_ties_count_by_score_equality():
    # two target candidates tie at the best score; only the second is in the top-5 -> kept (A2.8 score equality)
    r = G.v12_inputs(np.array([1]), np.array([1, 1, 1]), np.array([2.0, 2.0, 1.0]), np.array([False, False, True]),
                     np.array([False, True, True]))
    assert r["kept_t"].tolist() == [True] and r["kept_d"].tolist() == [True]


def test_noRT_v6_is_its_own_pipeline(clean):
    D, _, _ = clean
    assert set(D["v6"]) == {"v2", "v2_noRT"}
    assert D["v6"]["v2"]["score"].size == D["v6"]["v2_noRT"]["score"].size
    assert not np.array_equal(D["v6"]["v2"]["score"], D["v6"]["v2_noRT"]["score"])
