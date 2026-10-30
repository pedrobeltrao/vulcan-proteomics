"""Pure evaluation layer of the step-5 driver (scripts/a2/gates_a2.evaluate / overall) on synthetic pools."""
import numpy as np

import gates_a2 as G
from vulcan.gates import a2_config

CFG = a2_config.load()


def synth(seed=1, leak=0.0, nA=90_000, nB=30_000, nFS=40, nOFF=500):
    rng = np.random.default_rng(seed)
    n = nA + nB + nFS + nOFF
    pool = np.array(["A"] * nA + ["B"] * nB + ["FS"] * nFS + ["OFF"] * nOFF)
    proj = np.where(pool == "OFF", "host-offset", np.array([f"PXD{k:06d}" for k in rng.integers(0, 12, n)]))
    P = {}
    base_t = rng.random(n) < 0.5
    hs = rng.gumbel(size=n)
    for p in G.PIPES:
        wt = base_t.copy() if p == "fixed" else np.where(rng.random(n) < 0.95, base_t, rng.random(n) < 0.5)
        sc = hs if p == "fixed" else 0.8 * hs + 0.6 * rng.standard_normal(n) + (leak * wt if p != "fixed" else 0)
        tail = sc >= np.quantile(sc, 0.99)
        P[p] = {"win_t": wt, "score": sc, "tail1": tail, "logit": rng.standard_normal(n) + 0.74 * wt * 0}
    return {"pool": pool, "project": proj}, P


def test_clean_synthetic_passes_and_structure():
    R, P = synth()
    cfg = dict(CFG)
    cfg["gates"] = dict(CFG["gates"])
    cfg["gates"]["tail_g2"] = dict(CFG["gates"]["tail_g2"], bootstrap_B=2000)
    out = G.evaluate(R, P, cfg)
    g = out["learned"]["gbm"]
    assert out["n_by_pool"]["A"] == 90_000
    assert g["tail_G2_A"]["verdict"] == "pass", g["tail_G2_A"]
    assert g["G1_A"]["verdict"] == "pass" and g["G1_A"]["n_top"] >= 7200
    assert g["FS"]["tail_G2"]["verdict"] == "not-evaluable" and g["FS"]["G1"]["verdict"] == "not-evaluable"
    ov = G.overall(out, {"A": {"verdict": "pass"}, "AB": {"verdict": "pass"}})
    assert ov["fdr_claims_allowed"], ov
    assert not G.overall(out, None)["fdr_claims_allowed"]                       # G3 missing blocks


def test_leak_blocks():
    R, P = synth(leak=1.5)
    cfg = dict(CFG)
    cfg["gates"] = dict(CFG["gates"])
    cfg["gates"]["tail_g2"] = dict(CFG["gates"]["tail_g2"], bootstrap_B=1000)
    out = G.evaluate(R, P, cfg)
    ov = G.overall(out, {"A": {"verdict": "pass"}, "AB": {"verdict": "pass"}})
    assert not ov["fdr_claims_allowed"] and "tail_G2_A" in ov["blocked_by"]
