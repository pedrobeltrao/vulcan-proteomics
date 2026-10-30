#!/usr/bin/env python3
"""Monte-Carlo checks of the IMPLEMENTED A2 gate code (vulcan.gates.a2_stats) against the draft's checks
(reviews/checks/a2_power.json, seed 131).  Synthetic data only; reads no spectrum, PSM or label.

  S3  coverage of the one-sided 95% max(CP, cluster bootstrap) bound, through a2_stats.cluster_bootstrap_tail
      (draft design: J = 10, Dirichlet(3) sizes, n_tail 800, pooled mean 0.55, between-project SD 0 / 0.05 / 0.10)
  S5  G1 pass probability of the 99%-CI rule at n_top and rho (MC through vulcan.gates.stats.spearman_ci, the code
      a2_stats.g1_gate_a2 calls) vs the draft's analytic values; n_min (90% pass at rho = 0) <= 7,200
  S6  G3 per-bin TOST (a2_stats.tost_bin), all 10 bins must pass, vs the draft's table
  S7  FS delta: z * delta gaps to physical spacings (analytic), delta from config/a2.yaml
  T   full A2 tail-G2 rule (a2_stats.simulate_tail_g2_a2): verdict rates incl. the pool-validity control, and the
      structural check that A2 never passes where the A1 rule (same data) fails
Agreement criterion: |implemented - draft| <= 4 combined MC standard errors (both are MC estimates; analytic draft
values have SE 0).

    python scripts/a2/power_checks.py [--n-rep-scale 1.0] [--out results/2026-10-27_A2_implementation/power_checks.json]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates import a2_config  # noqa: E402
from vulcan.gates import a2_stats as S  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402

DRAFT = os.path.join(ROOT, "reviews", "checks", "a2_power.json")
SEED = 131


def _se(p, n):
    return math.sqrt(max(p * (1 - p), 1e-12) / n) if n else 0.0


def _cmp(mine, n_mine, draft, n_draft, k=4.0):
    """SE under H0 (same true rate): pooled rate if both are MC estimates, the draft's rate if it is analytic."""
    p0 = (mine * n_mine + draft * n_draft) / (n_mine + n_draft) if n_draft else draft
    se = math.sqrt(_se(p0, n_mine) ** 2 + (_se(p0, n_draft) ** 2 if n_draft else 0.0))
    se = max(se, 1.0 / max(n_mine, 1))
    return {"implemented": mine, "draft": draft, "combined_se": se, "z": (mine - draft) / se,
            "agree": bool(abs(mine - draft) <= k * se)}


# ------------------------------------------------------------------------------------------------ S3
def s3(n_rep=1000, B=2000, J=10, n_tail=800):
    g = np.random.default_rng(SEED + 3)
    draft = {r["between_project_sd"]: r["coverage_max_cp_cluster"]
             for r in json.load(open(DRAFT))["S3_heterogeneity_coverage"]}
    out = []
    for tau in (0.0, 0.05, 0.10):
        cov = 0
        for _ in range(n_rep):
            w = g.dirichlet(np.full(J, 3.0))
            n = np.maximum(g.multinomial(n_tail, w), 1)
            mu = 0.55
            if tau > 0:
                s = mu * (1 - mu) / tau ** 2 - 1
                p = g.beta(mu * s, (1 - mu) * s, size=J)
            else:
                p = np.full(J, mu)
            k = g.binomial(n, p)
            truth = float((n * p).sum() / n.sum())
            proj = np.repeat(np.arange(J), n)
            tl = np.ones(proj.size, bool)
            tg = np.concatenate([np.arange(nj) < kj for nj, kj in zip(n, k)])
            cb = S.cluster_bootstrap_tail(proj, tl, tg, tl, tg, B=B, seed=g)
            ucb = max(GS.cp_upper(int(k.sum()), int(n.sum())), cb["ucb_abs"])
            cov += ucb >= truth
        out.append({"between_project_sd": tau, **_cmp(cov / n_rep, n_rep, draft[tau], 1000)})
    return out


# ------------------------------------------------------------------------------------------------ S5
def s5(n_rep=2000):
    d = json.load(open(DRAFT))["S5_g1"]
    out = []
    for row in d["pass_prob_table"]:
        n = row["n_top"]
        if n not in (5000, 7200, 8139, 15000):
            continue
        for rho, key in ((0.0, "pass_rho0"), (0.02, "pass_rho0.02"), (0.05, "pass_rho0.05")):
            m = S.g1_pass_rate_mc(n, rho, n_rep, seed=SEED + n + int(rho * 1000))
            out.append({"n_top": n, "rho": rho, **_cmp(m, n_rep, row[key], 0)})
    nmin = {"draft_n_min_90pct": d["n_min_for_90pct_pass_at_rho0"], "config_n_min": S.G1_N_MIN,
            "config_ge_draft": bool(S.G1_N_MIN >= d["n_min_for_90pct_pass_at_rho0"])}
    return {"table": out, "n_min": nmin}


# ------------------------------------------------------------------------------------------------ S6
def s6(n_rep=2000):
    g = np.random.default_rng(SEED + 6)
    draft = json.load(open(DRAFT))["S6_g3"]
    out = []
    for row in draft:
        n, r0, dtop = row["n_per_bin_per_label"], row["recall"], row["top_bin_decoy_recall_deficit"]
        if n < 1000:
            continue
        ok = 0
        for _ in range(n_rep):
            good = True
            for b in range(10):
                kT = int(g.binomial(n, r0))
                kD = int(g.binomial(n, r0 - dtop if b == 9 else r0))
                good &= S.tost_bin(kT, n, kD, n)["pass"]
            ok += good
        out.append({"n": n, "recall": r0, "top_bin_decoy_deficit": dtop,
                    **_cmp(ok / n_rep, n_rep, row["pass_A2_tost"], 2000)})
    return out


# ------------------------------------------------------------------------------------------------ S7
def s7(cfg):
    d = cfg["fs"]["delta_th"]
    draft = [r for r in json.load(open(DRAFT))["S7_fs_delta"] if r["delta_th"] == round(2 / 7, 4)]
    phys = {"isotope_13C": 1.003355, "deamidation": 0.984016, "NH3": 17.026549, "H2O": 18.010565,
            "CO(a-ion)": 27.994915, "zero": 0.0}
    out = []
    for r in draft:
        m = r["z"] * d
        gap = min(min(abs(m - v), abs(m - (v % 1.0))) for v in phys.values())
        out.append({"z": r["z"], "gap_da": gap, "draft_gap_da": r["gap_da"], "agree": bool(abs(gap - r["gap_da"]) < 1e-3),
                    "ge_0.12": bool(gap >= 0.12)})
    return out


# ------------------------------------------------------------------------------------------------ full tail rule
def tail_rule(n_rep=500, B=2000):
    out = []
    g = np.random.default_rng(SEED + 1)
    never_looser = True
    for n_tail in (200, 400, 750, 1500):
        for pl in (0.50, 0.55, 0.60, 0.65):
            v = {"pass": 0, "fail": 0, S.NOT_EVALUABLE: 0, S.INVALID: 0}
            a1_pass = 0
            for _ in range(n_rep):
                J = 10
                w = g.dirichlet(np.full(J, 3.0))
                nL = g.multinomial(g.poisson(n_tail), w)
                nF = g.multinomial(g.poisson(n_tail), w)
                kL, kF = g.binomial(nL, pl), g.binomial(nF, 0.5)
                proj = np.concatenate([np.full(a + b, j) for j, (a, b) in enumerate(zip(nL, nF))])
                tl = np.concatenate([np.r_[np.ones(a, bool), np.zeros(b, bool)] for a, b in zip(nL, nF)])
                gl = np.concatenate([np.r_[np.arange(a) < k, np.zeros(b, bool)] for a, b, k in zip(nL, nF, kL)])
                tf = ~tl
                gf = np.concatenate([np.r_[np.zeros(a, bool), np.arange(b) < k] for a, b, k in zip(nL, nF, kF)])
                r = S.tail_g2_a2(proj, tl, gl, tf, gf, B=B, seed=g)
                v[r["verdict"]] += 1
                a1 = r["a1_statistics"]["verdict"] == "pass"
                a1_pass += a1
                if r["verdict"] == "pass" and not a1:
                    never_looser = False
            out.append({"n_tail_mean": n_tail, "p_learned": pl, "p_fixed": 0.5,
                        **{k: x / n_rep for k, x in v.items()}, "a1_rule_pass_same_data": a1_pass / n_rep})
    return {"grid": out, "a2_never_passes_where_a1_fails": never_looser,
            "note": "independent learned / hyperscore tails (conservative for the diff clause); homogeneous projects"}


def run(scale=1.0):
    cfg = a2_config.load()
    k = lambda n: max(50, int(n * scale))  # noqa: E731
    res = {"seed": SEED, "draft": os.path.relpath(DRAFT, ROOT), "n_rep_scale": scale,
           "S3": s3(k(1000)), "S5": s5(k(2000)), "S6": s6(k(2000)), "S7": s7(cfg), "tail_rule": tail_rule(k(500))}
    flat = [r["agree"] for r in res["S3"] + res["S5"]["table"] + res["S6"] + res["S7"]]
    res["all_agree"] = bool(all(flat)) and res["S5"]["n_min"]["config_ge_draft"] \
        and res["tail_rule"]["a2_never_passes_where_a1_fails"]
    res["n_comparisons"] = len(flat)
    return res


def main():  # pragma: no cover
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rep-scale", type=float, default=1.0)
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "2026-10-27_A2_implementation", "power_checks.json"))
    a = ap.parse_args()
    res = run(a.n_rep_scale)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    for sec in ("S3", "S6"):
        for r in res[sec]:
            print(sec, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    for r in res["S5"]["table"]:
        print("S5", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print("S5 n_min", res["S5"]["n_min"])
    for r in res["S7"]:
        print("S7", r)
    for r in res["tail_rule"]["grid"]:
        print("T", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print("A2 never passes where A1 fails:", res["tail_rule"]["a2_never_passes_where_a1_fails"])
    print("ALL AGREE:", res["all_agree"], "comparisons", res["n_comparisons"])


if __name__ == "__main__":
    main()
