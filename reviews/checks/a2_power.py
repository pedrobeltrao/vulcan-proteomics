#!/usr/bin/env python3
"""Amendment A2 draft (reviews/amendment_A2_draft.md), item 6: simulation and analytic checks of the proposed rules.

Synthetic data only.  No null spectra, no gate output, and no fresh-null data are read.  The only repository code
used is vulcan.gates.stats (the A1 gate statistics, unchanged) and vulcan.stats.g2.wilson_ci.

Sections (numbered as in the draft, §A2.9):
  S1  tail-G2 verdict rates vs. true tail target fraction and n_tail (A1 margins, unchanged), contamination-free null.
  S2  tail-G2 under contamination: verdict rate of a NON-leaking scorer vs. the share c of truly-human tail spectra.
  S3  between-project heterogeneity: coverage of the exact CP bound vs. the A2 rule max(CP, cluster bootstrap).
  S4  spectrum-level screens: asymmetric "best hit in human+own DB is human" vs. the label-symmetric family screen
      (A2.2(d)); target fraction in the learned tail with and without a leak.
  S5  G1: pass probability vs. n (rho = 0) and at the boundary; rho induced by contamination for a non-leaking
      scorer; rho of a target-favouring leak.
  S6  G3: per-bin A1 point rule at n ~ 100 vs. the A2 per-bin TOST at n >= 1000.
  S7  FS nulls: z * delta for fragment charges 1-4 vs. physical mass differences (analytic).

    python reviews/checks/a2_power.py            # writes reviews/checks/a2_power.json
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates import stats as GS  # noqa: E402
from vulcan.stats.g2 import wilson_ci  # noqa: E402

SEED = 131  # SEEDS[4]; synthetic checks only
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "a2_power.json")


def rng(k=0):
    return np.random.default_rng(SEED + k)


# ------------------------------------------------------------------------------------------------ S1
def s1_tail_power(n_rep=2000):
    out = []
    for n_tail in (200, 400, 750, 1500):
        for pl in (0.50, 0.55, 0.60, 0.65, 0.70):
            r = GS.simulate_tail_g2(n_rep, n_null=int(n_tail / 0.01), q_tail=0.01, p_fixed=0.5, p_learned=pl,
                                    overlap=0.6, seed=SEED)
            out.append({"n_tail": n_tail, "p_learned": pl, "p_fixed": 0.5, "pass": r["pass"], "fail": r["fail"],
                        "not_evaluable": r["not-evaluable"]})
    return out


# ------------------------------------------------------------------------------------------------ S2
def s2_contamination(n_rep=2000, n_tail=750):
    """Non-leaking scorer.  A share c of the learned tail is truly-human spectra (target wins w.p. ~1).  The learned
    scorer recognises them better than hyperscore: hyperscore's tail holds c_fixed = c / 2 such spectra (as in
    M3e's S. aureus source: 375 vs 131)."""
    g = rng(2)
    out = []
    for c in (0.0, 0.002, 0.01, 0.05, 0.10, 0.20, 0.40):
        v = {"pass": 0, "fail": 0}
        for _ in range(n_rep):
            nL = n_tail
            nF = n_tail // 2
            kL = g.binomial(int(round(nL * (1 - c))), 0.5) + int(round(nL * c))
            kF = g.binomial(int(round(nF * (1 - c / 2))), 0.5) + int(round(nF * c / 2))
            # independent tails (conservative for the diff clause: no overlap credit)
            a = np.r_[np.ones(nL, bool), np.zeros(nF, bool)]
            b = np.r_[np.arange(nL) < kL, np.zeros(nF, bool)]
            cc = np.r_[np.zeros(nL, bool), np.ones(nF, bool)]
            d = np.r_[np.zeros(nL, bool), np.arange(nF) < kF]
            r = GS.tail_g2(a, b, cc, d)
            v[r["verdict"] if r["verdict"] in v else "fail"] += 1
        out.append({"contaminated_share_of_learned_tail": c, "pass": v["pass"] / n_rep,
                    "expected_target_fraction": 0.5 + c / 2})
    return out


# ------------------------------------------------------------------------------------------------ S3
def _cluster_ucb(k, n, conf=0.95, B=2000, g=None):
    """Percentile cluster (project) bootstrap one-sided upper bound of sum(k)/sum(n)."""
    J = len(n)
    idx = g.integers(0, J, size=(B, J))
    est = k[idx].sum(1) / np.maximum(n[idx].sum(1), 1)
    return float(np.quantile(est, conf))


def s3_heterogeneity(n_rep=1000, J=10, n_tail=800):
    """Projects differ in their tail target fraction (Beta around a pooled mean of 0.55, SD tau).  Coverage of the
    one-sided 95% upper bound for the pooled (n-weighted) fraction."""
    g = rng(3)
    out = []
    for tau in (0.0, 0.05, 0.10):
        cov_cp = cov_max = 0
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
            truth = float((n * p).sum() / n.sum())  # pooled estimand given the project mix
            cp = GS.cp_upper(int(k.sum()), int(n.sum()))
            cl = _cluster_ucb(k, n, g=g)
            cov_cp += cp >= truth
            cov_max += max(cp, cl) >= truth
        out.append({"between_project_sd": tau, "J": J, "coverage_cp": cov_cp / n_rep,
                    "coverage_max_cp_cluster": cov_max / n_rep})
    return out


# ------------------------------------------------------------------------------------------------ S4
def s4_screens(n_spec=400_000, n_rep=20, leak=0.0):
    """Per spectrum: own-proteome best hyperscore o (the certified own peptide; large for certified spectra), best
    human target h_T and best human decoy h_D (exchangeable chance matches), learned scores l = 0.8 h + 0.6 e (+leak
    on the target side, in SD units of l).  Tail = learned winner >= its top-1% quantile (over the unscreened pool,
    standing in for the host-run threshold).  Screens:
      none        : all spectra
      asym        : drop if h_T > o            (human target DB only, no human decoys -> label-asymmetric)
      sym_family  : drop if max(h_T, h_D) > o  (human targets AND their paired decoys -> label-symmetric, A2.2(d))
    Returns the learned-tail target fraction under each screen."""
    g = rng(4)
    res = {k: [] for k in ("none", "asym", "sym_family")}
    ntail = {k: [] for k in res}
    for _ in range(n_rep):
        hT, hD = g.gumbel(0, 1, n_spec), g.gumbel(0, 1, n_spec)
        o = g.gumbel(0, 1, n_spec) + g.choice([0.0, 2.5], n_spec, p=[0.5, 0.5])  # half the pool well explained by own
        lT = 0.8 * hT + 0.6 * g.standard_normal(n_spec) + leak
        lD = 0.8 * hD + 0.6 * g.standard_normal(n_spec)
        win_t = lT > lD
        lw = np.maximum(lT, lD)
        thr = np.quantile(lw, 0.99)
        tail = lw >= thr
        keep = {"none": np.ones(n_spec, bool), "asym": ~(hT > o), "sym_family": ~(np.maximum(hT, hD) > o)}
        for k, m in keep.items():
            t = tail & m
            res[k].append(float(win_t[t].mean()))
            ntail[k].append(int(t.sum()))
    return {k: {"tail_target_fraction": float(np.mean(v)), "mc_se": float(np.std(v) / math.sqrt(len(v))),
                "mean_n_tail": float(np.mean(ntail[k]))} for k, v in res.items()} | {"leak_sd_units": leak}


# ------------------------------------------------------------------------------------------------ S5
def g1_pass_prob(n, rho, bound=0.05, conf=0.99):
    se = math.sqrt((1 + rho * rho / 2) / (n - 3))
    z = stats.norm.ppf(0.5 + conf / 2)
    zr = math.atanh(rho)
    lo, hi = math.atanh(-bound), math.atanh(bound)
    # pass iff z_hat - z*se > lo and z_hat + z*se < hi
    a, b = lo + z * se, hi - z * se
    if a >= b:
        return 0.0
    return float(stats.norm.cdf((b - zr) / se) - stats.norm.cdf((a - zr) / se))


def g1_n_min(power=0.9):
    n = 100
    while g1_pass_prob(n, 0.0) < power:
        n += 50
    return n


def s5_g1(n_mc=200_000):
    g = rng(5)
    tab = [{"n_top": n, "pass_rho0": g1_pass_prob(n, 0.0), "pass_rho0.02": g1_pass_prob(n, 0.02),
            "pass_rho0.05": g1_pass_prob(n, 0.05)} for n in (1000, 2500, 5000, 7200, 8139, 15000, 30000)]
    # contamination-induced rho for a NON-leaking scorer: top-decile winners; a share c are truly-human spectra
    # (target winners, score shifted up by 1 SD, logit drawn from the observed-peptide distribution, which the
    # classifier ranks like targets); the rest are exchangeable nulls.  Classifier: logit ~ N(0.74, 1) for targets,
    # N(0, 1) for decoys (AUC ~0.70 as M3e's held-out AUC).
    induced = []
    for c in (0.0, 0.05, 0.10, 0.20, 0.40):
        n = n_mc
        cont = g.random(n) < c
        is_t = np.where(cont, True, g.random(n) < 0.5)
        score = g.standard_normal(n) + np.where(cont, 1.0, 0.0)
        logit = g.standard_normal(n) + np.where(is_t, 0.74, 0.0)
        induced.append({"contaminated_share": c, "rho": float(stats.spearmanr(score, logit).statistic)})
    # a target-favouring leak on a clean null: learned score = s + beta * 1[target]
    leak = []
    for beta in (0.05, 0.1, 0.2, 0.3):
        n = n_mc
        is_t = g.random(n) < 0.5
        score = g.standard_normal(n) + beta * is_t
        logit = g.standard_normal(n) + 0.74 * is_t
        leak.append({"beta_sd": beta, "rho": float(stats.spearmanr(score, logit).statistic)})
    return {"pass_prob_table": tab, "n_min_for_90pct_pass_at_rho0": g1_n_min(),
            "contamination_induced_rho": induced, "label_leak_rho": leak}


# ------------------------------------------------------------------------------------------------ S6
def newcombe(k1, n1, k2, n2, conf):
    p1, l1, u1 = wilson_ci(k1, n1, conf)
    p2, l2, u2 = wilson_ci(k2, n2, conf)
    d = p1 - p2
    return d, d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)


def g3_rules(n_per_bin, rT, rD, n_bins, g, n_rep):
    a1 = a2 = 0
    for _ in range(n_rep):
        ok1 = ok2 = True
        for b in range(n_bins):
            rt = rT[b] if isinstance(rT, (list, tuple)) else rT
            kT = g.binomial(n_per_bin, rt)
            kD = g.binomial(n_per_bin, rD)
            d, lo, hi = newcombe(kT, n_per_bin, kD, n_per_bin, 0.90)  # 90% two-sided = TOST at 5% each side
            ok1 &= abs(d) <= 0.01
            ok2 &= (lo >= -0.01) and (hi <= 0.01)
        a1 += ok1
        a2 += ok2
    return a1 / n_rep, a2 / n_rep


def s6_g3(n_rep=2000):
    g = rng(6)
    out = []
    for n in (100, 1000, 2000, 5000):
        for r0 in (0.999, 0.995, 0.99):
            for delta_top in (0.0, 0.01, 0.02, 0.03):
                rT = [r0] * 9 + [min(1.0, r0 + delta_top)] if delta_top else r0
                rD = r0 if delta_top == 0 else r0
                if delta_top:  # put the asymmetry in the top bin, decoys lost (harmful direction)
                    rT = [r0] * 9 + [r0]
                    rD_top = r0 - delta_top
                    p1 = p2 = 0
                    for _ in range(n_rep):
                        ok1 = ok2 = True
                        for b in range(10):
                            kT = g.binomial(n, r0)
                            kD = g.binomial(n, rD_top if b == 9 else r0)
                            d, lo, hi = newcombe(kT, n, kD, n, 0.90)
                            ok1 &= abs(d) <= 0.01
                            ok2 &= (lo >= -0.01) and (hi <= 0.01)
                        p1 += ok1
                        p2 += ok2
                    pa1, pa2 = p1 / n_rep, p2 / n_rep
                else:
                    pa1, pa2 = g3_rules(n, rT, rD, 10, g, n_rep)
                out.append({"n_per_bin_per_label": n, "recall": r0, "top_bin_decoy_recall_deficit": delta_top,
                            "pass_A1_point_rule": float(pa1), "pass_A2_tost": float(pa2)})
    return out


# ------------------------------------------------------------------------------------------------ S7
def s7_fs_delta():
    phys = {"isotope_13C": 1.003355, "deamidation": 0.984016, "NH3": 17.026549, "H2O": 18.010565,
            "CO(a-ion)": 27.994915, "zero": 0.0}
    rows = []
    for dlt in (round(2 / 7, 4), 0.27, 0.31):
        for z in (1, 2, 3, 4):
            m = z * dlt
            near = min(phys.items(), key=lambda kv: min(abs(m - kv[1]), abs(m - (kv[1] % 1.0))))
            gap_da = min(abs(m - near[1]), abs(m - (near[1] % 1.0)))
            rows.append({"delta_th": dlt, "z": z, "mass_shift_da": round(m, 4), "nearest_physical": near[0],
                         "gap_da": round(gap_da, 4), "gap_ppm_at_mz_500": round(gap_da / z / 500 * 1e6, 1)})
    return rows


def main():
    res = {"seed": SEED, "note": "synthetic only; see module docstring",
           "S1_tail_g2_power": s1_tail_power(),
           "S2_contamination_nonleaking": s2_contamination(),
           "S3_heterogeneity_coverage": s3_heterogeneity(),
           "S4_screens_no_leak": s4_screens(leak=0.0),
           "S4_screens_leak_0.25sd": s4_screens(leak=0.25),
           "S5_g1": s5_g1(),
           "S6_g3": s6_g3(),
           "S7_fs_delta": s7_fs_delta()}
    json.dump(res, open(OUT, "w"), indent=1)
    for k, v in res.items():
        print("==", k)
        if isinstance(v, list):
            for r in v:
                print("  ", r)
        else:
            print("  ", v)


if __name__ == "__main__":
    main()
