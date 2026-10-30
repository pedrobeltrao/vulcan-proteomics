#!/usr/bin/env python3
"""Amendment A3.1.9: power check of the A3 gate rows (V1-V14) for v2_sep, before the tag.  SYNTHETIC DATA ONLY.

No DN2 file, no gate output and no E7b output is read.  Seed 131 (SEEDS[4]).  The rules are the implemented ones:
vulcan.gates.a2_stats (tail-G2, G1, G3 TOST), vulcan.gates.a3_stats (V6, V12, tiers), vulcan.stats.g2 (G2-original).

Projection of the DN2 sizes from DEV structure (A3.1.9: E7's G2 mixtures, 35% secondary units, D0 neighbour-filter
drop rate 11.8%), as a function of the admitted A count |A| and of B (all non-A spectra of the admitted runs):
  * |A| is bounded by A2.9.3's STOPPING RULE (independent review M3): the walk stops at the first admitted project
    count >= 10 at which |A| >= 120,000 (then also n_top >= 72,000 and the tail proxy >= 400 at the constants below),
    and at 16 projects otherwise.  Per-run A density d in {2k, 4k, 8k, 15k, 25k} (3 runs per project) gives the
    reachable points (J projects, |A|) = (16, 96k), (10, 120k), (10, 240k), (10, 450k), (10, 750k); 16 projects
    happen only when the criteria are not met (low yield).  B/|A| in {1, 3} (B is not projected by any committed dev
    number; both are shown);
  * fraction of A spectra (and of every null unit) with a candidate (a winner): 0.60 (A2.9.3 design 72,000 / 120,000;
    dev check: 0.64 of D7 nulls in results/2026-10-28_A3_preconditions/truncation_dev.json);
  * secondary units per B spectrum: 0.35 / 0.65 = 0.54 (E7 G2: 35% of null units are secondary); B spectra with >= 2
    group winners: 0.40 (so class-1 V6 winners = 0.40 |B| c, class-2 = 0.54 |B| c);
  * tail fraction of A units (and of A-syn units at the stratum-2/>=3 thresholds): 400 / 120,000 (A2.9.3 tail proxy);
  * G3 recall 0.998 (E2) for every pool; V6 drop rate 0.118 (E7 D0); V12 miss rate 0.0981 and 63.75% of A spectra
    with a candidate per label (dev max-variance decile on a D0-host mixture with production calibration,
    results/2026-10-28_A3_preconditions/truncation_dev.json); |FS| = |A| by construction (A2.4: one FS copy per A
    spectrum), 5% of FS spectra with a winner (A2.4: most none);
  * V13: the A-syn tail fraction at the stratum-2/>=3 thresholds is not measured on dev; it is set to A's and a
    sensitivity at half that fraction is reported (independent review minor 1).
Outputs (reviews/checks/a3_power.json + printed table):
  * n_min of V6 and V12 (the V6 rule: smallest n per label per bin at which a symmetric process at the dev rate
    passes one 90% Newcombe TOST bin with probability >= 0.90; 20,000 MC replicates per n) and V14 (A2.6 rule:
    7,200; its pass probability at 7,200 is reported); V11 is intentionally unused (A3.2.4);
  * per V row and |A| scenario: projected n, evaluable or not, pass probability under a clean scorer, and the
    pass ("firing-miss") rate at the margin (tail fraction 0.60 / diff 0.10, rho 0.05, recall gap 0.01, drop gap
    0.01, G2 excess 0.01, AUC 0.55);
  * the joint pass probability per tier: rows are simulated separately (a single synthetic DN2 per replicate was not
    built), so the joint is reported as BOUNDS that hold under any dependence between rows -- Frechet lower bound
    max(0, sum p - (k - 1)) and upper bound min p -- plus the product under independence (independent review M4;
    a disclosed departure from "simulated joint"); a not-evaluable row has pass probability 0;
  * flags, recorded in advance (A3.1.9, read strictly as "can never PASS"): a row is flagged when its clean-scorer
    pass probability is < 0.05 at EVERY reachable point (per B/|A| and over both); a tier is "expected not
    evaluable / not passable" when its joint upper bound is < 0.05 at every reachable point.
    python reviews/checks/a3_power.py      # writes reviews/checks/a3_power.json
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates import a2_stats as A2S  # noqa: E402
from vulcan.gates import a3_stats as A3S  # noqa: E402
from vulcan.stats import g2 as G2  # noqa: E402

SEED = 131
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "a3_power.json")
DENSITY = (2_000, 4_000, 8_000, 15_000, 25_000)        # A spectra per run; 3 runs per project
REACH = []
for _d in DENSITY:
    _j = 10
    while _j < 16 and _j * 3 * _d < 120_000:
        _j += 1
    REACH.append((_j, _j * 3 * _d))
A_GRID = tuple(a for _, a in REACH)
B_RATIO = (1.0, 3.0)
C_CAND = 0.60
SEC_PER_B = 0.35 / 0.65
B_MULTI = 0.40
TAIL_FRAC = 400 / 120_000
RECALL = 0.998
DROP = 0.118
V12_MISS = 0.0981
V12_CAND = 0.6375
FS_WIN = 0.05
N_REP_TAIL = 600
N_REP_G1 = 600
N_REP_BIN = 4000


def rng(k=0):
    return np.random.default_rng(SEED + k)


def projects(nA):
    return dict((a, j) for j, a in REACH)[nA]


# ------------------------------------------------------------------------------------------------ row simulators
def p_tail(n_tail, J, p_learned=0.5, k=0):
    if n_tail < 200:
        return {"evaluable": False, "pass": 0.0}
    r = A2S.simulate_tail_g2_a2(N_REP_TAIL, n_tail, J=J, p_fixed=0.5, p_learned=p_learned, B=2000, seed=SEED + k)
    return {"evaluable": True, "pass": r["pass"], "not_evaluable_invalid": r[A2S.INVALID], "fail": r["fail"]}


def p_g1(n_top, rho=0.0, k=0):
    if n_top < 7200:
        return {"evaluable": False, "pass": 0.0}
    p = A2S.g1_pass_rate_mc(int(n_top), rho, N_REP_G1, seed=SEED + k)
    ctl = 0.99                      # hyperscore control: two-sided 99% CI contains 0 under rho = 0
    return {"evaluable": True, "pass": p * ctl, "learned_ci_inside": p, "control_valid": ctl}


def p_bins(n_per_label, n_bins, p_t, p_d, n_min, k=0, n_rep=N_REP_BIN):
    """All bins must pass (TOST per bin); bins i.i.d.; the gap (p_t != p_d) is put in ONE bin (the top bin)."""
    if n_per_label < n_min:
        return {"evaluable": False, "pass": 0.0, "n_per_label_per_bin": int(n_per_label)}
    g = rng(k)
    n = int(n_per_label)
    ok = np.ones(n_rep, bool)
    for b in range(n_bins):
        pt = p_t if b == n_bins - 1 else p_d
        kT = g.binomial(n, pt, size=n_rep)
        kD = g.binomial(n, p_d, size=n_rep)
        ok &= np.array([A2S.tost_bin(int(a), n, int(c), n)["pass"] for a, c in zip(kT, kD)])
    return {"evaluable": True, "pass": float(ok.mean()), "n_per_label_per_bin": n}


def p_g2(n, diff=0.0):
    r = G2.g2_power(int(n), diff)
    return r["p_pass_joint_lb"]


def p_g5(nT, nD, auc=0.5):
    if nT == 0 or nD == 0:
        return 0.0
    se = math.sqrt((nT + nD + 1) / (12.0 * nT * nD))
    from scipy.stats import norm
    return float(norm.cdf((0.55 - auc) / se))


# ------------------------------------------------------------------------------------------------ main
def n_min_section():
    t0 = time.perf_counter()
    v6 = A3S.n_min_symmetric(1 - DROP, n_rep=20_000, seed=SEED)
    v12 = A3S.n_min_symmetric(1 - V12_MISS, n_rep=20_000, seed=SEED + 1)
    v14_at = A2S.g1_pass_rate_mc(7200, 0.0, 2000, seed=SEED + 2)
    return {"V6": {**v6, "rate": "keep rate 1 - 0.118 (dev drop 11.8%)"},
            "V11": "unused (A3.2.4: number intentionally unused)",
            "V12": {**v12, "rate": f"recall 1 - {V12_MISS} (dev max-variance decile)"},
            "V14": {"n_min_top": 7200, "rule": "A2.6 unchanged", "pass_prob_at_7200_rho0_learned_ci": v14_at},
            "V8": {"n_min": int(G2.g2_sample_size()), "rule": "g2_sample_size() (§6.4 budget)"},
            "wall_s": time.perf_counter() - t0}


def scenario(nA, b_ratio, nmin):
    nB = b_ratio * nA
    J = projects(nA)
    c = C_CAND
    out = {"n_A": nA, "B_over_A": b_ratio, "J_projects": J}
    rows, marg = {}, {}
    nt = TAIL_FRAC * nA
    rows["V1"] = {"n_tail": nt, "clean": p_tail(int(nt), J, 0.5, 1), "margin_p0.60": p_tail(int(nt), J, 0.60, 2)}
    ntop = 0.1 * c * nA
    rows["V2"] = {"n_top": ntop, "clean": p_g1(ntop, 0.0, 3), "margin_rho0.05": p_g1(ntop, 0.05, 4)}
    rows["V3"] = {"clean": p_bins(c * nA / 10, 10, RECALL, RECALL, 2000, 5),
                  "margin_gap0.01": p_bins(c * nA / 10, 10, RECALL - 0.01, RECALL, 2000, 6)}
    units = c * (nA + nB * (1 + SEC_PER_B))
    rows["V4"] = {"n_units": units, "clean": p_bins(units / 10, 10, RECALL, RECALL, 2000, 7),
                  "margin_gap0.01": p_bins(units / 10, 10, RECALL - 0.01, RECALL, 2000, 8)}
    sec = c * nB * SEC_PER_B
    rows["V5"] = {"n_units": sec, "clean": p_bins(sec / 10, 10, RECALL, RECALL, 2000, 9),
                  "margin_gap0.01": p_bins(sec / 10, 10, RECALL - 0.01, RECALL, 2000, 10)}
    for cls, size, k in (("V6.class1", B_MULTI * nB * c, 11), ("V6.class2", sec, 13)):
        top = 0.05 * size          # per label: top decile (10%) of the class, half targets
        rest = 0.45 * size
        nm = nmin["V6"]["n_min"]
        ev = top >= nm and rest >= nm
        if ev:
            a = p_bins(top, 1, 1 - DROP, 1 - DROP, nm, k)["pass"] * p_bins(rest, 1, 1 - DROP, 1 - DROP, nm, k + 100)["pass"]
            m = p_bins(top, 1, 1 - DROP - 0.01, 1 - DROP, nm, k + 1)["pass"] * \
                p_bins(rest, 1, 1 - DROP, 1 - DROP, nm, k + 101)["pass"]
        else:
            a = m = 0.0
        rows[cls] = {"n_winners": size, "n_per_label_top_bin": top, "clean": {"evaluable": ev, "pass": a},
                     "margin_drop_gap0.01": {"evaluable": ev, "pass": m}}
    rows["V7"] = {"n": units, "clean": {"evaluable": True, "pass": p_g2(units)},
                  "margin_excess0.01": {"evaluable": True, "pass": p_g2(units, 0.01)}}
    g2n = nmin["V8"]["n_min"]
    rows["V8"] = {"n": sec, "clean": {"evaluable": sec >= g2n, "pass": p_g2(sec) if sec >= g2n else 0.0},
                  "margin_excess0.01": {"evaluable": sec >= g2n, "pass": p_g2(sec, 0.01) if sec >= g2n else 0.0}}
    for k_, n_ in (("V9.a", units), ("V9.b", c * nA), ("V9.c", sec)):
        ev = (n_ / 2 >= 200) if k_ == "V9.c" else True
        rows[k_] = {"n": n_, "clean": {"evaluable": ev, "pass": p_g5(n_ / 2, n_ / 2) if ev else 0.0},
                    "margin_auc0.55": {"evaluable": ev, "pass": p_g5(n_ / 2, n_ / 2, 0.55) if ev else 0.0}}
    fsw = FS_WIN * nA
    fs_g5_ev = fsw / 2 >= 200
    fs_tail_ev = TAIL_FRAC * FS_WIN * nA >= 200
    fs_g1_ev = 0.1 * fsw >= 7200
    p_block = (1 - p_g5(fsw / 2, fsw / 2)) if fs_g5_ev else 0.0
    rows["V10"] = {"n_fs_winners": fsw, "evaluable_parts": {"tail_G2": fs_tail_ev, "G1": fs_g1_ev, "G5": fs_g5_ev},
                   "clean": {"evaluable": True, "pass": 1 - p_block, "note": "not evaluable never blocks (A2.4)"}}
    n12 = V12_CAND * nA / 10
    nm12 = nmin["V12"]["n_min"]
    rows["V12"] = {"n_per_label_per_bin": n12, "clean": p_bins(n12, 10, 1 - V12_MISS, 1 - V12_MISS, nm12, 15),
                   "margin_gap0.01": p_bins(n12, 10, 1 - V12_MISS - 0.01, 1 - V12_MISS, nm12, 16)}
    ns = TAIL_FRAC * 2 * nA
    rows["V13"] = {"n_tail": ns, "clean": p_tail(int(ns), J, 0.5, 17), "margin_p0.60": p_tail(int(ns), J, 0.60, 18),
                   "sensitivity_half_tail_fraction": p_tail(int(ns / 2), J, 0.5, 21)}
    ntop14 = 0.1 * c * 2 * nA
    rows["V14"] = {"n_top": ntop14, "clean": p_g1(ntop14, 0.0, 19), "margin_rho0.05": p_g1(ntop14, 0.05, 20)}
    for k, v in rows.items():
        marg[k] = v["clean"]["pass"] if v["clean"]["evaluable"] else 0.0
    t1 = ["V1", "V2", "V3", "V4", "V6.class1", "V7", "V9.a", "V9.b", "V10", "V12"]
    t2 = ["V5", "V6.class2", "V8", "V9.c", "V13", "V14"]
    def bounds(keys):
        ps = [marg[k] for k in keys]
        return {"independence": float(np.prod(ps)), "frechet_lower": float(max(0.0, sum(ps) - (len(ps) - 1))),
                "upper_min": float(min(ps))}
    b1, b2 = bounds(t1), bounds(t1 + t2)
    j1, j2 = b1["independence"], b2["independence"]
    out["joint_bounds_tier1"], out["joint_bounds_tier2"] = b1, b2
    out["rows"] = rows
    out["marginal_pass_clean"] = marg
    out["not_evaluable_rows"] = sorted(k for k, v in rows.items() if not v["clean"]["evaluable"])
    out["joint_pass_tier1"] = j1
    out["joint_pass_tier2"] = j2
    out["joint_assumption"] = ("rows simulated separately; joint_pass_* = product (independence); bounds valid under "
                               "any dependence in joint_bounds_*")
    return out


def main():
    t0 = time.perf_counter()
    nmin = n_min_section()
    print("n_min:", json.dumps({k: (v["n_min"] if isinstance(v, dict) and "n_min" in v else v) for k, v in nmin.items()
                                if k != "wall_s"}), flush=True)
    sc = []
    for nA in A_GRID:
        for br in B_RATIO:
            s = scenario(nA, br, nmin)
            sc.append(s)
            print(f"|A|={nA:>9,} J={s['J_projects']} B/A={br}: tier1 {s['joint_pass_tier1']:.3f} "
                  f"[{s['joint_bounds_tier1']['frechet_lower']:.3f}, {s['joint_bounds_tier1']['upper_min']:.3f}] "
                  f"tier2 {s['joint_pass_tier2']:.3f} "
                  f"not evaluable {s['not_evaluable_rows']}", flush=True)
    t1 = ["V1", "V2", "V3", "V4", "V6.class1", "V7", "V9.a", "V9.b", "V10", "V12"]
    t2 = ["V5", "V6.class2", "V8", "V9.c", "V13", "V14"]
    design = [s for s in sc if s["n_A"] == 120_000]

    def flagged_rows(ss):
        return sorted(k for k in ss[0]["marginal_pass_clean"] if all(s["marginal_pass_clean"][k] < 0.05 for s in ss))

    def tier_flag(ss, key):
        return bool(all(s[key]["upper_min"] < 0.05 for s in ss))
    flags = {"reachable_points": [{"J": j, "n_A": a} for j, a in REACH],
             "rule": "row: clean pass prob < 0.05 at every reachable point; tier: joint upper bound < 0.05 at every reachable point"}
    for lab, ss in (("all_B_ratios", sc), *[(f"B_over_A_{br}", [s for s in sc if s["B_over_A"] == br]) for br in B_RATIO]):
        flags[lab] = {"rows_never_passable": flagged_rows(ss),
                      "tier1_expected_not_evaluable": tier_flag(ss, "joint_bounds_tier1"),
                      "tier2_expected_not_evaluable": tier_flag(ss, "joint_bounds_tier2")}
    flags["design_point_A120k"] = {"not_evaluable_rows": sorted(set.intersection(*[set(s["not_evaluable_rows"]) for s in design])),
                                   "tier1_not_evaluable": bool(set.intersection(*[set(s["not_evaluable_rows"]) for s in design]) & set(t1))}
    flags["smallest_reachable_A_with_tier1_upper_bound_ge_0.05"] = min(
        [s["n_A"] for s in sc if s["joint_bounds_tier1"]["upper_min"] >= 0.05], default=None)
    out = {"seed": SEED, "spec": __doc__.split("    python")[0].strip(), "n_min": nmin, "scenarios": sc,
           "expected_not_evaluable": flags, "tiers": {"tier1": t1, "tier2": t2},
           "constants": {"REACH": REACH, "A_GRID": A_GRID, "B_RATIO": B_RATIO, "V12_CAND": V12_CAND, "C_CAND": C_CAND, "SEC_PER_B": SEC_PER_B,
                         "B_MULTI": B_MULTI, "TAIL_FRAC": TAIL_FRAC, "RECALL": RECALL, "DROP": DROP,
                         "V12_MISS": V12_MISS, "FS_WIN": FS_WIN},
           "wall_s": time.perf_counter() - t0}
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print("expected not evaluable:", json.dumps(flags))


if __name__ == "__main__":
    main()
