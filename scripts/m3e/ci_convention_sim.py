#!/usr/bin/env python3
"""A1.6 CI-convention simulation (DESIGN §6.0 A1.6, §6.7): which 95% CI for the arm-P statistic (mean paired FDP over
replicates) has coverage closest to nominal without going below 93%:
    (a) Poisson-bootstrap variance of the replicate mean + between-seed variance / J   (v0.3 rule)
    (b) Poisson-bootstrap variance only
Both CIs come from ONE call of the library estimator (vulcan.stats.entrapment.evaluate_arm, n_boot = 2000, the
production default) per simulated dataset.

Generative model (fixed before running; one simulated dataset = one replicate):
  * m original target peptides, a fraction pi present.  Present: target score ~ N(mu, 1), own decoy ~ N(0, 1).
    Absent: null competition, winner score = max of two N(0,1), target wins w.p. p_null.
  * The target side (spectra, targets, the one decoy seed) is drawn ONCE per dataset and is shared by the J = 3
    entrapment seeds (2001-2003), as in the dev design (pseudo-reverse decoys x 3 entrapment seeds).
  * Per entrapment seed: m paired shuffled entrapment peptides (r = 1), all absent: null competition with their own
    decoys (winner score max of two N(0,1), entrapment wins w.p. p_null).
  * Per seed: peptide-level TDC over T u E vs their decoys (vulcan.stats.tdc.tdc_qvalues), list at q <= alpha,
    paired counts (vulcan.stats.entrapment.paired_counts) -> Replicate.
Estimand of the CI: mu_P = E[mean paired FDP over the J seeds] over the whole generative process (dataset and
seeds), estimated by the grand mean over the R replicates (MC SE reported).  Coverage = fraction of replicates whose
CI contains mu_P.  Diagnostics: coverage of the expected true FDP, mean CI widths.
Scenarios (fixed):
  S1 D0-like 1%:  m = 30000, pi = 0.8, mu = 4.0, p_null = 0.5,  alpha = 0.01
  S2 D0-like 5%:  m = 30000, pi = 0.8, mu = 4.0, p_null = 0.5,  alpha = 0.05
  S3 D1-like 1%:  m = 7000,  pi = 0.85, mu = 4.0, p_null = 0.5, alpha = 0.01
  S4 biased 1%:   m = 30000, pi = 0.8, mu = 4.0, p_null = 0.55, alpha = 0.01
R = 10,000 per scenario, split over SEEDS (11, 23, 47, 83, 131); bootstrap seed per replicate = (seed, rep).
Decision rule (A1.6, applied mechanically): a convention is eligible iff its coverage >= 0.93 in EVERY scenario;
among eligible ones adopt the one with the smallest mean |coverage - 0.95| over scenarios; if none is eligible or
they tie, keep (a) (A1.6 interim rule: the more conservative decides).

    python scripts/m3e/ci_convention_sim.py [--reps 10000] [--workers 12] --out results/.../ci_convention_sim.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.stats._util import ENTRAP_SEEDS_DEV, SEEDS, make_rng  # noqa: E402

SCEN = {
    "S1_D0like_1pct": dict(m=30000, pi=0.8, mu=4.0, p_null=0.5, alpha=0.01),
    "S2_D0like_5pct": dict(m=30000, pi=0.8, mu=4.0, p_null=0.5, alpha=0.05),
    "S3_D1like_1pct": dict(m=7000, pi=0.85, mu=4.0, p_null=0.5, alpha=0.01),
    "S4_biased_1pct": dict(m=30000, pi=0.8, mu=4.0, p_null=0.55, alpha=0.01),
}


def one(args):
    name, seed, rep = args
    from vulcan.stats.entrapment import Replicate, evaluate_arm, paired_counts
    from vulcan.stats.tdc import tdc_qvalues
    c = SCEN[name]
    m, pi, mu, pn, a = c["m"], c["pi"], c["mu"], c["p_null"], c["alpha"]
    rng = make_rng(seed, 61, rep)
    npres = int(round(m * pi))
    pres = np.arange(m) < npres
    t_s = np.where(pres, rng.standard_normal(m) + mu, np.maximum(rng.standard_normal(m), rng.standard_normal(m)))
    d_s = rng.standard_normal(m)
    t_dec = np.where(pres, d_s > t_s, rng.random(m) >= pn)        # decoy wins the target pair
    wt = np.where(pres, np.maximum(t_s, d_s), t_s)
    reps, true = [], []
    for es in ENTRAP_SEEDS_DEV:
        r2 = make_rng(seed, 62, rep, es)
        we = np.maximum(r2.standard_normal(m), r2.standard_normal(m))
        e_dec = r2.random(m) >= pn
        W = np.r_[wt, we]
        D = np.r_[t_dec, e_dec]
        q = tdc_qvalues(W, D)
        acc = (q <= a) & ~D
        pair = np.r_[np.arange(m), np.arange(m)]
        ent = np.r_[np.zeros(m, bool), np.ones(m, bool)]
        ts = ~D
        reps.append(Replicate(paired_counts(pair[ts], ent[ts], W[ts], acc[ts]), es, 1001))
        false = np.r_[~pres, np.ones(m, bool)]
        true.append((acc & false).sum() / max(acc.sum(), 1))
    ev = evaluate_arm(reps, a, "paired", n_boot=2000, seed=int(seed) * 100003 + rep)
    return (name, ev["mean_fdp"], ev["ci95_a"][0], ev["ci95_a"][1], ev["ci95_b"][0], ev["ci95_b"][1],
            float(np.mean(true)), ev["var_boot"], ev["var_between_over_J"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10000)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--scenarios", nargs="+", default=list(SCEN))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    jobs = []
    for name in a.scenarios:
        per = a.reps // len(SEEDS)
        for s in SEEDS:
            jobs += [(name, s, r) for r in range(per)]
    t0 = time.time()
    with Pool(a.workers) as pool:
        res = pool.map(one, jobs, chunksize=20)
    out = {"rule": "eligible iff coverage >= 0.93 in every scenario; adopt eligible convention with smallest mean "
                   "|coverage - 0.95|; none/tie -> (a)", "scenarios": {}, "wall_s": time.time() - t0}
    arr = {}
    for r in res:
        arr.setdefault(r[0], []).append(r[1:])
    for name, rows in arr.items():
        x = np.array(rows)
        mean, la, ua, lb, ub, tr, vb, vbt = x.T
        mu = float(mean.mean())
        mu_se = float(mean.std(ddof=1) / np.sqrt(mean.size))
        mt = float(tr.mean())
        cov = lambda lo, hi, v: float(np.mean((lo <= v) & (v <= hi)))  # noqa: E731
        out["scenarios"][name] = {
            **SCEN[name], "R": int(mean.size), "mu_P": mu, "mu_P_mc_se": mu_se, "sd_mean_fdp_across_reps": float(mean.std(ddof=1)),
            "mean_true_fdp": mt,
            "coverage_a": cov(la, ua, mu), "coverage_b": cov(lb, ub, mu),
            "coverage_a_of_true_fdp": cov(la, ua, mt), "coverage_b_of_true_fdp": cov(lb, ub, mt),
            "mean_halfwidth_a": float(((ua - la) / 2).mean()), "mean_halfwidth_b": float(((ub - lb) / 2).mean()),
            "mean_se_a": float(np.sqrt(vb + vbt).mean()), "mean_se_b": float(np.sqrt(vb).mean()),
            "coverage_mc_se": float(np.sqrt(0.95 * 0.05 / mean.size)),
        }
        print(name, json.dumps(out["scenarios"][name]), flush=True)
    sc = out["scenarios"].values()
    elig = {k: all(s[f"coverage_{k}"] >= 0.93 for s in sc) for k in ("a", "b")}
    dev = {k: float(np.mean([abs(s[f"coverage_{k}"] - 0.95) for s in sc])) for k in ("a", "b")}
    ok = [k for k in ("a", "b") if elig[k]]
    if not ok or (len(ok) == 2 and dev["a"] == dev["b"]):
        choice = "a"
    else:
        choice = min(ok, key=lambda k: dev[k])
    out.update(eligible=elig, mean_abs_dev_from_095=dev, chosen=choice)
    print("chosen", choice, elig, dev)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
