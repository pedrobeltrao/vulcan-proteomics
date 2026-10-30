#!/usr/bin/env python
"""Run the §6.7 Monte-Carlo estimator verification at >= 1e4 replicates and write
results/2026-10-26_stats_sims/{metrics.json, config.yaml, env.txt, git_sha}.

Usage: [VULCAN_THREADS=4] python scripts/run_stats_sims.py [--R 10000] [--out DIR]
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vulcan import hw  # noqa: E402

hw.pin_thread_pools()  # affinity mask, or $VULCAN_THREADS (=4 for M1 continuity)

import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from vulcan.stats import sim  # noqa: E402
from vulcan.stats._util import SEEDS  # noqa: E402
from vulcan.stats.g2 import g2_power, g2_sample_size, simulate_g2  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=10_000)
    ap.add_argument("--R-group-general", type=int, default=10_000)
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "2026-10-26_stats_sims"))
    args = ap.parse_args()
    R = args.R
    os.makedirs(args.out, exist_ok=True)
    metrics = {"replicates": R, "seeds": list(SEEDS), "alphas": list(sim.ALPHAS), "timing_s": {}}

    def run(name, fn, *a, **k):
        t0 = time.perf_counter()
        r = fn(*a, **k)
        metrics["timing_s"][name] = round(time.perf_counter() - t0, 1)
        metrics[name] = r
        print(f"[{metrics['timing_s'][name]:7.1f}s] {name}", flush=True)
        return r

    run("tdc", sim.sim_tdc, R)
    run("tdc_tight", sim.sim_tdc, R, m=4000, mu=3.0, true_always_win=True)
    run("tdc_tight_biased_p055", sim.sim_tdc, R, m=4000, mu=3.0, true_always_win=True, p_null_target=0.55)
    run("double_competition", sim.sim_double_competition, R)
    run("group_walk", sim.sim_groupwalk, R, frac_true=(0.05, 0.3, 0.8), mus=(2.0, 3.0, 4.0))
    run("reset_single", sim.sim_reset, R, mode="single")
    run("reset_two", sim.sim_reset, R, mode="two")
    run("reset_single_then_groupwalk_factor2_exploratory", sim.sim_reset, R, mode="single", then_groupwalk=True)
    # weak-initial-score regime: training often cannot start (fallback = initial score)
    run("reset_single_weak_init", sim.sim_reset, R, mode="single", m=2000, pi0=0.6, mu=(2.5, 2.0, 1.5, 0.0, 0.0))
    run("reset_two_weak_init", sim.sim_reset, R, mode="two", m=2000, pi0=0.6, mu=(2.5, 2.0, 1.5, 0.0, 0.0))
    run("reset_naive_selftrained_control", sim.sim_reset, min(R, 2000), m=500, naive=True, n_noise=40)
    run("picked_protein", sim.sim_picked_protein, R)
    run("picked_protein_highnull", sim.sim_picked_protein, R, null_psm_rate=0.9)
    run("picked_group_general", sim.sim_picked_group_general, args.R_group_general)
    run("multi_decoy_d3", sim.sim_multi_decoy, R, d=3)
    run("entrapment_estimators", sim.sim_entrapment_estimators, R)

    # entrapment decision rule demonstrations
    verdicts = {}
    for name, kw in {"fair_tight_mu3": {}, "fair_conservative_qscale0.6": {"q_scale": 0.6},
                     "biased_T_and_E_p055": {"p_target_vs_decoy": 0.55},
                     "biased_T_only_p055": {"p_target_vs_decoy": 0.55, "p_entrap_vs_decoy": 0.5}}.items():
        for a in (0.01, 0.05):
            verdicts[f"{name}@{a}"] = sim.sim_entrapment_verdict(alpha=a, **kw)
    metrics["entrapment_verdict_demo"] = verdicts

    # G2 power / sample size
    g2 = {"n_required_power0.90_rho0": g2_sample_size(0.9), "n_required_power0.92_rho0": g2_sample_size(0.92),
          "n_required_power0.90_rho0.5": g2_sample_size(0.9, rho=0.5), "mc": {}}
    for n in (16_600, 44_500, 110_000, g2_sample_size(0.9), g2_sample_size(0.92)):
        for rho in (0.0, 0.5):
            for pl in (0.5, 0.52, 0.55):
                key = f"n={n},rho={rho},p_learned={pl}"
                g2["mc"][key] = simulate_g2(n, R, 0.5, pl, rho=rho, seed=11)
                g2["mc"][key]["analytic"] = g2_power(n, pl - 0.5, rho=rho)
    # v0.2 absolute equivalence test at its stated n: pass probability for a fair scorer
    rng = np.random.default_rng(11)
    from vulcan.stats.g2 import wilson_ci
    k = rng.binomial(16_600, 0.5, size=R)
    _, lo, hi = wilson_ci(k, 16_600, 0.99)
    g2["v02_absolute_test_n16600_pass_rate_fair"] = float(((lo >= 0.49) & (hi <= 0.51)).mean())
    metrics["g2"] = g2

    with open(os.path.join(args.out, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=1, default=float)
    with open(os.path.join(args.out, "env.txt"), "w") as f:
        f.write(f"python {sys.version}\nnumpy {np.__version__}\nplatform {platform.platform()}\n")
        f.write(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}\n")
        f.write(f"hwinfo: {json.dumps(hw.fingerprint())}\n")
    with open(os.path.join(args.out, "config.yaml"), "w") as f:
        f.write(f"script: scripts/run_stats_sims.py\nreplicates: {R}\nreplicates_picked_group_general: "
                f"{args.R_group_general}\nseeds: {list(SEEDS)}\nalphas: {list(sim.ALPHAS)}\n")
    try:
        sha = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "vulcan/stats"], capture_output=True,
                               text=True).stdout.strip()
    except Exception:
        sha, dirty = "", ""
    with open(os.path.join(args.out, "git_sha"), "w") as f:
        f.write((sha or "no-commit") + (" (vulcan/stats uncommitted)" if dirty else "") + "\n")
    print("done", json.dumps(metrics["timing_s"]))


if __name__ == "__main__":
    main()
