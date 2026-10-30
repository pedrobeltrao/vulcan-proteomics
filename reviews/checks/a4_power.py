#!/usr/bin/env python3
"""Amendment A4 (draft): power, cluster-robustness and feasibility check of the fresh null set DN3.

SYNTHETIC SIMULATION + LABEL-FREE DN2 YIELDS ONLY.  Seed 131 (SEEDS[4]).
Inputs that are read (all label-free; none is a per-label, tail-label or gate quantity):
  * results/2026-10-27_A2_dn2/admission_counts.json -- per-project |A|, A with an exhaustive candidate, tail proxy
    (label-collapsed: max over T u D hyperscore >= host 1% threshold), B candidates;
  * data/cache/a2/admission/<acc>/<run>.npz -- only len(codes) and len(A_spec) (per-run spectra and |A|);
  * results/2026-10-27_A2_dn2/walk_events.jsonl -- only download sizes and timestamps (feasibility);
  * results/2026-10-27_A2_dn2/selection_log.tsv is NOT needed (walk decision counts are copied from NOTES.md).
No DN2 mixture exists, no label array is opened, no sealed file and no E7b output is read.

Rules simulated are the A4 rules (reviews/amendment_A4_draft.md), built from the implemented A2/A3 statistics:
  * tail-G2 = vulcan.gates.a2_stats.tail_g2_a2 unchanged in form, with conf 0.975 (A4 alpha split) and the cluster
    unit of A4.3 (genus/lab component) instead of the project;
  * G1 = A2.6 statistic with a 99.5% Fisher-z CI, hulled with a cluster-jackknife 99.5% t-interval (A4.5);
  * TOST rows (G3, V4-V6, V12) = Newcombe interval at conf 0.95 (two one-sided 97.5%), hulled with a cluster
    percentile bootstrap 95% interval (A4.5);
  * the A2 rule (conf 0.95/0.99/0.90, no cluster term outside tail-G2) is simulated next to it for comparison.
Sections (json keys):
  S0 dn2_yields      label-free per-run / per-project yields and their variance components (design-effect inputs)
  S1 n_min           A4 n_min of G1, G3, V6, V12 (definitions of A2.6 / A3.1.9, re-derived under the A4 rule)
  S2 tost_cluster    TOST rows: clean pass and margin false-pass, with between-cluster heterogeneity of the gap
  S3 g1_cluster      G1: clean pass and margin false-pass with heterogeneity of rho
  S1b V12_row        per-bin n at which the whole 10-bin V12 row passes (sets the A_cand stopping criterion)
  S4 walk            the DN3 walk under the A4 stopping rule (cluster cap 20% of the tail proxy): P(full stop),
                     P(v1 minimum), runs, admitted projects, clusters, downloaded files, |A|, A_cand, B, over a grid
                     of unknowns (median runs per project, protocol pass rate, tiers walked)
  S5 tail_g2         tail-G2 at the v1-minimum state (smallest DN3 on which gates run) and at the full stop: clean
                     pass, margin false-pass (tau = 0 / 0.05 / 0.10), A2 rule vs A4 rule, leaks concentrated in a
                     fraction of clusters, with and without the cluster cap; kappa (n_tail / proxy) in {1, 0.5} and
                     tail-proxy rate x {1, 0.5}
  (rev 2) dn2_tier0_projects(): DN2's admitted projects as tier 0 of the DN3 walk (D26 option (b)); the realistic
                     walk that uses it is reviews/checks/a4_r1_walk.py, which supersedes S4 for evaluability
  S6 joint           joint clean pass per scorer / tier (independence product, Frechet lower, min upper) at the v1
                     minimum and at the full stop
    conda_envs/pe/bin/python reviews/checks/a4_power.py   # writes reviews/checks/a4_power.json (+ .out via tee)
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
import time

import numpy as np
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates import a2_stats as A2S  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402
from vulcan.stats import g2 as G2  # noqa: E402

SEED = 131
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "a4_power.json")
QUICK = os.environ.get("A4_QUICK") == "1"

# ---------------------------------------------------------------------------------------------- A4 rule constants
TAIL_CONF_A4 = 0.975          # A4.6: one-sided level halved for re-tested rows (A2: 0.95)
G1_CONF_A4 = 0.995            # A2: 0.99 (two-sided CI)
TOST_CONF_A4 = 0.95           # A2: 0.90 (two one-sided 95%) -> two one-sided 97.5%
TOL = 0.01
G1_BOUND = 0.05
KEEP_V6 = 1 - 0.118           # dev neighbour-filter drop 11.8% (E7 D0), as A3.1.9
KEEP_V12 = 1 - 0.0981         # dev top-5 truncation, max-variance decile (A3 preconditions)
RECALL_G3 = 0.998             # E2

# ---------------------------------------------------------------------------------------------- S0 DN2 yields
def dn2_yields():
    ac = json.load(open(os.path.join(ROOT, "results/2026-10-27_A2_dn2/admission_counts.json")))
    proj = {}
    for acc, d in ac["per_project"].items():
        c = d["counts"]
        proj[acc] = {"runs": len(d["runs"]), "n_spectra": c["n_spectra"], "n_A": c["n_A"],
                     "n_A_cand": c["n_A_with_candidate"], "tail_proxy": c["tail_proxy"],
                     "n_B": c["n_B_candidates"]}
    runs = []
    for f in sorted(glob.glob(os.path.join(ROOT, "data/cache/a2/admission/*/*.npz"))):
        z = np.load(f)
        runs.append({"project": f.split("/")[-2], "run": os.path.basename(f)[:-4],
                     "n_spectra": int(len(z["codes"])), "n_A": int(len(z["A_spec"]))})
    A = np.array([r["n_A"] for r in runs], float)
    # variance components of log(|A| per run); runs with |A| = 0 (1 of 18) are a separate zero-inflation term
    pos = A > 0
    la = np.log(A[pos])
    pj = np.array([r["project"] for r in runs])[pos]
    u = sorted(set(pj))
    pm = np.array([la[pj == p].mean() for p in u])
    within = np.concatenate([la[pj == p] - la[pj == p].mean() for p in u])
    s_w = float(np.sqrt((within ** 2).sum() / (len(la) - len(u))))
    s_b = float(np.sqrt(max(pm.var(ddof=1) - s_w ** 2 / 3, 0)))
    zero_frac = float((~pos).mean())
    nA = np.array([proj[p]["n_A"] for p in u], float)
    tp = np.array([proj[p]["tail_proxy"] for p in u], float)
    ca = np.array([proj[p]["n_A_cand"] for p in u], float)
    tail_rate = tp / nA
    cand = ca / nA
    tot = {k: int(sum(proj[p][k] for p in u)) for k in ("runs", "n_spectra", "n_A", "n_A_cand", "tail_proxy", "n_B")}
    # ICC-type summary for the label-free counts: share of the largest project, effective number of projects
    def jeff(x):
        w = x / x.sum()
        return float(1 / (w ** 2).sum())
    ev = [json.loads(l) for l in open(os.path.join(ROOT, "results/2026-10-27_A2_dn2/walk_events.jsonl"))]
    sizes = [e["expected_size"] for e in ev if e["ev"] == "downloaded"]
    return {"projects": proj, "runs": runs, "totals": tot,
            "per_run_A": {"mean": float(A.mean()), "median": float(np.median(A)), "sd": float(A.std(ddof=1)),
                          "log_sd_between_projects": s_b, "log_sd_within_project": s_w,
                          "zero_run_fraction": zero_frac, "geo_mean_positive": float(np.exp(la.mean()))},
            "per_run_spectra_mean": float(np.mean([r["n_spectra"] for r in runs])),
            "B_over_A": tot["n_B"] / tot["n_A"],
            "tail_proxy_per_A": {"pooled": tot["tail_proxy"] / tot["n_A"], "per_project": tail_rate.tolist(),
                                 "log_sd_between_projects": float(np.log(tail_rate).std(ddof=1))},
            "A_cand_per_A": {"pooled": tot["n_A_cand"] / tot["n_A"], "per_project": cand.tolist(),
                             "logit_sd_between_projects": float(np.log(cand / (1 - cand)).std(ddof=1))},
            "effective_projects": {"A": jeff(nA), "tail_proxy": jeff(tp), "A_cand": jeff(ca)},
            "download": {"n_files": len(sizes), "mean_GB": float(np.mean(sizes) / 1e9),
                         "median_GB": float(np.median(sizes) / 1e9)},
            "walk_decisions_dn2": {"candidates": 58, "protocol_text_or_sdrf": 29, "genus_slot": 10,
                                   "marker_screen_excluded": 6, "admitted": 6, "no_reference_proteome": 3,
                                   "too_few_raws": 4, "source": "results/2026-10-27_A2_dn2/NOTES.md"}}


def dn2_tier0_projects():
    """Revision 2 (D26 option (b)): DN2's 6 admitted projects as walk tier 0.  LABEL-FREE inputs only:
    per-project |A| per run, A_cand/|A| and tail proxy/|A| from admission_counts.json (proxy rate with +0.5
    continuity, since counts are 1-10), and n_raw_records from DN2's own file-listing events (walk_events.jsonl)."""
    ac = json.load(open(os.path.join(ROOT, "results/2026-10-27_A2_dn2/admission_counts.json")))
    nraw = {}
    for l in open(os.path.join(ROOT, "results/2026-10-27_A2_dn2/walk_events.jsonl")):
        e = json.loads(l)
        if e["ev"] == "listed" and e.get("accession") in ac["per_project"]:
            nraw.setdefault(e["accession"], e["n_raw_records"])
    out = []
    for acc in ac["admitted"]:
        c = ac["per_project"][acc]["counts"]
        R = len(ac["per_project"][acc]["runs"])
        out.append({"accession": acc, "n_raw": int(nraw[acc]), "A_per_run": c["n_A"] / R,
                    "cand_rate": c["n_A_with_candidate"] / max(c["n_A"], 1),
                    "proxy_rate": (c["tail_proxy"] + 0.5) / max(c["n_A"], 1)})
    return out


# ---------------------------------------------------------------------------------------------- cluster weights
def cluster_weights(J, rng, log_sd=1.0):
    """Unequal cluster sizes: lognormal weights (log-sd 1.0 ~ DN2's between-project spread of |A| per run, damped
    by the per-project run cap)."""
    w = rng.lognormal(0.0, log_sd, size=J)
    return w / w.sum()


# ---------------------------------------------------------------------------------------------- TOST + cluster hull
def tost_hull_pass(kT, nT, kD, nD, conf, B, rng, tol=TOL, use_cluster=True):
    """kT/nT/kD/nD: per-cluster arrays.  A4 rule: Newcombe(conf) on totals, hulled with the cluster percentile
    bootstrap interval at the same two-sided level; pass iff the hull is inside [-tol, tol]."""
    KT, NT, KD, ND = kT.sum(), nT.sum(), kD.sum(), nD.sum()
    if NT == 0 or ND == 0:
        return False
    _, (lo, hi) = GS.newcombe_diff(int(KT), int(NT), int(KD), int(ND), conf)
    if use_cluster:
        J = kT.size
        idx = rng.integers(0, J, size=(B, J))
        sNT, sND = nT[idx].sum(1), nD[idx].sum(1)
        d = kT[idx].sum(1) / np.maximum(sNT, 1) - kD[idx].sum(1) / np.maximum(sND, 1)
        d = np.where((sNT > 0) & (sND > 0), d, 1.0)
        a = (1 - conf) / 2
        lo, hi = min(lo, float(np.quantile(d, a))), max(hi, float(np.quantile(d, 1 - a)))
    return bool(lo >= -tol and hi <= tol)


def tost_pass_prob(n, p_d, gap, tau, J, conf, n_rep, rng, B=400, log_sd=1.0, use_cluster=True):
    """n per label (total over clusters); target rate p_d + gap_j with gap_j ~ N(gap, tau) per cluster."""
    ok = 0
    for _ in range(n_rep):
        w = cluster_weights(J, rng, log_sd)
        nT = rng.multinomial(n, w)
        nD = rng.multinomial(n, w)
        g = gap + (rng.normal(0, tau, J) if tau > 0 else 0.0)
        pt = np.clip(p_d + g, 0, 1)
        kT = rng.binomial(nT, pt)
        kD = rng.binomial(nD, p_d)
        ok += tost_hull_pass(kT.astype(float), nT.astype(float), kD.astype(float), nD.astype(float), conf, B, rng,
                             use_cluster=use_cluster)
    return ok / n_rep


def n_min_tost(p, J, conf, n_rep, seed, target=0.90, lo=500, hi=200_000, step=100, use_cluster=True):
    def prob(n):
        return tost_pass_prob(n, p, 0.0, 0.0, J, conf, n_rep, np.random.default_rng(seed), use_cluster=use_cluster)
    if prob(hi) < target:
        return {"n_min": None, "note": f"not reached at {hi}"}
    a, b = lo, hi
    while b - a > step:
        m = (a + b) // 2
        if prob(m) >= target:
            b = m
        else:
            a = m
    b = int(math.ceil(b / step) * step)
    while prob(b) < target:
        b += step
    return {"n_min": b, "pass_prob_at_n_min": prob(b), "p": p, "J": J, "conf": conf, "n_rep": n_rep,
            "cluster_hull": use_cluster}


# ---------------------------------------------------------------------------------------------- G1 + jackknife
def spearman_fz(x, y):
    r = stats.spearmanr(x, y).statistic
    return float(r), math.atanh(max(min(r, 0.999999), -0.999999))


def g1_a4_pass(x, y, cl, J, conf=G1_CONF_A4, bound=G1_BOUND, use_cluster=True):
    n = x.size
    r, z = spearman_fz(x, y)
    se = math.sqrt((1 + r * r / 2) / (n - 3))
    h = stats.norm.ppf(0.5 + conf / 2) * se
    lo, hi = z - h, z + h
    if use_cluster:
        zs = []
        for j in range(J):
            m = cl != j
            zs.append(spearman_fz(x[m], y[m])[1])
        zs = np.array(zs)
        se_j = math.sqrt((J - 1) / J * ((zs - zs.mean()) ** 2).sum())
        hj = stats.t.ppf(0.5 + conf / 2, J - 1) * se_j
        lo, hi = min(lo, z - hj), max(hi, z + hj)
    return bool(math.tanh(lo) > -bound and math.tanh(hi) < bound)


def g1_pass_prob(n_top, rho, tau, J, n_rep, rng, conf=G1_CONF_A4, use_cluster=True, log_sd=1.0):
    ok = 0
    for _ in range(n_rep):
        w = cluster_weights(J, rng, log_sd)
        nj = rng.multinomial(n_top, w)
        cl = np.repeat(np.arange(J), nj)
        rj = np.clip(rho + (rng.normal(0, tau, J) if tau > 0 else np.zeros(J)), -0.9, 0.9)
        rl = 2 * np.sin(np.pi * rj / 6)
        x = rng.standard_normal(n_top)
        e = rng.standard_normal(n_top)
        rr = rl[cl]
        y = rr * x + np.sqrt(1 - rr * rr) * e
        ok += g1_a4_pass(x, y, cl, J, conf, use_cluster=use_cluster)
    return ok / n_rep


def n_min_g1(J, n_rep, seed, target=0.90):
    def prob(n):
        return g1_pass_prob(n, 0.0, 0.0, J, n_rep, np.random.default_rng(seed))
    a, b = 4000, 40_000
    if prob(b) < target:
        return {"n_min": None}
    while b - a > 200:
        m = (a + b) // 2
        if prob(m) >= target:
            b = m
        else:
            a = m
    b = int(math.ceil(b / 100) * 100)
    while prob(b) < target:
        b += 200
    return {"n_min": b, "pass_prob_at_n_min": prob(b), "J": J, "conf": G1_CONF_A4, "n_rep": n_rep}


# ---------------------------------------------------------------------------------------------- S4 walk model
def cap_factors(c, share):
    """A4.3.9 cluster cap: per-cluster retention factors f_j in (0, 1] such that no cluster holds more than `share`
    of the capped total of c (here the label-collapsed tail proxy).  Applied as a seeded uniform subsample of that
    cluster's A spectra (and, with B's own shares, of its B spectra), so it is label-free (P1/P2)."""
    c = np.asarray(c, float)
    if c.sum() <= 0 or c.size * share < 1:
        return np.ones_like(c)
    x = c.copy()
    for _ in range(200):
        L = share * x.sum()
        if (x <= L * (1 + 1e-9)).all():
            break
        x = np.minimum(x, L)
    return np.where(c > 0, x / np.maximum(c, 1e-12), 1.0)


def capped_state(st_, share):
    J = len(st_["clusters"])
    pc = np.asarray(st_["proj_cluster"])
    cp = np.bincount(pc, weights=st_["proj_proxy"], minlength=J)
    ca = np.bincount(pc, weights=st_["proj_A_cand"], minlength=J)
    f = cap_factors(cp, share)
    return {"J": J, "proxy": float((cp * f).sum()), "A_cand": float((ca * f).sum()), "cluster_proxy": cp * f,
            "factors": f}


def _criteria_met(st_, c):
    cs = capped_state(st_, c["cap_share"])
    return cs["J"] >= c["J_min"] and cs["proxy"] >= c["proxy_min"] and cs["A_cand"] >= c["A_cand_min"]


def simulate_walk(rng, y, *, pools, med_runs, p_protocol, p_marker, p_files, run_cap, kappa, crit_full, crit_min,
                  max_runs, max_proj_per_cluster=3, link_prob=0.15, tail_scale=1.0):
    """One DN3 walk over the tiers in `pools` (candidate counts, walked in order).  Each candidate: metadata pass
    w.p. p_protocol, >= 1 usable RAW w.p. p_files, marker screen pass w.p. p_marker (stage 1 on <= 3 runs: an
    excluded project costs min(R, 3) downloads).  Admitted project: R = min(run_cap, max(1, round(lognormal(med_runs,
    0.8)))) runs; per-run |A| lognormal with the project mean drawn from DN2's between-project spread and DN2's
    within-project spread; per-project tail-proxy and candidate rates from DN2's between-project spread.  A project
    joins an existing cluster (genus/lab link) w.p. link_prob if that cluster has < max_proj_per_cluster projects.
    The walk stops when crit_full holds (A4.3.7 stopping rule), at the run cap, or at pool exhaustion.  The first
    state at which crit_min (v1 minimum, A4.3.8) holds is kept as a snapshot."""
    mA = math.log(y["per_run_A"]["geo_mean_positive"])
    sb, sw = y["per_run_A"]["log_sd_between_projects"], y["per_run_A"]["log_sd_within_project"]
    z0 = y["per_run_A"]["zero_run_fraction"]
    lt = math.log(y["tail_proxy_per_A"]["pooled"] * tail_scale)
    st = y["tail_proxy_per_A"]["log_sd_between_projects"]
    lc = math.log(y["A_cand_per_A"]["pooled"] / (1 - y["A_cand_per_A"]["pooled"]))
    sc = y["A_cand_per_A"]["logit_sd_between_projects"]
    spr = y["per_run_spectra_mean"]
    st_ = {"cand": 0, "admitted": 0, "runs": 0, "dl_files": 0, "A": 0.0, "A_cand": 0.0, "proxy": 0.0, "B": 0.0,
           "clusters": [], "proj_cluster": [], "proj_A_cand": [], "proj_proxy": [], "proj_B": [],
           "stopped": False, "met_min": False, "snap_min": None, "end": "pool exhausted"}
    n_pool = sum(pools)
    while st_["cand"] < n_pool:
        st_["cand"] += 1
        if rng.random() > p_protocol or rng.random() > p_files:
            continue
        R = int(min(run_cap, max(1, round(rng.lognormal(math.log(med_runs), 0.8)))))
        if rng.random() > p_marker:
            st_["dl_files"] += min(R, 3)
            continue
        if st_["runs"] + R > max_runs:
            st_["end"] = "run cap"
            break
        st_["dl_files"] += R
        mu = rng.normal(mA, sb)
        a = float((rng.lognormal(mu, sw, size=R) * (rng.random(R) > z0)).sum())
        tr = math.exp(rng.normal(lt, st))
        cr = 1 / (1 + math.exp(-rng.normal(lc, sc)))
        proxy = rng.poisson(a * tr)
        acand = rng.binomial(int(a), cr)
        b = R * spr - a
        open_cl = [k for k, c in enumerate(st_["clusters"]) if c < max_proj_per_cluster]
        if open_cl and rng.random() < link_prob:
            k = int(rng.choice(open_cl))
            st_["clusters"][k] += 1
        else:
            st_["clusters"].append(1)
            k = len(st_["clusters"]) - 1
        st_["admitted"] += 1
        st_["runs"] += R
        st_["A"] += a
        st_["A_cand"] += acand
        st_["proxy"] += proxy
        st_["B"] += b
        st_["proj_cluster"].append(k)
        st_["proj_A_cand"].append(acand)
        st_["proj_proxy"].append(proxy)
        st_["proj_B"].append(b)
        if not st_["met_min"] and _criteria_met(st_, crit_min):
            st_["met_min"] = True
            cs = capped_state(st_, crit_min["cap_share"])
            st_["snap_min"] = {"runs": st_["runs"], "J": cs["J"], "proxy": cs["proxy"], "A_cand": cs["A_cand"],
                               "proj_cluster": list(range(cs["J"])), "proj_proxy": list(cs["cluster_proxy"]),
                               "dl_files": st_["dl_files"],
                               "uncapped": {"J": cs["J"], "proj_cluster": list(range(cs["J"])),
                                            "proj_proxy": list(cs["cluster_proxy"] / cs["factors"])}}
        if _criteria_met(st_, crit_full):
            st_["stopped"] = True
            st_["end"] = "stopping rule met"
            break
    st_["J"] = len(st_["clusters"])
    st_["kappa"] = kappa
    return st_


def tail_g2_at(walk, rng, p_learned_mean, tau, kappa, conf, B, concentrated=None, cap_share=None):
    """One tail-G2 evaluation on a walk outcome: per-cluster tails ~ Poisson(kappa * proxy); learned and fixed
    tails independent (conservative, as A2.10 S1/S3)."""
    J = walk["J"]
    cp = np.bincount(walk["proj_cluster"], weights=walk["proj_proxy"], minlength=J)
    if cap_share is not None:
        cp = cp * cap_factors(cp, cap_share)
    nL = rng.poisson(kappa * cp)
    nF = rng.poisson(kappa * cp)
    if concentrated is not None:                    # leak only in a fraction of clusters (by count)
        frac, p_hi = concentrated
        pl = np.full(J, 0.5)
        pl[rng.permutation(J)[:max(1, int(round(frac * J)))]] = p_hi
    elif tau > 0:
        s = p_learned_mean * (1 - p_learned_mean) / tau ** 2 - 1
        pl = rng.beta(p_learned_mean * s, (1 - p_learned_mean) * s, size=J)
    else:
        pl = np.full(J, p_learned_mean)
    kL = rng.binomial(nL, pl)
    kF = rng.binomial(nF, 0.5)
    cl = np.concatenate([np.full(nL[j] + nF[j], j) for j in range(J)])
    tl = np.concatenate([np.r_[np.ones(nL[j], bool), np.zeros(nF[j], bool)] for j in range(J)])
    gl = np.concatenate([np.r_[np.arange(nL[j]) < kL[j], np.zeros(nF[j], bool)] for j in range(J)])
    tf = np.concatenate([np.r_[np.zeros(nL[j], bool), np.ones(nF[j], bool)] for j in range(J)])
    gf = np.concatenate([np.r_[np.zeros(nL[j], bool), np.arange(nF[j]) < kF[j]] for j in range(J)])
    r = A2S.tail_g2_a2(cl, tl, gl, tf, gf, B=B, seed=rng, conf=conf)
    return r["verdict"], r["n_tail"]


# ---------------------------------------------------------------------------------------------- main
def main():
    t0 = time.perf_counter()
    rng0 = np.random.default_rng(SEED)
    out = {"seed": SEED, "spec": __doc__.split("    conda_envs")[0].strip()}
    y = dn2_yields()
    out["S0_dn2_yields"] = {k: v for k, v in y.items() if k not in ("runs",)}
    out["S0_dn2_yields"]["runs"] = y["runs"]
    print("S0 per-run |A|:", json.dumps(y["per_run_A"]), "tail/A", round(y["tail_proxy_per_A"]["pooled"], 6),
          "cand/A", round(y["A_cand_per_A"]["pooled"], 3), "B/A", round(y["B_over_A"], 1),
          "J_eff (A, tail, cand):", json.dumps(y["effective_projects"]), flush=True)

    # ---------------- S1 n_min under the A4 rule (J = 20 clusters, lognormal weights) and A2 rule for comparison
    J_REF = 20
    nrep_nmin = 300 if QUICK else 1000
    nmin = {}
    for name, p in (("V6", KEEP_V6), ("V12", KEEP_V12)):
        nmin[name] = {"A4": n_min_tost(p, J_REF, TOST_CONF_A4, nrep_nmin, SEED + 1),
                      "A4_no_cluster_term": n_min_tost(p, J_REF, TOST_CONF_A4, nrep_nmin, SEED + 1,
                                                       use_cluster=False),
                      "A3_reference": {"V6": 22700, "V12": 19200}[name]}
        print("S1", name, json.dumps(nmin[name]), flush=True)
    g3 = {}
    for n in (2000, 3000):
        g3[n] = tost_pass_prob(n, RECALL_G3, 0.0, 0.0, J_REF, TOST_CONF_A4, nrep_nmin, np.random.default_rng(SEED + 2))
    nmin["G3"] = {"A2_n_min": 2000, "A4_pass_prob_one_bin_at_recall_0.998": g3,
                  "rule": "A4 keeps 2,000 iff one-bin pass >= 0.99 (10 bins must pass); else raise"}
    print("S1 G3", json.dumps(nmin["G3"]), flush=True)
    nmin["G1"] = n_min_g1(J_REF, 150 if QUICK else 400, SEED + 3)
    nmin["G1"]["A2_n_min"] = 7200
    print("S1 G1", json.dumps(nmin["G1"]), flush=True)
    out["S1_n_min"] = nmin

    # ---------------- S2 TOST rows: clean pass and false pass at the margin, with heterogeneity
    s2 = {}
    nrep = 300 if QUICK else 1000
    for name, p in (("V6", KEEP_V6), ("V12", KEEP_V12), ("G3", RECALL_G3)):
        n = nmin[name]["A4"]["n_min"] if name in ("V6", "V12") else 2000
        cell = {"n_per_label_per_bin": n}
        for J in (10, 20, 40):
            r = np.random.default_rng(SEED + 10 + J)
            cell[f"J{J}"] = {
                "clean": tost_pass_prob(n, p, 0.0, 0.0, J, TOST_CONF_A4, nrep, r),
                "gap0.01_tau0": tost_pass_prob(n, p, -TOL, 0.0, J, TOST_CONF_A4, nrep, r),
                "gap0.01_tau0.005": tost_pass_prob(n, p, -TOL, 0.005, J, TOST_CONF_A4, nrep, r),
                "A2_rule_gap0.01_tau0.005": tost_pass_prob(n, p, -TOL, 0.005, J, 0.90, nrep, r, use_cluster=False),
            }
        s2[name] = cell
        print("S2", name, json.dumps(cell), flush=True)
    out["S2_tost_cluster"] = s2

    # ---------------- S3 G1
    s3 = {}
    n_top = max(nmin["G1"]["n_min"] or 9000, 9000)
    nrep_g1 = 150 if QUICK else 300
    for J in (10, 20, 40):
        r = np.random.default_rng(SEED + 30 + J)
        s3[f"J{J}"] = {"n_top": n_top,
                       "clean": g1_pass_prob(n_top, 0.0, 0.0, J, nrep_g1, r),
                       "rho0.05_tau0": g1_pass_prob(n_top, 0.05, 0.0, J, nrep_g1, r),
                       "rho0.05_tau0.03": g1_pass_prob(n_top, 0.05, 0.03, J, nrep_g1, r),
                       "A2_rule_rho0.05_tau0.03": g1_pass_prob(n_top, 0.05, 0.03, J, nrep_g1, r, conf=0.99,
                                                               use_cluster=False)}
        print("S3", f"J{J}", json.dumps(s3[f"J{J}"]), flush=True)
    out["S3_g1_cluster"] = s3

    # ---------------- S1b V12 row-level size (all ~10 bins must pass) -> A_cand criterion
    r12 = {}
    for n in (30_000, 35_000, 40_000, 45_000):
        pb = tost_pass_prob(n, KEEP_V12, 0.0, 0.0, J_REF, TOST_CONF_A4, 500 if QUICK else 2000,
                            np.random.default_rng(SEED + 5))
        r12[n] = {"one_bin": pb, "row_10_bins_indep": pb ** 10}
    n_row = min([n for n, v in r12.items() if v["row_10_bins_indep"] >= 0.80], default=45_000)
    out["S1_n_min"]["V12_row"] = {"grid": r12, "n_per_label_per_bin_for_row_pass_0.80": n_row,
                                  "note": "all bins at the max-variance dev rate (conservative); bins independent"}
    print("S1b V12 row", json.dumps(out["S1_n_min"]["V12_row"]), flush=True)
    a_cand_full = int(math.ceil(10 * n_row / 0.9 / 10_000) * 10_000)
    g1_n = nmin["G1"]["n_min"] or 11_100
    a_cand_min = int(math.ceil(max(10 * g1_n, 120_000) / 10_000) * 10_000)

    # ---------------- S4 walk under the A4 stopping rule
    crit_full = {"J_min": 20, "proxy_min": 1200, "A_cand_min": a_cand_full, "cap_share": 0.20}
    crit_min = {"J_min": 20, "proxy_min": 1200, "A_cand_min": a_cand_min, "cap_share": 0.20}
    MAX_RUNS = 2000
    out["S4_stopping_rule"] = {"full_stop": crit_full, "v1_minimum": crit_min, "max_runs": MAX_RUNS,
                               "run_cap_per_project": 24, "max_projects_per_cluster": 3}
    print("S4 criteria", json.dumps(out["S4_stopping_rule"]), flush=True)
    POOLS = {"I": [198], "I+II": [198, 259], "I+II+III": [198, 259, 135]}
    grid = []
    nwalk = 200 if QUICK else 1000
    for med_runs in (6, 12, 24):
        for p_protocol in (0.45, 0.7):
            for pool_name in ("I", "I+II", "I+II+III"):
                r = np.random.default_rng(SEED + 40)
                ws = [simulate_walk(r, y, pools=POOLS[pool_name], med_runs=med_runs, p_protocol=p_protocol,
                                    p_marker=0.5, p_files=0.93, run_cap=24, kappa=1.0, crit_full=crit_full,
                                    crit_min=crit_min, max_runs=MAX_RUNS) for _ in range(nwalk)]
                def q(ws_, key):
                    return [float(np.percentile([w[key] for w in ws_], pp)) for pp in (10, 50, 90)] if ws_ else None
                full = [w for w in ws if w["stopped"]]
                cell = {"med_runs_per_project": med_runs, "p_protocol_pass": p_protocol, "p_marker_pass": 0.5,
                        "pool": pool_name, "pool_sizes": POOLS[pool_name],
                        "P_full_stop": len(full) / nwalk, "P_v1_minimum": float(np.mean([w["met_min"] for w in ws])),
                        "runs_at_end_p10_50_90": q(ws, "runs"), "admitted_p10_50_90": q(ws, "admitted"),
                        "clusters_p10_50_90": q(ws, "J"), "candidates_walked_p10_50_90": q(ws, "cand"),
                        "A_p10_50_90": q(ws, "A"), "A_cand_p10_50_90": q(ws, "A_cand"),
                        "proxy_p10_50_90": q(ws, "proxy"), "B_p10_50_90": q(ws, "B"),
                        "downloaded_files_p10_50_90": q(ws, "dl_files"),
                        "runs_when_v1_minimum_first_met_p10_50_90":
                            [float(np.percentile([w["snap_min"]["runs"] for w in ws if w["met_min"]], pp))
                             for pp in (10, 50, 90)] if any(w["met_min"] for w in ws) else None}
                grid.append(cell)
                print("S4", json.dumps({k: cell[k] for k in ("med_runs_per_project", "p_protocol_pass", "pool",
                                                             "P_full_stop", "P_v1_minimum", "runs_at_end_p10_50_90",
                                                             "clusters_p10_50_90", "downloaded_files_p10_50_90",
                                                             "runs_when_v1_minimum_first_met_p10_50_90")}),
                      flush=True)
    out["S4_walk"] = grid

    # ---------------- S5 tail-G2 at the v1-minimum point (the smallest DN3 on which gates run)
    s5 = {}
    nrep_t = 150 if QUICK else 400
    Bt = 1000 if QUICK else 2000
    for tail_scale in (1.0, 0.5):
        for kappa in (1.0, 0.5):
            r = np.random.default_rng(SEED + 50)
            walks, fulls = [], []
            while len(walks) < nrep_t:
                w = simulate_walk(r, y, pools=POOLS["I+II+III"], med_runs=12, p_protocol=0.45, p_marker=0.5,
                                  p_files=0.93, run_cap=24, kappa=kappa, crit_full=crit_full, crit_min=crit_min,
                                  max_runs=MAX_RUNS, tail_scale=tail_scale)
                if w["met_min"]:
                    walks.append(w["snap_min"])
                if w["stopped"]:
                    cs = capped_state(w, crit_full["cap_share"])
                    fulls.append({"J": cs["J"], "proj_cluster": list(range(cs["J"])),
                                  "proj_proxy": list(cs["cluster_proxy"])})
            cell = {}
            for lab, pm, tau, conf, conc in (("clean_A4", 0.5, 0.0, TAIL_CONF_A4, None),
                                              ("clean_A2rule", 0.5, 0.0, 0.95, None),
                                              ("margin0.60_tau0_A4", 0.60, 0.0, TAIL_CONF_A4, None),
                                              ("margin0.60_tau0.05_A4", 0.60, 0.05, TAIL_CONF_A4, None),
                                              ("margin0.60_tau0.10_A4", 0.60, 0.10, TAIL_CONF_A4, None),
                                              ("margin0.60_tau0.10_A2rule", 0.60, 0.10, 0.95, None),
                                              ("mean0.60_leak0.90_in_25pct_clusters_A4", None, None, TAIL_CONF_A4,
                                               (0.25, 0.90)),
                                              ("mean0.60_leak0.90_in_25pct_clusters_A2rule", None, None, 0.95,
                                               (0.25, 0.90)),
                                              ("mean0.60_leak1.00_in_20pct_clusters_A4", None, None, TAIL_CONF_A4,
                                               (0.20, 1.00))):
                v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
                nt = []
                for w in walks:
                    verdict, n_tail = tail_g2_at(w, r, pm if pm is not None else 0.5, tau or 0.0, kappa, conf, Bt,
                                                 concentrated=conc)
                    v[verdict] += 1
                    nt.append(n_tail)
                cell[lab] = {k: x / len(walks) for k, x in v.items()} | {"n_tail_median": float(np.median(nt))}
            for lab, conc in (("UNCAPPED_mean0.60_leak0.90_in_25pct_clusters_A4", (0.25, 0.90)),
                              ("UNCAPPED_clean_A4", None)):
                v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
                for w in walks:
                    verdict, _ = tail_g2_at(w["uncapped"], r, 0.5, 0.0, kappa, TAIL_CONF_A4, Bt, concentrated=conc)
                    v[verdict] += 1
                cell[lab] = {k: x / len(walks) for k, x in v.items()}
            v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
            for w in fulls:
                verdict, _ = tail_g2_at(w, r, 0.5, 0.0, kappa, TAIL_CONF_A4, Bt)
                v[verdict] += 1
            cell["FULL_STOP_clean_A4"] = {k: x / max(len(fulls), 1) for k, x in v.items()} | {"n_walks": len(fulls)}
            v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
            for w in fulls:
                verdict, _ = tail_g2_at(w, r, 0.60, 0.10, kappa, TAIL_CONF_A4, Bt)
                v[verdict] += 1
            cell["FULL_STOP_margin0.60_tau0.10_A4"] = {k: x / max(len(fulls), 1) for k, x in v.items()}
            cell["J_median"] = float(np.median([w["J"] for w in walks]))
            cell["runs_median"] = float(np.median([w["runs"] for w in walks]))
            s5[f"tail_scale{tail_scale}_kappa{kappa}"] = cell
            print("S5", f"tail_scale{tail_scale}_kappa{kappa}", json.dumps(cell), flush=True)
    out["S5_tail_g2"] = s5

    # ---------------- S6 joint bounds: at the v1-minimum DN3 (worst case on which gates run) and at the full stop
    rr = np.random.default_rng(SEED + 60)
    nrep6 = 150 if QUICK else 400
    g1_min = g1_pass_prob(max(a_cand_min // 10, nmin["G1"]["n_min"] or 0), 0.0, 0.0, 20, nrep6, rr) * 0.99
    g1_full = g1_pass_prob(a_cand_full // 10, 0.0, 0.0, 20, nrep6 // 2, rr) * 0.99
    g3c = s2["G3"]["J20"]["clean"] ** 10
    v6n = 4 * nmin["V6"]["A4"]["n_min"]
    v6c = tost_pass_prob(v6n, KEEP_V6, 0.0, 0.0, 20, TOST_CONF_A4, 500 if QUICK else 1000, rr) ** 2
    v12c = out["S1_n_min"]["V12_row"]["grid"][n_row]["row_10_bins_indep"]
    g2c = G2.g2_power(10 ** 6, 0.0)["p_pass_joint_lb"]
    k1 = "tail_scale1.0_kappa1.0"
    k5 = "tail_scale1.0_kappa0.5"

    def bounds(d):
        ps = list(d.values())
        return {"rows": d, "independence": float(np.prod(ps)), "frechet_lower": float(max(0, sum(ps) - (len(ps) - 1))),
                "upper_min": float(min(ps))}
    res = {}
    for lab, k in (("kappa1", k1), ("kappa0.5", k5)):
        tmin = s5[k]["clean_A4"]["pass"]
        tfull = s5[k]["FULL_STOP_clean_A4"]["pass"]
        res[lab] = {
            "v1_at_v1_minimum": bounds({"tail_G2": tmin, "G1": g1_min, "G3_A": g3c, "G3_AuB": g3c, "G2_original": g2c,
                                        "G5": 1.0}),
            "v1_at_full_stop": bounds({"tail_G2": tfull, "G1": g1_full, "G3_A": g3c, "G3_AuB": g3c,
                                       "G2_original": g2c, "G5": 1.0}),
            "v2_sep_tier1_at_full_stop": bounds({"V1": tfull, "V2": g1_full, "V3": g3c, "V4": g3c, "V6_class1": v6c,
                                                 "V7": g2c, "V9ab": 1.0, "V10": 1.0, "V12": v12c}),
            "v2_sep_tier2_at_full_stop": bounds({"V1": tfull, "V2": g1_full, "V3": g3c, "V4": g3c, "V6_class1": v6c,
                                                 "V7": g2c, "V9ab": 1.0, "V10": 1.0, "V12": v12c, "V5": g3c,
                                                 "V6_class2": v6c, "V8": g2c, "V9c": 1.0, "V13": tfull,
                                                 "V14": g1_full})}
    res["assumptions"] = ("G1 at n_top = A_cand/10 (J = 20); V6 at 4x n_min (B/A ~ 40 in DN2 makes V6 large); V12 row at "
                          "the A_cand criterion; G5 at n >= 1e5 passes ~1; FS rows never block when not evaluable; "
                          "V13 taken equal to V1 (A-syn tail fraction unknown); rows treated as independent for the "
                          "product; Frechet/min bounds hold under any dependence")
    out["S6_joint"] = res
    print("S6", json.dumps({k: {kk: {x: round(vv[x], 3) for x in ("independence", "frechet_lower", "upper_min")}
                                 for kk, vv in v.items()} for k, v in res.items() if k != "assumptions"}), flush=True)
    out["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print("wrote", OUT, "wall", round(out["wall_s"]), "s")


if __name__ == "__main__":
    main()
