#!/usr/bin/env python3
"""C9 pilot step 3: projected evaluability and gate power of a human removed-peptide null set, A4 rules.

SYNTHETIC SIMULATION + LABEL-FREE PILOT YIELDS + METADATA ONLY.  Seed 131 (SEEDS[4], as the A4 checks).
Inputs (all label-free):
  * results/2026-10-28_C9_pilot/per_run/*.json -- per run: certified count, certified with a single MS1 feature,
    A_cand (certified with >= 1 candidate after removal), tail proxy (label-collapsed max >= host 1% threshold);
  * results/2026-10-28_C9_pilot/pool_count.json -- the eligible human pool (seeded walk order over the real projects
    and their normalised lab keys) and E9's label-free usable-RAW-per-project listing counts.
Machinery reused unchanged: reviews/checks/a4_power.py (cap_factors, g1_pass_prob, tost_pass_prob) and
vulcan.gates.a2_stats.tail_g2_a2 (tail-G2 with CP + cluster bootstrap + Wilson validity control), under the A4 levels
(tail 0.975, G1 0.995, TOST 0.95) and the A4 intersection-union rule over capped and uncapped pools.
Rules simulated (A4.3.6/A4.3.7 carried over; nothing relaxed):
  * cluster = normalised crawl lab; <= 5 projects per cluster; 10% cap = min of the water-filling factors of the tail
    proxy, A_cand and |A| per cluster;
  * v1 minimum: J >= 20, J_eff >= 10 (capped proxy and capped A_cand), capped proxy >= 1,200, capped A_cand >= 120,000;
    full stop (v2_sep tier 1): the same with capped A_cand >= 450,000;
  * R_cap runs per project (2 or 3, C9 design), drawn from that project's usable RAWs (E9 listing counts).
Scenarios for tail-G2: clean; margin 0.60 (tau 0.10 between clusters); 25% of clusters leak at 0.90; and RESIDUAL
CONTAMINATION (unremoved co-isolates / homologs, pro-target for every scorer): hyperscore tail target fraction 0.5+d,
learned 0.5+1.5d, d in {0.02, 0.05}.
    conda_envs/pe/bin/python scripts/c9/c9_power.py   # -> results/2026-10-28_C9_pilot/c9_power.json
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import random
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "reviews", "checks"))
sys.path.insert(0, ROOT)
import a4_power as AP  # noqa: E402
from vulcan.gates import a2_stats as A2S  # noqa: E402

SEED = 131
RES = os.path.join(ROOT, "results", "2026-10-28_C9_pilot")
OUT = os.path.join(RES, "c9_power.json")
QUICK = os.environ.get("C9_QUICK") == "1"
CAP_SHARE, PER_CLUSTER = 0.10, 5
CRIT_MIN = {"J_min": 20, "Jeff_min": 10, "proxy_min": 1200, "A_cand_min": 120_000}
CRIT_FULL = dict(CRIT_MIN, A_cand_min=450_000)
P_DECODE_FAIL = 2 / 1000          # E9: 2 of 1,000 runs failed to decode
P_PROJECT_DEAD = 2 / 142          # E9: 2 of 142 projects unreadable / 0 PSMs


# ---------------------------------------------------------------------------------------------- pilot yields
def pilot_yields(variant: str, thr: str):
    """variant: 'f0.20' (A2's ladder threshold), 'f0.05' (strict ladder screen: a surviving candidate shares
    < 5% of p's b/y ions), 'f0.20_single_ms1' (also <= 1 MS1 feature in the isolation window).
    thr: 'own' (the run's own hyperscore 1% threshold) or 'D0max' (A2/A4 proxy convention; f0.20 only)."""
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(RES, "per_run", "*.json")))
            if not os.path.basename(f).startswith("_")]
    out = []
    for r in recs:
        A = r["n_certified"]
        if variant == "f0.20":
            ac = r["A_cand"]
            tp = r["tail_proxy_own_thr"] if thr == "own" else r["tail_proxy_D0max_thr"]
        elif variant == "f0.05":
            assert thr == "own"
            ac = r["diag_ladder_threshold"]["0.05"]["A_cand"]
            tp = r["diag_ladder_threshold"]["0.05"]["tail_proxy_own_thr"]
        elif variant == "f0.20_single_ms1":
            assert thr == "own"
            A, ac = r["certified_single_ms1_feature"], r["A_cand_single_ms1_feature"]
            tp = r["tail_proxy_own_thr_single_ms1_feature"]
        out.append({"project": r["project"], "run": r["run"], "A": A, "A_cand": ac, "proxy": tp})
    zero = [o for o in out if o["A"] == 0]
    out_pos = [o for o in out if o["A"] > 0]
    z0 = len(zero) / len(out)
    out_all, out = out, out_pos
    proj = sorted({o["project"] for o in out})
    la = np.log([max(o["A"], 1) for o in out])
    pj = np.array([o["project"] for o in out])
    pm = np.array([la[pj == p].mean() for p in proj])
    wi = np.concatenate([la[pj == p] - la[pj == p].mean() for p in proj])
    dfw = len(la) - len(proj)
    s_w = float(np.sqrt((wi ** 2).sum() / dfw)) if dfw > 0 else 0.32
    s_w = max(s_w, 0.32)                      # DN2's within-project log-sd as a floor (one 2-run project here)
    s_b = float(max(pm.std(ddof=1), 0.0))
    agg = {p: {k: sum(o[k] for o in out if o["project"] == p) for k in ("A", "A_cand", "proxy")} for p in proj}
    cr = np.array([agg[p]["A_cand"] / agg[p]["A"] for p in proj])
    tr = np.array([(agg[p]["proxy"] + 0.5) / max(agg[p]["A_cand"], 1) for p in proj])
    return {"variant": variant, "thr": thr, "runs": out_all, "n_runs": len(out_all), "n_projects": len(proj),
            "zero_run_fraction": z0,
            "log_A_mean_project": float(pm.mean()), "log_A_sd_between": s_b, "log_A_sd_within": s_w,
            "cand_rate_per_project": cr.tolist(), "cand_logit_mean": float(np.mean(np.log(cr / (1 - cr)))),
            "cand_logit_sd": float(np.std(np.log(cr / (1 - cr)), ddof=1)),
            "proxy_per_cand_per_project": tr.tolist(), "proxy_log_mean": float(np.mean(np.log(tr))),
            "proxy_log_sd": float(np.std(np.log(tr), ddof=1)),
            "totals": {k: int(sum(o[k] for o in out_all)) for k in ("A", "A_cand", "proxy")}}


# ---------------------------------------------------------------------------------------------- pool
def norm_lab(s):
    return " ".join(sorted(t for t in re.findall(r"[a-z0-9]+", (s or "").lower()) if len(t) > 1))


def pool():
    sys.path.insert(0, HERE)
    import pool_count as PCNT
    import select_e9 as S
    pc = json.load(open(os.path.join(RES, "pool_count.json")))
    frame = list(csv.DictReader(open(os.path.join(ROOT, "data", "manifests", "t4_frame.tsv")), delimiter="\t"))
    # rebuild the same pool (same exclusions) deterministically
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        PCNT.main()
    excl = set()
    for f in sorted(glob.glob(os.path.join(ROOT, "data", "manifests", "*.tsv"))):
        if os.path.basename(f) == "t4_frame.tsv":
            continue
        try:
            for r in csv.DictReader(open(f), delimiter="\t"):
                for k in ("accession", "Accession", "project"):
                    if (r.get(k) or "").startswith("PXD"):
                        excl.add(r[k])
        except Exception:  # noqa: BLE001
            pass
    for f in glob.glob(os.path.join(ROOT, "results", "**", "t4_selector", "**", "t4_selection.tsv"), recursive=True):
        excl |= {r["accession"] for r in csv.DictReader(open(f), delimiter="\t")}
    excl |= {t[0] for t in json.load(open(os.path.join(ROOT, "data", "manifests", "e9_d8.log.json")))["tried"]}
    cand = [r for r in frame if r["eligible"] == "1" and r["organisms"] == "Homo sapiens (human)"
            and r["acquisition"] == "DDA" and r["label"] == "LFQ" and r["enrichment"] == "none"
            and r["vendor"] == "thermo" and S.INSTR.match(r["instruments"])
            and not S.BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])
            and r["accession"] not in excl]
    assert len(cand) == pc["pool_projects"], (len(cand), pc["pool_projects"])
    cand.sort(key=lambda r: r["accession"])
    usable = [v for v in pc["e9_usable_raw_per_project"]["values"]]
    return [(r["accession"], norm_lab(r["lab"]) or "?" + r["accession"]) for r in cand], usable, pc


# ---------------------------------------------------------------------------------------------- walk
def jeff(x):
    x = np.asarray(x, float)
    return float(x.sum() ** 2 / (x ** 2).sum()) if x.sum() > 0 else 0.0


def capped3(cp, ca, cA):
    if len(cp) * CAP_SHARE < 1:
        return None
    f = np.minimum.reduce([AP.cap_factors(cp, CAP_SHARE), AP.cap_factors(ca, CAP_SHARE),
                           AP.cap_factors(cA, CAP_SHARE)])
    return cp * f, ca * f, cA * f, f


def crit(cp, ca, cA, c):
    r = capped3(np.asarray(cp, float), np.asarray(ca, float), np.asarray(cA, float))
    if r is None:
        return False, None
    xp, xa, xA, f = r
    ok = (len(cp) >= c["J_min"] and xp.sum() >= c["proxy_min"] and xa.sum() >= c["A_cand_min"]
          and jeff(xp) >= c["Jeff_min"] and jeff(xa) >= c["Jeff_min"])
    return ok, r


def walk(rng, order_seed, projects, usable, y, R_cap, max_runs, yield_scale=1.0):
    order = list(projects)
    random.Random(order_seed).shuffle(order)
    cl_index, cl_count = {}, []
    cp, ca, cA = [], [], []
    st = {"runs": 0, "projects_walked": 0, "admitted": 0, "dl_files": 0, "met_min": False, "met_full": False,
          "snap_min": None, "snap_full": None, "end": "pool exhausted"}
    for acc, lab in order:
        st["projects_walked"] += 1
        k = cl_index.get(lab)
        if k is not None and cl_count[k] >= PER_CLUSTER:
            continue                                   # cluster full: skipped before listing
        nu = usable[rng.integers(len(usable))]
        if nu == 0:
            continue
        R = int(min(R_cap, nu))
        if st["runs"] + R > max_runs:
            st["end"] = "run cap"
            break
        st["dl_files"] += R
        if rng.random() < P_PROJECT_DEAD:
            continue
        ok_runs = int((rng.random(R) > P_DECODE_FAIL + y["zero_run_fraction"]).sum())
        if ok_runs == 0:
            continue
        mu = rng.normal(y["log_A_mean_project"], y["log_A_sd_between"])
        a = float(rng.lognormal(mu, y["log_A_sd_within"], size=ok_runs).sum()) * yield_scale
        crate = 1 / (1 + math.exp(-rng.normal(y["cand_logit_mean"], y["cand_logit_sd"])))
        acand = rng.binomial(int(a), crate)
        trate = math.exp(rng.normal(y["proxy_log_mean"], y["proxy_log_sd"]))
        proxy = rng.poisson(acand * min(trate, 1.0))
        if k is None:
            k = len(cl_count)
            cl_index[lab] = k
            cl_count.append(0)
            cp.append(0.0)
            ca.append(0.0)
            cA.append(0.0)
        cl_count[k] += 1
        cp[k] += proxy
        ca[k] += acand
        cA[k] += a
        st["admitted"] += 1
        st["runs"] += ok_runs
        if not st["met_min"]:
            ok, r = crit(cp, ca, cA, CRIT_MIN)
            if ok:
                st["met_min"] = True
                st["snap_min"] = {"runs": st["runs"], "J": len(cp), "dl_files": st["dl_files"],
                                  "cp_uncapped": list(cp), "cp_capped": list(r[0]), "ca_capped": list(r[1]),
                                  "ca_uncapped": list(ca), "admitted": st["admitted"]}
        ok, r = crit(cp, ca, cA, CRIT_FULL)
        if ok:
            st["met_full"] = True
            st["snap_full"] = {"runs": st["runs"], "J": len(cp), "dl_files": st["dl_files"],
                               "cp_uncapped": list(cp), "cp_capped": list(r[0]), "ca_capped": list(r[1]),
                               "ca_uncapped": list(ca), "admitted": st["admitted"]}
            st["end"] = "full stop"
            break
    st["J"] = len(cp)
    st["A_cand_uncapped"] = float(sum(ca))
    st["proxy_uncapped"] = float(sum(cp))
    return st


# ---------------------------------------------------------------------------------------------- tail-G2 (A4 IUT)
def tail_iut(snap, rng, kappa, scen, conf=AP.TAIL_CONF_A4, B=2000):
    """As reviews/checks/a4_r1_walk.tail_iut, plus residual contamination: the FIXED (hyperscore) tail target
    fraction is p_fixed (0.5 on a valid null), the learned per-cluster fractions from the scenario."""
    cu = np.asarray(snap["cp_uncapped"], float)
    cc = np.asarray(snap["cp_capped"], float)
    J = cu.size
    kind, par = scen
    p_fixed = 0.5
    if kind == "clean":
        pl = np.full(J, 0.5)
    elif kind == "tau":
        mu, tau = par
        s = mu * (1 - mu) / tau ** 2 - 1
        pl = rng.beta(mu * s, (1 - mu) * s, size=J)
    elif kind == "concentrated":
        frac, phi = par
        pl = np.full(J, 0.5)
        pl[rng.permutation(J)[:max(1, int(round(frac * J)))]] = phi
    elif kind == "contam":
        d = par
        p_fixed = 0.5 + d
        pl = np.full(J, 0.5 + 1.5 * d)
    nL = rng.poisson(kappa * cu)
    nF = rng.poisson(kappa * cu)
    kL = rng.binomial(nL, pl)
    kF = rng.binomial(nF, p_fixed)
    f = np.where(cu > 0, cc / np.maximum(cu, 1e-12), 1.0)
    nLc = rng.binomial(nL, f)
    kLc = rng.hypergeometric(np.maximum(kL, 0), np.maximum(nL - kL, 0), nLc)
    nFc = rng.binomial(nF, f)
    kFc = rng.hypergeometric(np.maximum(kF, 0), np.maximum(nF - kF, 0), nFc)

    def gate(nl, kl, nf, kf):
        cl = np.concatenate([np.full(nl[j] + nf[j], j) for j in range(J)])
        tl = np.concatenate([np.r_[np.ones(nl[j], bool), np.zeros(nf[j], bool)] for j in range(J)])
        gl = np.concatenate([np.r_[np.arange(nl[j]) < kl[j], np.zeros(nf[j], bool)] for j in range(J)])
        tf = np.concatenate([np.r_[np.zeros(nl[j], bool), np.ones(nf[j], bool)] for j in range(J)])
        gf = np.concatenate([np.r_[np.zeros(nl[j], bool), np.arange(nf[j]) < kf[j]] for j in range(J)])
        return A2S.tail_g2_a2(cl, tl, gl, tf, gf, B=B, seed=rng, conf=conf)["verdict"]
    vu, vc = gate(nL, kL, nF, kF), gate(nLc, kLc, nFc, kFc)
    for o in (A2S.INVALID, "fail", A2S.NOT_EVALUABLE):     # A4.3.6 IUT precedence
        if o in (vu, vc):
            return o, int(nLc.sum())
    return "pass", int(nLc.sum())


def quant(xs, qs=(10, 50, 90)):
    return [float(np.percentile(xs, q)) for q in qs] if len(xs) else None


# ---------------------------------------------------------------------------------------------- main
def main():
    t0 = time.perf_counter()
    projects, usable, pc = pool()
    out = {"seed": SEED, "spec": __doc__.split("    conda_envs")[0].strip(), "pool_projects": len(projects),
           "pool_labs": pc["pool_labs"], "criteria": {"v1_minimum": CRIT_MIN, "full_stop": CRIT_FULL,
                                                      "cap_share": CAP_SHARE, "per_cluster": PER_CLUSTER}}
    nwalk = 100 if QUICK else 400
    ntail = 60 if QUICK else 200
    Bt = 1000 if QUICK else 2000
    budgets = (60, 100, 150, 200, 300, 500, 1000, 2000)
    cells = []
    yields = {}
    snaps = {}
    for variant, thr in (("f0.20", "own"), ("f0.20", "D0max"), ("f0.05", "own"), ("f0.20_single_ms1", "own")):
        if True:
            y = pilot_yields(variant, thr)
            yields[f"{variant}/{thr}"] = {k: v for k, v in y.items() if k != "runs"}
            for R_cap in (2, 3):
                for yscale in (1.0, 0.5):
                    rng = np.random.default_rng(SEED)
                    ws = [walk(rng, 47 + i, projects, usable, y, R_cap, max(budgets), yscale) for i in range(nwalk)]
                    mins = [w["snap_min"] for w in ws if w["met_min"]]
                    fulls = [w["snap_full"] for w in ws if w["met_full"]]
                    c = {"variant": variant, "thr": thr, "R_cap": R_cap, "yield_scale": yscale,
                         "P_v1_minimum_by_run_budget": {b: float(np.mean([s is not None and s["runs"] <= b
                                                                          for s in [w["snap_min"] for w in ws]]))
                                                        for b in budgets},
                         "P_full_stop_by_run_budget": {b: float(np.mean([s is not None and s["runs"] <= b
                                                                         for s in [w["snap_full"] for w in ws]]))
                                                       for b in budgets},
                         "runs_at_v1_minimum_p10_50_90": quant([s["runs"] for s in mins]),
                         "files_downloaded_at_v1_minimum_p10_50_90": quant([s["dl_files"] for s in mins]),
                         "projects_at_v1_minimum_p10_50_90": quant([s["admitted"] for s in mins]),
                         "J_at_v1_minimum_p10_50_90": quant([s["J"] for s in mins]),
                         "capped_proxy_at_v1_minimum_p50": quant([sum(s["cp_capped"]) for s in mins])[1]
                         if mins else None,
                         "capped_A_cand_at_v1_minimum_p50": quant([sum(s["ca_capped"]) for s in mins])[1]
                         if mins else None,
                         "runs_at_full_stop_p10_50_90": quant([s["runs"] for s in fulls]),
                         "binding_at_v1_minimum": None}
                    # which criterion binds: first-met order of proxy vs A_cand vs J along the walk is not tracked;
                    # report the capped excess ratios at the v1 minimum instead
                    if mins:
                        c["binding_at_v1_minimum"] = {
                            "proxy_over_min_p50": float(np.median([sum(s["cp_capped"]) / 1200 for s in mins])),
                            "A_cand_over_min_p50": float(np.median([sum(s["ca_capped"]) / 120_000 for s in mins])),
                            "J_over_min_p50": float(np.median([s["J"] / 20 for s in mins]))}
                    cells.append(c)
                    print("WALK", json.dumps({k: c[k] for k in ("variant", "thr", "R_cap", "yield_scale",
                                                                "runs_at_v1_minimum_p10_50_90",
                                                                "binding_at_v1_minimum")}),
                          json.dumps(c["P_v1_minimum_by_run_budget"]), flush=True)
                    if thr == "own" and yscale == 1.0 and R_cap == 3:
                        snaps[variant] = mins
    out["pilot_yields"] = yields
    out["walk_cells"] = cells

    # ---------------- gates at the v1-minimum state (reference cell, and the strict single-MS1 variant)
    gates = {}
    scen = [("clean", ("clean", None)), ("margin0.60_tau0.10", ("tau", (0.60, 0.10))),
            ("margin0.60_tau0", ("tau", (0.60, 1e-3))),
            ("leak0.90_in_25pct_clusters", ("concentrated", (0.25, 0.90))),
            ("contam_d0.02", ("contam", 0.02)), ("contam_d0.05", ("contam", 0.05))]
    for lab_, mins in ((f"{v}/own/R3", snaps[v]) for v in ("f0.20", "f0.05", "f0.20_single_ms1") if v in snaps):
        if not mins:
            continue
        rng = np.random.default_rng(SEED + 50)
        g = {}
        for kappa in (1.0, 0.5):
            for name, sc in scen:
                v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
                nt = []
                for i in range(ntail):
                    vv, n = tail_iut(mins[i % len(mins)], rng, kappa, sc, B=Bt)
                    v[vv] += 1
                    nt.append(n)
                g[f"kappa{kappa}/{name}"] = {k: x / ntail for k, x in v.items()} | {"n_tail_capped_p50":
                                                                                    float(np.median(nt))}
                print("TAIL", lab_, kappa, name, json.dumps(g[f"kappa{kappa}/{name}"]), flush=True)
        # G1 at n_top = capped A_cand / 10 with the walk's J (A4 jackknife hull), clean; G3 per bin per label
        Js = int(np.median([s["J"] for s in mins]))
        ntop = int(np.median([sum(s["ca_capped"]) for s in mins]) / 10)
        g["G1_clean_pass"] = AP.g1_pass_prob(ntop, 0.0, 0.0, min(Js, 60), 60 if QUICK else 200,
                                             np.random.default_rng(SEED + 60)) * 0.99   # x hyperscore control
        g["G1_rho0.05_pass"] = AP.g1_pass_prob(ntop, 0.05, 0.0, min(Js, 60), 60 if QUICK else 200,
                                               np.random.default_rng(SEED + 61))
        n_bin = ntop // 2
        g3 = AP.tost_pass_prob(n_bin, AP.RECALL_G3, 0.0, 0.0, min(Js, 60), AP.TOST_CONF_A4, 100 if QUICK else 300,
                               np.random.default_rng(SEED + 62))
        g["G3_clean_pass_all_bins"] = g3 ** 10
        g["G1_n_top"], g["G3_n_per_label_per_bin"], g["J_median"] = ntop, n_bin, Js
        gates[lab_] = g
        print("G1/G3", lab_, json.dumps({k: g[k] for k in ("G1_clean_pass", "G1_rho0.05_pass",
                                                           "G3_clean_pass_all_bins", "G1_n_top", "J_median")}),
              flush=True)
    out["gates_at_v1_minimum"] = gates

    # ---------------- summary probabilities for the D28 decision
    summ = {}
    for lab_, g in gates.items():
        cellkey = lab_.split("/")
        c = next(x for x in cells if x["variant"] == cellkey[0] and x["thr"] == cellkey[1] and x["R_cap"] == 3
                 and x["yield_scale"] == 1.0)
        for kappa in (1.0, 0.5):
            for budget in (100, 200, 300):
                pmin = c["P_v1_minimum_by_run_budget"][budget]
                jc = g[f"kappa{kappa}/clean"]["pass"] * g["G1_clean_pass"] * g["G3_clean_pass_all_bins"]
                caught = g[f"kappa{kappa}/margin0.60_tau0.10"]["fail"]
                summ[f"{lab_}/kappa{kappa}/budget{budget}"] = {
                    "P_v1_minimum": pmin, "joint_clean_pass_given_evaluable": jc,
                    "P_clean_v1_passes": pmin * jc, "P_margin_leak_caught_given_evaluable": caught,
                    "P_margin_leak_caught": pmin * caught,
                    "P_false_pass_margin_given_evaluable": g[f"kappa{kappa}/margin0.60_tau0.10"]["pass"]}
    out["summary"] = summ
    out["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print("wrote", OUT, round(out["wall_s"]), "s")


if __name__ == "__main__":
    main()
