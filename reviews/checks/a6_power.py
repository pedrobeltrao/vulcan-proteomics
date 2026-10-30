#!/usr/bin/env python3
"""Amendment A6 drafting check: c9_power re-run with the FINAL A6 rules (pilot condition 6).

SYNTHETIC SIMULATION + LABEL-FREE YIELDS + METADATA ONLY.  Seed 131.  Reuses scripts/c9/c9_power.py (tail_iut, quant)
and reviews/checks/a4_power.py / a4_r3_walk.py (cap factors, G1/TOST power, person/lab keys) unchanged.

What changes against scripts/c9/c9_power.py (the pilot):
  * YIELDS: reviews/checks/a6_cert_yield.json, i.e. the final admission rules measured on the pilot's 12 DEV runs:
    two-engine certification (hyperscore q<=1e-3 AND Sage rank-1 agreement at spectrum_q<=1e-3), ladder screen
    f_max = 0.05, contaminant pairs dropped.  'cert_hs' (one engine) is kept as a sensitivity cell.
  * POOL / EXCLUSIONS (A6.3.2): the pilot's exclusions PLUS every PXD string in any tracked file under data/manifests/,
    results/ or config/ (except t4_frame.tsv), every accession with a cached PRIDE file listing, and A4.3.3 item 5
    (any person/lab key shared with an A1 null project, a DNIT walk candidate or a DN2 walk candidate).
  * CLUSTERS (A6.3.3 = A4.3.4): connected components of ADMITTED projects linked by any shared order-insensitive
    person/lab key (crawl lab + PRIDE labPIs + submitters); a candidate whose component would exceed 5 projects is
    skipped before listing ("cluster full").  The pilot used the crawl lab string only.
  * WALK ORDER seeds 9001.. (the pilot's simulation used 47+i; A6's real walk seed is derived from the tag commit).
  * v1 is tested on the v1-MINIMUM PREFIX (A6.6.1); v2_sep / ce_s on the walk end (full stop or hard cap).
  * extra residual-contamination cell d = 0.01, and tail-G2 of the hyperscore-referenced v2_sep row V1 at the
    full-stop state (a proxy: v2_sep's tail is modelled as the hyperscore proxy times kappa).
    conda_envs/pe/bin/python reviews/checks/a6_power.py   # -> reviews/checks/a6_power.json (+ .out)
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import random
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "c9"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "e9"))
sys.path.insert(0, ROOT)
import a4_power as AP  # noqa: E402
import a4_r3_walk as R3  # noqa: E402
import c9_power as CP  # noqa: E402
from vulcan.gates import a2_stats as A2S  # noqa: E402

SEED = 131
QUICK = os.environ.get("A6_QUICK") == "1"
OUT = os.path.join(HERE, "a6_power.json")
CAP_SHARE, PER_CLUSTER, R_CAP, HARD_CAP = 0.10, 5, 3, 2000
CRIT_MIN = {"J_min": 20, "Jeff_min": 10, "proxy_min": 1200, "A_cand_min": 120_000}
CRIT_FULL = dict(CRIT_MIN, A_cand_min=450_000)


# ------------------------------------------------------------------------------------------- yields (final rules)
def yields(variant):
    d = json.load(open(os.path.join(HERE, "a6_cert_yield.json")))
    out = [{"project": r["project"], "run": r["run"], "A": r[variant]["n_certified"],
            "A_cand": r[variant]["A_cand"], "proxy": r[variant]["tail_proxy_own_thr"]} for r in d["per_run"].values()]
    z0 = sum(o["A"] == 0 for o in out) / len(out)
    pos = [o for o in out if o["A"] > 0]
    proj = sorted({o["project"] for o in pos})
    la = np.log([o["A"] for o in pos])
    pj = np.array([o["project"] for o in pos])
    pm = np.array([la[pj == p].mean() for p in proj])
    wi = np.concatenate([la[pj == p] - la[pj == p].mean() for p in proj])
    dfw = len(la) - len(proj)
    s_w = max(float(np.sqrt((wi ** 2).sum() / dfw)) if dfw > 0 else 0.32, 0.32)
    agg = {p: {k: sum(o[k] for o in pos if o["project"] == p) for k in ("A", "A_cand", "proxy")} for p in proj}
    cr = np.array([agg[p]["A_cand"] / agg[p]["A"] for p in proj])
    tr = np.array([(agg[p]["proxy"] + 0.5) / max(agg[p]["A_cand"], 1) for p in proj])
    return {"variant": variant, "n_runs": len(out), "zero_run_fraction": z0,
            "log_A_mean_project": float(pm.mean()), "log_A_sd_between": float(pm.std(ddof=1)), "log_A_sd_within": s_w,
            "cand_logit_mean": float(np.mean(np.log(cr / (1 - cr)))),
            "cand_logit_sd": float(np.std(np.log(cr / (1 - cr)), ddof=1)),
            "proxy_log_mean": float(np.mean(np.log(tr))), "proxy_log_sd": float(np.std(np.log(tr), ddof=1)),
            "totals": {k: int(sum(o[k] for o in out)) for k in ("A", "A_cand", "proxy")}}


# ------------------------------------------------------------------------------------------- pool with A6 exclusions
def tracked_pxd():
    files = subprocess.run(["git", "ls-files", "data/manifests", "results", "config"], cwd=ROOT, capture_output=True,
                           text=True, check=True).stdout.split()
    acc = set()
    for f in files:
        if os.path.basename(f) == "t4_frame.tsv":
            continue
        try:
            acc |= set(re.findall(r"PXD\d{6}", open(os.path.join(ROOT, f), errors="ignore").read()))
        except (OSError, UnicodeDecodeError):
            pass
    return acc


def pool():
    import io
    import contextlib
    import pool_count as PCNT
    PCNT.OUT = os.devnull                     # CP.pool() re-runs pool_count.main(); never rewrite the pilot's record
    with contextlib.redirect_stdout(io.StringIO()):
        projects_c9, usable, pc = CP.pool()                     # the pilot's pool (asserts its count)
    c9 = {a for a, _ in projects_c9}
    excl = {"tracked_pxd": tracked_pxd(),
            "cached_file_listing": {re.match(r"files_(PXD\d{6})", os.path.basename(f)).group(1)
                                    for f in glob.glob(os.path.join(ROOT, "data/cache/pride_meta/files_PXD*"))}}
    crawl = pd.read_csv(os.path.join(ROOT, "data/cache/pride_meta/all_projects_classified.tsv"), sep="\t",
                        dtype=str).fillna("").set_index("accession")
    prior = R3.prior_null_accessions()
    kp = R3.keysets(sorted(prior), crawl)
    prior_keys = set().union(*kp.values()) if kp else set()
    cand = sorted(c9)
    kc = R3.keysets(cand, crawl)
    excl["same_lab_item5"] = {a for a in cand if kc[a] & prior_keys}
    allx = set().union(*excl.values())
    keep = [a for a in cand if a not in allx]
    return keep, {a: kc[a] for a in keep}, usable, {"pilot_pool": len(c9), "a6_pool": len(keep),
                                                    "removed_by": {k: len(c9 & v) for k, v in excl.items()},
                                                    "removed_total": len(c9) - len(keep)}


# ------------------------------------------------------------------------------------------- walk (A4.3.4 clusters)
def walk(rng, order_seed, accs, keys, usable, y, max_runs=HARD_CAP):
    order = sorted(accs)
    random.Random(order_seed).shuffle(order)
    par = []
    tot = {}                                  # root -> [proxy, A_cand, A, n_projects]

    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i
    key_owner = {}
    S = {"runs": 0, "dl_files": 0, "walked": 0, "admitted": 0, "skip_full": 0, "snap_min": None, "snap_full": None,
         "end": "pool exhausted"}

    def arrays():
        v = np.array(list(tot.values()), float)
        return v[:, 0], v[:, 1], v[:, 2]

    def snap(cp, ca, cA, r):
        return {"runs": S["runs"], "J": len(cp), "dl_files": S["dl_files"], "admitted": S["admitted"],
                "cp_uncapped": list(cp), "cp_capped": list(r[0]), "ca_capped": list(r[1]), "ca_uncapped": list(ca)}

    for acc in order:
        S["walked"] += 1
        ks = keys[acc] or {"?" + acc}
        linked = {find(key_owner[k]) for k in ks if k in key_owner}
        if linked and sum(tot[r][3] for r in linked) + 1 > PER_CLUSTER:
            S["skip_full"] += 1
            continue
        nu = usable[rng.integers(len(usable))]
        if nu == 0:
            continue
        R = int(min(R_CAP, nu))
        if S["runs"] + R > max_runs:
            S["end"] = "hard cap"
            break
        S["dl_files"] += R
        if rng.random() < CP.P_PROJECT_DEAD:
            continue
        ok_runs = int((rng.random(R) > CP.P_DECODE_FAIL + y["zero_run_fraction"]).sum())
        if ok_runs == 0:
            continue
        mu = rng.normal(y["log_A_mean_project"], y["log_A_sd_between"])
        a = float(rng.lognormal(mu, y["log_A_sd_within"], size=ok_runs).sum())
        crate = 1 / (1 + math.exp(-rng.normal(y["cand_logit_mean"], y["cand_logit_sd"])))
        ac = float(rng.binomial(int(a), crate))
        trate = math.exp(rng.normal(y["proxy_log_mean"], y["proxy_log_sd"]))
        pr = float(rng.poisson(ac * min(trate, 1.0)))
        i = len(par)
        par.append(i)
        t = [pr, ac, a, 1]
        for r in linked:
            par[r] = i
            o = tot.pop(r)
            t = [t[j] + o[j] for j in range(4)]
        tot[i] = t
        for k in ks:
            key_owner[k] = i
        S["admitted"] += 1
        S["runs"] += ok_runs
        cp, ca, cA = arrays()
        if S["snap_min"] is None:
            ok, r = CP.crit(cp, ca, cA, CRIT_MIN)
            if ok:
                S["snap_min"] = snap(cp, ca, cA, r)
        ok, r = CP.crit(cp, ca, cA, CRIT_FULL)
        if ok:
            S["snap_full"] = snap(cp, ca, cA, r)
            S["end"] = "full stop"
            break
    S["J"] = len(tot)
    return S


# ------------------------------------------------------------------------------------------- main
def main():
    t0 = time.perf_counter()
    accs, keys, usable, pinfo = pool()
    print("POOL", json.dumps(pinfo), flush=True)
    out = {"seed": SEED, "spec": __doc__.split("    conda_envs")[0].strip(), "pool": pinfo,
           "criteria": {"v1_minimum": CRIT_MIN, "full_stop": CRIT_FULL, "cap_share": CAP_SHARE,
                        "per_cluster": PER_CLUSTER, "R": R_CAP, "hard_cap_runs": HARD_CAP}}
    nwalk = 100 if QUICK else 400
    ntail = 60 if QUICK else 200
    Bt = 1000 if QUICK else 2000
    budgets = (100, 150, 200, 300, 500, 1000, 2000)
    out["yields"], out["walk_cells"], snaps = {}, [], {}
    for variant in ("cert_agree", "cert_hs"):
        y = yields(variant)
        out["yields"][variant] = y
        for yscale in (1.0, 0.5):
            yy = dict(y, log_A_mean_project=y["log_A_mean_project"] + math.log(yscale))
            rng = np.random.default_rng(SEED)
            ws = [walk(rng, 9001 + i, accs, keys, usable, yy) for i in range(nwalk)]
            mins = [w["snap_min"] for w in ws if w["snap_min"]]
            fulls = [w["snap_full"] for w in ws if w["snap_full"]]
            c = {"variant": variant, "yield_scale": yscale,
                 "P_v1_minimum_by_run_budget": {b: float(np.mean([w["snap_min"] is not None and
                                                                   w["snap_min"]["runs"] <= b for w in ws]))
                                                for b in budgets},
                 "P_full_stop_by_run_budget": {b: float(np.mean([w["snap_full"] is not None and
                                                                  w["snap_full"]["runs"] <= b for w in ws]))
                                               for b in budgets},
                 "runs_at_v1_minimum_p10_50_90": CP.quant([s["runs"] for s in mins]),
                 "files_downloaded_at_v1_minimum_p10_50_90": CP.quant([s["dl_files"] for s in mins]),
                 "projects_at_v1_minimum_p10_50_90": CP.quant([s["admitted"] for s in mins]),
                 "J_at_v1_minimum_p10_50_90": CP.quant([s["J"] for s in mins]),
                 "runs_at_full_stop_p10_50_90": CP.quant([s["runs"] for s in fulls]),
                 "files_downloaded_at_full_stop_p10_50_90": CP.quant([s["dl_files"] for s in fulls]),
                 "J_at_full_stop_p10_50_90": CP.quant([s["J"] for s in fulls]),
                 "projects_walked_p50": float(np.median([w["walked"] for w in ws])),
                 "skipped_cluster_full_p50": float(np.median([w["skip_full"] for w in ws])),
                 "binding_at_v1_minimum": {
                     "proxy_over_min_p50": float(np.median([sum(s["cp_capped"]) / 1200 for s in mins])),
                     "A_cand_over_min_p50": float(np.median([sum(s["ca_capped"]) / 120_000 for s in mins])),
                     "J_over_min_p50": float(np.median([s["J"] / 20 for s in mins]))} if mins else None}
            out["walk_cells"].append(c)
            print("WALK", json.dumps({k: c[k] for k in ("variant", "yield_scale", "runs_at_v1_minimum_p10_50_90",
                                                        "runs_at_full_stop_p10_50_90", "J_at_v1_minimum_p10_50_90",
                                                        "binding_at_v1_minimum")}),
                  json.dumps(c["P_v1_minimum_by_run_budget"]), json.dumps(c["P_full_stop_by_run_budget"]), flush=True)
            snaps[(variant, yscale)] = (mins, fulls)

    scen = [("clean", ("clean", None)), ("margin0.60_tau0.10", ("tau", (0.60, 0.10))),
            ("margin0.60_tau0", ("tau", (0.60, 1e-3))), ("leak0.90_in_25pct_clusters", ("concentrated", (0.25, 0.90))),
            ("contam_d0.01", ("contam", 0.01)), ("contam_d0.02", ("contam", 0.02)), ("contam_d0.05", ("contam", 0.05))]
    gates = {}
    for (variant, yscale), (mins, fulls) in snaps.items():
        if variant != "cert_agree":
            continue
        for state, ss in (("v1_prefix", mins),):
            if not ss or (state == "full_stop" and yscale != 1.0):
                continue
            lab_ = f"{variant}/x{yscale}/{state}"
            rng = np.random.default_rng(SEED + 50)
            g = {}
            for kappa in (1.0, 0.5):
                for name, sc in scen:
                    v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
                    nt = []
                    for i in range(ntail):
                        vv, n = CP.tail_iut(ss[i % len(ss)], rng, kappa, sc, B=Bt)
                        v[vv] += 1
                        nt.append(n)
                    g[f"kappa{kappa}/{name}"] = {k: x / ntail for k, x in v.items()} | {
                        "n_tail_capped_p50": float(np.median(nt))}
                    print("TAIL", lab_, kappa, name, json.dumps(g[f"kappa{kappa}/{name}"]), flush=True)
            Js = int(np.median([s["J"] for s in ss]))
            ntop = int(np.median([sum(s["ca_capped"]) for s in ss]) / 10)
            g["G1_clean_pass"] = AP.g1_pass_prob(ntop, 0.0, 0.0, min(Js, 60), 60 if QUICK else 200,
                                                 np.random.default_rng(SEED + 60)) * 0.99
            g["G1_rho0.05_pass"] = AP.g1_pass_prob(ntop, 0.05, 0.0, min(Js, 60), 60 if QUICK else 200,
                                                   np.random.default_rng(SEED + 61))
            g3 = AP.tost_pass_prob(ntop // 2, AP.RECALL_G3, 0.0, 0.0, min(Js, 60), AP.TOST_CONF_A4,
                                   100 if QUICK else 300, np.random.default_rng(SEED + 62))
            g["G3_clean_pass_all_bins"] = g3 ** 10
            g["G1_n_top"], g["G3_n_per_label_per_bin"], g["J_median"] = ntop, ntop // 2, Js
            gates[lab_] = g
            print("G1/G3", lab_, json.dumps({k: g[k] for k in ("G1_clean_pass", "G1_rho0.05_pass",
                                                               "G3_clean_pass_all_bins", "G1_n_top", "J_median")}),
                  flush=True)
    out["gates"] = gates

    summ = {}
    for yscale in (1.0, 0.5):
        c = next(x for x in out["walk_cells"] if x["variant"] == "cert_agree" and x["yield_scale"] == yscale)
        g = gates.get(f"cert_agree/x{yscale}/v1_prefix")
        if g is None:
            continue
        for kappa in (1.0, 0.5):
            for budget in (200, 300, 500, 2000):
                pmin = c["P_v1_minimum_by_run_budget"][budget]
                jc = g[f"kappa{kappa}/clean"]["pass"] * g["G1_clean_pass"] * g["G3_clean_pass_all_bins"]
                summ[f"x{yscale}/kappa{kappa}/budget{budget}"] = {
                    "P_v1_minimum": pmin, "joint_clean_pass_given_evaluable": jc, "P_clean_v1_passes": pmin * jc,
                    "P_margin_leak_caught": pmin * g[f"kappa{kappa}/margin0.60_tau0.10"]["fail"],
                    "P_false_pass_margin_given_evaluable": g[f"kappa{kappa}/margin0.60_tau0.10"]["pass"],
                    "P_false_pass_concentrated_given_evaluable": g[f"kappa{kappa}/leak0.90_in_25pct_clusters"]["pass"],
                    "P_invalid_clean_given_evaluable": g[f"kappa{kappa}/clean"][A2S.INVALID],
                    "P_pass_contam_d0.01": g[f"kappa{kappa}/contam_d0.01"]["pass"],
                    "P_pass_contam_d0.02": g[f"kappa{kappa}/contam_d0.02"]["pass"]}
    out["summary"] = summ
    out["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print("wrote", OUT, round(out["wall_s"]), "s")


if __name__ == "__main__":
    main()
