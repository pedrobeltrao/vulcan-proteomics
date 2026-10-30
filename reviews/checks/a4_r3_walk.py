#!/usr/bin/env python3
"""Amendment A4 revision 3: walk and tail-G2 check implementing the FULL A4 selection spec (second review B1, M4,
M5, minor 1).  METADATA + SYNTHETIC ONLY, seed 131.  Supersedes the walk numbers of a4_r1_walk.py (which is kept for
the record).

Differences from a4_r1_walk.py:
  * person/lab keys (A4.3.4): normalised crawl `lab` plus every `labPIs` and `submitters` name from the cached PRIDE
    project record; normalisation is ORDER-INSENSITIVE (lower-case, alphanumeric tokens, single-letter tokens dropped,
    tokens sorted), so "zhang yongqian" == "Yongqian Zhang";
  * cluster links: shared genus OR any shared key;
  * A4.3.3 item 5 (same-lab exclusion) applied to the pool AND to tier 0: a project is excluded if any key is shared
    with an A1 null project (g2null.tsv + g2null "tried"), a DNIT walk candidate, or a DN2 walk candidate not in
    tier 0;
  * the cap takes the minimum of the water-filling factors of the tail proxy, A_cand AND |A| (A4.3.6);
  * tier 0 is shuffled correctly (random.Random(order_seed).shuffle) and, because 3 of its 5 protocol-eligible
    projects already showed one non-keratin panel member at 3 runs (DN2 NOTES), its stage-2 extra exclusion is 0.5;
  * a REFERENCE cell draws each project's raw-file count from the empirical distribution of DN2's 21 listing events
    (label-free; usable ~ Binomial(n_raw, 0.8)) instead of an assumed median;
  * IUT precedence (A4.3.6, second review minor 1): any invalid pool -> invalid; else any fail -> fail; else any
    not-evaluable -> not-evaluable; else pass.
    conda_envs/pe/bin/python reviews/checks/a4_r3_walk.py   # -> reviews/checks/a4_r3_walk.json (+ .out via tee)
"""
from __future__ import annotations

import glob
import gzip
import json
import math
import os
import random
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import a4_power as AP  # noqa: E402
import a4_r1_walk as W1  # noqa: E402
from vulcan.gates import a2_stats as A2S  # noqa: E402

SEED = 131
QUICK = os.environ.get("A4_QUICK") == "1"
OUT = os.path.join(HERE, "a4_r3_walk.json")
DN2_ADMITTED = ["PXD006818", "PXD018953", "PXD044078", "PXD011591", "PXD062500", "PXD038243"]


def norm_person(s: str) -> str:
    toks = [t for t in re.split(r"[^a-z0-9]+", str(s).lower()) if len(t) > 1]
    return " ".join(sorted(toks))


def load_records(accs):
    want, recs = set(accs), {}
    for f in sorted(glob.glob(os.path.join(ROOT, "data/cache/pride_meta/projects_p*.json.gz"))):
        for r in json.load(gzip.open(f)):
            if r.get("accession") in want:
                recs[r["accession"]] = r
    return recs


def person_names(r):
    out = []
    for fld in ("labPIs", "submitters"):
        for p in (r.get(fld) or []):
            if isinstance(p, dict):
                nm = p.get("name") or " ".join(str(p.get(x, "")) for x in ("firstName", "lastName"))
            else:
                nm = str(p)
            out.append(nm)
    return out


def keysets(accs, crawl):
    recs = load_records(accs)
    ks = {}
    for a in accs:
        k = set()
        if a in crawl.index and crawl.loc[a, "lab"]:
            k.add(norm_person(crawl.loc[a, "lab"]))
        if a in recs:
            k |= {norm_person(n) for n in person_names(recs[a])}
        ks[a] = {x for x in k if x}
    return ks


def prior_null_accessions():
    acc = set(pd.read_csv(os.path.join(ROOT, "data/manifests/g2null.tsv"), sep="\t", dtype=str).accession.dropna())
    acc |= {x[0] for x in json.load(open(os.path.join(ROOT, "data/manifests/g2null_select.log.json")))["tried"]}
    acc |= {x["accession"] for x in json.load(open(os.path.join(ROOT, "data/manifests/dnit_select.log.json")))["log"]}
    dn2 = {x["accession"] for x in json.load(open(os.path.join(ROOT, "data/manifests/dn2_select.log.json")))["log"]}
    acc |= dn2 - set(DN2_ADMITTED)
    return {a for a in acc if str(a).startswith("PXD")}


def capped3(cp, ca, cA, share):
    if len(cp) * share < 1:
        return None
    f = np.minimum.reduce([AP.cap_factors(cp, share), AP.cap_factors(ca, share), AP.cap_factors(cA, share)])
    return cp * f, ca * f, cA * f, f


def criteria3(cp, ca, cA, c):
    r = capped3(np.asarray(cp, float), np.asarray(ca, float), np.asarray(cA, float), c["cap_share"])
    if r is None:
        return False, None
    xp, xa, xA, f = r
    ok = (len(cp) >= c["J_min"] and xp.sum() >= c["proxy_min"] and xa.sum() >= c["A_cand_min"]
          and W1.jeff(xp) >= c["Jeff_min"] and W1.jeff(xa) >= c["Jeff_min"])
    return ok, r


def walk3(rng, y, tiers, flags, keys, *, order_seed, med_runs, run_cap, p_proteome, p_files, p_marker1, e_stage2,
          crit_full, crit_min, max_runs=2000, max_per_cluster=5, tier0=None, u_tier0=0.8, p_marker1_tier0=0.9,
          e_stage2_tier0=0.5, raw_empirical=None, u_emp=0.8):
    mA = math.log(y["per_run_A"]["geo_mean_positive"])
    sb, sw = y["per_run_A"]["log_sd_between_projects"], y["per_run_A"]["log_sd_within_project"]
    z0 = y["per_run_A"]["zero_run_fraction"]
    lt = math.log(y["tail_proxy_per_A"]["pooled"])
    st = y["tail_proxy_per_A"]["log_sd_between_projects"]
    lc = math.log(y["A_cand_per_A"]["pooled"] / (1 - y["A_cand_per_A"]["pooled"]))
    sc = y["A_cand_per_A"]["logit_sd_between_projects"]
    adm, proxy, acand, aa = [], [], [], []
    S = {"runs": 0, "dl": 0, "skip_full": 0, "walked": 0, "met_min": False, "snap": None, "full": False,
         "tier0_admitted": 0, "tier0_runs": 0}

    def comps():
        n = len(adm)
        par = list(range(n))

        def f(i):
            while par[i] != i:
                par[i] = par[par[i]]
                i = par[i]
            return i
        by = {}
        for i, (g, ks) in enumerate(adm):
            for key in [("g", g)] + [("k", k) for k in ks]:
                if key in by:
                    par[f(i)] = f(by[key])
                else:
                    by[key] = i
        return [f(i) for i in range(n)]

    def sums():
        c = comps()
        u = {r: j for j, r in enumerate(dict.fromkeys(c))}
        idx = np.array([u[r] for r in c])
        J = len(u)
        return (np.bincount(idx, weights=proxy, minlength=J), np.bincount(idx, weights=acand, minlength=J),
                np.bincount(idx, weights=aa, minlength=J))

    def full_slot(g, ks):
        linked = [i for i, (gg, kk) in enumerate(adm) if gg == g or (ks & kk)]
        if not linked:
            return False
        c = comps()
        roots = {c[i] for i in linked}
        return sum(1 for x in c if x in roots) + 1 > max_per_cluster

    def admit(g, ks, R, a, pr, cr):
        adm.append((g, ks))
        proxy.append(float(rng.poisson(a * pr)))
        acand.append(float(rng.binomial(int(a), cr)))
        aa.append(a)
        S["runs"] += R
        cp, ca, cA = sums()
        if not S["met_min"]:
            ok, cs = criteria3(cp, ca, cA, crit_min)
            if ok:
                S["met_min"] = True
                S["snap"] = {"runs": S["runs"], "J": len(cp), "cp_capped": cs[0], "cp_uncapped": cp}
        ok, _ = criteria3(cp, ca, cA, crit_full)
        if ok:
            S["full"] = True
            return True
        return False

    def draw_R():
        if raw_empirical is not None:
            n = int(rng.choice(raw_empirical))
            return int(min(run_cap, max(1, rng.binomial(n, u_emp))))
        return int(min(run_cap, max(1, round(rng.lognormal(math.log(med_runs), 0.8)))))

    t0 = list(tier0 or [])
    random.Random(order_seed).shuffle(t0)
    for p0 in t0:
        S["walked"] += 1
        g, ks = p0["genus"], p0["keys"]
        if full_slot(g, ks) or not p0["eligible"]:
            continue
        R = int(min(run_cap, max(1, rng.binomial(p0["n_raw"], u_tier0))))
        if rng.random() > p_marker1_tier0:
            S["dl"] += min(R, 3)
            continue
        S["dl"] += R
        if R > 3 and rng.random() < e_stage2_tier0:
            continue
        a = float((rng.lognormal(math.log(max(p0["A_per_run"], 1.0)), sw, size=R) * (rng.random(R) > z0)).sum())
        S["tier0_admitted"] += 1
        S["tier0_runs"] += R
        if admit(g, ks, R, a, p0["proxy_rate"], p0["cand_rate"]):
            break
    if not S["full"]:
        for t in ("I", "II", "III"):
            cand = list(tiers[t][["accession", "genus"]].itertuples(index=False))
            random.Random(order_seed).shuffle(cand)
            stop = False
            for acc, g in cand:
                S["walked"] += 1
                ks = keys[acc]
                if full_slot(g, ks):
                    S["skip_full"] += 1
                    continue
                if not flags[acc]["rev1"] or rng.random() > p_proteome or rng.random() > p_files:
                    continue
                R = draw_R()
                if rng.random() > p_marker1:
                    S["dl"] += min(R, 3)
                    continue
                S["dl"] += R
                if R > 3 and rng.random() < e_stage2:
                    continue
                if S["runs"] + R > max_runs:
                    stop = True
                    break
                mu = rng.normal(mA, sb)
                a = float((rng.lognormal(mu, sw, size=R) * (rng.random(R) > z0)).sum())
                tr = math.exp(rng.normal(lt, st))
                cr = 1 / (1 + math.exp(-rng.normal(lc, sc)))
                if admit(g, ks, R, a, tr, cr):
                    stop = True
                    break
            if stop:
                break
    S["J"] = len(set(comps())) if adm else 0
    S["admitted"] = len(adm)
    return S


def tail_iut3(snap, rng, kappa, scenario):
    v, vc, vu, n = W1.tail_iut(snap, rng, kappa, scenario)
    order = [A2S.INVALID, "fail", A2S.NOT_EVALUABLE]
    for o in order:
        if o in (vc, vu):
            return o, vc, vu, n
    return "pass", vc, vu, n


def main():
    t_start = time.perf_counter()
    out = {"seed": SEED, "spec": __doc__.split("    conda_envs")[0].strip()}
    k, exinfo = W1.pool_r1()
    flags = W1.protocol_flags(list(k.accession))
    crawl = pd.read_csv(os.path.join(ROOT, "data/cache/pride_meta/all_projects_classified.tsv"), sep="\t",
                        dtype=str).fillna("").set_index("accession")
    prior = prior_null_accessions()
    keys = keysets(list(k.accession) + DN2_ADMITTED + sorted(prior), crawl)
    prior_keys = set().union(*[keys.get(a, set()) for a in prior])
    k["item5"] = [bool(keys[a] & prior_keys) for a in k.accession]
    k["rev1"] = [flags[a]["rev1"] for a in k.accession]
    item5 = {t: {"removed": int(k[(k.tier == t) & k.item5].shape[0]),
                 "removed_rev1_passers": int(k[(k.tier == t) & k.item5 & k.rev1].shape[0])} for t in ("I", "II", "III")}
    k = k[~k.item5]
    tiers = {t: k[k.tier == t].sort_values("accession") for t in ("I", "II", "III")}
    y = AP.dn2_yields()
    t0 = AP.dn2_tier0_projects()
    f0 = W1.protocol_flags([p["accession"] for p in t0])
    for p in t0:
        p["genus"] = crawl.loc[p["accession"], "organisms"].split(" ")[0].capitalize()
        p["keys"] = keys[p["accession"]]
        p["protocol_rev1"] = bool(f0[p["accession"]]["rev1"])
        p["item5_hit"] = bool(p["keys"] & prior_keys)
        p["eligible"] = p["protocol_rev1"] and not p["item5_hit"]
    out["pool_r3"] = {"exclusions": exinfo, "item5_removed": item5,
                      "tiers": {t: {"projects": int(len(v)), "rev1_passers": int(v.rev1.sum())} for t, v in tiers.items()}}
    out["tier0"] = [{x: (sorted(v) if isinstance(v, set) else v) for x, v in p.items()} for p in t0]
    print("POOL", json.dumps(out["pool_r3"]), flush=True)
    print("TIER0", json.dumps([(p["accession"], p["protocol_rev1"], p["item5_hit"], p["eligible"]) for p in t0]),
          flush=True)
    emp = []
    ac = set(json.load(open(os.path.join(ROOT, "results/2026-10-27_A2_dn2/admission_counts.json")))["admitted"])
    seen = set()
    for l in open(os.path.join(ROOT, "results/2026-10-27_A2_dn2/walk_events.jsonl")):
        e = json.loads(l)
        if e["ev"] == "listed" and e["accession"] not in seen:
            seen.add(e["accession"])
            emp.append(int(e["n_raw_records"]))
    out["empirical_raw_counts_dn2_listings"] = sorted(emp)

    D = {"cap_share": 0.10, "Jeff_min": 10, "proxy_min": 1200, "max_per_cluster": 5, "run_cap": 48}
    crit = lambda a: {"J_min": 20, "proxy_min": D["proxy_min"], "A_cand_min": a, "cap_share": D["cap_share"],
                      "Jeff_min": D["Jeff_min"]}
    nw = 100 if QUICK else 300
    grid = []
    cells = [("empirical DN2 raw counts (REFERENCE)", None), ("median 6", 6), ("median 12", 12), ("median 24", 24)]
    for with_t0 in (True, False):
        for lab, med in cells:
            for e2 in (0.15, 0.35):
                r = np.random.default_rng(SEED + 70)
                ws = [walk3(r, y, tiers, flags, keys, order_seed=1000 + i, med_runs=med or 12, run_cap=D["run_cap"],
                            p_proteome=16 / 19, p_files=0.93, p_marker1=0.5, e_stage2=e2,
                            crit_full=crit(450_000), crit_min=crit(120_000), max_per_cluster=D["max_per_cluster"],
                            tier0=t0 if with_t0 else None, raw_empirical=emp if med is None else None)
                      for i in range(nw)]
                q = lambda key: [float(np.percentile([w[key] for w in ws], p)) for p in (10, 50, 90)]
                cell = {"tier0": with_t0, "files": lab, "stage2_extra_exclusion": e2,
                        "P_v1_min": float(np.mean([w["met_min"] for w in ws])),
                        "P_full": float(np.mean([w["full"] for w in ws])),
                        "runs": q("runs"), "J": q("J"), "admitted": q("admitted"), "downloaded": q("dl"),
                        "tier0_admitted": q("tier0_admitted")}
                grid.append(cell)
                print("WALK", json.dumps(cell), flush=True)
    out["walk_r3"] = grid

    tg = {}
    nrep = 120 if QUICK else 300
    for kappa in (1.0, 0.5):
        r = np.random.default_rng(SEED + 80)
        snaps, i = [], 0
        while len(snaps) < nrep and i < 30 * nrep:
            w = walk3(r, y, tiers, flags, keys, order_seed=5000 + i, med_runs=12, run_cap=D["run_cap"],
                      p_proteome=16 / 19, p_files=0.93, p_marker1=0.5, e_stage2=0.35, crit_full=crit(450_000),
                      crit_min=crit(120_000), max_per_cluster=D["max_per_cluster"], tier0=t0, raw_empirical=emp)
            i += 1
            if w["met_min"]:
                snaps.append(w["snap"])
        cell = {"n_snapshots": len(snaps), "J_median": float(np.median([s["J"] for s in snaps])),
                "Jeff_capped_proxy_min": float(np.min([W1.jeff(s["cp_capped"]) for s in snaps])),
                "max_capped_share_max": float(np.max([s["cp_capped"].max() / s["cp_capped"].sum() for s in snaps]))}
        for lab, scen in (("clean", ("clean", None)), ("margin0.60_tau0.10", ("tau", (0.60, 0.10))),
                          ("mean0.60_leak0.90_in_25pct", ("concentrated", (0.25, 0.90))),
                          ("mean0.60_leak1.00_in_20pct", ("concentrated", (0.20, 1.00))),
                          ("uncapped_mean0.60_leak_in_largest", ("largest", None))):
            v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
            vc_pass = vu_pass = 0
            for s in snaps:
                verdict, vc, vu, _ = tail_iut3(s, r, kappa, scen)
                v[verdict] += 1
                vc_pass += vc == "pass"
                vu_pass += vu == "pass"
            cell[lab] = {kk: x / len(snaps) for kk, x in v.items()} | {"capped_only_pass": vc_pass / len(snaps),
                                                                      "uncapped_only_pass": vu_pass / len(snaps)}
        tg[f"kappa{kappa}"] = cell
        print("TAIL", f"kappa{kappa}", json.dumps(cell), flush=True)
    out["tail_r3"] = tg
    out["wall_s"] = time.perf_counter() - t_start
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print("wrote", OUT, round(out["wall_s"]), "s")


if __name__ == "__main__":
    main()
