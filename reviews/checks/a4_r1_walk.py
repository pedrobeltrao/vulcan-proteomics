#!/usr/bin/env python3
"""Amendment A4 revisions 1-2: realistic DN3 walk (rev 2: DN2's admitted projects as tier 0, D26) (review B2, M2, M4, M5) and tail-G2 under the revision-1 cap rule
(review B3, M1).  METADATA + SYNTHETIC ONLY, seed 131.

What is new relative to a4_power.py (revision 0):
  * pool = revision-1 frame: lineage denies extended (Amoebozoa, Ciliophora, Mollicutes, Leptospira, dermatophytes),
    exclusions extended to every PXD accession in any tracked file under data/manifests/ or results/ and every
    accession with a cached PRIDE file listing (review M5);
  * the protocol screen is applied DETERMINISTICALLY per project to the cached 2026-10-26 crawl records
    (data/cache/pride_meta/projects_p*.json.gz: title, keywords, description, sample/data protocol), under the
    revision-1 rule (A4.3.5(a)); A2-verbatim and revision-0 rules are reported for comparison;
  * the real genus/lab component structure of the pool with <= 3 admitted projects per cluster ("cluster full"
    skips before listing), the reference-proteome step (DN2: 16/19), stage-2 marker exclusions that cost all R
    downloads, at an extra exclusion rate e for projects with R > 3;
  * the cluster cap is 10% and takes the minimum of the water-filling factors on the tail proxy and on A_cand; the
    stopping rule needs J_eff >= 15 on both capped quantities (review M2);
  * tail-G2 is evaluated with the revision-1 intersection-union rule: pass on the capped AND the uncapped pool
    (review B3), at the v1-minimum states of the realistic walk.
No PRIDE call, no file listing, no spectrum, no label.  The yield model is a4_power's (DN2 label-free yields).
    conda_envs/pe/bin/python reviews/checks/a4_r1_walk.py   # -> reviews/checks/a4_r1_walk.json (+ .out via tee)
"""
from __future__ import annotations

import glob
import gzip
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
sys.path.insert(0, ROOT)
import a4_pool_count as PC  # noqa: E402
import a4_power as AP  # noqa: E402
from vulcan.gates import a2_stats as A2S  # noqa: E402

SEED = 131
QUICK = os.environ.get("A4_QUICK") == "1"
OUT = os.path.join(HERE, "a4_r1_walk.json")

# ------------------------------------------------------------------------------------------- revision-1 frame
EXTRA_DENY_LINEAGE = ("Amoebozoa", "Ciliophora", "Mollicutes")
EXTRA_DENY_GENERA = {"Leptospira", "Trichophyton", "Microsporum", "Epidermophyton", "Mycoplasmopsis",
                     "Mesomycoplasma", "Malacoplasma"}


def tracked_pxd() -> set[str]:
    files = subprocess.run(["git", "ls-files", "data/manifests", "results"], cwd=ROOT, capture_output=True,
                           text=True).stdout.split()
    acc = set()
    for f in files:
        p = os.path.join(ROOT, f)
        if f.endswith("t4_frame.tsv") or os.path.getsize(p) > 50_000_000:
            continue
        try:
            raw = open(p, "rb").read()
            if f.endswith(".gz"):
                raw = gzip.decompress(raw)
            acc |= set(re.findall(rb"PXD\d{6}", raw))
        except Exception:
            continue
    acc = {a.decode() for a in acc}
    cached = {os.path.basename(x).split("_")[1] for x in glob.glob(os.path.join(ROOT, "data/cache/pride_meta/files_PXD*"))}
    return acc | cached


def pool_r1():
    a = pd.read_csv(os.path.join(ROOT, "data/cache/pride_meta/all_projects_classified.tsv"), sep="\t", dtype=str).fillna("")
    a = a[a.exclude_reasons.isin(["", "no_quant_or_phospho_evidence"])]
    a = a[(a.vendor == "thermo") & (a.acquisition == "DDA") & (a.enrichment == "none") & (a.label.isin(["LFQ", "UNK"]))]
    inst = a.instruments.str.split(";")
    a = a.assign(instr_tier=[("A1" if all(re.match(PC.A1_INSTR, x) for x in s) else
                              "EXT" if all(re.match(PC.A1_INSTR, x) or re.match(PC.EXT_INSTR, x) for x in s) else "")
                             if s != [""] else "" for s in inst])
    a = a[a.instr_tier != ""]
    ex0 = PC.exclusions()
    ex1 = ex0 | tracked_pxd()
    n_ex0, n_ex1 = int(a.accession.isin(ex0).sum()), int(a.accession.isin(ex1).sum())
    a = a[~a.accession.isin(ex1)]
    lin = json.load(open(os.path.join(HERE, "a4_genus_lineage.json")))
    a = a.assign(genus=a.organisms.str.split(" ").str[0].str.capitalize())
    a = a.assign(lineage=[lin.get(g, "") for g in a.genus])
    a = a.assign(cls=[PC.org_class(L) for L in a.lineage])
    deny = (a.genus.isin(PC.DENY_GENERA) | a.genus.isin(EXTRA_DENY_GENERA)
            | a.organisms.str.lower().str.contains(PC.DENY_SPECIES_RE)
            | a.lineage.apply(lambda L: any(k in L for k in EXTRA_DENY_LINEAGE)))
    k = a[~deny & a.cls.isin(["I_microbe", "II_land_plant"])].copy()
    k["tier"] = np.where(k.instr_tier == "EXT", "III", np.where(k.cls == "I_microbe", "I", "II"))
    k["labn"] = k.lab.str.lower().str.replace(r"[^a-z0-9]+", " ", regex=True).str.strip()
    return k, {"excluded_rev0": n_ex0, "excluded_rev1": n_ex1, "n_exclusion_accessions_rev1": len(ex1)}


# ------------------------------------------------------------------------------------------- protocol screens
A1_EXP = ("HLA|MHC|immunopept|cross-?link|XL-MS|glyco|SILAC|DIA\\b|SWATH|data.independent|phospho|single.cell|TMT|"
          "iTRAQ|ubiquitin|acetyl|N-?termin|top-down|crosslink|metaproteom|interactome|AP-MS|BioID|pull.?down|"
          "immunoprecip|labell?ing|dimethyl|PRM|SRM|targeted")
A1_CONT = "infect|host|co-?culture"
A2_ADD = ("serum|FBS|FCS|fetal|blood|plasma|patient|clinical|skin|wound|swab|sputum|tissue|biopsy|macrophage|"
          "cell line|HeLa|HEK|K562|spike|UPS|standard|benchmark|mixed proteome|entrapment|human")
A2_PROT = "brain heart|BHI|blood agar|peptone from (?:meat|casein)|Homo sapiens|human (?:database|proteome)|co-?culture|host"
A4_MULTI = "symbio|nodul|mycorrhiz|endophyt|aphid|herbivor|nematode|infest|microbiome|metagenom|allergen"
A4_MEDIA = "proteose peptone|meat extract|HL-?5\\b|peptone|tryptone|casein|yeast extract"
# revision-1 split: CONTAMINATION terms on every field; EXPERIMENT-TYPE / generic terms on sample-description fields
R1_CONT = ("serum|FBS|FCS|fetal|blood|plasma|patient|clinical|skin|wound|swab|sputum|biopsy|macrophage|"
           "cell line|HeLa|HEK|K562|brain heart|BHI|blood agar|Homo sapiens|human (?:database|proteome|sample|cell|"
           "tissue|plasma|serum)|co-?culture|host|infect|" + A4_MULTI + "|proteose peptone|meat extract|HL-?5\\b")
R1_SAMPLE = A1_EXP + "|tissue|spike|UPS|standard|benchmark|mixed proteome|entrapment|human|patient"
SAMPLE_FIELDS = ("title", "projectDescription", "sampleProcessingProtocol")
ALL_FIELDS = ("title", "projectDescription", "sampleProcessingProtocol", "dataProcessingProtocol")


def _txt(r, fields):
    t = " ".join(str(r.get(f) or "") for f in fields)
    return t + " " + " ".join(map(str, r.get("keywords") or []))


def protocol_flags(accs):
    recs = {}
    want = set(accs)
    for f in sorted(glob.glob(os.path.join(ROOT, "data/cache/pride_meta/projects_p*.json.gz"))):
        for r in json.load(gzip.open(f)):
            if r.get("accession") in want:
                recs[r["accession"]] = r
    out = {}
    for acc in accs:
        r = recs.get(acc)
        if r is None:
            out[acc] = {"record": False, "a2_verbatim": False, "rev0": False, "rev1": False}
            continue
        allt, samp = _txt(r, ALL_FIELDS), _txt(r, SAMPLE_FIELDS)
        a2 = not re.search("|".join([A1_EXP, A1_CONT, A2_ADD, A2_PROT]), allt, re.I)
        rev0 = a2 and not re.search(A4_MULTI + "|IgE|proteose peptone|meat extract|HL5", allt, re.I)
        rev1 = (not re.search(r"\b(?:%s)" % R1_CONT, allt, re.I) and not re.search(r"\bIgE\b", allt)
                and not re.search(r"\b(?:%s)" % R1_SAMPLE, samp, re.I))
        out[acc] = {"record": True, "a2_verbatim": bool(a2), "rev0": bool(rev0), "rev1": bool(rev1),
                    "media_peptone_mention": bool(re.search(r"\b(?:%s)" % A4_MEDIA, allt, re.I))}
    return out


# ------------------------------------------------------------------------------------------- cap and criteria
def jeff(x):
    x = np.asarray(x, float)
    return float(x.sum() ** 2 / (x ** 2).sum()) if x.sum() > 0 else 0.0


def capped(cp, ca, share):
    J = len(cp)
    if J * share < 1:
        return None
    f = np.minimum(AP.cap_factors(cp, share), AP.cap_factors(ca, share))
    return cp * f, ca * f, f


def criteria(cp, ca, c):
    r = capped(np.asarray(cp, float), np.asarray(ca, float), c["cap_share"])
    if r is None:
        return False, None
    xp, xa, f = r
    ok = (len(cp) >= c["J_min"] and xp.sum() >= c["proxy_min"] and xa.sum() >= c["A_cand_min"]
          and jeff(xp) >= c["Jeff_min"] and jeff(xa) >= c["Jeff_min"])
    return ok, (xp, xa, f)


# ------------------------------------------------------------------------------------------- realistic walk
def walk(rng, y, tiers, flags, *, order_seed, med_runs, run_cap, p_proteome, p_files, p_marker1, e_stage2,
         crit_full, crit_min, max_runs=2000, tail_scale=1.0, max_per_cluster=3, tier0=None, u_tier0=0.8,
         p_marker1_tier0=0.9):
    """tier0 (revision 2, D26 option (b)): DN2's admitted projects, walked FIRST under every A4 screen and cap.
    Their yields are DN2's label-free per-project values (per-run |A|, A_cand/|A|, tail proxy/|A|); their number of
    usable RAWs is Binomial(n_raw_records from DN2's listing, u_tier0), capped at run_cap; stage 1 of the marker
    screen passes w.p. p_marker1_tier0 (they passed A2's 3-run screen; A4 re-picks files), stage 2 as for others."""
    mA = math.log(y["per_run_A"]["geo_mean_positive"])
    sb, sw = y["per_run_A"]["log_sd_between_projects"], y["per_run_A"]["log_sd_within_project"]
    z0 = y["per_run_A"]["zero_run_fraction"]
    lt = math.log(y["tail_proxy_per_A"]["pooled"] * tail_scale)
    st = y["tail_proxy_per_A"]["log_sd_between_projects"]
    lc = math.log(y["A_cand_per_A"]["pooled"] / (1 - y["A_cand_per_A"]["pooled"]))
    sc = y["A_cand_per_A"]["logit_sd_between_projects"]
    adm, proxy, acand = [], [], []
    S = {"runs": 0, "dl": 0, "skip_full": 0, "walked": 0, "met_min": False, "snap": None, "full": False}

    def comps():
        n = len(adm)
        par = list(range(n))

        def f(i):
            while par[i] != i:
                par[i] = par[par[i]]
                i = par[i]
            return i
        by = {}
        for i, (g, l) in enumerate(adm):
            for key in (("g", g), ("l", l)):
                if not key[1]:
                    continue
                if key in by:
                    par[f(i)] = f(by[key])
                else:
                    by[key] = i
        return [f(i) for i in range(n)]

    def cluster_sums():
        c = comps()
        u = {r: j for j, r in enumerate(dict.fromkeys(c))}
        idx = np.array([u[r] for r in c])
        J = len(u)
        return (np.bincount(idx, weights=proxy, minlength=J), np.bincount(idx, weights=acand, minlength=J),
                np.array(proxy, float), np.array(acand, float), idx)

    def admit_and_check(g, l, R, a, pr, cr):
        adm.append((g, l))
        proxy.append(float(rng.poisson(a * pr)))
        acand.append(float(rng.binomial(int(a), cr)))
        S["runs"] += R
        cp, ca, pp, pa, idx = cluster_sums()
        if not S["met_min"]:
            ok, cs = criteria(cp, ca, crit_min)
            if ok:
                S["met_min"] = True
                S["snap"] = {"runs": S["runs"], "J": len(cp), "cp_capped": cs[0], "cp_uncapped": cp,
                             "proj_proxy": pp, "proj_idx": idx, "admitted": len(adm)}
        ok, cs = criteria(cp, ca, crit_full)
        if ok:
            S["full"] = True
            S["J"] = len(cp)
            S["admitted"] = len(adm)
            return True
        return False

    S["tier0_admitted"] = 0
    S["tier0_runs"] = 0
    for p0 in (sorted(tier0, key=lambda d: random.Random(order_seed).random()) if tier0 else []):
        S["walked"] += 1
        g, l = p0["genus"], p0["labn"]
        linked = [i for i, (gg, ll) in enumerate(adm) if gg == g or (l and ll == l)]
        if linked:
            c = comps()
            roots = {c[i] for i in linked}
            if sum(1 for x in c if x in roots) + 1 > max_per_cluster:
                S["skip_full"] += 1
                continue
        if not p0["protocol_rev1"]:
            continue
        R = int(min(run_cap, max(1, rng.binomial(p0["n_raw"], u_tier0))))
        if rng.random() > p_marker1_tier0:
            S["dl"] += min(R, 3)
            continue
        S["dl"] += R
        if R > 3 and rng.random() < e_stage2:
            continue
        a = float((rng.lognormal(math.log(max(p0["A_per_run"], 1.0)), sw, size=R) * (rng.random(R) > z0)).sum())
        S["tier0_admitted"] += 1
        S["tier0_runs"] += R
        if admit_and_check(g, l, R, a, p0["proxy_rate"] * tail_scale, p0["cand_rate"]):
            return S | {"end": "full stop"}

    for t in ("I", "II", "III"):
        cand = list(tiers[t][["accession", "genus", "labn"]].itertuples(index=False))
        random.Random(order_seed).shuffle(cand)
        for acc, g, l in cand:
            S["walked"] += 1
            linked = [i for i, (gg, ll) in enumerate(adm) if gg == g or (l and ll == l)]
            if linked:
                c = comps()
                roots = {c[i] for i in linked}
                if sum(1 for x in c if x in roots) + 1 > max_per_cluster:
                    S["skip_full"] += 1
                    continue
            if not flags[acc]["rev1"] or rng.random() > p_proteome or rng.random() > p_files:
                continue
            R = int(min(run_cap, max(1, round(rng.lognormal(math.log(med_runs), 0.8)))))
            if rng.random() > p_marker1:
                S["dl"] += min(R, 3)
                continue
            S["dl"] += R
            if R > 3 and rng.random() < e_stage2:
                continue
            if S["runs"] + R > max_runs:
                return S | {"end": "run cap"}
            mu = rng.normal(mA, sb)
            a = float((rng.lognormal(mu, sw, size=R) * (rng.random(R) > z0)).sum())
            tr = math.exp(rng.normal(lt, st))
            cr = 1 / (1 + math.exp(-rng.normal(lc, sc)))
            adm.append((g, l))
            proxy.append(float(rng.poisson(a * tr)))
            acand.append(float(rng.binomial(int(a), cr)))
            S["runs"] += R
            cp, ca, pp, pa, idx = cluster_sums()
            if not S["met_min"]:
                ok, cs = criteria(cp, ca, crit_min)
                if ok:
                    S["met_min"] = True
                    S["snap"] = {"runs": S["runs"], "J": len(cp), "cp_capped": cs[0], "cp_uncapped": cp,
                                 "proj_proxy": pp, "proj_idx": idx, "admitted": len(adm)}
            ok, cs = criteria(cp, ca, crit_full)
            if ok:
                S["full"] = True
                S["J"] = len(cp)
                S["admitted"] = len(adm)
                return S | {"end": "full stop"}
    cp, ca, *_ = cluster_sums() if adm else (np.zeros(0), np.zeros(0))
    S["J"] = len(cp)
    S["admitted"] = len(adm)
    return S | {"end": "pool exhausted"}


# ------------------------------------------------------------------------------------------- tail-G2 IUT
def tail_iut(snap, rng, kappa, scenario, conf=AP.TAIL_CONF_A4, B=2000):
    """Rev-1 rule: pass iff tail-G2 passes on the capped AND on the uncapped pool.  Clusters are the A4 units; the
    project-unit bound (review M1) coincides here because projects are drawn inside clusters with their own counts;
    the simulation uses the cluster unit for both (conservative only in expectation; see A4.5.2)."""
    cu = np.asarray(snap["cp_uncapped"], float)
    cc = np.asarray(snap["cp_capped"], float)
    J = cu.size
    kind, par = scenario
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
    elif kind == "largest":          # leak in the largest clusters so the UNCAPPED tail-weighted mean is 0.60
        pl = np.full(J, 0.5)
        order = np.argsort(-cu)
        need = 0.10 * cu.sum()
        for j in order:
            add = min(0.5 * cu[j], need)
            pl[j] = 0.5 + add / max(cu[j], 1e-9)
            need -= add
            if need <= 0:
                break
    # realise the uncapped tails, then subsample uniformly to the capped sizes (the cap subsamples spectra)
    nL = rng.poisson(kappa * cu)
    nF = rng.poisson(kappa * cu)
    kL = rng.binomial(nL, pl)
    kF = rng.binomial(nF, 0.5)
    f = np.where(cu > 0, cc / np.maximum(cu, 1e-12), 1.0)
    nLc = rng.binomial(nL, f)
    kLc = rng.hypergeometric(np.maximum(kL, 0), np.maximum(nL - kL, 0), nLc) if nL.sum() else nLc
    nFc = rng.binomial(nF, f)
    kFc = rng.hypergeometric(np.maximum(kF, 0), np.maximum(nF - kF, 0), nFc) if nF.sum() else nFc

    def gate(nl, kl, nf, kf):
        cl = np.concatenate([np.full(nl[j] + nf[j], j) for j in range(J)])
        tl = np.concatenate([np.r_[np.ones(nl[j], bool), np.zeros(nf[j], bool)] for j in range(J)])
        gl = np.concatenate([np.r_[np.arange(nl[j]) < kl[j], np.zeros(nf[j], bool)] for j in range(J)])
        tf = np.concatenate([np.r_[np.zeros(nl[j], bool), np.ones(nf[j], bool)] for j in range(J)])
        gf = np.concatenate([np.r_[np.zeros(nl[j], bool), np.arange(nf[j]) < kf[j]] for j in range(J)])
        return A2S.tail_g2_a2(cl, tl, gl, tf, gf, B=B, seed=rng, conf=conf)["verdict"]
    vu = gate(nL, kL, nF, kF)
    vc = gate(nLc, kLc, nFc, kFc)
    if vu == "pass" and vc == "pass":
        v = "pass"
    elif "fail" in (vu, vc):
        v = "fail"
    else:
        v = vu if vu != "pass" else vc
    return v, vc, vu, int(nLc.sum())


# ------------------------------------------------------------------------------------------- main
def main():
    t0 = time.perf_counter()
    out = {"seed": SEED, "spec": __doc__.split("    conda_envs")[0].strip()}
    k, exinfo = pool_r1()
    flags = protocol_flags(list(k.accession))
    k["rev1"] = [flags[a]["rev1"] for a in k.accession]
    k["a2v"] = [flags[a]["a2_verbatim"] for a in k.accession]
    k["rev0"] = [flags[a]["rev0"] for a in k.accession]
    tiers = {t: k[k.tier == t].sort_values("accession") for t in ("I", "II", "III")}
    pool = {"exclusions": exinfo}
    for t, v in tiers.items():
        gcount = v.genus.value_counts()
        pool[t] = {"projects": int(len(v)), "genera": int(v.genus.nunique()),
                   "admissible_by_genus_slots": int(gcount.clip(upper=3).sum()),
                   "protocol_pass": {"A2_verbatim": round(float(v.a2v.mean()), 3), "rev0": round(float(v.rev0.mean()), 3),
                                     "rev1": round(float(v.rev1.mean()), 3)},
                   "rev1_passers": int(v.rev1.sum()),
                   "rev1_passers_admissible_by_genus_slots": int(v[v.rev1].genus.value_counts().clip(upper=3).sum()),
                   "media_peptone_mention_among_rev1_passers": int(sum(flags[a].get("media_peptone_mention", False)
                                                                       for a in v[v.rev1].accession))}
    out["pool_r1"] = pool
    print("POOL", json.dumps(pool), flush=True)
    y = AP.dn2_yields()

    # Revision 2 (D26 option (b)): DN2's 6 admitted projects as tier 0, walked first, under every A4 screen.
    a_all = pd.read_csv(os.path.join(ROOT, "data/cache/pride_meta/all_projects_classified.tsv"), sep="\t",
                        dtype=str).fillna("").set_index("accession")
    tz = AP.dn2_tier0_projects()
    f0 = protocol_flags([p0["accession"] for p0 in tz])
    for p0 in tz:
        row = a_all.loc[p0["accession"]]
        p0["genus"] = row.organisms.split(" ")[0].capitalize()
        p0["labn"] = re.sub(r"[^a-z0-9]+", " ", row.lab.lower()).strip()
        p0["protocol_rev1"] = bool(f0[p0["accession"]]["rev1"])
        p0["instruments"] = row.instruments
    out["tier0_dn2"] = tz
    print("TIER0", json.dumps(tz), flush=True)

    # Revision-1/2 design (A4.3.6-A4.3.8): 10% cap on min(tail proxy, A_cand) factors, J_eff >= 10 on both capped
    # quantities, capped proxy >= 1,200, <= 5 projects per cluster, <= 48 files per project.
    DESIGN = {"cap_share": 0.10, "Jeff_min": 10, "proxy_min": 1200, "max_per_cluster": 5, "run_cap": 48}
    out["design_r2"] = DESIGN
    crit = lambda acand, d=DESIGN: {"J_min": 20, "proxy_min": d["proxy_min"], "A_cand_min": acand,
                                     "cap_share": d["cap_share"], "Jeff_min": d["Jeff_min"]}
    nw = 100 if QUICK else 300
    grid = []
    variants = [("rev2 design (DN2 tier 0, u=0.8)", dict(DESIGN, tier0=True, u0=0.8)),
                ("rev2 design (DN2 tier 0, u=0.5)", dict(DESIGN, tier0=True, u0=0.5)),
                ("rev1 design (no tier 0)", dict(DESIGN, tier0=False, u0=0.8))]
    for vname, d in variants:
        for med_runs in (6, 12, 24):
            for e2 in (0.15, 0.35):
                r = np.random.default_rng(SEED + 70)
                ws = [walk(r, y, tiers, flags, order_seed=1000 + i, med_runs=med_runs, run_cap=d["run_cap"],
                           p_proteome=16 / 19, p_files=0.93, p_marker1=0.5, e_stage2=e2,
                           crit_full=crit(450_000, d), crit_min=crit(120_000, d),
                           max_per_cluster=d["max_per_cluster"], tier0=tz if d["tier0"] else None,
                           u_tier0=d["u0"]) for i in range(nw)]
                q = lambda key: [float(np.percentile([w.get(key, 0) for w in ws], p)) for p in (10, 50, 90)]
                cell = {"variant": vname, "run_cap": d["run_cap"], "med_runs": med_runs, "stage2_extra_exclusion": e2,
                        "P_full": float(np.mean([w["full"] for w in ws])),
                        "P_v1_min": float(np.mean([w["met_min"] for w in ws])),
                        "runs": q("runs"), "J": q("J"), "admitted": q("admitted"), "downloaded": q("dl"),
                        "walked": q("walked"), "cluster_full_skips": q("skip_full"),
                        "tier0_admitted": q("tier0_admitted"), "tier0_runs": q("tier0_runs"),
                        "runs_at_v1_min": [float(np.percentile([w["snap"]["runs"] for w in ws if w["met_min"]], p))
                                           for p in (10, 50, 90)] if any(w["met_min"] for w in ws) else None}
                grid.append(cell)
                print("WALK", json.dumps({x: cell[x] for x in ("variant", "run_cap", "med_runs", "stage2_extra_exclusion",
                                                                "P_full", "P_v1_min", "runs", "J", "downloaded",
                                                                "tier0_admitted", "tier0_runs")}), flush=True)
    out["walk_r1"] = grid

    # tail-G2 rev-1 rule at the v1-minimum states (run cap 48, 12 runs/project, e2 0.35)
    tg = {}
    nrep = 120 if QUICK else 300
    for kappa in (1.0, 0.5):
        r = np.random.default_rng(SEED + 80)
        snaps = []
        i = 0
        while len(snaps) < nrep and i < 20 * nrep:
            w = walk(r, y, tiers, flags, order_seed=5000 + i, med_runs=12, run_cap=DESIGN["run_cap"],
                     p_proteome=16 / 19, p_files=0.93, p_marker1=0.5, e_stage2=0.35, crit_full=crit(450_000),
                     crit_min=crit(120_000), max_per_cluster=DESIGN["max_per_cluster"], tier0=tz)
            i += 1
            if w["met_min"]:
                snaps.append(w["snap"])
        cell = {"n_snapshots": len(snaps), "J_median": float(np.median([s["J"] for s in snaps])),
                "Jeff_capped_proxy_median": float(np.median([jeff(s["cp_capped"]) for s in snaps])),
                "Jeff_uncapped_proxy_median": float(np.median([jeff(s["cp_uncapped"]) for s in snaps]))}
        for lab, scen in (("clean", ("clean", None)), ("margin0.60_tau0.10", ("tau", (0.60, 0.10))),
                          ("mean0.60_leak0.90_in_25pct", ("concentrated", (0.25, 0.90))),
                          ("mean0.60_leak1.00_in_20pct", ("concentrated", (0.20, 1.00))),
                          ("uncapped_mean0.60_leak_in_largest", ("largest", None))):
            v = {"pass": 0, "fail": 0, A2S.NOT_EVALUABLE: 0, A2S.INVALID: 0}
            vc_pass = vu_pass = 0
            nt = []
            for s in snaps:
                verdict, vc, vu, n = tail_iut(s, r, kappa, scen)
                v[verdict] = v.get(verdict, 0) + 1
                vc_pass += vc == "pass"
                vu_pass += vu == "pass"
                nt.append(n)
            cell[lab] = {kk: x / len(snaps) for kk, x in v.items()} | {
                "capped_only_pass": vc_pass / len(snaps), "uncapped_only_pass": vu_pass / len(snaps),
                "n_tail_capped_median": float(np.median(nt))}
        tg[f"kappa{kappa}"] = cell
        print("TAIL", f"kappa{kappa}", json.dumps(cell), flush=True)
    out["tail_r1"] = tg
    out["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print("wrote", OUT, round(out["wall_s"]), "s")


if __name__ == "__main__":
    main()
