#!/usr/bin/env python3
"""A3.1.8(iii): DESIGN §6.7 Monte-Carlo verification of the FULL v2_sep competition and of pooled picked protein groups.
Synthetic data only.  >= 1e4 replicates per scenario; a list passes iff mean FDP - 2.326 * SE <= alpha.

The estimator under test is the committed code, called unchanged:
  vulcan.chimeric.compete.chimeric_competition (stage-1 winner per query unit -> CONGA neighbour filter over ALL group
  winners of a spectrum in descending score order -> paired double competition (peptidoform units, pair = tdc_pair)
  -> strata = chimeric rank of the PAIR WINNER pooled {1, 2, >=3} -> label-blind merge_adjacent (>= 200 pairs)
  -> q_sep = TDC within each merged stratum), and for proteins bench/t1/evaluate's rule: vulcan.stats.protein.
  picked_group_fdr (subset grouping, rescue, seed 11) over sequence-level best scores of the pair-winner table, POOLED
  across strata.

Synthetic replicate (scenario parameters below; every random quantity from numpy SeedSequence(131).spawn):
  * a peptide universe of n_pairs target/decoy pairs (peptidoform id = 2*pair + is_decoy; tdc_pair = pair), grouped
    into proteins (lengths 1..L, a share of peptides shared by two proteins; decoy proteins mirror the targets, rev_);
    present proteins are drawn, and true peptides come from present proteins only.
  * spectra with 1..4 precursor groups (query units; chimeric rank = group + 1) with the D0-like distribution p_groups.
    A unit is TRUE with probability pi[rank]: its true target gets a signal score.  Every unit gets K = 5 candidates
    (top-K analogue): the true peptide if any, plus null candidates = random pairs, side T or D with probability 1/2
    (exchangeable nulls), scores i.i.d. from the rank-specific null law Gumbel(loc[rank], scale[rank]) -- the
    stratum-heterogeneous nulls.
  * shadows: a secondary unit of a spectrum whose rank-1 unit is true carries, with probability p_shadow, the rank-1
    true peptide as a candidate (score = rank-1 true score + N(0, sd_shadow)): the same peptide seen through another
    precursor hypothesis, which the neighbour filter must remove (shared ions = 100%); if the shadow outscores the
    rank-1 winner, the RANK-1 winner is the one dropped.
  * ion sets for the neighbour filter (label-blind): each peptidoform has n_ions random bins; peptidoforms are in
    families of size fam (targets and decoys alike, separately) whose members share half their bins, so random null
    winners of one spectrum also collide (null drops, symmetric in label).
  * small universes make the same pair appear in units of different rank, so null pairs whose target and decoy best
    PSMs come from units of different rank (and different strata) are frequent; the pair's stratum is that of its
    winner, as in production.
  Truth: a target peptidoform is false iff its peptide is not a true peptide of any unit of the replicate; a target
  protein group is false iff none of its proteins is present.  FDP of an empty list = 0.
Lists (each at alpha = 1% and 5%): each MERGED stratum list (the lists A3.3.6 gates; q_sep within the merged block),
the union (all strata, common q_sep cutoff, as bench/t1 uses q_sep), and pooled picked protein groups.  Raw strata
{1, 2, >=3} under merged q are reported.  Verdict: pass iff every gated list passes in every scenario.

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/a3/dev_mc_v2sep.py --reps 10000 --jobs 12 --out <json>
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

SEED = 131
ALPHAS = (0.01, 0.05)
Z = 2.326

BASE = dict(n_spec=6000, p_groups=(0.41, 0.35, 0.16, 0.08), pi={1: 0.6, 2: 0.4, 3: 0.3}, K=5,
            null_loc={1: 0.0, 2: 0.4, 3: 0.8}, null_scale={1: 1.0, 2: 1.15, 3: 1.3},
            sig_mu={1: 9.0, 2: 9.0, 3: 9.0}, sig_sd={1: 2.0, 2: 2.2, 3: 2.4},
            p_shadow=0.5, sd_shadow=0.7, n_pairs=40000, prot_len_max=20, p_shared=0.05, frac_present=0.25,
            n_ions=24, n_bins=6000, fam=4, protein=True)
SCENARIOS = {
    "S1_dev_like": {},
    "S2_small_merge": dict(n_spec=700),
    "S3_strong_heterogeneity": dict(null_loc={1: 0.0, 2: 0.9, 3: 1.6}, null_scale={1: 1.0, 2: 1.35, 3: 1.7},
                                    sig_mu={1: 9.0, 2: 8.0, 3: 7.5}, pi={1: 0.6, 2: 0.3, 3: 0.2}),
    "S4_cross_rank_pairs": dict(n_pairs=4000, frac_present=0.4),
}


class FakeIons:
    """IonSets stand-in: .sets(rows) -> frozensets of ion bins, rows = peptidoform ids (built lazily from the bin
    matrix, so only winners are materialised)."""

    def __init__(self, bins):
        self.bins = bins
        self.cache = {}

    def sets(self, rows):
        out = []
        for r in rows:
            r = int(r)
            s = self.cache.get(r)
            if s is None:
                s = frozenset(self.bins[r].tolist())
                self.cache[r] = s
            out.append(s)
        return out


def make_ions(rng, n_pf, P):
    """Bin matrix (n_pf, n_ions): half the bins from the peptidoform's family core (families of size fam, built
    separately for targets and decoys), half its own."""
    n, h = P["n_ions"], P["n_ions"] // 2
    fam_id = np.empty(n_pf, np.int64)
    n_f = 0
    for side in (0, 1):
        ids = np.arange(side, n_pf, 2)
        perm = rng.permutation(ids)
        fam_id[perm] = n_f + np.arange(perm.size) // P["fam"]
        n_f += -(-perm.size // P["fam"])
    core = rng.integers(P["n_bins"], size=(n_f, h))
    own = rng.integers(P["n_bins"], size=(n_pf, n - h))
    return np.concatenate([core[fam_id], own], axis=1)


def make_proteins(rng, P):
    n = P["n_pairs"]
    lens = np.minimum(rng.geometric(0.15, size=n), P["prot_len_max"])
    starts = np.r_[0, np.cumsum(lens)]
    n_prot = int(np.searchsorted(starts, n, side="right"))
    owner = np.searchsorted(starts, np.arange(n), side="right") - 1
    pep_prot = [[int(o)] for o in owner]
    shared = np.flatnonzero(rng.random(n) < P["p_shared"])
    other = rng.integers(n_prot, size=shared.size)
    for k, q in zip(shared, other):
        if int(q) != pep_prot[k][0]:
            pep_prot[k].append(int(q))
    return n_prot, pep_prot


def replicate(seed_seq, P):
    from vulcan.chimeric.compete import chimeric_competition
    rng = np.random.default_rng(seed_seq)
    n_pf = 2 * P["n_pairs"]
    ions = FakeIons(make_ions(rng, n_pf, P))
    n_prot, pep_prot = make_proteins(rng, P)
    present_prot = rng.random(n_prot) < P["frac_present"]
    pres_pep = np.array([any(present_prot[q] for q in pp) for pp in pep_prot])
    present_peps = np.flatnonzero(pres_pep)
    ng = rng.choice(np.arange(1, 5), size=P["n_spec"], p=P["p_groups"])
    U = int(ng.sum())
    spec_u = np.repeat(np.arange(P["n_spec"]), ng)
    first_u = np.r_[0, np.cumsum(ng)[:-1]]
    grp_u = np.arange(U) - np.repeat(first_u, ng)
    rank_u = grp_u + 1
    rr_u = np.minimum(rank_u, 3)
    lut = lambda d: np.array([0.0, d[1], d[2], d[3]])  # noqa: E731
    true_u = rng.random(U) < lut(P["pi"])[rr_u]
    tp_u = present_peps[rng.integers(present_peps.size, size=U)]
    ts_u = rng.normal(lut(P["sig_mu"])[rr_u], lut(P["sig_sd"])[rr_u])
    prim = first_u[spec_u]                                   # rank-1 unit of each unit's spectrum
    shadow_u = (rank_u > 1) & true_u[prim] & (rng.random(U) < P["p_shadow"])
    sh_score = ts_u[prim] + rng.normal(0, P["sd_shadow"], size=U)
    K = P["K"]
    n_null = K - true_u.astype(int) - shadow_u.astype(int)
    nu = np.repeat(np.arange(U), n_null)
    pairs = rng.integers(P["n_pairs"], size=nu.size)
    side = rng.random(nu.size) < 0.5
    nscore = rng.gumbel(lut(P["null_loc"])[rr_u[nu]], lut(P["null_scale"])[rr_u[nu]])
    tu, su = np.flatnonzero(true_u), np.flatnonzero(shadow_u)
    unit = np.r_[tu, su, nu]
    pf = np.r_[2 * tp_u[tu], 2 * tp_u[prim[su]], 2 * pairs + side].astype(np.int64)
    score = np.r_[ts_u[tu], sh_score[su], nscore]
    spec, grp, rank = spec_u[unit], grp_u[unit], rank_u[unit]
    is_dec = (pf % 2).astype(bool)
    pair = pf // 2
    cc = chimeric_competition(spec, grp, rank, pf, pair, score, is_dec, pf, ions, seed=11)
    res = cc["res"]
    true_set = np.zeros(P["n_pairs"], bool)
    true_set[tp_u[tu]] = True
    tgt = ~res.is_decoy
    false_t = tgt & ~true_set[res.pair]
    merged = cc["merged"]
    out = {"merge": "|".join(sorted(set(merged.tolist()))), "n_dropped_nf": cc["n_dropped_by_neighbour_filter"],
           "n_group_winners": cc["n_group_winners"],
           "cross_rank_pairs": cross_rank_pairs(spec, grp, rank, pf, score, cc)}
    fdp = {}
    for a in ALPHAS:
        acc = tgt & (cc["q_sep"] <= a)
        lists = {"union": acc}
        for m in np.unique(merged):
            lists["merged:" + m] = acc & (merged == m)
        for k in (1, 2, 3):
            lists[f"raw:{k}"] = acc & (cc["stratum"] == k)
        for name, L in lists.items():
            n = int(L.sum())
            fdp[(name, a)] = (int((L & false_t).sum()) / n if n else 0.0, n)
    if P["protein"]:
        fdp.update(protein_fdp(res, true_set, pep_prot, present_prot))
    out["fdp"] = fdp
    return out


def cross_rank_pairs(spec, grp, rank, pf, score, cc):
    """Number of competed pairs whose best target PSM and best decoy PSM (among neighbour-filter survivors) come from
    units of different chimeric stratum."""
    sv = cc["survivors"]
    pr, d, st, sc = pf[sv] // 2, pf[sv] % 2, np.minimum(rank[sv], 3), score[sv]
    o = np.lexsort((-sc, d, pr))
    first = np.r_[True, (pr[o][1:] != pr[o][:-1]) | (d[o][1:] != d[o][:-1])]
    b = o[first]
    bt = {int(p): int(s) for p, dd, s in zip(pr[b], d[b], st[b]) if dd == 0}
    bd = {int(p): int(s) for p, dd, s in zip(pr[b], d[b], st[b]) if dd == 1}
    return int(sum(1 for p, s in bt.items() if p in bd and bd[p] != s))


def protein_fdp(res, true_set, pep_prot, present_prot):
    from vulcan.stats.protein import picked_group_fdr
    pep2prot, pep_score = {}, {}
    for p, d, s in zip(res.pair, res.is_decoy, res.score):
        k = ("D:" if d else "T:") + str(int(p))
        pep2prot[k] = [("rev_P" if d else "P") + str(q) for q in pep_prot[int(p)]]
        pep_score[k] = float(s)
    allp = {q for v in pep2prot.values() for q in v}
    base = lambda q: q[4:] if q.startswith("rev_") else q  # noqa: E731
    r = picked_group_fdr(pep2prot, pep_score, {q: base(q) for q in allp}, {q: q.startswith("rev_") for q in allp},
                         seed=11)
    out = {}
    for a in ALPHAS:
        n = f = 0
        for g, q, d, o in zip(r.groups, r.q, r.is_decoy, r.is_obsolete):
            if q <= a and not d and not o:
                n += 1
                f += not any(present_prot[int(x[1:])] for x in g if not x.startswith("rev_"))
        out[("protein", a)] = (f / n if n else 0.0, n)
    return out


def _work(args):
    name, seeds = args
    P = {**BASE, **SCENARIOS[name]}
    return [replicate(s, P) for s in seeds]


def run(reps, jobs, scen):
    from concurrent.futures import ProcessPoolExecutor
    out = {}
    root = np.random.SeedSequence(SEED)
    children = root.spawn(len(SCENARIOS))
    for i, name in enumerate(SCENARIOS):
        if scen and name not in scen:
            continue
        t0 = time.perf_counter()
        ss = children[i].spawn(reps)
        chunks = [(name, ss[k::jobs * 4]) for k in range(jobs * 4)]
        rows = []
        with ProcessPoolExecutor(jobs) as ex:
            for part in ex.map(_work, chunks):
                rows.extend(part)
        out[name] = summarize(rows, {**BASE, **SCENARIOS[name]}, time.perf_counter() - t0)
        print(name, json.dumps({k: v for k, v in out[name].items() if k in ("pass", "worst", "wall_s")}), flush=True)
    return out


def summarize(rows, P, wall):
    keys = sorted({k for r in rows for k in r["fdp"]}, key=lambda x: (x[1], x[0]))
    lists = {}
    ok = True
    worst = None
    for k in keys:
        v = np.array([r["fdp"][k][0] for r in rows if k in r["fdp"]])
        nd = np.array([r["fdp"][k][1] for r in rows if k in r["fdp"]])
        m = float(v.mean())
        se = float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan")
        gated = not k[0].startswith("raw:")
        p = bool(m - Z * se <= k[1]) if v.size > 1 else False
        rec = {"alpha": k[1], "n_reps": int(v.size), "mean_fdp": m, "se": se, "mean_minus_2.326se": m - Z * se,
               "mean_discoveries": float(nd.mean()), "gated": gated, "pass": p}
        lists[f"{k[0]}@{k[1]}"] = rec
        if gated:
            ok &= p
            marg = (m - Z * se) - k[1]
            if worst is None or marg > worst[1]:
                worst = (f"{k[0]}@{k[1]}", marg)
    merges = {}
    for r in rows:
        merges[r["merge"]] = merges.get(r["merge"], 0) + 1
    gw = np.array([r["n_group_winners"] for r in rows], float)
    dr = np.array([r["n_dropped_nf"] for r in rows], float)
    return {"params": {k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()}) for k, v in P.items()},
            "n_reps": len(rows), "lists": lists, "merge_patterns": merges,
            "neighbour_filter_drop_rate": float(dr.sum() / gw.sum()),
            "mean_cross_rank_pairs": float(np.mean([r["cross_rank_pairs"] for r in rows])),
            "pass": bool(ok), "worst": worst, "wall_s": wall}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10000)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--scenario", nargs="*")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    res = run(a.reps, a.jobs, a.scenario)
    out = {"spec": __doc__.split("    /tmp")[0].strip(), "seed": SEED, "reps_per_scenario": a.reps,
           "scenarios": res, "pass": bool(all(v["pass"] for v in res.values()))}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print("MC v2_sep pass:", out["pass"])


if __name__ == "__main__":
    main()
