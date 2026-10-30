#!/usr/bin/env python3
"""E7 rank-stratified read-out (pre-registered, config.yaml `e7_gain.by_rank`).  Dev only.

For a T1-driver work dir and a chimeric method (vulcan_chim_v1 | vulcan_v2):
  * per chimeric-rank stratum s in {1, 2, >=3} (stratum of the pair winner): arm P and arm C FDP estimates and §6.3
    verdicts at alpha in {1%, 5%} (convention (a) decides; (b) reported), counts restricted to the stratum:
      P: N_T = accepted targets in s; N_E = accepted entrapments in s; E_only / E>=T for entrapments in s, with the
         paired target's acceptance and score taken from the FULL list (any stratum);
      C: N_T, N_E = accepted forward / foreign peptidoforms in s; r as the driver (effective search space).
  * production IDs (fam T, dataset level, mean over decoy seeds) by stratum at own q (Group-walk) and at the method's
    matched arm-P FDP 1% / 5% threshold (from the driver's matched_fdp_curves in metrics.json), plus the pooled-TDC q.
  * closed v1 (`vulcan`) at the same matched FDP for the gain.
    python scripts/e7/strata.py --dataset D0 --method vulcan_chim_v1 --metrics results/.../t1/D0_e7/metrics.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from bench.t1 import evaluate as EV  # noqa: E402
from bench.t1 import ours_chimeric as OC  # noqa: E402
from bench.t1.dbs import parse_rep, replicate_names  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402
from vulcan.stats._util import DECOY_SEEDS, ENTRAP_SEEDS_DEV  # noqa: E402

COL = {"vulcan_chim_v1": "score_v1", "vulcan_v2": "score_v2", "vulcan_chim_v1_sep": "score_v1", "vulcan_v2_sep": "score_v2"}
STRATA = (1, 2, 3)
ALPHAS = (0.01, 0.05)


def paired_stratum(tab, alpha, k):
    t = tab[~tab.is_decoy & (tab.fam != -1)]
    key = EV._pair_key(t, "peptidoform").values
    ent = (t.fam == 1).values
    acc = (t.q <= alpha).values
    s = t.score.values
    st = t.stratum.values
    keys, inv = np.unique(key, return_inverse=True)
    m = len(keys)
    t_acc = np.zeros(m, np.int64)
    np.maximum.at(t_acc, inv[~ent], acc[~ent].astype(np.int64))
    t_s = np.full(m, -np.inf)
    np.maximum.at(t_s, inv[~ent & acc], s[~ent & acc])
    t_k = np.zeros(m, np.int64)
    np.maximum.at(t_k, inv[~ent], (acc & (st == k))[~ent].astype(np.int64))
    e_k = np.zeros(m, np.int64)
    np.maximum.at(e_k, inv[ent], (acc & (st == k))[ent].astype(np.int64))
    e_s = np.full(m, -np.inf)
    sel = ent & acc & (st == k)
    np.maximum.at(e_s, inv[sel], s[sel])
    e_only = (e_k == 1) & (t_acc == 0)
    e_gt = (e_k == 1) & (t_acc == 1) & (e_s >= t_s)
    keep = (t_k + e_k) > 0
    return ENT.EntrapmentCounts(key=keys[keep], t=t_k[keep], e=e_k[keep], e_only=e_only[keep].astype(np.int64),
                                e_gt_t=e_gt[keep].astype(np.int64))


def combined_stratum(tab, alpha, k):
    t = tab[~tab.is_decoy & (tab.fam != 1) & (tab.stratum == k)]
    key = (t.seq.astype(str) + "|" + t["mask"].astype(str)).values
    return ENT.combined_counts(key, (t.fam == -1).values, (t.q <= alpha).values)


def q_at_match(curve, f):
    fdp = np.asarray(curve["fdp_P_mean"])
    ids = np.asarray(curve["ids_run_summed_mean"])
    return EV.matched_ids(fdp, ids, f)[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--method", required=True, choices=sorted(COL))
    ap.add_argument("--metrics", required=True, help="driver metrics.json with matched_fdp_curves for method + vulcan")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=None)
    a = ap.parse_args()
    work = a.work or f"/tmp/vulcan_t1_{os.getuid()}/{a.dataset}"
    met = json.load(open(a.metrics))
    runs = met["runs"]
    curves = met["matched_fdp_curves"]
    reps = replicate_names(DECOY_SEEDS, ENTRAP_SEEDS_DEV)
    out = {"dataset": a.dataset, "method": a.method, "runs": runs, "by_stratum": {}, "ids": {}}
    qm = {m: {str(f): q_at_match(curves[m]["peptidoform"], f) for f in (0.01, 0.05)} for m in (a.method, "vulcan")}
    out["q_at_matched_fdp"] = qm
    P = {k: {al: [] for al in ALPHAS} for k in STRATA}
    C = {k: {al: [] for al in ALPHAS} for k in STRATA}
    ids = {k: {lab: [] for lab in ("own_q_0.01_gw", "own_q_0.01_pooled", "matched_0.01", "matched_0.05")} for k in STRATA}
    dec = {k: {str(al): [] for al in ALPHAS} for k in STRATA}     # decoy-estimated FDR in the stratum: D_k / T_k
    ids_closed = {lab: [] for lab in ("own_q_0.01", "matched_0.01", "matched_0.05")}
    diag = {}
    rC = None
    for rep in reps:
        kind, j, s = parse_rep(rep)
        rl = EV.RepLabels(os.path.join(work, "index", rep))
        tab, cc = OC.chim_units(os.path.join(work, "ours_chim", rep), rl, runs, COL[a.method], return_detail=True,
                                qcol="q_sep" if a.method.endswith("_sep") else "q_gw")
        if kind == "prod":
            fw = ~tab.is_decoy & (tab.fam == 0)
            for k in STRATA:
                sk = fw & (tab.stratum == k)
                ids[k]["own_q_0.01_gw"].append(int((sk & (tab.q <= 0.01)).sum()))
                ids[k]["own_q_0.01_pooled"].append(int((sk & (tab.q_pooled <= 0.01)).sum()))
                for f in ("0.01", "0.05"):
                    ids[k][f"matched_{f}"].append(int((sk & (tab.q <= qm[a.method][f])).sum()))
            for k in STRATA:
                for al in ALPHAS:
                    acc = (tab.stratum == k) & (tab.q <= al)
                    dec[k][str(al)].append((int((acc & tab.is_decoy).sum()), int((acc & ~tab.is_decoy).sum())))
            cl = EV.ours_units(os.path.join(work, "ours", rep), rl, runs)
            fc = ~cl.is_decoy & (cl.fam == 0)
            ids_closed["own_q_0.01"].append(int((fc & (cl.q <= 0.01)).sum()))
            for f in ("0.01", "0.05"):
                ids_closed[f"matched_{f}"].append(int((fc & (cl.q <= qm["vulcan"][f])).sum()))
            diag[rep] = {k: cc[k] for k in ("pair_counts", "merge", "n_group_winners", "n_after_neighbour_filter",
                                             "n_dropped_by_neighbour_filter")}
        elif kind == "P":
            for k in STRATA:
                for al in ALPHAS:
                    P[k][al].append(ENT.Replicate(paired_stratum(tab, al, k), s, j))
        else:
            rC = rl.n_fwd_pf[-1] / rl.n_fwd_pf[0]
            for k in STRATA:
                for al in ALPHAS:
                    C[k][al].append(ENT.Replicate(combined_stratum(tab, al, k), -1, j))
        print(rep, flush=True)
    for k in STRATA:
        rec = {"ids_prod_dataset_level_mean_over_decoy_seeds": {lab: float(np.mean(v)) for lab, v in ids[k].items()},
               "verdicts": {},
               "decoy_estimated_fdr_in_stratum_prod": {al: {"D": float(np.mean([d for d, _ in v])),
                                                            "T": float(np.mean([t for _, t in v])),
                                                            "D_over_T": float(np.sum([d for d, _ in v]) / max(1, np.sum([t for _, t in v])))}
                                                       for al, v in dec[k].items()}}
        for al in ALPHAS:
            v = {}
            for arm, reps_, est, r in (("P", P[k][al], "paired", 1.0), ("C", C[k][al], "combined", rC)):
                v[arm] = {}
                for cv in ("a", "b"):
                    e = ENT.evaluate_arm(reps_, al, est, r=r, n_boot=EV.N_BOOT, seed=EV.SEED, ci_convention=cv)
                    v[arm][cv] = {x: e[x] for x in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")}
            v["overall_a"] = ENT.combine_arms({"P": v["P"]["a"]["verdict"], "C": v["C"]["a"]["verdict"]},
                                              "static_low_capacity")
            v["overall_b"] = ENT.combine_arms({"P": v["P"]["b"]["verdict"], "C": v["C"]["b"]["verdict"]},
                                              "static_low_capacity")
            rec["verdicts"][str(al)] = v
        out["by_stratum"][str(k) if k < 3 else ">=3"] = rec
    out["closed_v1_ids_prod_mean"] = {lab: float(np.mean(v)) for lab, v in ids_closed.items()}
    out["competition_diag_prod"] = diag
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    for k, v in out["by_stratum"].items():
        print(k, v["ids_prod_dataset_level_mean_over_decoy_seeds"], v["decoy_estimated_fdr_in_stratum_prod"],
              {al: (round(x["P"]["a"]["mean_fdp"], 4), x["P"]["a"]["verdict"], round(x["C"]["a"]["mean_fdp"], 4),
                    x["C"]["a"]["verdict"]) for al, x in v["verdicts"].items()})
    print("closed", out["closed_v1_ids_prod_mean"])


if __name__ == "__main__":
    main()
