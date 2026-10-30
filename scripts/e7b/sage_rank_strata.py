#!/usr/bin/env python3
"""E7b descriptive read-out (pre-registered in config.yaml `additional_calibration_readout`): arm-P entrapment FDP of
the Sage-chimera comparators at their own q, split by the Sage rank of the peptidoform's best PSM (1 vs >= 2).

  sage_chim<R>               results.sage.tsv, ALL ranks; q = peptide_q, score = sage_discriminant_score; stratum = rank
                             of the highest-scoring PSM of the peptidoform
  sage_chim<R>+pred+mokapot  mokapot.{,decoy.}peptides.txt; stratum = rank of the peptide's representative PSM
                             (ScanNr = scan * 16 + rank, config.yaml)
Counts per stratum as scripts/e7/strata.py (paired estimator; the paired target's acceptance from the full list);
convention (a) CIs.  Also prod IDs by stratum (mean over decoy seeds).  Descriptive only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "e7"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from bench.t1 import comparators as CMP  # noqa: E402
from bench.t1 import evaluate as EV  # noqa: E402
from bench.t1.dbs import parse_rep, replicate_names  # noqa: E402
from strata import paired_stratum  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402
from vulcan.stats._util import DECOY_SEEDS, ENTRAP_SEEDS_DEV  # noqa: E402

ALPHAS = (0.01, 0.05)


def units(rows):
    """rows: (peptide, is_decoy, q, score, rank) -> peptidoform DataFrame with stratum (1, or 2 for rank >= 2)."""
    d = {}
    for pep, dec, q, s, rk in rows:
        k = CMP.parse_peptide(pep) + (bool(dec),)
        cur = d.get(k)
        if cur is None:
            d[k] = [q, s, rk]
        else:
            cur[0] = min(cur[0], q)
            if s > cur[1]:
                cur[1], cur[2] = s, rk
    df = pd.DataFrame([(k[0], k[1], k[2], v[0], v[1], 1 if v[2] == 1 else 2) for k, v in d.items()],
                      columns=["seq", "mask", "is_decoy", "q", "score", "stratum"])
    df["mask"] = df["mask"].astype(np.uint64)
    df["runbits"] = 0
    return df


def load(method, cmp_root, rep):
    R = int(method.split("sage_chim")[1].split("+")[0])
    if method.endswith("+pred+mokapot"):
        d = os.path.join(cmp_root, f"sage_chim{R}_pred_mokapot", rep)
        rows = []
        for fn, dec in (("mokapot.peptides.txt", False), ("mokapot.decoy.peptides.txt", True)):
            p = os.path.join(d, fn)
            if os.path.exists(p):
                for x in CMP.read_ragged(p):
                    rows.append((x["Peptide"], dec, float(x["mokapot q-value"]), float(x["mokapot score"]),
                                 int(x["ScanNr"]) % 16))
        return rows
    t = pd.read_csv(os.path.join(cmp_root, f"sage_chim{R}", rep, "results.sage.tsv"), sep="\t",
                    usecols=["peptide", "label", "peptide_q", "sage_discriminant_score", "rank"])
    return list(zip(t.peptide.values, t.label.values == -1, t.peptide_q.values, t.sage_discriminant_score.values,
                    t["rank"].values))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--methods", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    work = f"/tmp/vulcan_t1_{os.getuid()}/{a.dataset}"
    runs = None
    reps = replicate_names(DECOY_SEEDS, ENTRAP_SEEDS_DEV)
    out = {"dataset": a.dataset, "methods": {}}
    for m in a.methods:
        P = {k: {al: [] for al in ALPHAS} for k in (1, 2)}
        ids = {k: {al: [] for al in ALPHAS} for k in (1, 2)}
        for rep in reps:
            kind, j, s = parse_rep(rep)
            if kind == "C":
                continue
            rl = EV.RepLabels(os.path.join(work, "index", rep))
            tab = EV._finish(units(load(m, os.path.join(work, "cmp"), rep)), rl, runs or [])
            if kind == "prod":
                fw = ~tab.is_decoy & (tab.fam == 0)
                for k in (1, 2):
                    for al in ALPHAS:
                        ids[k][al].append(int((fw & (tab.stratum == k) & (tab.q <= al)).sum()))
            else:
                for k in (1, 2):
                    for al in ALPHAS:
                        P[k][al].append(ENT.Replicate(paired_stratum(tab, al, k), s, j))
        rec = {}
        for k in (1, 2):
            rk = {}
            for al in ALPHAS:
                e = ENT.evaluate_arm(P[k][al], al, "paired", r=1.0, n_boot=EV.N_BOOT, seed=EV.SEED, ci_convention="a")
                rk[str(al)] = {"prod_ids_mean": float(np.mean(ids[k][al])),
                               **{x: e[x] for x in ("mean_fdp", "ci95", "verdict", "mean_lb", "mean_n_t")}}
            rec["rank_1" if k == 1 else "rank_ge2"] = rk
        out["methods"][m] = rec
        print(m, json.dumps({kk: {al: (round(v["mean_fdp"], 4), v["verdict"], round(v["prod_ids_mean"]))
                                  for al, v in vv.items()} for kk, vv in rec.items()}), flush=True)
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
