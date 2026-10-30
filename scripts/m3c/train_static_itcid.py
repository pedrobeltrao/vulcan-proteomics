#!/usr/bin/env python3
"""Train the FROZEN IT-CID static-first scorers (tree 3c) on the project-disjoint IT-CID training runs
(data/manifests/train_itcid.tsv, fold map train_itcid_folds.tsv).  Recipe identical to the OT-HCD M3 scorer
(scripts/m3/train_static.py), only the instrument class, search tolerance and predictor differ:

Labels (Percolator-style, training runs only): per run, PSM-level TDC on rank-1 PSMs with the current score;
positives = rank-1 targets at q <= 1%; negatives = rank-1 decoys + rank>=2 candidates of positive spectra.
3 iterations (hyperscore -> model -> model).  Hyper-parameters by leave-one-PROJECT-fold-out (3 folds), metric =
rank-1 target PSMs at q <= 1% on held-out projects; final model refit on all folds and frozen (sha256 in .json).
Variants {linear (standardised logistic regression, C in {0.01, 0.1, 1}), gbm (HistGBM depth 3 / 8 leaves /
min leaf 200, n_iter in {50, 150, 300})} x {base (K2 label-blind), pred (+ Prosit CID/TMT-CID features)}.
Model choice: the pre-stated rule in results/2026-10-27_3c_itcid/NOTES.md (written before any dev evaluation).

    python scripts/m3c/train_static_itcid.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m3"))
import common as C  # noqa: E402
import m3_data as MD  # noqa: E402
import train_static as TS  # noqa: E402  (M3 recipe: run_targets, fit, train, evaluate)

import numpy as np  # noqa: E402

from vulcan.rescore.itcid import ITCID_PRED_FEATURES  # noqa: E402
from vulcan.rescore.static import BASE_FEATURES  # noqa: E402


def load_training():
    rows = C.train_rows()
    runs = []
    for r in rows:
        fam = C.train_family(r)
        run = os.path.splitext(r["filename"])[0]
        d = C.psm_dir("train")
        lab = MD.labels(C.index_dir(fam, "prod"))
        import pyarrow.parquet as pq
        t = pq.read_table(os.path.join(d, run + ".psms.parquet"))
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        ft = pq.read_table(os.path.join(d, run + ".feat.parquet"))
        feat = {c: ft.column(c).to_numpy() for c in ft.column_names}
        psm["is_decoy"] = lab["is_decoy"][psm["pf_row"]]
        runs.append({"run": run, "accession": r["accession"], "fold": int(r["fold"]), "species": r["species"],
                     "stratum": r["stratum"], "psm": psm, "feat": feat})
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(C.RES, "models"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.perf_counter()
    runs = load_training()
    log = {"n_runs": len(runs), "runs": [{k: r[k] for k in ("run", "accession", "fold", "species", "stratum")}
                                        for r in runs],
           "n_psms": int(sum(len(r["psm"]["rank"]) for r in runs)), "variants": {}}
    folds = sorted({r["fold"] for r in runs})
    grids = {"linear": [{"C": c} for c in (0.01, 0.1, 1.0)], "gbm": [{"n": n} for n in (50, 150, 300)]}
    hs_only = sum(TS.run_targets(r, r["psm"]["hyperscore"].astype(float))[2] for r in runs)
    hs_fold = {f: sum(TS.run_targets(r, r["psm"]["hyperscore"].astype(float))[2] for r in runs if r["fold"] == f)
               for f in folds}
    log["hyperscore_psms_at_1pct"] = {"total": hs_only, "per_fold": hs_fold}
    for fs_name, names in (("base", BASE_FEATURES), ("pred", BASE_FEATURES + ITCID_PRED_FEATURES)):
        for kind in ("linear", "gbm"):
            cv, cv_fold = {}, {}
            for hp in grids[kind]:
                tot, per = 0, {}
                for f in folds:
                    tr = [r for r in runs if r["fold"] != f]
                    te = [r for r in runs if r["fold"] == f]
                    m, _ = TS.train(tr, names, kind, hp)
                    per[f] = TS.evaluate(m, te)
                    tot += per[f]
                cv[json.dumps(hp)] = tot
                cv_fold[json.dumps(hp)] = per
            best = json.loads(max(cv, key=cv.get))
            model, npos = TS.train(runs, names, kind, best)
            name = f"itcid_{kind}_{fs_name}"
            model.meta = {"name": name, "instrument_class": "IT-CID", "hyperparams": best,
                          "cv_heldout_psms_at_1pct": cv, "cv_heldout_per_fold": cv_fold,
                          "train_positives_per_iter": npos,
                          "train_psms_at_1pct_final": TS.evaluate(model, runs),
                          "train_psms_at_1pct_hyperscore": hs_only,
                          "training_runs": [r["run"] for r in runs],
                          "search": {k: list(v) if isinstance(v, tuple) else v for k, v in C.SEARCH.items()},
                          "predictor": ("Prosit_2020_intensity_CID / Prosit_2020_intensity_TMT (CID, CE 35) + "
                                        "Prosit_2019_irt / Prosit_2020_irt_TMT via Koina") if fs_name == "pred" else None}
            sha = model.save(os.path.join(a.out, name + ".pkl"))
            log["variants"][name] = {"hyperparams": best, "cv": cv, "cv_per_fold": cv_fold, "sha256": sha,
                                     "npos": npos, "train_psms_at_1pct_final": model.meta["train_psms_at_1pct_final"]}
            print(name, json.dumps(log["variants"][name]), flush=True)
    log["train_wall_s"] = time.perf_counter() - t0
    json.dump(log, open(os.path.join(a.out, "training_log.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
