#!/usr/bin/env python3
"""Train scorer v2 by the pre-registered rule (results/2026-10-28_E7_chimeric_v2/config.yaml `training`), on the frozen
project-disjoint training set (train_v2_folds.tsv) prepared by scripts/e7/prep_train.py.  No dev data are read.

Labels (as M3, query unit = (spectrum, precursor group)): per run, PSM-level TDC on rank-1 PSMs of each query unit with
the current score; positives = rank-1 targets at q <= 1%; negatives = rank-1 decoys + rank >= 2 candidates of positive
units; 3 iterations (hyperscore -> model -> model).
Candidates: {linear, gbm} x {v1 features, v2 features}, grids linear C {0.01, 0.1, 1}, gbm max_iter {50, 150, 300}.
Model choice: leave-one-fold-out over project folds; metric = sum over held-out runs of target peptidoforms at q <= 1%
after the full chimeric pipeline (group winners -> CONGA neighbour filter -> paired double competition, pooled TDC,
per run).  Highest wins; a candidate within 0.5% of the best is preferred if it uses the v1 feature set, then linear.
Refit on all folds -> models/static_v2.pkl (sha256 in models/training_log.json).

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/e7/train_v2.py [--jobs N]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from vulcan.chimeric.compete import IonSets, chimeric_competition  # noqa: E402
from vulcan.rescore.features_v2 import V2_EXTRA, model_psm  # noqa: E402
from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, StaticModel, feature_matrix  # noqa: E402
from vulcan.stats.tdc import tdc_qvalues  # noqa: E402

W = os.environ.get("VULCAN_E7_WORK", "/tmp/vulcan_e7")
RES = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2")
INDEX = {"human": "/tmp/vulcan_e2/index/human_contam_prod", "mouse": "/tmp/vulcan_m3/index/mouse_contam_prod",
         "ecoli": "/tmp/vulcan_m3/index/ecoli_contam_prod"}
ALPHA = 0.01
N_ITER = 3
FS = {"v1": BASE_FEATURES + PRED_FEATURES, "v2": BASE_FEATURES + PRED_FEATURES + V2_EXTRA}
GRIDS = {"linear": [{"C": c} for c in (0.01, 0.1, 1.0)], "gbm": [{"n": n} for n in (50, 150, 300)]}
TOL = 0.005


def load_labels(index_dir):
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    return {k: t.column(k).to_numpy(zero_copy_only=False) for k in ("is_decoy", "tdc_pair")}


def load():
    from vulcan.index.build import L1Index
    rows = list(csv.DictReader(open(os.path.join(RES, "train_v2_folds.tsv")), delimiter="\t"))
    plog = json.load(open(os.path.join(W, "train", "prep_log.json")))
    labs, ions = {}, {}
    runs, excluded = [], {}
    for r in rows:
        run = r["filename"][:-4]
        if "excluded" in plog.get(run, {"excluded": "not prepared"}):
            excluded[run] = plog.get(run, {}).get("excluded", "not prepared")
            continue
        sp = r["species"]
        if sp not in labs:
            labs[sp] = load_labels(INDEX[sp])
            ions[sp] = IonSets(L1Index(INDEX[sp]))
        t = pq.read_table(os.path.join(W, "train", run + ".chim.psms.parquet"))
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in t.column_names}
        ft = pq.read_table(os.path.join(W, "train", run + ".chim.feat.parquet"))
        feat = {c: ft.column(c).to_numpy() for c in ft.column_names}
        psm["is_decoy"] = labs[sp]["is_decoy"][psm["pf_row"]]
        psm["pair"] = labs[sp]["tdc_pair"][psm["pf_row"]]
        runs.append({"run": run, "accession": r["accession"], "fold": int(r["fold"]), "species": sp,
                     "source": r["source"], "psm": psm, "feat": feat, "ions": ions[sp]})
    return runs, excluded


def run_targets(run, score):
    psm = run["psm"]
    r1 = psm["rank"] == 1
    dec = psm["is_decoy"]
    q = np.full(len(r1), np.inf)
    q[r1] = tdc_qvalues(score[r1], dec[r1])
    pos = r1 & ~dec & (q <= ALPHA)
    unit = psm["spec_idx"].astype(np.int64) * 64 + psm["group"].astype(np.int64)
    pu = np.unique(unit[pos])
    neg = (r1 & dec) | (~r1 & np.isin(unit, pu))
    return pos, neg


def fit(kind, X, y, hp, seed=11):
    if kind == "linear":
        from sklearn.linear_model import LogisticRegression
        mu, sd = X.mean(0), X.std(0)
        sd[sd == 0] = 1
        m = LogisticRegression(C=hp["C"], max_iter=2000, class_weight="balanced").fit((X - mu) / sd, y)
        return StaticModel("linear", None, {"w": m.coef_[0].copy(), "b": float(m.intercept_[0])}, mu, sd)
    from sklearn.ensemble import HistGradientBoostingClassifier
    m = HistGradientBoostingClassifier(max_depth=3, max_leaf_nodes=8, learning_rate=0.1, max_iter=hp["n"],
                                       l2_regularization=1.0, min_samples_leaf=200, early_stopping=False,
                                       class_weight="balanced", random_state=seed).fit(X, y)
    return StaticModel("gbm", None, m)


def train(runs, names, kind, hp, Xcache):
    Xs = [Xcache[(r["run"], tuple(names))] for r in runs]
    scores = [r["psm"]["hyperscore"].astype(float) for r in runs]
    model, npos = None, []
    for _ in range(N_ITER):
        P, N = [], []
        for r, X, s in zip(runs, Xs, scores):
            pos, neg = run_targets(r, s)
            P.append(X[pos])
            N.append(X[neg])
        Xp, Xn = np.vstack(P), np.vstack(N)
        npos.append(len(Xp))
        model = fit(kind, np.vstack([Xp, Xn]), np.r_[np.ones(len(Xp)), np.zeros(len(Xn))], hp)
        model.names = list(names)
        scores = [model.score(X_) for X_ in Xs]
    return model, npos


def heldout_ids(model, run, X):
    return pipeline_ids(run, model.score(X))


def pipeline_ids(run, s):
    """Target peptidoforms at q <= 1% after the full chimeric pipeline (pooled TDC) for one run."""
    psm = run["psm"]
    cc = chimeric_competition(psm["spec_idx"].astype(np.int64), psm["group"], psm["chim_rank"], psm["pf_id"],
                              psm["pair"], s, psm["is_decoy"], psm["pf_row"], run["ions"], seed=11)
    r = cc["res"]
    return int(((r.q <= ALPHA) & ~r.is_decoy).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(RES, "models"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.perf_counter()
    runs, excluded = load()
    Xc = {}
    for r in runs:
        for names in FS.values():
            Xc[(r["run"], tuple(names))] = feature_matrix(model_psm(r["psm"]), r["feat"], names)
    folds = sorted({r["fold"] for r in runs})
    log = {"rule": __doc__.split("Model choice:")[1].split("Refit")[0].strip(), "excluded_runs": excluded,
           "runs": [{k: r[k] for k in ("run", "accession", "fold", "species", "source")} for r in runs],
           "n_runs": len(runs), "n_psms": int(sum(len(r["psm"]["rank"]) for r in runs)), "cv": {}}
    print(f"{len(runs)} runs, {log['n_psms']} PSMs, folds {folds}; excluded {excluded}", flush=True)
    hs_ids = sum(pipeline_ids(r, r["psm"]["hyperscore"].astype(float)) for r in runs)
    log["hyperscore_ids_all_runs"] = hs_ids
    cand = []
    for fs_name, names in FS.items():
        for kind, grid in GRIDS.items():
            for hp in grid:
                tot, per_fold = 0, {}
                for f in folds:
                    tr = [r for r in runs if r["fold"] != f]
                    te = [r for r in runs if r["fold"] == f]
                    m, _ = train(tr, names, kind, hp, Xc)
                    k = sum(heldout_ids(m, r, Xc[(r["run"], tuple(names))]) for r in te)
                    per_fold[f] = k
                    tot += k
                key = f"{fs_name}|{kind}|{json.dumps(hp)}"
                log["cv"][key] = {"heldout_ids": tot, "per_fold": per_fold}
                cand.append((tot, fs_name, kind, hp))
                print(key, tot, per_fold, f"{time.perf_counter() - t0:.0f}s", flush=True)
                json.dump(log, open(os.path.join(a.out, "training_log.json"), "w"), indent=1, default=float)
    best = max(c[0] for c in cand)
    near = [c for c in cand if c[0] >= (1 - TOL) * best]
    pref = sorted(near, key=lambda c: (c[1] != "v1", c[2] != "linear", -c[0]))
    tot, fs_name, kind, hp = pref[0]
    log["best_cv"] = best
    log["within_tolerance"] = [f"{c[1]}|{c[2]}|{json.dumps(c[3])}: {c[0]}" for c in near]
    log["chosen"] = {"features": fs_name, "kind": kind, "hyperparams": hp, "cv_heldout_ids": tot}
    model, npos = train(runs, FS[fs_name], kind, hp, Xc)
    model.meta = {"name": "static_v2", "instrument_class": "OT-HCD", "features_set": fs_name, "hyperparams": hp,
                  "train_positives_per_iter": npos, "training_runs": [r["run"] for r in runs],
                  "chimeric": True, "predictor": "AlphaPeptDeep generic (peptdeep 1.5.1, pretrained_models_v3)",
                  "rule": log["rule"]}
    sha = model.save(os.path.join(a.out, "static_v2.pkl"))
    log["chosen"]["sha256"] = sha
    log["chosen"]["train_positives_per_iter"] = npos
    log["train_wall_s"] = time.perf_counter() - t0
    json.dump(log, open(os.path.join(a.out, "training_log.json"), "w"), indent=1, default=float)
    print("chosen", log["chosen"], flush=True)


if __name__ == "__main__":
    main()
