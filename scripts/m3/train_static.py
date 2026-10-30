#!/usr/bin/env python3
"""Train the FROZEN static-first scorers (DESIGN §4 L3) on project-disjoint training runs (data/manifests/train.tsv,
role=train; fold map data/manifests/train_folds.tsv).  Instrument class OT-HCD only (runs whose MS2 are not
FTMS+HCD are excluded and logged).

Training labels (Percolator-style, documented in NOTES): per run, PSM-level TDC q-values on rank-1 PSMs with the
current score; positives = rank-1 targets with q <= 1%; negatives = all rank-1 decoys PLUS the rank>=2 candidates
(either label) of spectra whose rank-1 is a positive (incorrect candidates; teaches within-spectrum re-ranking).
3 iterations (hyperscore -> model -> model).  Hyper-parameters chosen by leave-one-fold-out over PROJECT folds,
metric = rank-1 target PSMs at q <= 1% on held-out projects.  The final model is refit on all folds and frozen.
Variants: {linear, gbm} x {base (label-blind K2 features only), pred (+ AlphaPeptDeep fragment/RT features)}.

    conda_envs/pred/bin/python scripts/m3/train_static.py --out results/2026-10-27_M3_scoring/models
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402

from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, StaticModel, feature_matrix  # noqa: E402
from vulcan.stats.tdc import tdc_qvalues  # noqa: E402

ALPHA = 0.01
N_ITER = 3


def is_ot_hcd(store):
    t = ipc.open_file(store).read_all().select(["ms_level", "activation", "analyzer"])
    m = t.column("ms_level").to_numpy() == 2
    act = t.column("activation").to_numpy(zero_copy_only=False)[m]
    ana = t.column("analyzer").to_numpy(zero_copy_only=False)[m]
    return float(np.mean((act == 2) & (ana == 1)))


def load_training():
    rows = C.train_rows()
    fold = {r["accession"]: int(r["fold"]) for r in rows}
    by_run = {r["filename"][:-4]: r for r in rows}
    runs, excluded = [], {}
    for sp in ("human", "mouse", "ecoli"):
        d = C.psm_dir(f"train_{sp}")
        lab = MD.labels(C.index_dir(f"{sp}_contam", "prod"))
        for run, psm, feat in MD.load_runs(d):
            r = by_run[run]
            frac = is_ot_hcd(C.train_store(r))
            if frac < 0.95:
                excluded[run] = f"not OT-HCD (FTMS+HCD fraction of MS2 = {frac:.2f})"
                continue
            psm["is_decoy"] = lab["is_decoy"][psm["pf_row"]]
            runs.append({"run": run, "accession": r["accession"], "fold": fold[r["accession"]], "species": sp,
                         "psm": psm, "feat": feat})
    return runs, excluded


def run_targets(run, score):
    """Positives/negatives for one run given the current score (rank-1 TDC at PSM level)."""
    psm = run["psm"]
    r1 = psm["rank"] == 1
    dec = psm["is_decoy"]
    q = np.full(len(r1), np.inf)
    q[r1] = tdc_qvalues(score[r1], dec[r1])
    pos = r1 & ~dec & (q <= ALPHA)
    pos_spec = np.zeros(int(psm["spec_idx"].max()) + 1, bool)
    pos_spec[psm["spec_idx"][pos]] = True
    neg = (r1 & dec) | (~r1 & pos_spec[psm["spec_idx"]])
    return pos, neg, int(pos.sum())


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


def train(runs, names, kind, hp):
    Xs = [feature_matrix(r["psm"], r["feat"], names) for r in runs]
    scores = [r["psm"]["hyperscore"].astype(float) for r in runs]
    model, npos = None, []
    for it in range(N_ITER):
        P, N = [], []
        tot = 0
        for r, X, s in zip(runs, Xs, scores):
            pos, neg, n = run_targets(r, s)
            tot += n
            P.append(X[pos]); N.append(X[neg])
        npos.append(tot)
        Xp, Xn = np.vstack(P), np.vstack(N)
        X = np.vstack([Xp, Xn])
        y = np.r_[np.ones(len(Xp)), np.zeros(len(Xn))]
        model = fit(kind, X, y, hp)
        model.names = list(names)
        scores = [model.score(X_) for X_ in Xs]
    return model, npos


def evaluate(model, runs):
    n = 0
    for r in runs:
        s = model.apply(r["psm"], r["feat"])
        _, _, k = run_targets(r, s)
        n += k
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(C.RES, "models"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.perf_counter()
    runs, excluded = load_training()
    log = {"excluded_runs": excluded, "n_runs": len(runs),
           "runs": [{k: r[k] for k in ("run", "accession", "fold", "species")} for r in runs],
           "n_psms": int(sum(len(r["psm"]["rank"]) for r in runs)), "variants": {}}
    folds = sorted({r["fold"] for r in runs})
    grids = {"linear": [{"C": c} for c in (0.01, 0.1, 1.0)], "gbm": [{"n": n} for n in (50, 150, 300)]}
    for fs_name, names in (("base", BASE_FEATURES), ("pred", BASE_FEATURES + PRED_FEATURES)):
        for kind in ("linear", "gbm"):
            cv = {}
            for hp in grids[kind]:
                tot = 0
                for f in folds:
                    tr = [r for r in runs if r["fold"] != f]
                    te = [r for r in runs if r["fold"] == f]
                    m, _ = train(tr, names, kind, hp)
                    tot += evaluate(m, te)
                cv[json.dumps(hp)] = tot
            best = json.loads(max(cv, key=cv.get))
            model, npos = train(runs, names, kind, best)
            hs_only = sum(run_targets(r, r["psm"]["hyperscore"].astype(float))[2] for r in runs)
            name = f"static_{kind}_{fs_name}"
            model.meta = {"name": name, "instrument_class": "OT-HCD", "hyperparams": best,
                          "cv_heldout_psms_at_1pct": cv, "train_positives_per_iter": npos,
                          "train_psms_at_1pct_final": evaluate(model, runs), "train_psms_at_1pct_hyperscore": hs_only,
                          "training_runs": [r["run"] for r in runs],
                          "predictor": "AlphaPeptDeep generic (peptdeep 1.5.1, pretrained_models_v3)" if fs_name == "pred" else None}
            sha = model.save(os.path.join(a.out, name + ".pkl"))
            log["variants"][name] = {"hyperparams": best, "cv": cv, "sha256": sha, "npos": npos,
                                     "train_psms_at_1pct_final": model.meta["train_psms_at_1pct_final"],
                                     "train_psms_at_1pct_hyperscore": hs_only}
            print(name, log["variants"][name], flush=True)
    log["train_wall_s"] = time.perf_counter() - t0
    with open(os.path.join(a.out, "training_log.json"), "w") as fh:
        json.dump(log, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
