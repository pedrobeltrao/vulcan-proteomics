#!/usr/bin/env python3
"""G1 sequence-only target-vs-decoy classifier (DESIGN §6.4 G1).  Diagnostic tool, never part of the pipeline.

Features (pure sequence properties; pseudo-reverse decoys keep composition and length, so only order-sensitive
features can separate them): one-hot residues at positions 1-3 and -2..-4 (position -1 = K/R is shared), dipeptide
frequencies (400), length, and AlphaPeptDeep's predicted RT (rt_pred: a sequence property, excluded from the scorer
for exactly this reason).  Learner: sklearn HistGradientBoosting (max_iter 300).
Training data: candidate peptidoforms (all ranks, both labels) of the 9 human TRAINING runs (/tmp/vulcan_m3/psms/
train_human; human_contam prod index), unique per peptidoform; pairs split by tdc_pair parity -> train on even pairs,
held-out AUC on odd pairs; plus AUC on D0 production rank-1 PSMs (diagnostic AUC reported by G1).
    conda_envs/pred/bin/python scripts/m3e/seqclf.py --train
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
PATH = os.path.join(W, "seqclf.pkl")
AA = "ACDEFGHIKLMNPQRSTVWY"
_IDX = np.full(256, -1, np.int64)
for _k, _a in enumerate(AA):
    _IDX[ord(_a)] = _k


def featurize(seqs, rt_pred) -> np.ndarray:
    n = len(seqs)
    X = np.zeros((n, 3 * 20 + 3 * 20 + 400 + 2), np.float32)
    for i, s in enumerate(seqs):
        c = _IDX[np.frombuffer(s.encode(), np.uint8)]
        L = c.size
        for k in range(3):
            if k < L and c[k] >= 0:
                X[i, 20 * k + c[k]] = 1
            j = L - 2 - k
            if j >= 0 and c[j] >= 0:
                X[i, 60 + 20 * k + c[j]] = 1
        dp = c[:-1] * 20 + c[1:]
        dp = dp[(c[:-1] >= 0) & (c[1:] >= 0)]
        np.add.at(X[i], 120 + dp, 1.0 / max(L - 1, 1))
        X[i, 520] = L
    X[:, 521] = np.asarray(rt_pred, np.float32)
    return X


class SeqClf:
    def __init__(self, model, info):
        self.model, self.info = model, info

    def logit(self, seqs, rt_pred):
        return self.model.decision_function(featurize(seqs, rt_pred))

    def logit_rows(self, l1, rows, rt_pred):
        return self.logit([l1.sequence(int(r)) for r in rows], rt_pred)


def load():
    with open(PATH, "rb") as fh:
        return pickle.load(fh)


def train():
    import m3_data as MD
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from vulcan.index.build import L1Index
    idx = C.index_dir("human_contam", "prod")
    lab = MD.labels(idx)
    l1 = L1Index(idx)
    rows, rt = [], []
    for _, psm, feat in MD.load_runs(os.path.join(C.WORK, "psms", "train_human")):
        rows.append(psm["pf_row"])
        rt.append(feat["rt_pred"])
    rows, rt = np.concatenate(rows), np.concatenate(rt)
    u, first = np.unique(rows, return_index=True)
    rows, rt = u, rt[first]
    rng = np.random.default_rng(11)
    dec = lab["is_decoy"][rows]
    pair = lab["tdc_pair"][rows]
    keep = np.r_[rng.permutation(np.flatnonzero(~dec))[:300000], rng.permutation(np.flatnonzero(dec))[:300000]]
    rows, rt, dec, pair = rows[keep], rt[keep], dec[keep], pair[keep]
    X = featurize([l1.sequence(int(r)) for r in rows], rt)
    tr = (pair % 2) == 0
    m = HistGradientBoostingClassifier(max_iter=300, random_state=11).fit(X[tr], ~dec[tr])
    auc_ho = float(roc_auc_score(~dec[~tr], m.decision_function(X[~tr])))
    # D0 production rank-1 PSMs (diagnostic)
    aucs = {}
    for _, psm, feat in MD.load_runs(C.psm_dir("D0_prod"))[:1]:
        r1 = psm["rank"] == psm["rank"].min()
        rr = psm["pf_row"][r1]
        aucs["D0_rank1_first_run"] = float(roc_auc_score(~lab["is_decoy"][rr], m.decision_function(
            featurize([l1.sequence(int(r)) for r in rr], feat["rt_pred"][r1]))))
    m_no_rt = HistGradientBoostingClassifier(max_iter=300, random_state=11).fit(X[tr][:, :-1], ~dec[tr])
    auc_ho_no_rt = float(roc_auc_score(~dec[~tr], m_no_rt.decision_function(X[~tr][:, :-1])))
    info = {"n_train": int(tr.sum()), "n_heldout": int((~tr).sum()), "auc_heldout_pairs": auc_ho,
            "auc_heldout_pairs_without_rt_pred": auc_ho_no_rt, **aucs,
            "features": "pos1-3 one-hot, pos -2..-4 one-hot, dipeptide freq, length, rt_pred"}
    print(json.dumps(info), flush=True)
    sys.path.insert(0, HERE)
    import seqclf as _mod          # pickle under module name 'seqclf' (importable by gates_a1), not __main__
    with open(PATH, "wb") as fh:
        pickle.dump(_mod.SeqClf(m, info), fh)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", action="store_true")
    if ap.parse_args().train:
        train()
