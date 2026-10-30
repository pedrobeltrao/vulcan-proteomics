"""RESET: REScoring via Estimating and Training (Freestone, Käll, Noble, Keich).

References
----------
* J. Freestone, L. Käll, W.S. Noble, U. Keich. "How to train a post-processor for
  tandem mass spectrometry proteomics database search while maintaining control
  of the false discovery rate." bioRxiv 10.1101/2023.10.26.564068 (RECOMB 2024,
  LNCS 14758:448; J Proteome Res 2025 24(5):2266). Percolator-RESET.
* Barber & Candès 2015 (SeqStep+ / knockoff+), for the "+1" finite-sample form.

Exact estimator implemented here (DESIGN.md §6.4; written explicitly per the
round-2 review [m1])
--------------------------------------------------------------------------------
Input: competition winners *after* PSM-level and paired peptide-level competition
using the original (search-engine) score. Competition is never redone with the
learned score.

1. Split. Every winning decoy is independently assigned, with probability 1/2
   (seeded), to the *training* set or the *estimating* set.
   Pseudo-targets = winning targets + estimating decoys (labels hidden from the
   learner; they enter training as unlabeled pseudo-targets).

2. Training (``n_iter`` iterations, default 5, initial score = search-engine score).
   At each iteration rank pseudo-targets and training decoys by the current
   score and compute, at every threshold t,

       FDR_train(t) = c_train * (D_train(t) + 1) / max(PT(t), 1)

   with c_train = 3 (single-decoy: a null winner is a pseudo-target with
   probability 3/4 vs. a training decoy 1/4) or c_train = 2 (two-decoy: 2/3 vs
   1/3). Positives = pseudo-targets with monotone FDR_train q-value <= alpha_train
   (default 0.01); negatives = all training decoys. Fit the learner, rescore all
   winners.

3. Reporting. Discard training decoys. Rank targets and estimating decoys by the
   final learned score; at every threshold t

       FDR_final(t) = c_final * (D_est(t) + 1) / max(T(t), 1)

   with c_final = 2 (single-decoy; remaining nulls are targets w.p. 2/3) or
   c_final = 1 (two-decoy; target and estimating decoy equally likely). The "+1"
   is added to the estimating-decoy count *before* multiplying by c (SeqStep+ with
   c/(1-c) = c_final). q-values are the monotone version; report targets with
   q <= alpha. Consequence: at alpha = 0.1 % single-decoy RESET cannot report
   anything with fewer than 2000 targets.

Two-decoy mode: each target competes with two decoys from *independent* shuffles
(exchangeable; not reverse + shuffle). The winner is kept; winning decoys are then
split at random as above (not by decoy database).

Learners are pluggable: any object with ``fit(X_pos, X_neg) -> self`` and
``decision_function(X) -> scores``. Default: :class:`LDALearner` (closed-form,
fast, deterministic). :class:`LinearSVMLearner` (sklearn LinearSVC, Percolator-like
class-weight grid C+, C- in {0.1, 1, 10} chosen by 3-fold CV on the training data
only) is available.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._util import make_rng
from .tdc import compete, tdc_qvalues

__all__ = ["LDALearner", "LinearSVMLearner", "reset", "ResetResult", "RESET_CONSTANTS",
           "two_decoy_competition"]

RESET_CONSTANTS = {
    "single": {"c_train": 3.0, "c_final": 2.0},
    "two": {"c_train": 2.0, "c_final": 1.0},
}


class LDALearner:
    """Fisher LDA with ridge shrinkage on standardized features."""

    def __init__(self, shrinkage: float = 1e-2):
        self.shrinkage = shrinkage
        self.w = None
        self.b = 0.0

    def fit(self, X_pos, X_neg):
        mp, mn = X_pos.mean(axis=0), X_neg.mean(axis=0)
        Xc = np.vstack([X_pos - mp, X_neg - mn])
        S = Xc.T @ Xc / max(Xc.shape[0] - 2, 1)
        S += self.shrinkage * np.eye(S.shape[0]) * max(np.trace(S) / S.shape[0], 1e-12)
        self.w = np.linalg.solve(S, mp - mn)
        self.b = -0.5 * float((mp + mn) @ self.w)
        return self

    def decision_function(self, X):
        return X @ self.w + self.b


class LinearSVMLearner:
    """Linear SVM (sklearn LinearSVC) with Percolator-style class-weight grid.

    The grid (C_pos, C_neg) in {0.1, 1, 10}^2 is chosen by 3-fold CV on the
    *training data only* (selected positives vs training decoys), maximising the
    held-out ROC AUC. No estimating-decoy label is ever used.
    """

    def __init__(self, grid=(0.1, 1.0, 10.0), cv: int = 3, seed: int = 11, max_iter: int = 5000):
        self.grid = grid
        self.cv = cv
        self.seed = seed
        self.max_iter = max_iter
        self.model = None

    def _fit_one(self, X, y, cp, cn):
        from sklearn.svm import LinearSVC

        m = LinearSVC(C=1.0, class_weight={1: cp, 0: cn}, dual="auto", max_iter=self.max_iter)
        return m.fit(X, y)

    def fit(self, X_pos, X_neg):
        X = np.vstack([X_pos, X_neg])
        y = np.r_[np.ones(len(X_pos), dtype=int), np.zeros(len(X_neg), dtype=int)]
        rng = np.random.default_rng(self.seed)
        fold = rng.integers(0, self.cv, size=len(y))
        best, best_auc = None, -np.inf
        for cp in self.grid:
            for cn in self.grid:
                aucs = []
                for f in range(self.cv):
                    tr, te = fold != f, fold == f
                    if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                        continue
                    m = self._fit_one(X[tr], y[tr], cp, cn)
                    aucs.append(_auc(m.decision_function(X[te]), y[te]))
                auc = np.mean(aucs) if aucs else -np.inf
                if auc > best_auc:
                    best, best_auc = (cp, cn), auc
        cp, cn = best if best is not None else (1.0, 1.0)
        self.model = self._fit_one(X, y, cp, cn)
        return self

    def decision_function(self, X):
        return self.model.decision_function(X)


def _auc(score, y):
    from scipy.stats import rankdata

    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)


@dataclass
class ResetResult:
    score: np.ndarray             # final learned score for every input row
    is_training_decoy: np.ndarray
    q: np.ndarray                 # final q-value; NaN for training decoys
    mode: str
    c_train: float
    c_final: float
    n_positives: list = field(default_factory=list)  # positives used per iteration

    def accepted(self, alpha: float, is_decoy) -> np.ndarray:
        """Boolean mask of reported targets."""
        is_decoy = np.asarray(is_decoy, dtype=bool)
        return (~self.is_training_decoy) & (~is_decoy) & (self.q <= alpha)


def _standardize(X):
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def reset(X, is_decoy, init_score, mode: str = "single", seed=None, alpha_train: float = 0.01,
          n_iter: int = 5, learner=None, min_positives: int = 10, standardize: bool = True) -> ResetResult:
    """Run RESET on peptide-level competition winners.

    Parameters
    ----------
    X : (n, p) features of the winning PSM of each competition (label-blind).
    is_decoy : winner label (any decoy database in two-decoy mode).
    init_score : search-engine score used for competition and the first iteration
        (as in the RESET paper). ``None``: Percolator-style start, the single
        standardized feature (either sign) giving the most training positives;
        this choice uses only training-decoy labels, so it is RESET-compatible.
    mode : "single" or "two".
    seed : seed for the training/estimating split.
    """
    if mode not in RESET_CONSTANTS:
        raise ValueError("mode must be 'single' or 'two'")
    c_train = RESET_CONSTANTS[mode]["c_train"]
    c_final = RESET_CONSTANTS[mode]["c_final"]
    X = np.asarray(X, dtype=np.float64)
    if X.ndim == 1:
        X = X[:, None]
    is_decoy = np.asarray(is_decoy, dtype=bool)
    n = X.shape[0]
    rng = make_rng(seed)
    is_train = is_decoy & (rng.random(n) < 0.5)
    pseudo = ~is_train
    if standardize:
        X = _standardize(X)
    learner = LDALearner() if learner is None else learner
    if init_score is None:
        best, score = -1, X[:, 0].copy()
        for j in range(X.shape[1]):
            for sgn in (1.0, -1.0):
                cand = sgn * X[:, j]
                npos_j = int((pseudo & (tdc_qvalues(cand, is_train, factor=c_train) <= alpha_train)).sum())
                if npos_j > best:
                    best, score = npos_j, cand.copy()
    else:
        score = np.asarray(init_score, dtype=np.float64).copy()
    npos = []
    for _ in range(n_iter):
        # training competition: pseudo-targets vs training decoys (labels of
        # estimating decoys are never consulted)
        qt = tdc_qvalues(score, is_train, factor=c_train)
        pos = pseudo & (qt <= alpha_train)
        npos.append(int(pos.sum()))
        if pos.sum() < min_positives or is_train.sum() < 2:
            break
        learner.fit(X[pos], X[is_train])
        score = np.asarray(learner.decision_function(X), dtype=np.float64)
    keep = ~is_train
    q = np.full(n, np.nan)
    q[keep] = tdc_qvalues(score[keep], is_decoy[keep], factor=c_final)
    return ResetResult(score=score, is_training_decoy=is_train, q=q, mode=mode, c_train=c_train,
                       c_final=c_final, n_positives=npos)


def two_decoy_competition(target_score, decoy1_score, decoy2_score, seed=None):
    """Competition of each target with its two decoys (independent shuffles).

    Returns (winner index into {0: target, 1: decoy1, 2: decoy2}, winning score).
    Ties broken at random.
    """
    t = np.asarray(target_score, dtype=np.float64)
    S = np.column_stack([t, np.asarray(decoy1_score, float), np.asarray(decoy2_score, float)])
    n = len(t)
    grp = np.repeat(np.arange(n), 3)
    w = compete(grp, S.ravel(), seed=seed)
    which = (w % 3).astype(np.int64)
    return which, S.ravel()[w]
