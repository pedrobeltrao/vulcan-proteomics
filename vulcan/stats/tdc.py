"""Target-decoy competition (TDC), paired double competition, multi-decoy TDC.

Conventions (DESIGN.md §4 L4, §6):

* Scores: higher is better. NaN scores are treated as -inf.
* Competition (spectrum level, peptide level, pair level): the highest score in a
  group wins; exact ties are broken uniformly at random with a seeded generator,
  independently of the label.
* FDR estimate at a score threshold ``s`` (all items with score >= s):
  ``FDR(s) = factor * (D(s) + 1) / max(T(s), 1)`` ("+1" inside, multiplied by
  ``factor``; factor = 1 for plain TDC). Items with identical scores share one
  threshold, so every decoy tied with a target at score s is counted when that
  target is accepted ("decoys first at equal scores").
* q-values are the monotone version ``q_i = min_{s <= s_i} FDR(s)``, capped at 1.
  Reporting ``q <= alpha`` is identical to the SeqStep+/TDC threshold (the
  largest accepted list with ``FDR <= alpha``).

References: Elias & Gygi 2007; He et al. 2015 (arXiv:1501.00537); Barber &
Candès 2015 (knockoff+); Lin, Short, Noble, Keich 2022 (JPR 21:2023, double
competition); Emery, Hasam, Noble, Keich 2020 (RECOMB, multiple decoys).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._util import as_codes, make_rng, segment_cumsum, segment_reverse_cummin, segment_starts

__all__ = [
    "compete",
    "psm_competition",
    "tdc_qvalues",
    "tdc_fdr_curve",
    "count_discoveries",
    "double_competition",
    "DoubleCompetitionResult",
    "multi_decoy_competition",
    "multi_decoy_qvalues",
]


def _clean_scores(score) -> np.ndarray:
    s = np.asarray(score, dtype=np.float64)
    if np.isnan(s).any():
        s = np.where(np.isnan(s), -np.inf, s)
    return s


def compete(group, score, seed=None, tiebreak=None) -> np.ndarray:
    """Winner (row index) of each group: the max score; ties broken at random.

    Parameters
    ----------
    group : array of keys (any dtype). Rows sharing a key compete.
    score : float array, higher is better.
    seed  : int / Generator for tie-breaking (ignored if ``tiebreak`` is given).
    tiebreak : optional float/int array; among tied scores the *smallest*
        tiebreak wins (use e.g. hash(sequence, seed) for storage-order-free,
        reproducible tie-breaking, DESIGN §4 L1).

    Returns
    -------
    int64 array of winning row indices, one per distinct group, ordered by group key.
    """
    s = _clean_scores(score)
    n = s.shape[0]
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    g = as_codes(group)
    if tiebreak is None:
        tb = make_rng(seed).random(n)
    else:
        tb = np.asarray(tiebreak)
    order = np.lexsort((tb, -s, g))
    first = segment_starts(g[order])
    return order[first].astype(np.int64)


def psm_competition(spectrum_id, score, seed=None, tiebreak=None) -> np.ndarray:
    """PSM-level TDC: all hypotheses (charges, isotope errors, co-isolated
    features) and all candidates of a spectrum compete; one winner per spectrum."""
    return compete(spectrum_id, score, seed=seed, tiebreak=tiebreak)


def tdc_fdr_curve(score, is_decoy, factor: float = 1.0, block=None, plus_one: bool = True):
    """Raw (non-monotone) FDR estimate at each item's own score threshold.

    Returns ``(fdr, n_targets_at_threshold, n_decoys_at_threshold)`` in the
    original order. ``block`` (optional) computes independent curves per block
    (used to vectorise Monte-Carlo replicates).
    """
    s = _clean_scores(score)
    d = np.asarray(is_decoy, dtype=bool)
    n = s.shape[0]
    b = np.zeros(n, dtype=np.int64) if block is None else as_codes(block)
    order = np.lexsort((-s, b))
    ss, ds, bs = s[order], d[order], b[order]
    bstart = segment_starts(bs)
    cD = segment_cumsum(ds.astype(np.int64), bstart)
    cT = segment_cumsum((~ds).astype(np.int64), bstart)
    # counts at the END of each tie run (all items with score >= s)
    tstart = bstart.copy()
    if n > 1:
        tstart[1:] |= ss[1:] != ss[:-1]
    tid = np.cumsum(tstart) - 1
    ends = np.concatenate([np.flatnonzero(tstart)[1:] - 1, [n - 1]])
    Dt = cD[ends[tid]]
    Tt = cT[ends[tid]]
    fdr = factor * (Dt + (1.0 if plus_one else 0.0)) / np.maximum(Tt, 1)
    out_f = np.empty(n)
    out_T = np.empty(n, dtype=np.int64)
    out_D = np.empty(n, dtype=np.int64)
    out_f[order] = fdr
    out_T[order] = Tt
    out_D[order] = Dt
    return out_f, out_T, out_D


def tdc_qvalues(score, is_decoy, factor: float = 1.0, block=None, plus_one: bool = True) -> np.ndarray:
    """Monotone TDC q-values ``q_i = min_{s<=s_i} factor*(D(s)+1)/max(T(s),1)``.

    Decoys receive q-values too (useful for diagnostics); report targets only.
    ``block`` computes independent q-values per block.
    """
    s = _clean_scores(score)
    n = s.shape[0]
    if n == 0:
        return np.zeros(0)
    fdr, _, _ = tdc_fdr_curve(s, is_decoy, factor=factor, block=block, plus_one=plus_one)
    b = np.zeros(n, dtype=np.int64) if block is None else as_codes(block)
    order = np.lexsort((-s, b))
    q_sorted = segment_reverse_cummin(fdr[order], segment_starts(b[order]), cap=1.0)
    q = np.empty(n)
    q[order] = q_sorted
    return q


def count_discoveries(q, is_decoy, alpha: float) -> int:
    """Number of target discoveries at q <= alpha."""
    q = np.asarray(q)
    return int(np.count_nonzero((q <= alpha) & ~np.asarray(is_decoy, dtype=bool)))


@dataclass
class DoubleCompetitionResult:
    """Output of :func:`double_competition` (one row per competed pair unit)."""

    row: np.ndarray        # index into the input PSM table of the winning PSM
    pair: np.ndarray       # pair key of the unit
    score: np.ndarray
    is_decoy: np.ndarray
    q: np.ndarray
    psm_winners: np.ndarray  # rows that won spectrum-level competition

    def accepted(self, alpha: float) -> np.ndarray:
        return self.row[(self.q <= alpha) & ~self.is_decoy]


def double_competition(spectrum_id, unit_id, pair_id, score, is_decoy, seed=None,
                       factor: float = 1.0, tiebreak=None) -> DoubleCompetitionResult:
    """Paired double competition (Lin et al. 2022 "PSM-and-peptide").

    1. All hypotheses/candidates of a spectrum compete (PSM-level TDC).
    2. Among PSM winners, the best PSM per unit (peptide for DDA, precursor for
       DIA; e.g. ``unit_id`` = peptide or (peptide, charge)).
    3. Each target unit competes with its *own* paired decoy unit: units sharing
       ``pair_id`` compete (for precursor level pass pair key = (pair, charge)).
    4. (D+1)/T monotone q-values over the pair winners.

    Each stage uses an independent seeded tie-breaking stream, unless a
    storage-order-free ``tiebreak`` key per row is given (recommended:
    ``_util.hash_tiebreak(sequence, seed)``), which is then used at every stage.
    ``pair_id`` must be identical for a target unit and its decoy unit.
    """
    score = _clean_scores(score)
    is_decoy = np.asarray(is_decoy, dtype=bool)
    unit_id = np.asarray(unit_id)
    pair_id = np.asarray(pair_id)
    # unit ids must not collide between target and decoy side
    tb = None if tiebreak is None else np.asarray(tiebreak)
    sd = (lambda k: make_rng(seed, k) if seed is not None else None)
    w1 = psm_competition(spectrum_id, score, seed=sd(1), tiebreak=tb)
    unit_key = as_codes(unit_id[w1]) * 2 + is_decoy[w1].astype(np.int64)
    w2 = w1[compete(unit_key, score[w1], seed=sd(2), tiebreak=None if tb is None else tb[w1])]
    w3 = w2[compete(pair_id[w2], score[w2], seed=sd(3), tiebreak=None if tb is None else tb[w2])]
    q = tdc_qvalues(score[w3], is_decoy[w3], factor=factor)
    return DoubleCompetitionResult(row=w3, pair=pair_id[w3], score=score[w3], is_decoy=is_decoy[w3], q=q,
                                   psm_winners=w1)


# --------------------------------------------------------------------------------------
# Multiple decoys (Emery, Hasam, Noble, Keich 2020; arXiv:1907.01458)
# --------------------------------------------------------------------------------------

def multi_decoy_competition(target_score, decoy_scores, k_target: int | None = None,
                            k_decoy: int | None = None, seed=None):
    """Competition of each target with its d exchangeable decoys.

    For hypothesis i with original score Z_i0 and decoy scores Z_i1..Z_id the
    winning score is ``W_i = max_j Z_ij`` and the label is derived from the rank
    ``r_i`` (1 = best) of the original among the d+1 scores (random ties):

    * target-win if ``r_i <= k_target``
    * decoy-win  if ``r_i >  d + 1 - k_decoy``
    * ignored otherwise.

    Under the null (exchangeable original and decoys) ``P(target)=k_t/(d+1)``,
    ``P(decoy)=k_d/(d+1)`` and the label is independent of W, so
    ``FDR = (k_t/k_d) * (D+1)/T`` controls FDR (Emery et al. 2020, Thm 1;
    c = k_t/(d+1), lambda = 1 - k_d/(d+1)). d = 1, k_t = k_d = 1 is plain TDC.
    Default: ``k_t = k_d = floor((d+1)/2)``.

    Returns ``(winning_score, label)`` with label +1 target, -1 decoy, 0 ignored.
    """
    t = _clean_scores(target_score)
    D = np.atleast_2d(_clean_scores(decoy_scores))
    if D.shape[0] != t.shape[0]:
        D = D.T
    d = D.shape[1]
    kt = (d + 1) // 2 if k_target is None else int(k_target)
    kd = (d + 1) // 2 if k_decoy is None else int(k_decoy)
    if kt + kd > d + 1 or kt < 1 or kd < 1:
        raise ValueError("need 1 <= k_target, k_decoy and k_target + k_decoy <= d + 1")
    rng = make_rng(seed)
    allz = np.column_stack([t, D])
    jitter = rng.random(allz.shape)
    # rank of the original: 1 + #decoys strictly better (+ random share of exact ties)
    better = (D > t[:, None]) | ((D == t[:, None]) & (jitter[:, 1:] < jitter[:, :1]))
    rank = 1 + better.sum(axis=1)
    label = np.zeros(t.shape[0], dtype=np.int8)
    label[rank <= kt] = 1
    label[rank > d + 1 - kd] = -1
    return allz.max(axis=1), label


def multi_decoy_qvalues(winning_score, label, k_target: int, k_decoy: int) -> np.ndarray:
    """q-values for :func:`multi_decoy_competition` output (ignored rows get q=1)."""
    label = np.asarray(label)
    keep = label != 0
    q = np.ones(label.shape[0])
    q[keep] = tdc_qvalues(np.asarray(winning_score)[keep], label[keep] < 0, factor=k_target / k_decoy)
    return q
