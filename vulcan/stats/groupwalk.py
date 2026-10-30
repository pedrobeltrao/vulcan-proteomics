"""Group-walk (Freestone, Short, Noble, Keich 2022, Bioinformatics 38:i82) and the
pre-registered stratum merge rules (DESIGN.md §4 L4).

Algorithm (follows the reference R implementation, CRAN ``groupwalk::group_walk``):

* Within each group, hypotheses (competition winners) are sorted by winning score,
  ascending. The walk starts with *all* hypotheses in the candidate set and peels
  one hypothesis at a time from the bottom of one group (the "front").
* Before each peel the estimate ``FDR = factor * (D + 1) / max(T, 1)`` is computed
  over all not-yet-peeled hypotheses of all groups; the running minimum of this
  estimate is the q-value assigned to the hypothesis peeled at this step.
  Reporting targets with ``q <= alpha`` equals stopping at the first front with
  ``FDR <= alpha``.
* Group choice: while any non-exhausted group has had fewer than ``K`` peels,
  peel from a group with the fewest peels (phase 1). Afterwards peel from the group
  with the most decoys among its last ``K`` peeled hypotheses (phase 2). Ties
  between groups: seeded random (default, as in the R package) or lowest group
  index (``tie_break="first"``).
* Within-group ties in score: targets are peeled before decoys, so a target is
  never accepted without the decoys tied with it ("decoys first").

Validity needs (i) a group partition fixed independently of the target/decoy
labels (label-blind counts or dev data only), and (ii) the TDC null coin flips
(fair, independent) within each group. Group membership must be identical within a
target/decoy pair (see :func:`check_pair_invariant`).

With ``factor=2`` this is the "RESET + Group-walk" adaptation in DESIGN §6.4
(c = 2/3); it is exploratory and reported as unproven.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ._util import as_codes, make_rng

__all__ = [
    "group_walk",
    "group_walk_batch",
    "merge_adjacent",
    "merge_rare_to_other",
    "pool_chimeric_rank",
    "label_blind_pair_counts",
    "check_pair_invariant",
    "MIN_STRATUM_SIZE",
]

MIN_STRATUM_SIZE = 200  # pair competitions; a power rule, not a validity rule (DESIGN v0.3)


# --------------------------------------------------------------------------------------
# Group-walk
# --------------------------------------------------------------------------------------

def _group_orders(score, is_decoy, gcode, rng):
    """Per-group peel order (ascending score, targets before decoys at ties, random)."""
    n = score.shape[0]
    tb = rng.random(n)
    order = np.lexsort((tb, is_decoy.astype(np.int8), score, gcode))
    gs = gcode[order]
    bounds = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1], True])
    return order, bounds


def group_walk(score, is_decoy, group, K: int = 40, factor: float = 1.0, seed=None,
               tie_break: str = "random", return_frontier: bool = False):
    """Group-walk q-values for competition winners.

    Parameters
    ----------
    score : winning scores (one per competition), higher is better.
    is_decoy : winner label.
    group : stratum key per competition (fixed before estimation).
    K : window size (default 40, as in the paper/R package).
    factor : multiplies the decoy estimate (1 for TDC; 2 for the RESET adaptation).
    seed : seed for within-group score ties and between-group ties.
    tie_break : "random" or "first" (lowest group code wins ties).

    Returns
    -------
    q : float array (same order as input; decoys get values too).
    If ``return_frontier``: also the (n_steps+1, n_groups) array of peel counts.
    """
    score = np.asarray(score, dtype=np.float64)
    score = np.where(np.isnan(score), -np.inf, score)
    is_decoy = np.asarray(is_decoy, dtype=bool)
    n = score.shape[0]
    if n == 0:
        return (np.zeros(0), np.zeros((1, 0), dtype=np.int64)) if return_frontier else np.zeros(0)
    gcode = as_codes(group)
    _, gcode = np.unique(gcode, return_inverse=True)
    G = int(gcode.max()) + 1
    rng = make_rng(seed)
    order, bounds = _group_orders(score, is_decoy, gcode, rng)
    lists = [order[bounds[g]:bounds[g + 1]].tolist() for g in range(G)]
    labs = [is_decoy[np.asarray(lst, dtype=np.int64)].tolist() for lst in lists]
    totals = [len(lst) for lst in lists]
    peeled = [0] * G
    window = [0] * G  # decoys among the last K peeled of each group
    D = int(is_decoy.sum())
    T = n - D
    qrun = 1.0
    q = np.ones(n)
    active = [g for g in range(G) if totals[g] > 0]
    frontier = [list(peeled)] if return_frontier else None
    use_random = tie_break == "random"
    while active:
        fdr = factor * (D + 1) / max(T, 1)
        if fdr < qrun:
            qrun = fdr
        # choose group
        if any(peeled[g] < K for g in active):
            best = min(peeled[g] for g in active)
            cands = [g for g in active if peeled[g] == best]
        else:
            best = max(window[g] for g in active)
            cands = [g for g in active if window[g] == best]
        g = cands[int(rng.integers(len(cands)))] if (use_random and len(cands) > 1) else cands[0]
        k = peeled[g]
        idx = lists[g][k]
        lab = labs[g][k]
        q[idx] = qrun
        if lab:
            D -= 1
        else:
            T -= 1
        window[g] += lab
        if k >= K:
            window[g] -= labs[g][k - K]
        peeled[g] = k + 1
        if peeled[g] >= totals[g]:
            active.remove(g)
        if return_frontier:
            frontier.append(list(peeled))
    if return_frontier:
        return q, np.asarray(frontier, dtype=np.int64)
    return q


def group_walk_batch(labels, null, sizes, alphas: Sequence[float], K: int = 40, factor: float = 1.0):
    """Lock-step Group-walk over many Monte-Carlo replicates (simulation helper).

    Parameters
    ----------
    labels : bool array (R, G, L): decoy flag of each group's hypotheses in
        *ascending* score order (peel order), padded beyond ``sizes``.
    null : bool array (R, G, L): True for null *targets* (false discoveries).
    sizes : int array (G,) group sizes (same for every replicate).
    alphas : thresholds.

    Returns dict alpha -> (n_target_discoveries (R,), n_false (R,)).
    Ties between groups go to the lowest group index (``tie_break="first"``),
    so results equal :func:`group_walk` with that option.
    """
    labels = np.asarray(labels, dtype=bool)
    null = np.asarray(null, dtype=bool)
    R, G, L = labels.shape
    sizes = np.asarray(sizes, dtype=np.int64)
    valid = np.arange(L)[None, :] < sizes[:, None]  # (G, L)
    lab = labels & valid[None]
    nul = null & valid[None] & ~labels
    D = lab.sum(axis=(1, 2)).astype(np.int64)
    T = (valid.sum() - D).astype(np.int64)
    F = nul.sum(axis=(1, 2)).astype(np.int64)
    peeled = np.zeros((R, G), dtype=np.int64)
    window = np.zeros((R, G), dtype=np.int64)
    alphas = list(alphas)
    res = {a: [np.zeros(R, dtype=np.int64), np.zeros(R, dtype=np.int64), np.zeros(R, dtype=bool)]
           for a in alphas}
    n_total = int(sizes.sum())
    rr = np.arange(R)
    BIG = np.int64(1 << 40)
    for _ in range(n_total + 1):
        fdr = factor * (D + 1) / np.maximum(T, 1)
        all_done = True
        for a in alphas:
            T_a, F_a, done = res[a]
            hit = (~done) & (fdr <= a)
            T_a[hit] = T[hit]
            F_a[hit] = F[hit]
            done |= hit
            all_done &= bool(done.all())
        if all_done:
            break
        act = peeled < sizes[None, :]
        if not act.any():
            break
        # phase 1 rows: some active group with peeled < K
        p1 = (act & (peeled < K)).any(axis=1)
        key1 = np.where(act, peeled, BIG)
        g1 = key1.argmin(axis=1)
        key2 = np.where(act, window, -1)
        g2 = key2.argmax(axis=1)
        g = np.where(p1, g1, g2)
        k = peeled[rr, g]
        live = act.any(axis=1)
        kk = np.minimum(k, L - 1)
        lb = lab[rr, g, kk] & live
        nl = nul[rr, g, kk] & live
        D -= lb
        T -= (~lb) & live
        F -= nl
        window[rr, g] += lb
        drop = (k >= K) & live
        kd = np.clip(k - K, 0, L - 1)
        window[rr, g] -= np.where(drop, lab[rr, g, kd], False)
        peeled[rr, g] += live
    out = {}
    for a in alphas:
        T_a, F_a, done = res[a]
        out[a] = (np.where(done, T_a, 0), np.where(done, F_a, 0))
    return out


# --------------------------------------------------------------------------------------
# Strata construction and pre-registered merge rules (DESIGN v0.3 §4 L4)
# --------------------------------------------------------------------------------------

def label_blind_pair_counts(pair_id, stratum) -> dict:
    """Number of distinct pair competitions per stratum.

    Counts pairs, not winners or targets, so the partition is label-blind: it is the
    same whichever side of each pair wins.
    """
    pair_id = np.asarray(pair_id)
    stratum = np.asarray(stratum)
    check_pair_invariant(pair_id, stratum)
    keys, first = np.unique(as_codes(pair_id), return_index=True)
    st = stratum[first]
    u, c = np.unique(st, return_counts=True)
    return dict(zip(u.tolist(), c.tolist()))


def check_pair_invariant(pair_id, stratum) -> None:
    """Raise ValueError if a target/decoy pair spans more than one stratum."""
    p = as_codes(pair_id)
    s = as_codes(stratum)
    order = np.lexsort((s, p))
    ps, ss = p[order], s[order]
    same_pair = ps[1:] == ps[:-1]
    if np.any(same_pair & (ss[1:] != ss[:-1])):
        raise ValueError("stratum membership differs within a target/decoy pair")


def merge_adjacent(order: Sequence, counts: dict, min_size: int = MIN_STRATUM_SIZE) -> dict:
    """Merge ordered strata (e.g. mass-offset bins sorted by mass) until all have
    >= ``min_size`` pair competitions.

    Deterministic rule: repeatedly take the smallest block (ties: the earliest in
    ``order``) and merge it with its smaller adjacent neighbour (ties: the preceding
    one). Stops when every block has >= min_size or one block remains (pooled TDC).
    Strata absent from ``counts`` count as 0.

    Returns mapping stratum -> merged label (the tuple of member strata).
    """
    blocks = [[s] for s in order]
    sizes = [int(counts.get(s, 0)) for s in order]
    while len(blocks) > 1 and min(sizes) < min_size:
        i = int(np.argmin(sizes))
        if i == 0:
            j = 1
        elif i == len(blocks) - 1:
            j = i - 1
        else:
            j = i - 1 if sizes[i - 1] <= sizes[i + 1] else i + 1
        lo, hi = min(i, j), max(i, j)
        blocks[lo] = blocks[lo] + blocks[hi]
        sizes[lo] += sizes[hi]
        del blocks[hi], sizes[hi]
    mapping = {}
    for b in blocks:
        lab = tuple(b)
        for s in b:
            mapping[s] = lab
    return mapping


def merge_rare_to_other(counts: dict, min_size: int = MIN_STRATUM_SIZE, other="other") -> dict:
    """Unordered classes (PTM composition): classes with < min_size pair
    competitions go to ``other``. ``other`` is kept even if it stays below
    min_size (min size is a power rule, not a validity rule)."""
    return {k: (k if v >= min_size else other) for k, v in counts.items()}


def pool_chimeric_rank(rank, pool_from: int = 3) -> np.ndarray:
    """Label-free chimeric rank with ranks >= ``pool_from`` pooled into one stratum."""
    r = np.asarray(rank)
    return np.minimum(r, pool_from)
