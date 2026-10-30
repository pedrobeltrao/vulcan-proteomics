"""Faster, output-identical chimeric competition for production (E7b, results/2026-10-28_E7b_fair_speed).

Same steps and estimators as vulcan.chimeric.compete.chimeric_competition (stage-1 group winners, CONGA neighbour
filter, paired double competition, rank strata merged label-blind, per-stratum TDC q_sep).  What changes:

* ion-bin sets (IonSets.get) are computed in numba with numpy's exact arithmetic (sequential cumsum, numpy's pairwise
  summation for the total mass, round-half-even), as sorted unique int64 arrays instead of Python frozensets, and only
  for spectra with >= 2 group winners (a lone winner is always kept, its ion set is never looked at);
* the neighbour filter runs in numba over those arrays (|A & B| by sorted merge; same visiting order: spectra by key,
  winners by descending score with numpy's lexsort tie order);
* `group_walk=False` skips the Group-walk q-values (q_gw), which q_sep does not depend on.
`chimeric_competition(..., fast=False)` (vulcan.chimeric.compete) stays the reference.
"""
from __future__ import annotations

import numba as nb
import numpy as np

from vulcan.chimeric.compete import BIN_DA, FRAC, H2O, MIN_ION_INDEX, PROTON, IonSets
from vulcan.stats import groupwalk as GW
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.tdc import compete, double_competition, tdc_qvalues


@nb.njit(cache=True)
def _np_pairwise_sum(a, n):
    """numpy's pairwise_sum for float64 (numpy/_core/src/umath/loops_utils.h.src) for n <= PW_BLOCKSIZE (128):
    sequential below 8, else 8 interleaved partial sums combined pairwise, then the remainder.  Peptides are <= 50
    residues; longer inputs raise (numpy would split recursively; not replicated, non-recursive on purpose because
    numba's on-disk cache mishandles recursive functions)."""
    if n < 8:
        res = 0.0
        for i in range(n):
            res += a[i]
        return res
    if n > 128:
        raise ValueError("peptide longer than 128 residues")
    r0 = a[0]
    r1 = a[1]
    r2 = a[2]
    r3 = a[3]
    r4 = a[4]
    r5 = a[5]
    r6 = a[6]
    r7 = a[7]
    i = 8
    while i < n - (n % 8):
        r0 += a[i]
        r1 += a[i + 1]
        r2 += a[i + 2]
        r3 += a[i + 3]
        r4 += a[i + 4]
        r5 += a[i + 5]
        r6 += a[i + 6]
        r7 += a[i + 7]
        i += 8
    res = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
    while i < n:
        res += a[i]
        i += 1
    return res


@nb.njit(cache=True)
def _ion_bins(rows, seq_off, length, residues, mask, aa, var, bin_da, kmin, proton, h2o):
    """CSR of sorted unique ion bins per row (IonSets.get semantics)."""
    n_r = rows.shape[0]
    cnt = np.zeros(n_r + 1, np.int64)
    for t in range(n_r):
        n = np.int64(length[rows[t]])
        nb_ = max(0, n - 1 - (kmin - 1)) + max(0, n - kmin)
        cnt[t + 1] = cnt[t] + nb_
    out = np.empty(cnt[n_r], np.int64)
    uoff = np.zeros(n_r + 1, np.int64)
    w = 0
    for t in range(n_r):
        row = rows[t]
        o = np.int64(seq_off[row])
        n = np.int64(length[row])
        mk = np.uint64(mask[row])
        m = np.empty(n)
        for k in range(n):
            r = np.int64(residues[o + k])
            bit = (mk >> np.uint64(k)) & np.uint64(1)
            m[k] = aa[r] + (var[r] if bit else 0.0)
        tot = _np_pairwise_sum(m, n) + h2o
        pre = np.empty(max(n - 1, 0))
        acc = 0.0
        for k in range(n - 1):
            acc += m[k]
            pre[k] = acc
        nb_ = cnt[t + 1] - cnt[t]
        ions = np.empty(nb_, np.int64)
        j = 0
        for k in range(kmin - 1, n - 1):                    # b_i, i >= kmin
            ions[j] = np.int64(np.rint((pre[k] + proton) / bin_da))
            j += 1
        for k in range(max(0, n - kmin)):                   # y
            ions[j] = np.int64(np.rint(((tot - pre[k]) + proton) / bin_da))
            j += 1
        ions.sort()
        start = w
        for k in range(nb_):
            if k == 0 or ions[k] != ions[k - 1]:
                out[w] = ions[k]
                w += 1
        uoff[t + 1] = uoff[t] + (w - start)
    return uoff, out[:w]


@nb.njit(cache=True)
def _inter(a, a0, a1, b, b0, b1):
    i, j, c = a0, b0, 0
    while i < a1 and j < b1:
        if a[i] == b[j]:
            c += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return c


@nb.njit(cache=True)
def _filter(gstart, order, slot, uoff, bins, frac):
    """Greedy neighbour filter over spectrum groups [gstart[g], gstart[g+1]) of `order` (descending score)."""
    n = order.shape[0]
    keep = np.zeros(n, np.bool_)
    kept = np.empty(64, np.int64)
    for g in range(gstart.shape[0] - 1):
        a0, a1 = gstart[g], gstart[g + 1]
        if a1 - a0 == 1:
            keep[order[a0]] = True
            continue
        if a1 - a0 > kept.shape[0]:
            kept = np.empty(a1 - a0, np.int64)
        nk = 0
        for p in range(a0, a1):
            t = order[p]
            s = slot[t]
            la = uoff[s + 1] - uoff[s]
            ok = True
            for u in range(nk):
                su = slot[kept[u]]
                lb = uoff[su + 1] - uoff[su]
                c = _inter(bins, uoff[s], uoff[s + 1], bins, uoff[su], uoff[su + 1])
                if c > frac * min(la, lb):
                    ok = False
                    break
            if ok:
                kept[nk] = t
                nk += 1
                keep[t] = True
    return keep


def neighbour_filter_fast(spec_key, score, rows, ion_sets: IonSets, frac=FRAC):
    """Same keep mask as vulcan.ptm.conga.neighbour_filter(spec_key, score, ion_sets.sets(rows), frac)."""
    spec_key = np.asarray(spec_key)
    score = np.asarray(score, np.float64)
    rows = np.asarray(rows, np.int64)
    n = len(score)
    if n == 0:
        return np.zeros(0, bool)
    order = np.lexsort((-score, spec_key))
    ks = spec_key[order]
    brk = np.flatnonzero(ks[1:] != ks[:-1]) + 1
    gstart = np.concatenate([[0], brk, [n]]).astype(np.int64)
    gsize = np.diff(gstart)
    need = np.zeros(n, bool)
    need[order[np.repeat(gsize, gsize) >= 2]] = True        # winners of spectra with >= 2 winners
    nidx = np.flatnonzero(need)
    slot = np.full(n, -1, np.int64)
    slot[nidx] = np.arange(len(nidx))
    uoff, bins = _ion_bins(rows[nidx], *_l1_arrays(ion_sets), ion_sets.aa, ion_sets.var, float(ion_sets.bin),
                           int(ion_sets.min_ion), PROTON, H2O)
    return _filter(gstart, order.astype(np.int64), slot, uoff, bins, float(frac))


def _l1_arrays(ion_sets: IonSets):
    """(seq_off, length, residues, mask as uint64) of the IonSets' index, converted once and memoized."""
    a = getattr(ion_sets, "_fast_arrays", None)
    if a is None:
        l1 = ion_sets.l1
        mask = np.asarray(l1.mask)
        mask = mask.view(np.uint64) if mask.dtype == np.int64 else mask.astype(np.uint64)
        a = (np.asarray(l1.seq_off), np.asarray(l1.length), np.asarray(l1.residues), mask)
        ion_sets._fast_arrays = a
    return a


def ion_bins_of(rows, ion_sets: IonSets):
    """(off, bins) for explicit rows (used by the equality tests)."""
    return _ion_bins(np.asarray(rows, np.int64), *_l1_arrays(ion_sets), ion_sets.aa, ion_sets.var,
                     float(ion_sets.bin), int(ion_sets.min_ion), PROTON, H2O)


def chimeric_competition_fast(spec_key, group, chim_rank, unit, pair, score, is_decoy, rows, ion_sets: IonSets,
                              seed: int = 11, K: int = 40, pool_from: int = 3, min_size: int = GW.MIN_STRATUM_SIZE,
                              neighbour: bool = True, group_walk: bool = False) -> dict:
    """vulcan.chimeric.compete.chimeric_competition with the fast neighbour filter; q_gw only if group_walk."""
    spec_key = np.asarray(spec_key)
    unit = np.asarray(unit)
    score = np.asarray(score, np.float64)
    tb = hash_tiebreak(unit, seed)
    assert spec_key.dtype.kind not in "UO", "integer spectrum keys only"
    qkey = spec_key.astype(np.int64) * 64 + np.asarray(group, np.int64)
    w1 = compete(qkey, score, tiebreak=tb)
    if neighbour:
        keep = neighbour_filter_fast(spec_key[w1], score[w1], np.asarray(rows)[w1], ion_sets)
        w = w1[keep]
    else:
        keep = np.ones(len(w1), bool)
        w = w1
    res = double_competition(qkey[w], unit[w], np.asarray(pair)[w], score[w], np.asarray(is_decoy)[w], seed=seed,
                             tiebreak=tb[w])
    win = w[res.row]
    stratum = GW.pool_chimeric_rank(np.asarray(chim_rank)[win].astype(np.int64), pool_from)
    counts = GW.label_blind_pair_counts(res.pair, stratum)
    mp = GW.merge_adjacent(list(range(1, pool_from + 1)), counts, min_size)
    merged = np.array(["_".join(map(str, mp[int(s)])) for s in stratum])
    q_gw = GW.group_walk(res.score, res.is_decoy, merged, K=K, seed=seed) if group_walk else None
    q_sep = tdc_qvalues(res.score, res.is_decoy, block=merged)
    return {"res": res, "win": win, "survivors": w, "stratum": stratum, "merged": merged, "q_gw": q_gw,
            "q_pooled": res.q, "q_sep": q_sep, "pair_counts": {int(k): int(v) for k, v in counts.items()},
            "merge": {int(k): list(map(int, v)) for k, v in mp.items()},
            "n_group_winners": int(len(w1)), "n_after_neighbour_filter": int(len(w)),
            "n_dropped_by_neighbour_filter": int((~keep).sum())}
