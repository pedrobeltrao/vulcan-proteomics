"""Isotope-pattern grouping of hills into features (averagine, charge 1..6), biosaur2-style.

For every hill taken as a monoisotopic candidate and every charge z, isotopes k = 1, 2, ... are looked up at
mz + k*NEUTRON/z within `itol_ppm` among hills whose scan range overlaps the mono hill. The best match is the hill
whose intensity trace has the highest cosine with the mono trace. The trace cosine uses the full-hill norms
(`full_norm=True`, as biosaur2: scans outside the overlap count against the match) or only the overlap range.
The chain stops at the first missing isotope or a trace cosine below `min_trace_cos`.

The chain is then truncated to the longest prefix whose apex intensities have an averagine cosine >= `min_avg_cos`
and explain >= `min_explained` of the averagine envelope's total abundance (biosaur2's "averagineExplained": a
2-peak pattern of a heavy peptide, whose envelope is mostly in higher isotopes, is not accepted).

Charge harmonics (biosaur2 charge ban): a (mono, z) candidate is dropped if the same mono hill has a candidate at a
multiple charge k*z with at least as many isotopes (a z=2 envelope also looks like z=1 on every second peak).

Candidates (>= 2 isotopes) are accepted greedily, ordered by (number of isotopes, averagine cosine,
mono apex intensity) descending; an accepted candidate claims its hills. A later candidate whose mono hill
is taken is dropped; one whose isotope k is taken is truncated before k and re-scored.

Hot loops are numba (parallel over hills, visited in m/z order for memory locality); the greedy pass is sequential.

Fast exact path (E7b, results/2026-10-28_E7b_fair_speed; default when every hill covers consecutive MS1 scans, which
`detect_hills(gap=0)` guarantees): the same floating-point operations in the same order, organised cheaper --
  * full-norm trace cosines take each hill's sum of squares from a per-hill table filled by the same sequential loop
    (bitwise equal to the per-call loop), and the overlap is aligned by index arithmetic instead of a merge walk
    (contiguous hills: every overlap scan is present in both, so the merge visits the same pairs in the same order);
  * isotope lookups visit only hills in the mono hill's MS1-scan blocks (`_build_rt_bins`, 32 scans per block,
    m/z-sorted within a block) instead of every hill in the m/z window (D0: 98 window hills per lookup, 0.16% of them
    overlap in scans); the same candidates are scored and the "best cosine, first in m/z order on ties" rule is kept
    explicitly (global m/z rank as tie-break);
  * the first-isotope pre-screen stops at the first isotope that passes (it only needs existence).
`group_isotopes(..., exact_reference=True)` runs the original kernels (reference for the equality checks).
"""
from __future__ import annotations

import numba as nb
import numpy as np

from .averagine import STEP, get_table


@nb.njit(cache=True)
def _trace_cos(off, p_scan, p_int, a, b, full_norm):
    """Cosine of the intensity traces of hills a and b (missing scans = 0). Returns (cos, overlap length)."""
    s0 = max(p_scan[off[a]], p_scan[off[b]])
    s1 = min(p_scan[off[a + 1] - 1], p_scan[off[b + 1] - 1])
    if s1 < s0:
        return -1.0, 0
    i, j = off[a], off[b]
    ea, eb = off[a + 1], off[b + 1]
    saa = 0.0
    sbb = 0.0
    sab = 0.0
    if full_norm:
        for q in range(i, ea):
            saa += p_int[q] * p_int[q]
        for q in range(j, eb):
            sbb += p_int[q] * p_int[q]
    while i < ea and p_scan[i] < s0:
        i += 1
    while j < eb and p_scan[j] < s0:
        j += 1
    while (i < ea and p_scan[i] <= s1) or (j < eb and p_scan[j] <= s1):
        si = p_scan[i] if (i < ea and p_scan[i] <= s1) else 1 << 30
        sj = p_scan[j] if (j < eb and p_scan[j] <= s1) else 1 << 30
        if si == sj:
            x = p_int[i]
            y = p_int[j]
            sab += x * y
            if not full_norm:
                saa += x * x
                sbb += y * y
            i += 1
            j += 1
        elif si < sj:
            if not full_norm:
                saa += p_int[i] * p_int[i]
            i += 1
        else:
            if not full_norm:
                sbb += p_int[j] * p_int[j]
            j += 1
    if saa <= 0 or sbb <= 0:
        return 0.0, s1 - s0 + 1
    return sab / np.sqrt(saa * sbb), s1 - s0 + 1


@nb.njit(cache=True)
def _find_iso(mono, target, tol_ppm, smz, sidx, off, p_scan, p_int, min_trace_cos, min_overlap_for_cos, s_lo, s_hi,
              full_norm):
    lo = target * (1.0 - tol_ppm * 1e-6)
    hi = target * (1.0 + tol_ppm * 1e-6)
    p = np.searchsorted(smz, lo)
    best = -1
    bestc = -2.0
    m0 = p_scan[off[mono]]
    m1 = p_scan[off[mono + 1] - 1]
    while p < smz.shape[0] and smz[p] <= hi:
        if s_lo[p] > m1 or s_hi[p] < m0:   # no scan overlap: skip (contiguous, m/z-sorted arrays: cache friendly)
            p += 1
            continue
        h = sidx[p]
        c, ov = _trace_cos(off, p_scan, p_int, mono, h, full_norm)
        if ov >= 1:
            if ov < min_overlap_for_cos:
                c = max(c, min_trace_cos)   # too short to judge shape; accept on m/z + overlap
            if c >= min_trace_cos and c > bestc:
                bestc = c
                best = h
        p += 1
    return best, bestc


@nb.njit(cache=True, parallel=True)
def _hill_norm2(off, p_int):
    """Per-hill sum of squared intensities, accumulated exactly as _trace_cos(full_norm=True) does per call."""
    n_h = off.shape[0] - 1
    out = np.empty(n_h)
    for h in nb.prange(n_h):
        s = 0.0
        for q in range(off[h], off[h + 1]):
            s += p_int[q] * p_int[q]
        out[h] = s
    return out


@nb.njit(cache=True)
def _trace_cos_contig(off, h_start, h_end, p_int, norm2, a, b):
    """_trace_cos(full_norm=True) for hills covering consecutive scans (bitwise equal: same products, same order)."""
    s0 = max(h_start[a], h_start[b])
    s1 = min(h_end[a], h_end[b])
    if s1 < s0:
        return -1.0, 0
    i = off[a] + (s0 - h_start[a])
    j = off[b] + (s0 - h_start[b])
    sab = 0.0
    for k in range(s1 - s0 + 1):
        sab += p_int[i + k] * p_int[j + k]
    saa = norm2[a]
    sbb = norm2[b]
    if saa <= 0 or sbb <= 0:
        return 0.0, s1 - s0 + 1
    return sab / np.sqrt(saa * sbb), s1 - s0 + 1


@nb.njit(cache=True)
def _find_iso_fast(mono, target, tol_ppm, smz, sidx, off, h_start, h_end, p_int, norm2, min_trace_cos,
                   min_overlap_for_cos, s_lo, s_hi, first_only, p0=-1):
    """_find_iso (full_norm=True) on contiguous hills.  first_only: return the first passing isotope (existence).
    p0 >= 0: a position <= searchsorted(smz, lo) (monotone sweeps); advanced linearly to the same insertion point."""
    lo = target * (1.0 - tol_ppm * 1e-6)
    hi = target * (1.0 + tol_ppm * 1e-6)
    if p0 < 0:
        p = np.searchsorted(smz, lo)
    else:
        p = p0
        while p < smz.shape[0] and smz[p] < lo:
            p += 1
    best = -1
    bestc = -2.0
    m0 = h_start[mono]
    m1 = h_end[mono]
    while p < smz.shape[0] and smz[p] <= hi:
        if s_lo[p] > m1 or s_hi[p] < m0:
            p += 1
            continue
        h = sidx[p]
        c, ov = _trace_cos_contig(off, h_start, h_end, p_int, norm2, mono, h)
        if ov >= 1:
            if ov < min_overlap_for_cos:
                c = max(c, min_trace_cos)
            if c >= min_trace_cos and c > bestc:
                bestc = c
                best = h
                if first_only:
                    return best, bestc
        p += 1
    return best, bestc


SCAN_BIN = 32


@nb.njit(cache=True)
def _build_rt_bins(sidx, smz, h_start, h_end, n_scans, bs):
    """Hills bucketed by blocks of `bs` MS1 scans (a hill is in every block it touches), m/z ascending within a block
    (iterating the global m/z order, so ties keep the global order).  Entries carry the global m/z rank `pos`."""
    nbin = n_scans // bs + 1
    cnt = np.zeros(nbin + 1, np.int64)
    for q in range(sidx.shape[0]):
        h = sidx[q]
        for b in range(h_start[h] // bs, h_end[h] // bs + 1):
            cnt[b + 1] += 1
    for b in range(nbin):
        cnt[b + 1] += cnt[b]
    fill = cnt[:-1].copy()
    m = cnt[-1]
    b_mz = np.empty(m)
    b_pos = np.empty(m, np.int64)
    b_h = np.empty(m, np.int64)
    b_lo = np.empty(m, np.int64)
    b_hi = np.empty(m, np.int64)
    for q in range(sidx.shape[0]):
        h = sidx[q]
        for b in range(h_start[h] // bs, h_end[h] // bs + 1):
            k = fill[b]
            b_mz[k] = smz[q]
            b_pos[k] = q
            b_h[k] = h
            b_lo[k] = h_start[h]
            b_hi[k] = h_end[h]
            fill[b] += 1
    return cnt, b_mz, b_pos, b_h, b_lo, b_hi


@nb.njit(cache=True)
def _block_order(sidx, h_start, bs, nblk):
    """Stable counting sort of the m/z order by start block (== sidx[argsort(h_start[sidx] // bs, kind='stable')])."""
    cnt = np.zeros(nblk + 1, np.int64)
    for q in range(sidx.shape[0]):
        cnt[h_start[sidx[q]] // bs + 1] += 1
    for b in range(nblk):
        cnt[b + 1] += cnt[b]
    out = np.empty(sidx.shape[0], np.int64)
    for q in range(sidx.shape[0]):
        h = sidx[q]
        b = h_start[h] // bs
        out[cnt[b]] = h
        cnt[b] += 1
    return out


@nb.njit(cache=True)
def _candidates_in_visit_order(has1, visit, charges):
    """(ci, cz) in np.nonzero(has1) order, plus the processing permutation that visits them in `visit` order."""
    n_h, nz = has1.shape
    m = 0
    first = np.empty(n_h + 1, np.int64)
    for i in range(n_h):
        first[i] = m
        for zi in range(nz):
            if has1[i, zi]:
                m += 1
    first[n_h] = m
    ci = np.empty(m, np.int64)
    cz = np.empty(m, np.int64)
    k = 0
    for i in range(n_h):
        for zi in range(nz):
            if has1[i, zi]:
                ci[k] = i
                cz[k] = charges[zi]
                k += 1
    perm = np.empty(m, np.int64)
    k = 0
    for q in range(visit.shape[0]):
        i = visit[q]
        for t in range(first[i], first[i + 1]):
            perm[k] = t
            k += 1
    return ci, cz, perm


@nb.njit(cache=True)
def _find_iso_binned(mono, target, tol_ppm, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off, h_start, h_end, p_int,
                     norm2, min_trace_cos, min_overlap_for_cos, first_only, p_first=-1):
    """Same result as _find_iso(full_norm=True): among hills with lo <= mz <= hi that overlap the mono hill in scans,
    the highest trace cosine >= min_trace_cos, ties -> lowest global m/z rank (= first in the original sweep).
    Only the scan blocks of the mono hill are visited; a hill in several blocks counts in its first common block."""
    lo = target * (1.0 - tol_ppm * 1e-6)
    hi = target * (1.0 + tol_ppm * 1e-6)
    m0 = h_start[mono]
    m1 = h_end[mono]
    b0 = m0 // bs
    best = -1
    bestc = -2.0
    bestp = 1 << 62
    for b in range(b0, m1 // bs + 1):
        a = boff[b]
        e = boff[b + 1]
        if b == b0 and p_first >= 0:                 # caller's pointer: already the insertion point of lo
            p = p_first
        else:
            p = a + np.searchsorted(b_mz[a:e], lo)
        while p < e and b_mz[p] <= hi:
            if b_lo[p] > m1 or b_hi[p] < m0:
                p += 1
                continue
            if max(b0, b_lo[p] // bs) != b:          # already seen in an earlier common block
                p += 1
                continue
            h = b_h[p]
            c, ov = _trace_cos_contig(off, h_start, h_end, p_int, norm2, mono, h)
            if ov >= 1:
                if ov < min_overlap_for_cos:
                    c = max(c, min_trace_cos)
                if c >= min_trace_cos and (c > bestc or (c == bestc and b_pos[p] < bestp)):
                    bestc = c
                    best = h
                    bestp = b_pos[p]
                    if first_only:
                        return best, bestc
            p += 1
    return best, bestc


@nb.njit(cache=True)
def _avg_cos(table, mass, I, L):
    """(cosine of I[:L] with the averagine envelope, fraction of the envelope's abundance in its first L isotopes)."""
    r = mass / 10.0  # STEP
    r0 = int(r)
    if r0 >= table.shape[0] - 1:
        r0 = table.shape[0] - 2
        fr = 1.0
    else:
        fr = r - r0
    sab = 0.0
    saa = 0.0
    sbb = 0.0
    expl = 0.0
    for k in range(L):
        t = (1 - fr) * table[r0, k] + fr * table[r0 + 1, k] if k < table.shape[1] else 0.0
        sab += t * I[k]
        saa += I[k] * I[k]
        sbb += t * t
        expl += t
    if saa <= 0 or sbb <= 0:
        return 0.0, expl
    return sab / np.sqrt(saa * sbb), expl


@nb.njit(cache=True)
def _best_prefix(table, mass, I, n, min_avg_cos, min_explained):
    for L in range(n, 1, -1):
        c, e = _avg_cos(table, mass, I, L)
        if c >= min_avg_cos and e >= min_explained:
            return L, c
    return 1, 0.0


@nb.njit(cache=True)
def _chain(i, z, h_mz, h_apex_int, smz, sidx, off, p_scan, p_int, itol, max_iso, min_trace_cos,
           min_overlap_for_cos, table, min_avg_cos, min_explained, hills_out, cos_out, s_lo, s_hi, full_norm):
    """Fill hills_out[0..L) with the chain for (mono i, charge z); returns (L, averagine cos)."""
    hills_out[0] = i
    cos_out[0] = 1.0
    n = 1
    for k in range(1, max_iso):
        j, c = _find_iso(i, h_mz[i] + k * 1.0033548378 / z, itol, smz, sidx, off, p_scan, p_int,
                         min_trace_cos, min_overlap_for_cos, s_lo, s_hi, full_norm)
        if j < 0:
            break
        hills_out[n] = j
        cos_out[n] = c
        n += 1
    if n < 2:
        return 1, 0.0
    I = np.empty(n)
    for k in range(n):
        I[k] = h_apex_int[hills_out[k]]
    mass = (h_mz[i] - 1.007276466812) * z
    return _best_prefix(table, mass, I, n, min_avg_cos, min_explained)


@nb.njit(cache=True, parallel=True)
def _count_pass(h_mz, smz, sidx, off, p_scan, p_int, charges, itol, min_trace_cos, min_overlap_for_cos, s_lo, s_hi,
                full_norm):
    """Cheap pre-screen: (hill, charge) is a candidate iff its first isotope exists (necessary for >= 2 isotopes).
    Hills are visited in m/z order so consecutive lookups touch neighbouring memory."""
    n_h = h_mz.shape[0]
    nz = charges.shape[0]
    has1 = np.zeros((n_h, nz), np.bool_)
    for q in nb.prange(n_h):
        i = sidx[q]
        for zi in range(nz):
            j, c = _find_iso(i, h_mz[i] + 1.0033548378 / charges[zi], itol, smz, sidx, off, p_scan, p_int,
                             min_trace_cos, min_overlap_for_cos, s_lo, s_hi, full_norm)
            has1[i, zi] = j >= 0
    return has1


@nb.njit(cache=True, parallel=True)
def _fill_pass(ci, cz, h_mz, h_apex_int, smz, sidx, off, p_scan, p_int, itol, max_iso, min_trace_cos,
               min_overlap_for_cos, table, min_avg_cos, min_explained, s_lo, s_hi, full_norm):
    m = ci.shape[0]
    hills = np.full((m, max_iso), -1, np.int64)
    tcos = np.zeros((m, max_iso), np.float32)
    acos = np.zeros(m, np.float32)
    nn = np.zeros(m, np.int32)
    for q in nb.prange(m):
        hb = np.empty(max_iso, np.int64)
        cb = np.empty(max_iso)
        L, c = _chain(ci[q], cz[q], h_mz, h_apex_int, smz, sidx, off, p_scan, p_int, itol, max_iso,
                      min_trace_cos, min_overlap_for_cos, table, min_avg_cos, min_explained, hb, cb, s_lo, s_hi,
                      full_norm)
        nn[q] = L
        acos[q] = c
        for k in range(L):
            hills[q, k] = hb[k]
            tcos[q, k] = cb[k]
    return hills, tcos, acos, nn


@nb.njit(cache=True)
def _chain_fast(i, z, h_mz, h_apex_int, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off, h_start, h_end, p_int, norm2,
                itol, max_iso, min_trace_cos, min_overlap_for_cos, table, min_avg_cos, min_explained, hills_out,
                cos_out):
    hills_out[0] = i
    cos_out[0] = 1.0
    n = 1
    for k in range(1, max_iso):
        j, c = _find_iso_binned(i, h_mz[i] + k * 1.0033548378 / z, itol, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off,
                                h_start, h_end, p_int, norm2, min_trace_cos, min_overlap_for_cos, False)
        if j < 0:
            break
        hills_out[n] = j
        cos_out[n] = c
        n += 1
    if n < 2:
        return 1, 0.0
    I = np.empty(n)
    for k in range(n):
        I[k] = h_apex_int[hills_out[k]]
    mass = (h_mz[i] - 1.007276466812) * z
    return _best_prefix(table, mass, I, n, min_avg_cos, min_explained)


@nb.njit(cache=True, parallel=True)
def _count_pass_fast(h_mz, visit, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off, h_start, h_end, p_int, norm2, charges,
                     itol, min_trace_cos, min_overlap_for_cos):
    """visit: hill order (scan block, then m/z) so consecutive lookups search the same, cache-resident block; the
    result of each (hill, charge) does not depend on the visiting order."""
    n_h = h_mz.shape[0]
    nz = charges.shape[0]
    has1 = np.zeros((n_h, nz), np.bool_)
    # runs of `visit` with the same start block (m/z ascending inside): a forward pointer per charge into that block
    nblk = boff.shape[0] - 1
    vstart = np.zeros(nblk + 1, np.int64)
    for q in range(n_h):
        vstart[h_start[visit[q]] // bs + 1] += 1
    for b in range(nblk):
        vstart[b + 1] += vstart[b]
    for b in nb.prange(nblk):
        a = boff[b]
        e = boff[b + 1]
        for zi in range(nz):
            p = a
            for q in range(vstart[b], vstart[b + 1]):
                i = visit[q]
                target = h_mz[i] + 1.0033548378 / charges[zi]
                lo = target * (1.0 - itol * 1e-6)
                while p < e and b_mz[p] < lo:
                    p += 1
                j, c = _find_iso_binned(i, target, itol, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off, h_start, h_end,
                                        p_int, norm2, min_trace_cos, min_overlap_for_cos, True, p)
                has1[i, zi] = j >= 0
    return has1


@nb.njit(cache=True, parallel=True)
def _fill_pass_fast(ci, cz, perm, h_mz, h_apex_int, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off, h_start, h_end,
                    p_int, norm2, itol, max_iso, min_trace_cos, min_overlap_for_cos, table, min_avg_cos,
                    min_explained):
    """Candidates processed in `perm` order (cache locality), results written at their own index."""
    m = ci.shape[0]
    hills = np.full((m, max_iso), -1, np.int64)
    tcos = np.zeros((m, max_iso), np.float32)
    acos = np.zeros(m, np.float32)
    nn = np.zeros(m, np.int32)
    for qq in nb.prange(m):
        q = perm[qq]
        hb = np.empty(max_iso, np.int64)
        cb = np.empty(max_iso)
        L, c = _chain_fast(ci[q], cz[q], h_mz, h_apex_int, boff, b_mz, b_pos, b_h, b_lo, b_hi, bs, off, h_start,
                           h_end, p_int, norm2, itol, max_iso, min_trace_cos, min_overlap_for_cos, table, min_avg_cos,
                           min_explained, hb, cb)
        nn[q] = L
        acos[q] = c
        for k in range(L):
            hills[q, k] = hb[k]
            tcos[q, k] = cb[k]
    return hills, tcos, acos, nn


@nb.njit(cache=True)
def _greedy(order, hills, nn, cz, h_mz, h_apex_int, table, min_avg_cos, min_explained, n_h):
    used = np.zeros(n_h, np.bool_)
    m = order.shape[0]
    acc = np.zeros(m, np.bool_)
    acc_n = np.zeros(m, np.int32)
    acc_cos = np.zeros(m, np.float32)
    I = np.empty(hills.shape[1])
    for oi in range(m):
        q = order[oi]
        if used[hills[q, 0]]:
            continue
        L = 0
        for k in range(nn[q]):
            if used[hills[q, k]]:
                break
            L += 1
        if L < 2:
            continue
        for k in range(L):
            I[k] = h_apex_int[hills[q, k]]
        mass = (h_mz[hills[q, 0]] - 1.007276466812) * cz[q]
        best, bc = _best_prefix(table, mass, I, L, min_avg_cos, min_explained)
        if best < 2:
            continue
        for k in range(best):
            used[hills[q, k]] = True
        acc[q] = True
        acc_n[q] = best
        acc_cos[q] = bc
    return acc, acc_n, acc_cos


def charge_ban(ci: np.ndarray, cz: np.ndarray, nn: np.ndarray) -> np.ndarray:
    """Keep (mono, z) unless the same mono has a candidate at a multiple charge k*z (k >= 2) with >= as many isotopes."""
    keep = np.ones(len(ci), bool)
    if len(ci) == 0:
        return keep
    zmax = int(cz.max())
    n_h = int(ci.max()) + 1
    grid = np.zeros((n_h, zmax + 1), np.int32)
    np.maximum.at(grid, (ci, cz), nn)
    for z in range(1, zmax + 1):
        mult = [k * z for k in range(2, zmax // z + 1)]
        if not mult:
            continue
        m = cz == z
        hi = grid[ci[m]][:, mult].max(axis=1)
        keep[np.flatnonzero(m)[hi >= nn[m]]] = False
    return keep


def group_isotopes(h, *, charges=(1, 2, 3, 4, 5, 6), itol_ppm: float = 8.0, max_iso: int = 10,
                   min_trace_cos: float = 0.6, min_overlap_for_cos: int = 1, min_avg_cos: float = 0.8,
                   min_explained: float = 0.5, full_norm: bool = True, use_charge_ban: bool = True,
                   exact_reference: bool = False):
    """Return dict of arrays describing accepted features (see features.py for the table).  exact_reference: run the
    original kernels instead of the (bitwise-equal) fast path."""
    table = get_table()
    assert STEP == 10.0
    sidx = np.argsort(h.h_mz, kind="stable").astype(np.int64)
    smz = h.h_mz[sidx]
    charges = np.asarray(charges, np.int64)
    fast = (not exact_reference) and full_norm and bool(np.all(
        h.h_len.astype(np.int64) == h.h_end.astype(np.int64) - h.h_start.astype(np.int64) + 1))
    if fast:
        hs = np.ascontiguousarray(h.h_start, np.int64)
        he = np.ascontiguousarray(h.h_end, np.int64)
        norm2 = _hill_norm2(h.off, h.p_int)
        bs = SCAN_BIN
        rt_idx = _build_rt_bins(sidx, smz, hs, he, int(he.max()) + 1 if len(he) else 1, bs)
        visit = _block_order(sidx, hs, bs, len(rt_idx[0]) - 1)          # (start block, m/z rank)
        has1 = _count_pass_fast(h.h_mz, visit, *rt_idx, bs, h.off, hs, he, h.p_int, norm2, charges, itol_ppm,
                                min_trace_cos, min_overlap_for_cos)
    else:
        s_lo = np.ascontiguousarray(h.h_start[sidx])
        s_hi = np.ascontiguousarray(h.h_end[sidx])
        has1 = _count_pass(h.h_mz, smz, sidx, h.off, h.p_scan, h.p_int, charges, itol_ppm, min_trace_cos,
                           min_overlap_for_cos, s_lo, s_hi, full_norm)
    if fast:
        ci, cz, perm = _candidates_in_visit_order(has1, visit, charges)
        hills, tcos, acos, nn = _fill_pass_fast(ci, cz, perm, h.h_mz, h.h_apex_int, *rt_idx, bs, h.off, hs, he,
                                                h.p_int, norm2, itol_ppm, max_iso, min_trace_cos, min_overlap_for_cos,
                                                table, min_avg_cos, min_explained)
    else:
        ci, zi = np.nonzero(has1)
        cz = charges[zi]
        om = np.argsort(h.h_mz[ci], kind="stable")          # m/z order for memory locality, then restore
        hills, tcos, acos, nn = _fill_pass(ci[om].astype(np.int64), cz[om], h.h_mz, h.h_apex_int, smz, sidx, h.off,
                                           h.p_scan, h.p_int, itol_ppm, max_iso, min_trace_cos,
                                           min_overlap_for_cos, table, min_avg_cos, min_explained, s_lo, s_hi,
                                           full_norm)
        inv = np.empty_like(om)
        inv[om] = np.arange(len(om))
        hills, tcos, acos, nn = hills[inv], tcos[inv], acos[inv], nn[inv]
    n_cand = int(len(ci))
    ok = nn >= 2
    ci, cz, hills, tcos, acos, nn = ci[ok], cz[ok], hills[ok], tcos[ok], acos[ok], nn[ok]
    if use_charge_ban:
        k = charge_ban(ci, cz, nn)
        ci, cz, hills, tcos, acos, nn = ci[k], cz[k], hills[k], tcos[k], acos[k], nn[k]
    # priority: n isotopes desc, averagine cos desc, mono apex intensity desc, then hill index, charge
    order = np.lexsort((cz, ci, -h.h_apex_int[ci], -acos, -nn))
    acc, acc_n, acc_cos = _greedy(order, hills, nn, cz, h.h_mz, h.h_apex_int, table, min_avg_cos, min_explained, h.n)
    sel = np.flatnonzero(acc)
    return {
        "mono_hill": ci[sel].astype(np.int64),
        "charge": cz[sel].astype(np.int8),
        "n_isotopes": acc_n[sel],
        "iso_cos": acc_cos[sel],
        "hills": hills[sel],            # (n_feat, max_iso), -1 padded; only the first n_isotopes used
        "trace_cos": tcos[sel],
        "n_candidates": n_cand,
        "fast_path": fast,
    }
