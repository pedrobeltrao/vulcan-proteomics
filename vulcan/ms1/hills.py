"""Hill detection: link MS1 centroids across consecutive MS1 scans (biosaur2/Dinosaur style).

Pipeline (all batched over the whole run; no per-spectrum Python):

1. **Link** (`link_numpy` on CPU or `link_torch` on GPU). Every peak in MS1 scan s looks up the
   nearest peak (in m/z) of scan s-1 with one global searchsorted over the key scan*W + mz. It links
   to it if the ppm distance is within `htol_ppm`. A predecessor claimed by several peaks keeps only the
   closest (one-to-one). Optional gap linking (`gap` >= 1) joins hill starts to hill ends up to `gap`
   missing scans back.
2. **Label**: pointer jumping (root = root[root]) gives each peak its chain root -> hill id.
3. **Order**: counting sort of peaks by hill (stable, so scans stay ascending within a hill).
4. **Split** at valleys (numba): a local minimum of the 3-point smoothed trace splits the hill if
   min * hvf < both the left and right maxima of the remaining segments (biosaur2 `-hvf`, default 1.3),
   and both pieces keep >= minlh points.
5. **Filter** hills shorter than `minlh` scans and **summarise** (intensity-weighted m/z, apex, area).
"""
from __future__ import annotations

from dataclasses import dataclass

import numba as nb
import numpy as np

W_KEY = 16384.0  # key = scan_ordinal * W_KEY + mz; requires mz < W_KEY


@dataclass
class Hills:
    """Hills in CSR layout. Peak arrays are ordered by (hill, scan)."""
    off: np.ndarray        # int64 (n_hills+1)
    p_scan: np.ndarray     # int32 MS1 ordinal of each peak
    p_mz: np.ndarray       # float64
    p_int: np.ndarray      # float32
    p_src: np.ndarray      # int64 index into the run's flat MS1 peak arrays
    h_mz: np.ndarray       # float64 intensity-weighted mean m/z
    h_start: np.ndarray    # int32 first MS1 ordinal
    h_end: np.ndarray      # int32 last MS1 ordinal (inclusive)
    h_apex_scan: np.ndarray  # int32
    h_apex_int: np.ndarray   # float32
    h_area: np.ndarray       # float32, trapezoid over rt (intensity x minutes)
    h_len: np.ndarray        # int32 number of peaks

    @property
    def n(self) -> int:
        return len(self.off) - 1


# ----------------------------------------------------------------------------- linking (CPU)
def _best_prev(key, scan, mz, q_idx, q_key, q_mz, q_scan_target, htol_ppm):
    """Nearest key (same mz, target scan) among sorted `key`; returns (cand index or -1, ppm)."""
    n = len(key)
    pos = np.searchsorted(key, q_key)
    best = np.full(len(q_idx), -1, np.int64)
    bestd = np.full(len(q_idx), np.inf)
    for c in (pos - 1, pos):
        cc = np.clip(c, 0, n - 1)
        ok = (c >= 0) & (c < n) & (scan[cc] == q_scan_target)
        d = np.abs(mz[cc] - q_mz) / q_mz * 1e6
        d = np.where(ok, d, np.inf)
        better = d < bestd
        best = np.where(better, cc, best)
        bestd = np.where(better, d, bestd)
    best = np.where(bestd <= htol_ppm, best, -1)
    return best, bestd


def _one_to_one(q_idx, cand, dist):
    """Keep, for each claimed candidate, only the closest query. Returns (q, cand) pairs."""
    m = cand >= 0
    q, c, d = q_idx[m], cand[m], dist[m]
    if len(q) == 0:
        return q, c
    o = np.lexsort((q, d, c))   # by candidate, then distance, then query index (deterministic)
    cs = c[o]
    first = np.ones(len(cs), bool)
    first[1:] = cs[1:] != cs[:-1]
    return q[o][first], cs[first]


def link_numpy(scan: np.ndarray, mz: np.ndarray, htol_ppm: float, gap: int = 0) -> np.ndarray:
    """Predecessor index per peak (-1 = hill start). `scan` ascending, mz ascending within a scan."""
    n = len(mz)
    mz = mz.astype(np.float64, copy=False)
    key = scan.astype(np.float64) * W_KEY + mz
    pred = np.full(n, -1, np.int64)
    idx = np.arange(n, dtype=np.int64)
    cand, d = _best_prev(key, scan, mz, idx, key - W_KEY, mz, scan - 1, htol_ppm)
    q, c = _one_to_one(idx, cand, d)
    pred[q] = c
    for g in range(1, gap + 1):
        has_succ = np.zeros(n, bool)
        has_succ[pred[pred >= 0]] = True
        ends = np.flatnonzero(~has_succ)
        starts = np.flatnonzero(pred < 0)
        e_key = key[ends]
        cand, d = _best_prev(e_key, scan[ends], mz[ends], starts, key[starts] - (g + 1) * W_KEY,
                             mz[starts], scan[starts] - (g + 1), htol_ppm)
        cand = np.where(cand >= 0, ends[np.clip(cand, 0, None)], -1)
        q, c = _one_to_one(starts, cand, d)
        pred[q] = c
    return pred


def roots_numpy(pred: np.ndarray) -> np.ndarray:
    root = np.where(pred >= 0, pred, np.arange(len(pred)))
    while True:
        nxt = root[root]
        if np.array_equal(nxt, root):
            return root
        root = nxt


# ----------------------------------------------------------------------------- linking (GPU)
def link_torch(scan: np.ndarray, mz: np.ndarray, htol_ppm: float, gap: int = 0, device="cuda"):
    """Same as link_numpy + roots_numpy on a torch device; returns (pred, root) as numpy int64."""
    import torch
    t_scan = torch.from_numpy(np.ascontiguousarray(scan, dtype=np.int32)).to(device, non_blocking=True)
    t_mz = torch.from_numpy(np.ascontiguousarray(mz, dtype=np.float64)).to(device, non_blocking=True)
    n = t_mz.shape[0]
    key = t_scan.double() * W_KEY + t_mz
    idx = torch.arange(n, device=device)

    def best_prev(k_sorted, k_scan, k_mz, q_key, q_mz, q_target):
        pos = torch.searchsorted(k_sorted, q_key)
        m = k_sorted.shape[0]
        best = torch.full_like(pos, -1)
        bestd = torch.full_like(q_mz, float("inf"))
        for c in (pos - 1, pos):
            cc = c.clamp(0, max(m - 1, 0))
            ok = (c >= 0) & (c < m) & (k_scan[cc] == q_target)
            d = (k_mz[cc] - q_mz).abs() / q_mz * 1e6
            d = torch.where(ok, d, torch.full_like(d, float("inf")))
            better = d < bestd
            best = torch.where(better, cc, best)
            bestd = torch.where(better, d, bestd)
        best = torch.where(bestd <= htol_ppm, best, torch.full_like(best, -1))
        return best, bestd

    def one_to_one(q, c, d):
        m = c >= 0
        q, c, d = q[m], c[m], d[m]
        if q.numel() == 0:
            return q, c
        o = torch.sort(q, stable=True).indices          # lexsort (c, d, q) via stable passes
        o = o[torch.sort(d[o], stable=True).indices]
        o = o[torch.sort(c[o], stable=True).indices]
        cs = c[o]
        first = torch.ones_like(cs, dtype=torch.bool)
        first[1:] = cs[1:] != cs[:-1]
        return q[o][first], cs[first]

    pred = torch.full((n,), -1, dtype=torch.int64, device=device)
    cand, d = best_prev(key, t_scan, t_mz, key - W_KEY, t_mz, t_scan - 1)
    q, c = one_to_one(idx, cand, d)
    pred[q] = c
    for g in range(1, gap + 1):
        has_succ = torch.zeros(n, dtype=torch.bool, device=device)
        has_succ[pred[pred >= 0]] = True
        ends = torch.nonzero(~has_succ).squeeze(1)
        starts = torch.nonzero(pred < 0).squeeze(1)
        cand, d = best_prev(key[ends], t_scan[ends], t_mz[ends], key[starts] - (g + 1) * W_KEY,
                            t_mz[starts], t_scan[starts] - (g + 1))
        cand = torch.where(cand >= 0, ends[cand.clamp(0)], torch.full_like(cand, -1))
        q, c = one_to_one(starts, cand, d)
        pred[q] = c
    root = torch.where(pred >= 0, pred, idx)
    while True:
        nxt = root[root]
        if torch.equal(nxt, root):
            break
        root = nxt
    return pred.cpu().numpy(), root.cpu().numpy()


# ----------------------------------------------------------------------------- numba kernels
@nb.njit(cache=True, parallel=True)
def _link_adjacent(soff, mz, htol_ppm, pred):
    """pred[i] = nearest peak (m/z) of the previous MS1 scan within htol_ppm, one-to-one (the closest claimer wins,
    ties -> lower index). Same result as link_numpy(gap=0); O(n) two-pointer merge, parallel over scans."""
    n_s = soff.shape[0] - 1
    for s in nb.prange(1, n_s):
        a, b = soff[s], soff[s + 1]
        pa_, pb = soff[s - 1], soff[s]
        if b <= a or pb <= pa_:
            continue
        cand = np.full(b - a, -1, np.int64)
        dist = np.full(b - a, np.inf)
        j = pa_
        for i in range(a, b):
            x = mz[i]
            while j + 1 < pb and mz[j + 1] <= x:
                j += 1
            bd = np.inf
            bc = -1
            for c in (j, j + 1):
                if c < pb:
                    d = abs(mz[c] - x) / x * 1e6
                    if d < bd:
                        bd = d
                        bc = c
            if bd <= htol_ppm:
                cand[i - a] = bc
                dist[i - a] = bd
        claim_d = np.full(pb - pa_, np.inf)
        claim_q = np.full(pb - pa_, -1, np.int64)
        for i in range(a, b):
            c = cand[i - a]
            if c >= 0 and dist[i - a] < claim_d[c - pa_]:
                claim_d[c - pa_] = dist[i - a]
                claim_q[c - pa_] = i
        for c in range(pb - pa_):
            if claim_q[c] >= 0:
                pred[claim_q[c]] = pa_ + c


@nb.njit(cache=True)
def _roots_seq(pred):
    """pred[i] < i always (predecessors are in earlier scans), so one ascending pass resolves every root."""
    n = pred.shape[0]
    root = np.empty(n, np.int64)
    for i in range(n):
        p = pred[i]
        root[i] = i if p < 0 else root[p]
    return root


@nb.njit(cache=True)
def _hill_ids(pred):
    """Hill id per peak = rank (by index) of its chain root; same numbering as roots + cumsum(is_root) - 1."""
    n = pred.shape[0]
    hill = np.empty(n, np.int64)
    k = 0
    for i in range(n):
        p = pred[i]
        if p < 0:
            hill[i] = k
            k += 1
        else:
            hill[i] = hill[p]
    return hill, k


@nb.njit(cache=True)
def _sort_gather(hill, n_hills, scan, mz, inten, src):
    """Stable counting sort by hill fused with the gathers (same arrays as _counting_sort + fancy indexing)."""
    n = hill.shape[0]
    cnt = np.zeros(n_hills + 1, np.int64)
    for i in range(n):
        cnt[hill[i] + 1] += 1
    for h in range(n_hills):
        cnt[h + 1] += cnt[h]
    off = cnt.copy()
    pos = cnt[:-1].copy()
    p_scan = np.empty(n, np.int32)
    p_mz = np.empty(n, np.float64)
    p_int = np.empty(n, np.float32)
    p_src = np.empty(n, np.int64)
    for i in range(n):
        h = hill[i]
        q = pos[h]
        p_scan[q] = scan[i]
        p_mz[q] = mz[i]
        p_int[q] = inten[i]
        p_src[q] = src[i]
        pos[h] = q + 1
    return off, p_scan, p_mz, p_int, p_src


@nb.njit(cache=True)
def _compact(lab, minlh, p_scan, p_mz, p_int, p_src):
    """Drop hills with < minlh peaks (same result as bincount / keep mask / boolean compression / cumsum)."""
    n = lab.shape[0]
    n_h = lab[n - 1] + 1 if n else 0
    cnt = np.zeros(n_h, np.int64)
    for i in range(n):
        cnt[lab[i]] += 1
    nk = 0
    m = 0
    for h in range(n_h):
        if cnt[h] >= minlh:
            nk += 1
            m += cnt[h]
    off = np.zeros(nk + 1, np.int64)
    k = 0
    for h in range(n_h):
        if cnt[h] >= minlh:
            off[k + 1] = off[k] + cnt[h]
            k += 1
    o_scan = np.empty(m, np.int32)
    o_mz = np.empty(m, np.float64)
    o_int = np.empty(m, np.float32)
    o_src = np.empty(m, np.int64)
    w = 0
    for i in range(n):
        if cnt[lab[i]] >= minlh:
            o_scan[w] = p_scan[i]
            o_mz[w] = p_mz[i]
            o_int[w] = p_int[i]
            o_src[w] = p_src[i]
            w += 1
    return off, o_scan, o_mz, o_int, o_src


@nb.njit(cache=True)
def _counting_sort(hill, n_hills):
    cnt = np.zeros(n_hills + 1, np.int64)
    for i in range(hill.shape[0]):
        cnt[hill[i] + 1] += 1
    for h in range(n_hills):
        cnt[h + 1] += cnt[h]
    off = cnt.copy()
    pos = cnt[:-1].copy()
    order = np.empty(hill.shape[0], np.int64)
    for i in range(hill.shape[0]):
        h = hill[i]
        order[pos[h]] = i
        pos[h] += 1
    return order, off


@nb.njit(cache=True, parallel=True)
def _split_labels(off, inten, hvf, minlh):
    """New hill label per (hill-sorted) peak after valley splitting; labels are 0..m-1 ascending.

    A local minimum of the 3-point-smoothed trace starts a new segment if min * hvf < both the maximum of the current
    segment (left) and of the remainder (right), and both sides keep >= minlh points. Parallel over hills: segment
    starts are flagged per peak, labels are their cumulative sum."""
    n_h = off.shape[0] - 1
    start = np.zeros(off[-1], np.int64)
    for h in nb.prange(n_h):
        a, b = off[h], off[h + 1]
        L = b - a
        if L == 0:
            continue
        start[a] = 1
        if L < 2 * minlh + 1 or L < 5:
            continue
        s = np.empty(L)
        for i in range(L):
            lo = max(0, i - 1)
            hi = min(L - 1, i + 1)
            acc = 0.0
            for j in range(lo, hi + 1):
                acc += inten[a + j]
            s[i] = acc / (hi - lo + 1)
        rmax = np.empty(L + 1)          # suffix maxima: rmax[i] = max(s[i:])
        rmax[L] = 0.0
        for i in range(L - 1, -1, -1):
            rmax[i] = max(rmax[i + 1], s[i])
        seg = 0
        lmax = s[0]
        for i in range(1, L):
            if i - seg >= minlh and L - i >= minlh and i < L - 1:
                if s[i] < s[i - 1] and s[i] <= s[i + 1]:
                    if s[i] * hvf < lmax and s[i] * hvf < rmax[i + 1]:
                        start[a + i] = 1
                        seg = i
                        lmax = s[i]
                        continue
            if s[i] > lmax:
                lmax = s[i]
    lab = np.cumsum(start) - 1
    return lab


@nb.njit(cache=True)
def _summaries(off, p_scan, p_mz, p_int, ms1_rt):
    n_h = off.shape[0] - 1
    h_mz = np.empty(n_h)
    h_start = np.empty(n_h, np.int32)
    h_end = np.empty(n_h, np.int32)
    h_apex_scan = np.empty(n_h, np.int32)
    h_apex_int = np.empty(n_h, np.float32)
    h_area = np.empty(n_h, np.float32)
    h_len = np.empty(n_h, np.int32)
    for h in range(n_h):
        a, b = off[h], off[h + 1]
        sw = 0.0
        smz = 0.0
        best = -1.0
        bi = a
        area = 0.0
        for i in range(a, b):
            w = p_int[i]
            sw += w
            smz += w * p_mz[i]
            if w > best:
                best = w
                bi = i
            if i > a:
                area += 0.5 * (p_int[i] + p_int[i - 1]) * (ms1_rt[p_scan[i]] - ms1_rt[p_scan[i - 1]])
        h_mz[h] = smz / sw if sw > 0 else p_mz[a]
        h_start[h] = p_scan[a]
        h_end[h] = p_scan[b - 1]
        h_apex_scan[h] = p_scan[bi]
        h_apex_int[h] = best
        h_area[h] = area
        h_len[h] = b - a
    return h_mz, h_start, h_end, h_apex_scan, h_apex_int, h_area, h_len


# ----------------------------------------------------------------------------- driver
def detect_hills(scan: np.ndarray, mz: np.ndarray, inten: np.ndarray, ms1_rt: np.ndarray, *,
                 htol_ppm: float = 8.0, minlh: int = 2, hvf: float = 1.3, gap: int = 0,
                 min_intensity: float = 0.0, device: str = "cpu", timings: dict | None = None,
                 reference: bool = False) -> Hills:
    """Hills from flat MS1 peaks (scan = MS1 ordinal, ascending; mz ascending within a scan).
    The CPU gap-0 path uses fused numba passes (E7b; identical arrays); reference=True keeps the original steps."""
    import time
    t = timings if timings is not None else {}
    t0 = time.perf_counter()
    src = np.arange(len(mz), dtype=np.int64)
    if min_intensity > 0:
        keep = inten >= min_intensity
        scan, mz, inten, src = scan[keep], mz[keep], inten[keep], src[keep]
    if device.startswith("cuda"):
        _, root = link_torch(scan, mz, htol_ppm, gap, device)
    elif gap == 0:
        soff = np.searchsorted(scan, np.arange(int(ms1_rt.shape[0]) + 1)).astype(np.int64)
        pred = np.full(len(mz), -1, np.int64)
        _link_adjacent(soff, np.ascontiguousarray(mz, dtype=np.float64), float(htol_ppm), pred)
        if not reference:
            hill, n_h0 = _hill_ids(pred)
            t1 = time.perf_counter()
            t["link"] = t.get("link", 0.0) + t1 - t0
            off, p_scan, p_mz, p_int, p_src = _sort_gather(hill, n_h0, scan, mz, inten, src)
            del hill
            lab = _split_labels(off, p_int, hvf, minlh)
            off, p_scan, p_mz, p_int, p_src = _compact(lab, minlh, p_scan, p_mz, p_int, p_src)
            t2 = time.perf_counter()
            t["split_filter"] = t.get("split_filter", 0.0) + t2 - t1
            s = _summaries(off, p_scan, p_mz, p_int, ms1_rt.astype(np.float64))
            t["summaries"] = t.get("summaries", 0.0) + time.perf_counter() - t2
            return Hills(off, p_scan, p_mz, p_int, p_src, *s)
        root = _roots_seq(pred)
    else:
        root = _roots_seq(link_numpy(scan, mz, htol_ppm, gap))
    t1 = time.perf_counter()
    t["link"] = t.get("link", 0.0) + t1 - t0
    is_root = root == np.arange(len(root))
    hill_of_root = np.cumsum(is_root) - 1
    hill = hill_of_root[root]
    order, off = _counting_sort(hill, int(is_root.sum()))
    p_scan = scan[order].astype(np.int32)
    p_mz = mz[order].astype(np.float64)
    p_int = inten[order].astype(np.float32)
    p_src = src[order]
    lab = _split_labels(off, p_int, hvf, minlh)
    n_h = int(lab[-1]) + 1 if len(lab) else 0
    cnt = np.bincount(lab, minlength=n_h)
    keep_h = cnt >= minlh
    keep_p = keep_h[lab]
    p_scan, p_mz, p_int, p_src = p_scan[keep_p], p_mz[keep_p], p_int[keep_p], p_src[keep_p]
    cnt = cnt[keep_h]
    off = np.zeros(len(cnt) + 1, np.int64)
    np.cumsum(cnt, out=off[1:])
    t2 = time.perf_counter()
    t["split_filter"] = t.get("split_filter", 0.0) + t2 - t1
    s = _summaries(off, p_scan, p_mz, p_int, ms1_rt.astype(np.float64))
    t["summaries"] = t.get("summaries", 0.0) + time.perf_counter() - t2
    return Hills(off, p_scan, p_mz, p_int, p_src, *s)
