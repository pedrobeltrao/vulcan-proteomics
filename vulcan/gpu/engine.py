"""GPU closed search: K0 (peak view) -> hypotheses -> K1 (tiled inverted index, top-K, pair closure)
-> K2 (exact rescoring, label-blind features).  DESIGN.md §4 L1-L3, §10 E2.

Scoring semantics follow sage-core 99407db6 (Sage = oracle):
  * K0: peaks with m/z in [fragment_min_mz, fragment_max_mz], mass = m/z - PROTON, top `max_peaks` by intensity
        (global, Sage) or top-k per `window_th` Th window (DESIGN default view), then mass-sorted.
  * K1: preliminary count = #(peak, theoretical fragment) pairs within the fragment tolerance, fragments with
        ion index > min_ion_index and mass in [fragment_min_mz, fragment_max_mz] (Sage's fragment index).
        Thresholded append (count >= k1_threshold, identical for every candidate), then top-K per spectrum
        over all its hypotheses (ties: higher row first, i.e. heavier/hash order; never a label).
  * pair closure: the partner (target<->decoy, label-symmetric) of every shortlisted candidate is scored too.
  * K2: Sage hyperscore and features over all b/y ions; PSM rank by hyperscore; delta_next to the next
        candidate regardless of label.
Nothing in this module reads labels; vulcan/index/build.labels.arrow is joined only in L4.
"""
from __future__ import annotations

import math
import resource
import time
from dataclasses import dataclass, field

import numpy as np
import pyarrow as pa
import torch

from . import _lib as K

from .host import NEUTRON, PROTON, SearchParams, Spectra, make_hypotheses, pre_hypotheses  # noqa: F401  (torch-free)
SAGE_MONO = np.array([71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891, 113.08406,
                      0.0, 128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276, 128.05858, 156.1011,
                      87.03203, 101.04768, 150.95363, 99.06841, 186.07932, 0.0, 163.06332, 0.0], np.float32)


ROW_BITS = 28   # peptidoform rows < 2^28 (select_topk key packing); chunk_size 32768 x 2^13 chunks = 2^28


def _u32_to_i32(x: torch.Tensor) -> torch.Tensor:
    x = x.to(torch.int64)
    return torch.where(x >= 2 ** 31, x - 2 ** 32, x).to(torch.int32)


@dataclass
class Slice:
    chunk0: int
    chunk1: int
    offsets: torch.Tensor   # int32 bit pattern of u32, (chunk1-chunk0)*n_bins + 1
    entries: torch.Tensor   # int32 bit pattern of u32 (local << 16 | fine)


@dataclass
class Stats:
    h2d_bytes: int = 0
    gpu_ms: dict = field(default_factory=dict)
    wall_s: dict = field(default_factory=dict)
    counts: dict = field(default_factory=dict)


class GpuIndex:
    """Device-resident peptidoform table + K1 tiled inverted index (CSR over (chunk, log-ppm bin)).

    Slicing path: rows are cut at chunk boundaries into slices of <= max_slice_entries postings (u32 offsets);
    all slices are resident here (80 GB covers closed/phospho/semi human). For larger indexes the same Slice
    objects can live in pinned host memory and be streamed per batch (TODO when open/non-specific modes land).
    """

    def __init__(self, l1, params: SearchParams, static_mods: dict, var_mods: dict, device="cuda"):
        t0 = time.perf_counter()
        self.p = params
        self.dev = torch.device(device)
        L = K.lib()
        stat = np.zeros(26, np.float32)
        var = np.zeros(26, np.float32)
        for r, m in static_mods.items():
            stat[ord(r) - 65] = np.float32(m)
        nterm = 0.0
        for r, m in var_mods.items():
            if r == "^":          # M4: variable peptide-N-term mod, mask bit 63 (vulcan_core NTERM_BIT)
                nterm = float(m)
                continue
            var[ord(r) - 65] = np.float32(m)
        K.check(L.vk_set_params(SAGE_MONO.ctypes.data, stat.ctypes.data, var.ctypes.data), "set_params")
        K.check(L.vk_set_nterm(nterm), "set_nterm")
        assert params.fragment_unit in ("ppm", "da")
        K.check(L.vk_set_tol_da(int(params.fragment_unit == "da")), "set_tol_da")   # 3c (process-wide, like set_params)
        self.n = len(l1)
        d = self.dev
        self.h2d = 0

        def up(a):
            t = torch.from_numpy(np.array(a, copy=True)).to(d)
            self.h2d += t.numel() * t.element_size()
            return t
        self.mass32_np = l1.mass32
        self.seq_off = up(l1.seq_off.astype(np.int64))
        self.length = up(l1.length.astype(np.uint8))
        self.mask = up(l1.mask.view(np.int64))
        self.mass32 = up(l1.mass32.astype(np.float32))
        self.residues = up(l1.residues.astype(np.uint8))
        self.partner = up(l1.partner.astype(np.int64))
        self.inv_logw = 1.0 / math.log1p(params.bin_ppm * 1e-6)
        self.n_bins = int(math.floor(math.log(params.fragment_max_mz / params.fragment_min_mz) * self.inv_logw)) + 2
        assert self.n_bins < 2 ** 18, "bin index must fit 18 bits; increase bin_ppm"
        C = params.chunk_size
        assert C <= 32768 and C % params.tile == 0
        self.n_chunks = (self.n + C - 1) // C
        assert self.n_chunks < 2 ** 13 and self.n < 2 ** ROW_BITS
        stream = K.stream()
        cnt = torch.empty(self.n, dtype=torch.int32, device=d)
        K.check(L.vk_frag_count(self.n, K.ptr(self.seq_off), K.ptr(self.length), K.ptr(self.mask),
                                K.ptr(self.mass32), K.ptr(self.residues), params.min_ion_index,
                                params.fragment_min_mz, params.fragment_max_mz, K.ptr(cnt), stream), "frag_count")
        # slice boundaries (whole chunks)
        per_chunk = torch.zeros(self.n_chunks, dtype=torch.int64, device=d)
        per_chunk.index_add_(0, torch.arange(self.n, device=d) // C, cnt.to(torch.int64))
        pcs = per_chunk.cpu().numpy()
        bounds, acc, c0 = [], 0, 0
        for c, v in enumerate(pcs):
            if acc + v > params.max_slice_entries and c > c0:
                bounds.append((c0, c))
                c0, acc = c, 0
            acc += int(v)
        bounds.append((c0, self.n_chunks))
        self.slices: list[Slice] = []
        self.n_postings = int(pcs.sum())
        for c0, c1 in bounds:
            r0, r1 = c0 * C, min(c1 * C, self.n)
            cc = cnt[r0:r1].to(torch.int64)
            off = torch.zeros(r1 - r0 + 1, dtype=torch.int64, device=d)
            torch.cumsum(cc, 0, out=off[1:])
            tot = int(off[-1])
            keys = torch.empty(tot, dtype=torch.int64, device=d)
            K.check(L.vk_frag_emit(r1 - r0, K.ptr(self.seq_off[r0:r1]), K.ptr(self.length[r0:r1]),
                                   K.ptr(self.mask[r0:r1]), K.ptr(self.mass32[r0:r1]), K.ptr(self.residues),
                                   params.min_ion_index, params.fragment_min_mz, params.fragment_max_mz,
                                   self.inv_logw, C, r0, K.ptr(off), K.ptr(keys), stream), "frag_emit")
            keys = torch.sort(keys).values
            entries = _u32_to_i32(keys & 0xFFFFFFFF)
            lin = ((keys >> 50) - c0) * self.n_bins + ((keys >> 32) & 0x3FFFF)
            del keys
            nb = (c1 - c0) * self.n_bins
            hist = torch.bincount(lin, minlength=nb)
            del lin
            offs = torch.zeros(nb + 1, dtype=torch.int64, device=d)
            torch.cumsum(hist, 0, out=offs[1:])
            del hist
            assert int(offs[-1]) < 2 ** 32
            self.slices.append(Slice(c0, c1, _u32_to_i32(offs), entries))
        torch.cuda.synchronize()
        self.build_s = time.perf_counter() - t0
        self.device_bytes = sum(t.numel() * t.element_size() for t in
                                [self.seq_off, self.length, self.mask, self.mass32, self.residues, self.partner]) + \
            sum(s.offsets.numel() * 4 + s.entries.numel() * 4 for s in self.slices)


# ------------------------------------------------------------------------------------------------ spectra
def load_store(path: str) -> Spectra:
    from vulcan.io import store
    t = store.read_table(path)
    lvl = t.column("ms_level").to_numpy()
    t = t.filter(pa.array(lvl == 2))
    pk = store.peaks_numpy(t)
    z = t.column("precursor_charge").to_numpy(zero_copy_only=False)
    z = np.nan_to_num(np.asarray(z, dtype=np.float64), nan=0).astype(np.int8)
    pmz = np.asarray(t.column("precursor_mz").to_numpy(zero_copy_only=False), dtype=np.float64)
    name = store.run_meta(path).get("run_name", path)
    return Spectra(name, t.column("scan").to_numpy(), t.column("index").to_numpy(),
                   np.asarray(t.column("rt_min").to_numpy(), np.float64), pmz, z, pk.mz, pk.intensity, pk.offsets)


def k0_preprocess(sp: Spectra, p: SearchParams, dev, stats: Stats):
    """K0 on the GPU: range filter, top-k view, mass sort. Returns device peaks + host validity mask."""
    n = len(sp)
    mz = torch.from_numpy(sp.mz).to(dev)
    it = torch.from_numpy(sp.intensity).to(dev)
    lens = torch.from_numpy(np.diff(sp.offsets)).to(dev)
    stats.h2d_bytes += mz.numel() * 4 + it.numel() * 4 + lens.numel() * 8
    seg = torch.repeat_interleave(torch.arange(n, device=dev), lens)
    keep = (mz >= p.fragment_min_mz) & (mz <= p.fragment_max_mz) & ~torch.isnan(mz)
    mass = mz - torch.tensor(PROTON, device=dev)
    seg, mass, it, mzk = seg[keep], mass[keep], it[keep], mz[keep]
    grp = seg if p.peak_mode == "global" else seg * 1000 + torch.floor(mzk / p.window_th).to(torch.int64)
    kk = p.max_peaks if p.peak_mode == "global" else p.window_k
    o = torch.argsort(mass, descending=True, stable=True)
    o = o[torch.argsort(it[o], descending=True, stable=True)]
    o = o[torch.argsort(grp[o], stable=True)]
    g = grp[o]
    first = torch.searchsorted(g, g, side="left")
    rank = torch.arange(len(g), device=dev) - first
    o = o[rank < kk]
    seg, mass, it = seg[o], mass[o], it[o]
    o = torch.argsort(mass, stable=True)
    o = o[torch.argsort(seg[o], stable=True)]
    seg, mass, it = seg[o].contiguous(), mass[o].contiguous(), it[o].contiguous()
    if p.peak_mode == "window" and p.max_peaks:   # cap total peaks by intensity as well
        pass
    cnt = torch.bincount(seg, minlength=n)
    off = torch.zeros(n + 1, dtype=torch.int64, device=dev)
    torch.cumsum(cnt, 0, out=off[1:])

    # TIC per spectrum, deterministic: f64 row sums of a zero-padded (spectrum x peak) matrix (fixed reduction order).
    # M2b: the previous f64 global prefix-sum difference was not bit-stable run to run (CUDA scan association
    # varies with scheduling; seen as last-bit changes of matched_intensity_pct under concurrent streams).
    n_max = int(cnt.max()) if n else 0
    pad = torch.zeros((n, max(n_max, 1)), dtype=torch.float64, device=dev)
    pad[seg, torch.arange(len(seg), device=dev) - off[:-1][seg]] = it.to(torch.float64)
    tic = pad.sum(1).to(torch.float32)
    del pad
    valid = (cnt >= p.min_peaks).cpu().numpy()
    return {"off": off, "mass": mass, "int": it, "tic": tic, "valid": valid}


def work_items(hyp, p: SearchParams):
    """Tiles aligned to multiples of `tile` (which divides chunk_size): never cross a chunk boundary."""
    lo, hi = hyp["lo"], hyp["hi"]
    T = p.tile
    nz = hi > lo
    h = np.nonzero(nz)[0]
    t0 = lo[h] // T
    t1 = (hi[h] - 1) // T
    nt = (t1 - t0 + 1)
    hh = np.repeat(h, nt)
    tt = np.repeat(t0, nt) + (np.arange(nt.sum()) - np.repeat(np.cumsum(nt) - nt, nt))
    a = np.maximum(lo[hh], tt * T)
    b = np.minimum(hi[hh], (tt + 1) * T)
    chunk = (tt * T) // p.chunk_size
    base = chunk * p.chunk_size
    return {"hyp": hh.astype(np.int32), "chunk": chunk.astype(np.int32), "lo": (a - base).astype(np.int32),
            "hi": (b - base).astype(np.int32), "n_window": int((hi - lo)[nz].sum())}


class Timer:
    def __init__(self, stats: Stats, name: str):
        self.stats, self.name = stats, name

    def __enter__(self):
        self.e0 = torch.cuda.Event(enable_timing=True)
        self.e1 = torch.cuda.Event(enable_timing=True)
        self.e0.record()
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, *a):
        self.e1.record()
        self.e1.synchronize()
        self.stats.gpu_ms[self.name] = self.stats.gpu_ms.get(self.name, 0.0) + self.e0.elapsed_time(self.e1)
        self.stats.wall_s[self.name] = self.stats.wall_s.get(self.name, 0.0) + time.perf_counter() - self.t0


def run_k1(gi: GpuIndex, hyp, peaks, p: SearchParams, stats: Stats):
    L = K.lib()
    dev = gi.dev
    wi = work_items(hyp, p)
    n_h = len(hyp["spec"])
    d_spec = torch.from_numpy(hyp["spec"]).to(dev)
    d_maxfz = torch.from_numpy(hyp["maxfz"]).to(dev)
    scored = torch.zeros(n_h, dtype=torch.int32, device=dev)
    matched = torch.zeros(n_h, dtype=torch.int64, device=dev)
    cap = max(1, wi["n_window"])
    n_out = torch.zeros(1, dtype=torch.int64, device=dev)
    o_h = torch.empty(cap, dtype=torch.int32, device=dev)
    o_r = torch.empty(cap, dtype=torch.int32, device=dev)
    o_c = torch.empty(cap, dtype=torch.int32, device=dev)
    for sl in gi.slices:
        m = (wi["chunk"] >= sl.chunk0) & (wi["chunk"] < sl.chunk1)
        if not m.any():
            continue
        it = {k: torch.from_numpy(np.ascontiguousarray(wi[k][m])).to(dev) for k in ("hyp", "chunk", "lo", "hi")}
        stats.h2d_bytes += sum(v.numel() * 4 for v in it.values())
        K.check(L.vk_k1(int(m.sum()), K.ptr(it["hyp"]), K.ptr(it["chunk"]), K.ptr(it["lo"]), K.ptr(it["hi"]),
                        K.ptr(d_spec), K.ptr(d_maxfz), K.ptr(peaks["off"]), K.ptr(peaks["mass"]),
                        K.ptr(sl.offsets), K.ptr(sl.entries), gi.n_bins, sl.chunk0, p.chunk_size, p.tile,
                        p.fragment_ppm[0], p.fragment_ppm[1], p.fragment_min_mz, gi.inv_logw, p.k1_threshold,
                        cap, K.ptr(n_out), K.ptr(o_h), K.ptr(o_r), K.ptr(o_c), K.ptr(scored), K.ptr(matched),
                        K.stream()), "k1")
    n = int(n_out.item())
    assert n <= cap
    stats.counts["k1_work_items"] = stats.counts.get("k1_work_items", 0) + len(wi["hyp"])
    stats.counts["k1_window_candidates"] = stats.counts.get("k1_window_candidates", 0) + wi["n_window"]
    stats.counts["k1_appended"] = stats.counts.get("k1_appended", 0) + n
    return o_h[:n], o_r[:n], o_c[:n], scored, matched, d_spec, d_maxfz


def select_topk(o_h, o_r, o_c, d_spec, gi: GpuIndex, p: SearchParams):
    """Top-K per spectrum (count desc, then higher row), then pair closure. Returns (hyp, row, k1count) sorted
    by (spectrum, row)."""
    dev = gi.dev
    spec = d_spec[o_h.long()].long()
    # M4: 28 row bits (phospho/var-mod indexes up to 268M peptidoforms with chunk_size 32768); K1 counts are far
    # below 2^15, so the (spectrum, count desc, row desc) order -- and every E2/M2b output -- is unchanged.
    key = (spec << 43) | ((32767 - o_c.long().clamp(max=32767)) << ROW_BITS) | ((2 ** ROW_BITS - 1) - o_r.long())
    o = torch.argsort(key)
    ss = spec[o]
    first = torch.searchsorted(ss, ss, side="left")
    rank = torch.arange(len(ss), device=dev) - first
    sel = o[rank < p.topk]
    h, r, c = o_h[sel].long(), o_r[sel].long(), o_c[sel].long()
    if p.pair_closure:
        pr = gi.partner[r]
        h = torch.cat([h, h])
        r = torch.cat([r, pr])
        c = torch.cat([c, torch.full_like(c, -1)])   # partner's own count is looked up below if appended
    sp = d_spec[h].long()
    k = (sp << ROW_BITS) | r
    uk, inv = torch.unique(k, sorted=True, return_inverse=True)
    # lowest hypothesis index if a (spectrum, row) came from several hypotheses (index_put with duplicate indices
    # picks an arbitrary writer on CUDA); 0 such duplicates on D0/D1, so E2 outputs are unchanged
    hh = torch.full((len(uk),), 2 ** 62, dtype=torch.int64, device=dev).scatter_reduce_(0, inv, h, reduce="amin")
    cc = torch.full((len(uk),), -1, dtype=torch.int64, device=dev)
    cc.scatter_reduce_(0, inv, c, reduce="amax")
    rr = uk & (2 ** ROW_BITS - 1)
    return hh, rr, cc


def run_k2(gi: GpuIndex, h, r, d_spec, d_maxfz, peaks, p: SearchParams, shift: torch.Tensor | None = None):
    """K2 over candidates sorted by spectrum. Returns dict of device tensors (one row per candidate).
    shift (M4, optional): f32 per candidate, observed - candidate mass; scores shifted + unshifted ions (open and
    mass-offset modes). None = closed search (the E2 kernel path, bit-identical)."""
    L = K.lib()
    dev = gi.dev
    spec = d_spec[h].long()
    us, cnts = torch.unique_consecutive(spec, return_counts=True)
    off = torch.zeros(len(us) + 1, dtype=torch.int64, device=dev)
    torch.cumsum(cnts, 0, out=off[1:])
    n = len(r)
    row = r.to(torch.int32).contiguous()
    mfz = d_maxfz[h].contiguous()
    out = {"hyperscore": torch.empty(n, dtype=torch.float64, device=dev)}
    for k in ("summed_b", "summed_y", "ppm"):
        out[k] = torch.empty(n, dtype=torch.float32, device=dev)
    for k in ("matched_b", "matched_y", "longest_b", "longest_y"):
        out[k] = torch.empty(n, dtype=torch.int16, device=dev)
    sid = us.to(torch.int32).contiguous()
    if shift is None:
        K.check(L.vk_k2(len(us), K.ptr(off), K.ptr(row), K.ptr(mfz), K.ptr(sid), K.ptr(peaks["off"]),
                        K.ptr(peaks["mass"]), K.ptr(peaks["int"]), K.ptr(gi.seq_off), K.ptr(gi.length),
                        K.ptr(gi.mask), K.ptr(gi.mass32), K.ptr(gi.residues), p.fragment_ppm[0], p.fragment_ppm[1],
                        K.ptr(out["hyperscore"]), K.ptr(out["summed_b"]), K.ptr(out["summed_y"]), K.ptr(out["ppm"]),
                        K.ptr(out["matched_b"]), K.ptr(out["matched_y"]), K.ptr(out["longest_b"]),
                        K.ptr(out["longest_y"]), p.k2_threads, K.stream()), "k2")
    else:
        sh = shift.to(torch.float32).contiguous()
        K.check(L.vk_k2_shift(len(us), K.ptr(off), K.ptr(row), K.ptr(mfz), K.ptr(sid), K.ptr(peaks["off"]),
                              K.ptr(peaks["mass"]), K.ptr(peaks["int"]), K.ptr(gi.seq_off), K.ptr(gi.length),
                              K.ptr(gi.mask), K.ptr(gi.mass32), K.ptr(gi.residues), p.fragment_ppm[0],
                              p.fragment_ppm[1], K.ptr(sh), K.ptr(out["hyperscore"]), K.ptr(out["summed_b"]),
                              K.ptr(out["summed_y"]), K.ptr(out["ppm"]), K.ptr(out["matched_b"]),
                              K.ptr(out["matched_y"]), K.ptr(out["longest_b"]), K.ptr(out["longest_y"]),
                              p.k2_threads, K.stream()), "k2_shift")
    out["spec"] = spec
    return out


def lnfact(n):
    n = np.asarray(n, np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = n * np.log(n) - n + 0.5 * np.log(n) + 0.5 * np.log(2 * np.pi * n)
    return np.where(n == 0, 1.0, v)


def psm_table(sp: Spectra, gi: GpuIndex, l1, hyp, h, r, k1c, k2, scored, matched, peaks, p: SearchParams,
              report: int | None = None) -> pa.Table:
    """Rank candidates by hyperscore per spectrum (label-blind), compute Sage-style features, return top `report`."""
    report = p.report_psms if report is None else report
    dev = gi.dev
    mb = k2["matched_b"].long()
    my = k2["matched_y"].long()
    ok = (mb + my) >= p.min_matched_peaks
    idx = torch.nonzero(ok).squeeze(1)
    spec = k2["spec"][idx]
    hs = k2["hyperscore"][idx]
    rr = r[idx]
    # sort: spectrum asc, hyperscore desc, row asc (deterministic, label-free)
    o = torch.argsort(rr, stable=True)
    o = o[torch.argsort(hs[o], descending=True, stable=True)]
    o = o[torch.argsort(spec[o], stable=True)]
    idx, spec, hs = idx[o], spec[o], hs[o]
    first = torch.searchsorted(spec, spec, side="left")
    rank = torch.arange(len(spec), device=dev) - first
    nxt_same = torch.zeros_like(hs)
    if len(hs) > 1:
        same = spec[1:] == spec[:-1]
        nxt_same[:-1] = torch.where(same, hs[1:], torch.zeros_like(hs[1:]))
    best = hs[first]
    keep = rank < report
    sel = idx[keep]
    out = {k: v[sel].cpu().numpy() for k, v in k2.items() if k != "spec"}
    s = spec[keep].cpu().numpy()
    rk = rank[keep].cpu().numpy() + 1
    hsk = hs[keep].cpu().numpy()
    dn = hsk - nxt_same[keep].cpu().numpy()
    db = best[keep].cpu().numpy() - hsk
    hh = h[sel].cpu().numpy()
    rows = r[sel].cpu().numpy()
    k1 = k1c[sel].cpu().numpy()
    # per-spectrum K1 stats summed over hypotheses (Sage: scored_candidates, matched_peaks -> poisson lambda)
    n_spec = len(sp)
    sc = np.zeros(n_spec, np.int64)
    mt = np.zeros(n_spec, np.int64)
    np.add.at(sc, hyp["spec"], scored.cpu().numpy())
    np.add.at(mt, hyp["spec"], matched.cpu().numpy())
    tic = peaks["tic"].cpu().numpy()
    mbk, myk = out["matched_b"].astype(np.int64), out["matched_y"].astype(np.int64)
    kk = mbk + myk
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        lam = mt[s] / sc[s]
        pois = lam ** kk * np.exp(-lam) / np.exp(lnfact(kk))
        pois = np.where(np.isinf(pois), 1e-325, pois)
        log_pois = np.log10(pois)
    z = hyp["z"][hh].astype(np.int64)
    iso = hyp["iso"][hh].astype(np.int64)
    pm = hyp["pm"][hh]
    calc = l1.mass32[rows]
    ie = (iso.astype(np.float32) * NEUTRON).astype(np.float32)
    dm = (np.abs(pm - calc - ie) * np.float32(2e6) / (pm - ie + calc)).astype(np.float32)
    plen = l1.length[rows].astype(np.int64)
    cols = {
        "run": pa.array([sp.run] * len(s)),
        "spec_idx": s.astype(np.int64), "scan": sp.scan[s], "rt": sp.rt[s],
        "rank": rk.astype(np.int32), "pf_row": rows.astype(np.int64), "pf_id": l1.pf_id[rows],
        "charge": z.astype(np.int8), "isotope_error": iso.astype(np.int8),
        "expmass": pm.astype(np.float32), "calcmass": calc.astype(np.float32), "delta_mass_ppm": dm,
        "hyperscore": hsk, "delta_next": dn, "delta_best": db,
        "matched_b": mbk.astype(np.int16), "matched_y": myk.astype(np.int16), "matched_peaks": kk.astype(np.int16),
        "summed_intensity": (out["summed_b"] + out["summed_y"]).astype(np.float32),
        "matched_intensity_pct": (100.0 * (out["summed_b"] + out["summed_y"]) / tic[s]).astype(np.float32),
        "average_ppm": out["ppm"], "longest_b": out["longest_b"], "longest_y": out["longest_y"],
        "longest_y_pct": (out["longest_y"] / plen).astype(np.float32),
        "peptide_len": plen.astype(np.int16), "missed_cleavages": l1.missed[rows], "n_var": l1.n_var[rows],
        "k1_count": k1.astype(np.int32), "scored_candidates": sc[s].astype(np.int32), "poisson": log_pois,
    }
    return pa.table(cols)


def search(sp: Spectra, gi: GpuIndex, l1, p: SearchParams, report: int | None = None, pre_hyp=None):
    """Full closed search of one run. Returns (PSM table, stats). pre_hyp: optional (valid, hyp) from
    pre_hypotheses(); used only if `valid` equals K0's validity mask (so results are identical either way)."""
    st = Stats()
    ru0 = resource.getrusage(resource.RUSAGE_SELF)
    t0 = time.perf_counter()
    dev = gi.dev
    with Timer(st, "k0"):
        peaks = k0_preprocess(sp, p, dev, st)
    with Timer(st, "hyp"):
        if pre_hyp is not None and np.array_equal(pre_hyp[0], peaks["valid"]):
            hyp = pre_hyp[1]
            st.counts["pre_hyp_used"] = 1
        else:
            hyp = make_hypotheses(sp, peaks["valid"], gi.mass32_np, p)
            st.counts["pre_hyp_used"] = 0 if pre_hyp is None else -1   # -1: given but mismatched
    with Timer(st, "k1"):
        o_h, o_r, o_c, scored, matched, d_spec, d_maxfz = run_k1(gi, hyp, peaks, p, st)
    with Timer(st, "topk"):
        h, r, c = select_topk(o_h, o_r, o_c, d_spec, gi, p)
    with Timer(st, "k2"):
        k2 = run_k2(gi, h, r, d_spec, d_maxfz, peaks, p)
    with Timer(st, "features"):
        tbl = psm_table(sp, gi, l1, hyp, h, r, c, k2, scored, matched, peaks, p, report)
    ru1 = resource.getrusage(resource.RUSAGE_SELF)
    st.wall_s["total"] = time.perf_counter() - t0
    st.counts.update(n_spectra=len(sp), n_searched=int(peaks["valid"].sum()), n_hypotheses=len(hyp["spec"]),
                     n_rescored=int(len(r)), n_psms=tbl.num_rows)
    st.counts["cpu_core_s"] = (ru1.ru_utime - ru0.ru_utime) + (ru1.ru_stime - ru0.ru_stime)
    return tbl, st, {"peaks": peaks, "hyp": hyp, "h": h, "r": r, "k1c": c, "k2": k2}


# ------------------------------------------------------------------------------------------------ references
def window_candidates(hyp, dev, spec_subset: np.ndarray | None = None):
    """Every (hypothesis, row) in the precursor windows (the exhaustive candidate set), sorted by
    (spectrum, row). spec_subset restricts to those spectra."""
    lo, hi = hyp["lo"], hyp["hi"]
    hsel = np.arange(len(lo))
    if spec_subset is not None:
        hsel = hsel[np.isin(hyp["spec"], spec_subset)]
    n = (hi - lo)[hsel]
    hh = np.repeat(hsel, n)
    rr = np.repeat(lo[hsel], n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    s = hyp["spec"][hh]
    o = np.lexsort((rr, s))
    return torch.from_numpy(hh[o].astype(np.int64)).to(dev), torch.from_numpy(rr[o].astype(np.int64)).to(dev)


def k1_brute(gi: GpuIndex, hyp, peaks, h: torch.Tensor, r: torch.Tensor, p: SearchParams) -> torch.Tensor:
    """Brute-force preliminary counts for (hyp, row) pairs (no index): reference for K1."""
    L = K.lib()
    dev = gi.dev
    d_spec = torch.from_numpy(hyp["spec"]).to(dev)
    d_maxfz = torch.from_numpy(hyp["maxfz"]).to(dev)
    out = torch.empty(len(h), dtype=torch.int32, device=dev)
    hh = h.to(torch.int32).contiguous()
    rr = r.to(torch.int32).contiguous()
    K.check(L.vk_k1_brute(len(h), K.ptr(hh), K.ptr(rr), K.ptr(d_spec), K.ptr(d_maxfz), K.ptr(peaks["off"]),
                          K.ptr(peaks["mass"]), K.ptr(gi.seq_off), K.ptr(gi.length), K.ptr(gi.mask),
                          K.ptr(gi.mass32), K.ptr(gi.residues), p.min_ion_index, p.fragment_min_mz,
                          p.fragment_max_mz, p.fragment_ppm[0], p.fragment_ppm[1], gi.inv_logw, K.ptr(out),
                          K.stream()), "k1_brute")
    return out


def search_exhaustive(sp: Spectra, gi: GpuIndex, l1, p: SearchParams, report: int | None = None,
                      spec_subset: np.ndarray | None = None, batch: int = 20_000_000):
    """K2 over EVERY candidate in the precursor windows (no K1 shortlist): the brute-force reference."""
    st = Stats()
    dev = gi.dev
    peaks = k0_preprocess(sp, p, dev, st)
    hyp = make_hypotheses(sp, peaks["valid"], gi.mass32_np, p)
    o_h, o_r, o_c, scored, matched, d_spec, d_maxfz = run_k1(gi, hyp, peaks, p, st)   # only for poisson stats
    h, r = window_candidates(hyp, dev, spec_subset)
    # batch by whole spectra
    spec = d_spec[h].long()
    parts = []
    i = 0
    n = len(h)
    while i < n:
        j = min(n, i + batch)
        if j < n:
            j = int(torch.searchsorted(spec, spec[j - 1], side="right"))
        parts.append(run_k2(gi, h[i:j], r[i:j], d_spec, d_maxfz, peaks, p))
        i = j
    k2 = {k: torch.cat([q[k] for q in parts]) for k in parts[0]} if parts else None
    c = torch.full((len(h),), -1, dtype=torch.int64, device=dev)
    tbl = psm_table(sp, gi, l1, hyp, h, r, c, k2, scored, matched, peaks, p, report)
    return tbl, st, {"peaks": peaks, "hyp": hyp, "h": h, "r": r, "k2": k2}
