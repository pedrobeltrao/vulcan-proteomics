"""M4 search modes on the existing L1 index + K0/K1/K2 (DESIGN.md §4 L2 "Modes"): closed (incl. variable-mod and
phospho indexes), mass-offset (Unimod list + pilot histogram) and open (wide precursor window), with shifted +
unshifted ions in K2.

Hypotheses
----------
Closed hypotheses are exactly `vulcan.gpu.host.make_hypotheses` (spectrum x charge x isotope error), listed first.
Mass-offset mode adds, per spectrum x charge, one hypothesis per offset o (isotope error 0): candidate rows with
Sage f32 mass in M - o +- precursor ppm. Open mode adds one hypothesis per spectrum x charge whose window is
calc in [M - open_da[1], M - open_da[0]] (delta = M - calc in open_da). Every hypothesis of a spectrum competes in
one top-K (K1 count) and one ranking (K2 hyperscore), label-blind; a (spectrum, row) reached by several hypotheses
keeps the lowest hypothesis index (closed first, then offsets in list order, then open).

K2 gets shift = M[h] - calc for offset/open hypotheses (0 for closed) and scores shifted + unshifted b/y.

Open-mode K1 windows are huge, so K1 runs in spectrum batches with a bounded append buffer (retried with a
smaller batch on overflow); the inverted-index kernel itself is unchanged.

Added PSM columns (PSM_SCHEMA.md, M4): `mode` (u8: 0 closed, 1 mass-offset, 2 open), `hyp_offset` (f32, nominal
offset of the hypothesis, Da; 0 closed, NaN open), `mass_shift` (f32, Da: expmass - isotope_error*NEUTRON - calcmass).
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import pyarrow as pa
import torch

from vulcan.gpu import _lib as K
from vulcan.gpu import engine as E
from vulcan.gpu.host import NEUTRON, PROTON, SearchParams, make_hypotheses

MODE_CLOSED, MODE_OFFSET, MODE_OPEN = 0, 1, 2


@dataclass
class ModeParams:
    offsets: tuple = ()                 # mass-offset mode: Da (delta = observed - peptide), each its own hypothesis
    open_da: tuple | None = None        # open mode: (lo, hi) for delta = observed - peptide, e.g. (-150, 500)
    open_tile: int = 4096               # K1 tile for open hypotheses (shared mem 4 warps x 16 KB)
    k1_threshold_open: int = 1
    spectra_per_batch: int = 1000      # open mode only (closed/offset: one batch)
    append_cap: int = 600_000_000       # K1 append buffer (entries) per batch


def make_mode_hypotheses(sp, valid, idx_mass32, p: SearchParams, mp: ModeParams):
    hyp = make_hypotheses(sp, valid, idx_mass32, p)
    n0 = len(hyp["spec"])
    hyp["mode"] = np.zeros(n0, np.uint8)
    hyp["off"] = np.zeros(n0, np.float32)
    parts = [hyp]
    # base (spectrum, charge) rows at isotope 0
    base = hyp["iso"] == 0
    s, z = hyp["spec"][base], hyp["z"][base]
    M0 = hyp["M"][base]
    pm0 = hyp["pm"][base]
    mfz = hyp["maxfz"][base]
    plo, phi = np.float32(p.precursor_ppm[0]), np.float32(p.precursor_ppm[1])
    for o in mp.offsets:
        Mo = (M0 - np.float32(o)).astype(np.float32)
        lo = (Mo + (Mo * plo) / np.float32(1e6)).astype(np.float32)
        hi = (Mo + (Mo * phi) / np.float32(1e6)).astype(np.float32)
        parts.append({"spec": s, "z": z, "iso": np.zeros_like(z), "pm": pm0, "M": M0,
                      "lo": np.searchsorted(idx_mass32, lo, side="left").astype(np.int64),
                      "hi": np.searchsorted(idx_mass32, hi, side="right").astype(np.int64), "maxfz": mfz,
                      "mode": np.full(len(s), MODE_OFFSET, np.uint8), "off": np.full(len(s), o, np.float32)})
    if mp.open_da is not None:
        lo = (M0 - np.float32(mp.open_da[1])).astype(np.float32)
        hi = (M0 - np.float32(mp.open_da[0])).astype(np.float32)
        parts.append({"spec": s, "z": z, "iso": np.zeros_like(z), "pm": pm0, "M": M0,
                      "lo": np.searchsorted(idx_mass32, lo, side="left").astype(np.int64),
                      "hi": np.searchsorted(idx_mass32, hi, side="right").astype(np.int64), "maxfz": mfz,
                      "mode": np.full(len(s), MODE_OPEN, np.uint8), "off": np.full(len(s), np.nan, np.float32)})
    out = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    return out


def _sub(hyp, idx):
    return {k: v[idx] for k, v in hyp.items()}


def _k1_launch(gi, hyp_lo, hyp_hi, hh, tile, thr, d_spec, d_maxfz, peaks, p, cap, scored, matched, st):
    """One K1 pass over hypotheses `hh` (global indices). Returns (o_h, o_r, o_c, n) or None on overflow."""
    L = K.lib()
    dev = gi.dev
    pt = SearchParams(**{**p.to_dict(), "tile": tile})
    wi = E.work_items({"lo": hyp_lo[hh], "hi": hyp_hi[hh]}, pt)
    cap = int(min(max(1, wi["n_window"]), cap))
    n_out = torch.zeros(1, dtype=torch.int64, device=dev)
    o_h = torch.empty(cap, dtype=torch.int32, device=dev)
    o_r = torch.empty(cap, dtype=torch.int32, device=dev)
    o_c = torch.empty(cap, dtype=torch.int32, device=dev)
    ghyp = hh[wi["hyp"]].astype(np.int32)                   # work items -> global hypothesis index
    for sl in gi.slices:
        mm = (wi["chunk"] >= sl.chunk0) & (wi["chunk"] < sl.chunk1)
        if not mm.any():
            continue
        it = {"hyp": torch.from_numpy(np.ascontiguousarray(ghyp[mm])).to(dev)}
        for k in ("chunk", "lo", "hi"):
            it[k] = torch.from_numpy(np.ascontiguousarray(wi[k][mm])).to(dev)
        st.h2d_bytes += sum(v.numel() * 4 for v in it.values())
        K.check(L.vk_k1(int(mm.sum()), K.ptr(it["hyp"]), K.ptr(it["chunk"]), K.ptr(it["lo"]), K.ptr(it["hi"]),
                        K.ptr(d_spec), K.ptr(d_maxfz), K.ptr(peaks["off"]), K.ptr(peaks["mass"]),
                        K.ptr(sl.offsets), K.ptr(sl.entries), gi.n_bins, sl.chunk0, p.chunk_size, tile,
                        p.fragment_ppm[0], p.fragment_ppm[1], p.fragment_min_mz, gi.inv_logw, thr, cap,
                        K.ptr(n_out), K.ptr(o_h), K.ptr(o_r), K.ptr(o_c), K.ptr(scored), K.ptr(matched),
                        K.stream()), "k1")
    n = int(n_out.item())
    st.counts["k1_work_items"] = st.counts.get("k1_work_items", 0) + len(wi["hyp"])
    if n > cap:
        return None
    st.counts["k1_window_candidates"] = st.counts.get("k1_window_candidates", 0) + wi["n_window"]
    st.counts["k1_appended"] = st.counts.get("k1_appended", 0) + n
    return o_h[:n], o_r[:n], o_c[:n]


def k1_topk_modes(gi: E.GpuIndex, hyp, peaks, p: SearchParams, mp: ModeParams, st: E.Stats):
    """K1 + per-spectrum top-K + pair closure over all hypotheses, in spectrum batches (all hypotheses of a spectrum
    are in the same batch, so the top-K is exactly the one-shot top-K). Closed/offset hypotheses use p.tile, open
    ones mp.open_tile. On append-buffer overflow the batch is halved and redone (its counters reset first).
    Returns (h, r, c) like engine.select_topk plus per-hypothesis scored/matched counters and device arrays."""
    dev = gi.dev
    n_h = len(hyp["spec"])
    d_spec = torch.from_numpy(hyp["spec"]).to(dev)
    d_maxfz = torch.from_numpy(hyp["maxfz"]).to(dev)
    scored = torch.zeros(n_h, dtype=torch.int32, device=dev)
    matched = torch.zeros(n_h, dtype=torch.int64, device=dev)
    is_open = hyp["mode"] == MODE_OPEN
    order = np.argsort(hyp["spec"], kind="stable")
    sp_sorted = hyp["spec"][order]
    us = np.unique(sp_sorted)
    nb = mp.spectra_per_batch if is_open.any() else len(us)
    outs = []
    i = 0
    while i < len(us):
        a = np.searchsorted(sp_sorted, us[i], side="left")
        b = np.searchsorted(sp_sorted, us[min(i + nb, len(us)) - 1], side="right")
        hh_all = order[a:b]
        res = []
        ok = True
        for sel, tile, thr in ((hh_all[~is_open[hh_all]], p.tile, p.k1_threshold),
                               (hh_all[is_open[hh_all]], mp.open_tile, mp.k1_threshold_open)):
            if len(sel) == 0:
                continue
            r = _k1_launch(gi, hyp["lo"], hyp["hi"], np.sort(sel), tile, thr, d_spec, d_maxfz, peaks, p,
                           mp.append_cap, scored, matched, st)
            if r is None:
                ok = False
                break
            res.append(r)
        if not ok:
            th = torch.from_numpy(hh_all.astype(np.int64)).to(dev)
            scored[th] = 0
            matched[th] = 0
            if nb == 1:
                raise RuntimeError("K1 append buffer overflow for a single spectrum")
            nb = max(1, nb // 2)
            st.counts["k1_overflow_retries"] = st.counts.get("k1_overflow_retries", 0) + 1
            continue
        if res:
            o_h = torch.cat([x[0] for x in res]); o_r = torch.cat([x[1] for x in res])
            o_c = torch.cat([x[2] for x in res])
            outs.append(E.select_topk(o_h, o_r, o_c, d_spec, gi, p))
            del o_h, o_r, o_c, res
        i += nb
    if not outs:
        z = torch.zeros(0, dtype=torch.int64, device=dev)
        return z, z, z, scored, matched, d_spec, d_maxfz
    h = torch.cat([x[0] for x in outs]); r = torch.cat([x[1] for x in outs]); c = torch.cat([x[2] for x in outs])
    # global order by (spectrum, row), as engine.select_topk returns
    key = (d_spec[h].long() << E.ROW_BITS) | r
    o = torch.argsort(key)
    return h[o], r[o], c[o], scored, matched, d_spec, d_maxfz


def search_modes(sp, gi: E.GpuIndex, l1, p: SearchParams, mp: ModeParams, report: int | None = None):
    """Closed + mass-offset + open search of one run (modes per `mp`). Returns (PSM table, stats)."""
    st = E.Stats()
    t0 = time.perf_counter()
    dev = gi.dev
    with E.Timer(st, "k0"):
        peaks = E.k0_preprocess(sp, p, dev, st)
    with E.Timer(st, "hyp"):
        hyp = make_mode_hypotheses(sp, peaks["valid"], gi.mass32_np, p, mp)
    with E.Timer(st, "k1_topk"):
        h, r, c, scored, matched, d_spec, d_maxfz = k1_topk_modes(gi, hyp, peaks, p, mp, st)
    with E.Timer(st, "k2"):
        mode = torch.from_numpy(hyp["mode"]).to(dev)[h]
        M = torch.from_numpy(hyp["M"]).to(dev)[h]
        shift = torch.where(mode == MODE_CLOSED, torch.zeros_like(M), M - gi.mass32[r])
        parts, i, n = [], 0, len(h)
        spec = d_spec[h].long()
        while i < n:                                       # K2 in spectrum-contiguous batches
            j = min(n, i + 50_000_000)
            if j < n:
                j = int(torch.searchsorted(spec, spec[j - 1], side="right"))
            parts.append(E.run_k2(gi, h[i:j], r[i:j], d_spec, d_maxfz, peaks, p, shift=shift[i:j]))
            i = j
        k2 = {k: torch.cat([q[k] for q in parts]) for k in parts[0]}
    with E.Timer(st, "features"):
        tbl = E.psm_table(sp, gi, l1, hyp, h, r, c, k2, scored, matched, peaks, p, report)
        hh = _psm_hyp_index(tbl, sp, hyp, h, r, k2, p, report)
        rows = tbl.column("pf_row").to_numpy()
        iso = tbl.column("isotope_error").to_numpy().astype(np.float32)
        pm = tbl.column("expmass").to_numpy()
        ms = (pm - (iso * NEUTRON).astype(np.float32) - l1.mass32[rows]).astype(np.float32)
        tbl = tbl.append_column("mode", pa.array(hyp["mode"][hh], pa.uint8()))
        tbl = tbl.append_column("hyp_offset", pa.array(hyp["off"][hh], pa.float32()))
        tbl = tbl.append_column("mass_shift", pa.array(ms, pa.float32()))
    st.wall_s["total"] = time.perf_counter() - t0
    st.counts.update(n_spectra=len(sp), n_searched=int(peaks["valid"].sum()), n_hypotheses=len(hyp["spec"]),
                     n_rescored=int(len(r)), n_psms=tbl.num_rows)
    return tbl, st


def _psm_hyp_index(tbl, sp, hyp, h, r, k2, p, report):
    """Hypothesis index of each reported PSM, recomputed with psm_table's own ordering (spectrum asc, hyperscore
    desc, row asc; rank < report) so it aligns row by row with `tbl`."""
    report = p.report_psms if report is None else report
    mb, my = k2["matched_b"].long(), k2["matched_y"].long()
    idx = torch.nonzero((mb + my) >= p.min_matched_peaks).squeeze(1)
    spec, hs, rr = k2["spec"][idx], k2["hyperscore"][idx], r[idx]
    o = torch.argsort(rr, stable=True)
    o = o[torch.argsort(hs[o], descending=True, stable=True)]
    o = o[torch.argsort(spec[o], stable=True)]
    idx, spec = idx[o], spec[o]
    first = torch.searchsorted(spec, spec, side="left")
    rank = torch.arange(len(spec), device=spec.device) - first
    sel = idx[rank < report]
    hh = h[sel].cpu().numpy()
    assert len(hh) == tbl.num_rows and np.array_equal(r[sel].cpu().numpy(), tbl.column("pf_row").to_numpy())
    return hh
