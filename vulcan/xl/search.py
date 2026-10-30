"""XL-MS search on the L1 index with K0/K1/K2 (DESIGN.md §4 L2-L3, XL-MS bullets of L4).

Two candidate-generation strategies feed one pairing/scoring back end:

* **cleavable** (DSSO): MS2 doublets (vulcan.xl.doublets) give peptide masses P_i; their precursor complements
  M_k - linker - P_i (isotope error k) are added. Two linear-like searches are run over every candidate mass window:
  each candidate row is scored by K2's localization-aware shifted scoring with the light stub (alkene) and with the
  heavy stub (thiol) as the shift; the peptide's evidence is the better of the two.
* **non-cleavable** (DSS/BS3): alpha-first open search. One open hypothesis per spectrum (window: every row with
  mass in [min_pep_mass, M - linker - min_pep_mass + iso_max*NEUTRON]); K1 counts unshifted fragments and keeps the
  top-K alpha rows per spectrum, pair-closed (engine.select_topk). For each alpha and isotope error k the beta mass
  M_k - linker - m(alpha) is looked up in the index (the window contains a target and its paired decoy, which have
  equal masses, so beta lookup is pair-closed by construction). Alpha and beta are both scored by K2 with
  shift = M_k - m(peptide) (the partner + linker rides on one residue).

Pairs: every (a, b) with m(a) + m(b) + linker within the precursor tolerance of M_k. CSM score = min(e_a, e_b)
(fragment evidence for BOTH peptides, DESIGN L4). Within a spectrum the best pair is the CSM (ties: smaller
(pf_id_lo, pf_id_hi), a content hash -- never a label).

Competition (DESIGN L4: linear, mono-link and loop-link explanations compete in the same search): the same spectra
are searched closed (linear) and with mass offsets linker+H2O, linker+NH3 (mono-links) and linker (loop-link)
via vulcan.ptm.search.search_modes; the best of those is reported next to the XL explanation and the spectrum's
winner is the higher score (`xl_wins`).

Nothing here reads labels.arrow / proteins.arrow (label blindness, DESIGN §3; tests/gpu/test_xl.py G0).
"""
from __future__ import annotations

import time
from dataclasses import asdict, dataclass

import numpy as np
import pyarrow as pa
import torch

from vulcan.gpu import engine as E
from vulcan.gpu.host import NEUTRON, PROTON, SearchParams
from vulcan.ptm import search as PS

from . import doublets as DBL
from .linkers import Linker

KIND_LINEAR, KIND_MONO, KIND_LOOP = 0, 1, 2


@dataclass
class XLParams:
    z_min: int = 3                       # XL hypotheses for precursor charge >= z_min
    isotope_errors: tuple = (0, 3)       # precursor isotope error range for the XL mass constraint
    precursor_ppm: float = 10.0          # pair mass constraint (ppm of M)
    window_ppm: float = 20.0             # peptide-mass windows of the cleavable search (fragment-derived masses)
    min_pep_mass: float = 400.0
    alpha_topk: int = 20                 # non-cleavable: K1 top-K alpha rows per spectrum (before pair closure)
    max_doublets: int = 10               # cleavable: doublets kept per spectrum
    doublet_ppm: float = 15.0
    doublet_isotopes: int = 1            # also try P - j*NEUTRON, j <= doublet_isotopes, for every doublet mass
    min_matched_per_peptide: int = 3     # a peptide without >= this many matched b/y ions gets evidence 0
    competitors: bool = True
    report: int = 3                      # XL explanations kept per spectrum (rank 1 = CSM)
    window_topn: int = 5                 # rows kept per peptide-mass window by evidence (+ their pair partners)
    max_candidates_per_chunk: int = 8_000_000   # bounds GPU memory of the window scoring

    def to_dict(self):
        return asdict(self)


def _f32(x):
    return np.asarray(x, np.float32)


def _precursor(sp, valid, xp: XLParams):
    """XL-eligible spectra: K0-valid, finite precursor, instrument charge >= z_min. Returns spec, z, pm (f64)."""
    z = sp.charge.astype(np.int64)
    s = np.flatnonzero(valid & np.isfinite(sp.precursor_mz) & (z >= xp.z_min))
    pm = (sp.precursor_mz[s] - float(PROTON)) * z[s]
    return s, z[s], pm


def _expand(lo, hi):
    """(hyp, row) for every row in [lo, hi) of every hypothesis (torch, in hypothesis order)."""
    n = (hi - lo).clamp(min=0)
    hh = torch.repeat_interleave(torch.arange(len(n), device=n.device), n)
    starts = torch.cumsum(n, 0) - n
    rr = lo[hh] + (torch.arange(len(hh), device=n.device) - starts[hh])
    return hh, rr


def score_windows(gi: E.GpuIndex, peaks, p: SearchParams, h_spec, h_maxfz, h_lo_mass, h_hi_mass, h_shift,
                  st: E.Stats, batch: int = 40_000_000, h_obs=None):
    """Score EVERY index row with mass32 in [lo, hi] of each hypothesis with K2. Shift per candidate: h_shift[h]
    (constant per hypothesis, e.g. a linker stub) or, if h_obs is given, h_obs[h] - mass32[row] (the exact
    unexplained mass, e.g. partner + linker). Returns torch (hyp, row, hyperscore, matched), matched = b + y, on
    the device. Exhaustive within the (narrow) windows: no K1 shortlist, so retrieval is trivially symmetric between
    a target and its decoy (equal masses)."""
    dev = gi.dev
    order = np.argsort(h_spec, kind="stable")             # K2 needs candidates grouped by spectrum
    m32 = gi.mass32
    lo = torch.searchsorted(m32, torch.from_numpy(_f32(h_lo_mass)[order]).to(dev), side="left")
    hi = torch.searchsorted(m32, torch.from_numpy(_f32(h_hi_mass)[order]).to(dev), side="right")
    hh_o, rr = _expand(lo, hi)
    d_order = torch.from_numpy(order).to(dev)
    hh = d_order[hh_o]
    st.counts["k2_window_candidates"] = st.counts.get("k2_window_candidates", 0) + int(len(rr))
    if len(rr) == 0:
        z = torch.zeros(0, device=dev)
        return z.long(), z.long(), z.double(), z.long()
    d_spec = torch.from_numpy(np.ascontiguousarray(h_spec, np.int32)).to(dev)
    d_maxfz = torch.from_numpy(np.ascontiguousarray(h_maxfz, np.int8)).to(dev)
    if h_obs is not None:
        cshift = (torch.from_numpy(np.asarray(h_obs, np.float64)).to(dev)[hh] - m32[rr].double()).float()
    else:
        cshift = torch.from_numpy(_f32(h_shift)).to(dev)[hh]
    spec_sorted = d_spec[hh].long()
    out_s, out_m = [], []
    i = 0
    while i < len(rr):
        j = min(len(rr), i + batch)
        if j < len(rr):
            j = int(torch.searchsorted(spec_sorted, spec_sorted[j - 1], side="right"))
        k2 = E.run_k2(gi, hh[i:j], rr[i:j], d_spec, d_maxfz, peaks, p, shift=cshift[i:j].contiguous())
        out_s.append(k2["hyperscore"])
        out_m.append(k2["matched_b"].long() + k2["matched_y"].long())
        i = j
    return hh, rr, torch.cat(out_s), torch.cat(out_m)


def _maxfz(z, p: SearchParams):
    mfc = p.max_fragment_charge
    return np.maximum(np.minimum(z, mfc + 1) if mfc else z, 2).astype(np.int8)


def _best_per_key(key, ev):
    """Index of the max-evidence entry per key (ties: first in input order). torch."""
    o = torch.argsort(-ev, stable=True)
    o = o[torch.argsort(key[o], stable=True)]
    ks = key[o]
    if len(ks) == 0:
        return o
    first = torch.ones(len(ks), dtype=torch.bool, device=ks.device)
    first[1:] = ks[1:] != ks[:-1]
    return o[first]


def window_topn(win, row, ev, partner, n, ids):
    """Boolean mask (torch) over candidates: per window the top-n rows by evidence (ties: smaller pf_id, a content
    hash), closed under target<->decoy pairing (a kept row's partner in the same window is kept too).
    Label-symmetric: `partner` is the index's symmetric pair map (device tensor), never a label."""
    dev = win.device
    if len(win) == 0:
        return torch.zeros(0, dtype=torch.bool, device=dev)
    o = torch.argsort(ids[row], stable=True)
    o = o[torch.argsort(-ev[o], stable=True)]
    o = o[torch.argsort(win[o], stable=True)]
    ws = win[o]
    first = torch.searchsorted(ws, ws, side="left")
    top = torch.zeros(len(win), dtype=torch.bool, device=dev)
    top[o[(torch.arange(len(o), device=dev) - first) < n]] = True
    nr = len(partner) + 1
    key = win * nr + row
    sk = torch.argsort(key, stable=True)
    ks = key[sk]
    want = win[top] * nr + partner[row[top]]
    pos = torch.searchsorted(ks, want).clamp(max=len(sk) - 1)
    hit = ks[pos] == want
    top[sk[pos[hit]]] = True
    return top


# ------------------------------------------------------------------------------------------------ pairing
def pair_candidates(spec, row, mass, evid, extra, pm_of_spec, linker_mass, xp: XLParams, ids, device="cuda"):
    """All (a, b) pairs within one spectrum with mass_a + mass_b + linker = pm - k*NEUTRON (+- precursor_ppm),
    k in isotope_errors. Inputs (numpy) are per scored (spectrum, row) evidence. Returns dict: inputs re-sorted by
    (spec, mass) plus pair arrays a, b (indices into them) and iso. Each unordered pair appears once
    (int64 view of pf_id(a) <= pf_id(b); a == b allowed). Computed with torch on `device`."""
    o = np.lexsort((mass, spec))
    spec, row, mass, evid = spec[o], row[o], mass[o], evid[o]
    extra = {k: v[o] for k, v in extra.items()}
    dev = torch.device(device)
    d_spec = torch.from_numpy(spec.astype(np.int64)).to(dev)
    d_mass = torch.from_numpy(mass.astype(np.float64)).to(dev)
    key = d_spec.double() * 1.0e5 + d_mass
    d_id = torch.from_numpy(np.array(ids, copy=True).view(np.int64)).to(dev)[torch.from_numpy(row.astype(np.int64)).to(dev)]
    pm = torch.from_numpy(np.asarray(pm_of_spec, np.float64)).to(dev)[d_spec]
    pa_, pb_, pk_ = [], [], []
    for k in range(xp.isotope_errors[0], xp.isotope_errors[1] + 1):
        M = pm - k * float(NEUTRON)
        need = M - linker_mass - d_mass
        tol = M * xp.precursor_ppm * 1e-6
        lo = torch.searchsorted(key, d_spec.double() * 1.0e5 + need - tol, side="left")
        hi = torch.searchsorted(key, d_spec.double() * 1.0e5 + need + tol, side="right")
        a, b = _expand(lo, hi)
        keep = d_id[a] <= d_id[b]
        pa_.append(a[keep]); pb_.append(b[keep]); pk_.append(torch.full((int(keep.sum()),), k, dtype=torch.int8, device=dev))
    a, b, k = torch.cat(pa_), torch.cat(pb_), torch.cat(pk_)
    return {"spec": spec, "row": row, "mass": mass, "evid": evid, "extra": extra, "a": a, "b": b, "iso": k}


def best_pairs(P, ids, report: int):
    """Rank pairs per spectrum by min(evid_a, evid_b) desc, then max(evid) desc, then (pf_id a, pf_id b) asc
    (int64 view; a content hash). One entry per (spectrum, unordered row pair) -- its best isotope interpretation.
    Torch on the device of P["a"]; returns numpy."""
    a, b = P["a"], P["b"]
    dev = a.device
    ev = torch.from_numpy(np.asarray(P["evid"], np.float64)).to(dev)
    row = torch.from_numpy(P["row"].astype(np.int64)).to(dev)
    spec = torch.from_numpy(P["spec"].astype(np.int64)).to(dev)
    did = torch.from_numpy(np.array(ids, copy=True).view(np.int64)).to(dev)
    ea, eb = ev[a], ev[b]
    score = torch.minimum(ea, eb)
    smax = torch.maximum(ea, eb)
    ok = score > 0
    a, b, score, smax, iso = a[ok], b[ok], score[ok], smax[ok], P["iso"][ok]
    ra, rb = row[a], row[b]
    sp = spec[a]
    o = torch.argsort(did[rb], stable=True)
    for kk in (did[ra], -smax, -score, sp):
        o = o[torch.argsort(kk[o], stable=True)]
    a, b, score, smax, sp, iso, ra, rb = (x[o] for x in (a, b, score, smax, sp, iso, ra, rb))
    nr = int(len(ids)) + 1
    assert (int(sp.max()) + 1 if len(sp) else 1) * nr * nr < 2 ** 62, "pair key overflow"
    ukey = (sp * nr + ra) * nr + rb
    ou = torch.argsort(ukey, stable=True)            # stable: the best (earliest) entry of each key comes first
    ks = ukey[ou]
    first = torch.ones(len(ks), dtype=torch.bool, device=dev)
    first[1:] = ks[1:] != ks[:-1]
    keep = torch.zeros(len(a), dtype=torch.bool, device=dev)
    keep[ou[first]] = True
    a, b, score, smax, sp, iso = (x[keep] for x in (a, b, score, smax, sp, iso))
    fs = torch.searchsorted(sp, sp, side="left")
    rank = torch.arange(len(sp), device=dev) - fs + 1
    nxt = torch.zeros_like(score)
    if len(sp) > 1:
        same = sp[1:] == sp[:-1]
        nxt[:-1] = torch.where(same, score[1:], torch.zeros_like(score[1:]))
    sel = rank <= report
    out = {"a": a[sel], "b": b[sel], "score": score[sel], "score_max": smax[sel], "spec": sp[sel], "iso": iso[sel],
           "rank": rank[sel], "delta_next": (score - nxt)[sel]}
    return {k: v.cpu().numpy() for k, v in out.items()}


# ------------------------------------------------------------------------------------------------ strategies
def _ids(l1, dev):
    """pf_id (content hash) as an int64 device tensor: the label-free tie-break order."""
    return torch.from_numpy(np.array(l1.pf_id, copy=True).view(np.int64)).to(dev)


def candidates_cleavable(sp, gi, l1, peaks, p, xp: XLParams, linker: Linker, st: E.Stats):
    """Doublet masses + complements -> two stub-shifted window searches. Returns per (spec, row) evidence."""
    with E.Timer(st, "doublets"):
        d = DBL.detect(sp, linker, tol_ppm=xp.doublet_ppm, max_per_spec=xp.max_doublets,
                       min_pep_mass=xp.min_pep_mass, z_min=xp.z_min, device=gi.dev,
                       precursor_ppm=xp.precursor_ppm, iso_range=xp.isotope_errors)
    st.counts["n_doublets"] = int(len(d["spec"]))
    st.counts["n_spectra_with_doublet"] = int(len(np.unique(d["spec"])))
    st.counts["n_spectra_with_complementary_doublets"] = int(len(np.unique(d["spec"][d["paired"]])))
    with E.Timer(st, "xl_hyp"):
        s_ok, z_ok, pm_ok = _precursor(sp, peaks["valid"], xp)
        pm_of = np.full(len(sp), np.nan)
        z_of = np.zeros(len(sp), np.int64)
        pm_of[s_ok], z_of[s_ok] = pm_ok, z_ok
        sel = np.isfinite(pm_of[d["spec"]])
        ds, dm = d["spec"][sel], d["mass"][sel]
        # doublet isotope variants: the light doublet peak may be the +j isotope (mono peak missing/weak)
        ds = np.concatenate([ds] * (xp.doublet_isotopes + 1))
        dm = np.concatenate([dm - j * float(NEUTRON) for j in range(xp.doublet_isotopes + 1)])
        # windows: doublet mass +- window_ppm of P; complement M_k - linker - P +- (precursor_ppm of M + window_ppm
        # of P). The exact pair constraint (precursor_ppm of M) is applied when pairing.
        specs, masses, widths = [ds], [dm], [dm * xp.window_ppm * 1e-6]
        for k in range(xp.isotope_errors[0], xp.isotope_errors[1] + 1):
            comp = pm_of[ds] - k * float(NEUTRON) - linker.mass - dm
            ok = comp >= xp.min_pep_mass
            specs.append(ds[ok]); masses.append(comp[ok])
            widths.append((pm_of[ds[ok]] * xp.precursor_ppm + dm[ok] * xp.window_ppm) * 1e-6)
        hs, hm, w = np.concatenate(specs), np.concatenate(masses), np.concatenate(widths)
        lo_m, hi_m = hm - w, hm + w
        light, heavy = linker.doublet
        H = len(hs)
        h_spec = np.concatenate([hs, hs]).astype(np.int64)
        h_lo = np.concatenate([lo_m, lo_m]); h_hi = np.concatenate([hi_m, hi_m])
        h_shift = np.concatenate([np.full(H, light), np.full(H, heavy)])
        h_stub = np.concatenate([np.zeros(H, np.int8), np.ones(H, np.int8)])
        h_maxfz = _maxfz(z_of[h_spec], p)
    st.counts["n_mass_windows"] = int(H)
    dev = gi.dev
    ids = _ids(l1, dev)
    win_all = np.concatenate([np.arange(H), np.arange(H)])     # the two stub hypotheses of a window share it
    parts = []
    for hsel in _spectrum_chunks(h_spec, h_lo, h_hi, gi.mass32_np, xp.max_candidates_per_chunk):
        with E.Timer(st, "k2_linear_like"):
            hh, rr, hs_, mt = score_windows(gi, peaks, p, h_spec[hsel], h_maxfz[hsel], h_lo[hsel], h_hi[hsel],
                                            h_shift[hsel], st)
        with E.Timer(st, "evidence"):
            ev = torch.where(mt >= xp.min_matched_per_peptide, hs_, torch.zeros_like(hs_))
            win = torch.from_numpy(win_all[hsel]).to(dev)[hh]
            keep = window_topn(win, rr, ev, gi.partner, xp.window_topn, ids=ids)
            hh, rr, ev, mt = hh[keep], rr[keep], ev[keep], mt[keep]
            spec = torch.from_numpy(h_spec[hsel]).to(dev)[hh]
            stub = torch.from_numpy(h_stub[hsel]).to(dev)[hh]
            # evidence per (spec, row) = max over stubs (and over overlapping windows)
            w = _best_per_key(spec * (len(gi.mass32_np) + 1) + rr, ev)
            parts.append({"spec": spec[w], "row": rr[w], "evid": ev[w], "matched": mt[w], "stub": stub[w]})
            st.counts["rows_after_window_topn"] = st.counts.get("rows_after_window_topn", 0) + int(keep.sum())
            del hh, rr, hs_, mt, ev, win, keep
    out = {k: torch.cat([q[k] for q in parts]).cpu().numpy() for k in parts[0]}
    return {"spec": out["spec"], "row": out["row"], "evid": out["evid"],
            "extra": {"matched": out["matched"], "stub": out["stub"].astype(np.int8)}, "pm_of": pm_of, "z_of": z_of}


def _spectrum_chunks(h_spec, h_lo, h_hi, mass32, max_cand):
    """Hypothesis index arrays grouped by whole spectra with <= ~max_cand window candidates each."""
    n = (np.searchsorted(mass32, _f32(h_hi), side="right") - np.searchsorted(mass32, _f32(h_lo), side="left"))
    n = np.maximum(n, 0)
    o = np.argsort(h_spec, kind="stable")
    cs = np.cumsum(n[o])
    sp_sorted = h_spec[o]
    out, i = [], 0
    while i < len(o):
        j = int(np.searchsorted(cs, (cs[i - 1] if i else 0) + max_cand, side="right"))
        j = max(j, i + 1)
        if j < len(o):
            j = int(np.searchsorted(sp_sorted, sp_sorted[j - 1], side="right"))
        out.append(np.sort(o[i:j]))
        i = j
    return out


def candidates_noncleavable(sp, gi, l1, peaks, p, xp: XLParams, linker: Linker, st: E.Stats):
    """Alpha-first open search (K1 top-K, pair-closed) + beta mass lookup; K2 scoring of both with the partner
    shift. Returns per (spec, row, shift-interpretation) evidence: rows are scored once per isotope error k."""
    dev = gi.dev
    with E.Timer(st, "xl_hyp"):
        s_ok, z_ok, pm_ok = _precursor(sp, peaks["valid"], xp)
        pm_of = np.full(len(sp), np.nan)
        z_of = np.zeros(len(sp), np.int64)
        pm_of[s_ok], z_of[s_ok] = pm_ok, z_ok
        kmax = xp.isotope_errors[1]
        lo_m = np.full(len(s_ok), xp.min_pep_mass)
        hi_m = pm_ok - linker.mass - xp.min_pep_mass + 0.02
        hyp = {"spec": s_ok.astype(np.int32), "z": z_ok.astype(np.int8), "iso": np.zeros(len(s_ok), np.int8),
               "pm": _f32(pm_ok), "M": _f32(pm_ok),
               "lo": np.searchsorted(gi.mass32_np, _f32(lo_m), side="left").astype(np.int64),
               "hi": np.searchsorted(gi.mass32_np, _f32(hi_m), side="right").astype(np.int64),
               "maxfz": _maxfz(z_ok, p), "mode": np.full(len(s_ok), PS.MODE_OPEN, np.uint8),
               "off": np.full(len(s_ok), np.nan, np.float32)}
    st.counts["alpha_window_candidates"] = int((hyp["hi"] - hyp["lo"]).sum())
    with E.Timer(st, "k1_alpha"):
        pa_ = SearchParams(**{**p.to_dict(), "topk": xp.alpha_topk})
        mp = PS.ModeParams(open_da=(0.0, 1.0), spectra_per_batch=2000)
        h, r, c, *_ = PS.k1_topk_modes(gi, hyp, peaks, pa_, mp, st)
        a_spec = hyp["spec"][h.cpu().numpy()].astype(np.int64)
        a_row = r.cpu().numpy()
    st.counts["alpha_shortlist"] = int(len(a_row))
    with E.Timer(st, "beta_lookup"):
        ks = np.arange(xp.isotope_errors[0], kmax + 1)
        ma = gi.mass32_np[a_row].astype(np.float64)
        b_spec = np.repeat(a_spec, len(ks))
        b_k = np.tile(ks, len(a_spec))
        Mk = pm_of[b_spec] - b_k * float(NEUTRON)
        need = Mk - linker.mass - np.repeat(ma, len(ks))
        tol = Mk * xp.precursor_ppm * 1e-6
        okb = need >= xp.min_pep_mass - 1.0
        # hypotheses: alpha rescoring (exact shift for each k) + beta windows
        h_spec = np.concatenate([b_spec[okb], b_spec[okb]])
        a_rep = np.repeat(a_row, len(ks))[okb]
        a_m32 = gi.mass32_np[a_rep]
        h_lo = np.concatenate([a_m32, need[okb] - tol[okb]])
        h_hi = np.concatenate([a_m32, need[okb] + tol[okb]])
        Mk_ok = Mk[okb]
        n_ok = int(okb.sum())
    with E.Timer(st, "k2_alpha_beta"):
        # alpha: window = exactly its own mass32 (also holds its paired decoy and any isobaric row: all scored)
        # shift per candidate = Mk - mass32[row] (partner + linker on one residue)
        h_maxfz = _maxfz(z_of[h_spec], p)
        hh, rr, hs_, mt = score_windows(gi, peaks, p, h_spec, h_maxfz, h_lo, h_hi, None, st,
                                        h_obs=np.concatenate([Mk_ok, Mk_ok]))
    with E.Timer(st, "evidence"):
        ev = torch.where(mt >= xp.min_matched_per_peptide, hs_, torch.zeros_like(hs_))
        # beta windows: keep the top-N rows (+ partners); alpha rescoring hypotheses (first n_ok) are kept whole
        is_beta = hh >= n_ok
        kb = window_topn(hh[is_beta], rr[is_beta], ev[is_beta], gi.partner, xp.window_topn, ids=_ids(l1, dev))
        keep = ~is_beta
        keep[torch.nonzero(is_beta).squeeze(1)[kb]] = True
        hh, rr, ev, mt = hh[keep], rr[keep], ev[keep], mt[keep]
        spec = torch.from_numpy(h_spec.astype(np.int64)).to(dev)[hh]
        hk = torch.from_numpy(np.concatenate([b_k[okb], b_k[okb]]).astype(np.int64)).to(dev)[hh]
        w = _best_per_key((spec * (len(gi.mass32_np) + 1) + rr) * 8 + hk, ev)
        out = {k: v[w].cpu().numpy() for k, v in (("spec", spec), ("row", rr), ("evid", ev), ("matched", mt),
                                                     ("iso", hk))}
    st.counts["beta_hypotheses"] = n_ok
    st.counts["rows_after_window_topn"] = int(keep.sum())
    return {"spec": out["spec"], "row": out["row"], "evid": out["evid"],
            "extra": {"matched": out["matched"], "iso": out["iso"]}, "pm_of": pm_of, "z_of": z_of}


def _pairs_noncleavable(C, gi, linker, xp, ids):
    """Pairs for the alpha-first strategy: the evidence of a row depends on the isotope interpretation k, so pair
    per k with the k-specific evidence."""
    parts = []
    for k in range(xp.isotope_errors[0], xp.isotope_errors[1] + 1):
        m = C["extra"]["iso"] == k
        if not m.any():
            continue
        xk = XLParams(**{**xp.to_dict(), "isotope_errors": (k, k)})
        Pk = pair_candidates(C["spec"][m], C["row"][m], gi.mass32_np[C["row"][m]], C["evid"][m],
                             {"matched": C["extra"]["matched"][m]}, C["pm_of"], linker.mass, xk, ids, device=gi.dev)
        parts.append(Pk)
    # merge: offset a/b indices
    off = 0
    out = {"spec": [], "row": [], "mass": [], "evid": [], "extra": {"matched": []}, "a": [], "b": [], "iso": []}
    for q in parts:
        for kk in ("spec", "row", "mass", "evid"):
            out[kk].append(q[kk])
        out["extra"]["matched"].append(q["extra"]["matched"])
        out["a"].append(q["a"] + off); out["b"].append(q["b"] + off); out["iso"].append(q["iso"])
        off += len(q["spec"])
    res = {}
    for k, v in out.items():
        if k == "extra":
            res[k] = {kk: np.concatenate(vv) for kk, vv in v.items()}
        elif k in ("a", "b", "iso"):
            res[k] = torch.cat(v) if v else torch.zeros(0, dtype=torch.int64, device=gi.dev)
        else:
            res[k] = np.concatenate(v) if v else np.zeros(0)
    return res


def search(sp, gi: E.GpuIndex, l1, p: SearchParams, xp: XLParams, linker: Linker):
    """XL search of one run. Returns (CSM table, linear-competitor PSM table or None, Stats)."""
    st = E.Stats()
    t0 = time.perf_counter()
    with E.Timer(st, "k0"):
        peaks = E.k0_preprocess(sp, p, gi.dev, st)
    ids = l1.pf_id
    if linker.cleavable:
        C = candidates_cleavable(sp, gi, l1, peaks, p, xp, linker, st)
        with E.Timer(st, "pairing"):
            P = pair_candidates(C["spec"], C["row"], gi.mass32_np[C["row"]], C["evid"], C["extra"], C["pm_of"],
                                linker.mass, xp, ids, device=gi.dev)
    else:
        C = candidates_noncleavable(sp, gi, l1, peaks, p, xp, linker, st)
        with E.Timer(st, "pairing"):
            P = _pairs_noncleavable(C, gi, linker, xp, ids)
    with E.Timer(st, "rank_pairs"):
        B = best_pairs(P, ids, xp.report)
    st.counts["n_pairs_considered"] = int(len(P["a"]))
    lin = None
    if xp.competitors:
        with E.Timer(st, "competitors"):
            mp = PS.ModeParams(offsets=tuple(linker.mono_masses) + (linker.loop_mass,))
            lin, lst = PS.search_modes(sp, gi, l1, p, mp, report=1)
        st.counts["competitor_psms"] = lin.num_rows
        for k, v in lst.wall_s.items():
            st.wall_s[f"competitors.{k}"] = v
    with E.Timer(st, "csm_table"):
        tbl = csm_table(sp, l1, P, B, C, lin, linker)
    st.wall_s["total"] = time.perf_counter() - t0
    st.counts.update(n_spectra=len(sp), n_xl_eligible=int(np.isfinite(C["pm_of"]).sum()), n_csm_rows=tbl.num_rows)
    return tbl, lin, st


def csm_table(sp, l1, P, B, C, lin, linker) -> pa.Table:
    a, b = B["a"], B["b"]
    ra, rb = P["row"][a], P["row"][b]
    s = B["spec"]
    pm = C["pm_of"][s]
    cols = {
        "run": pa.array([sp.run] * len(s)), "spec_idx": s.astype(np.int64), "scan": sp.scan[s], "rt": sp.rt[s],
        "rank": B["rank"].astype(np.int32), "charge": C["z_of"][s].astype(np.int8), "isotope_error": B["iso"],
        "expmass": pm.astype(np.float64),
        "pf_row_a": ra.astype(np.int64), "pf_row_b": rb.astype(np.int64),
        "pf_id_a": l1.pf_id[ra], "pf_id_b": l1.pf_id[rb],
        "mass_a": l1.mass32[ra], "mass_b": l1.mass32[rb],
        "score": B["score"], "score_a": P["evid"][a], "score_b": P["evid"][b],
        "matched_a": P["extra"]["matched"][a].astype(np.int16), "matched_b": P["extra"]["matched"][b].astype(np.int16),
        "delta_next": B["delta_next"],
        "delta_mass_ppm": ((pm - B["iso"] * float(NEUTRON) - linker.mass - l1.mass32[ra].astype(np.float64)
                            - l1.mass32[rb].astype(np.float64)) / pm * 1e6).astype(np.float32),
    }
    # homo-family pair: the two sides are the same peptide sequence up to the target<->decoy pairing (row a == row b,
    # or b is a's pair partner, or same sequence with different var mods). Label-free (sequences + partner map).
    # Such a pair can only win as TT or DD (min(e_t, e_d) never beats max), so TT/TD/DD cannot estimate its error
    # rate: L4 reports these separately (DESIGN L4 note, M5).
    seq_a = [l1.sequence(int(x)) for x in ra]
    seq_b = [l1.sequence(int(x)) for x in rb]
    seq_pb = [l1.sequence(int(x)) for x in l1.partner[rb]]
    cols["homo_family"] = np.array([x == y or x == z for x, y, z in zip(seq_a, seq_b, seq_pb)], bool)
    if "stub" in P["extra"]:
        cols["stub_a"] = P["extra"]["stub"][a].astype(np.int8)
        cols["stub_b"] = P["extra"]["stub"][b].astype(np.int8)
    tbl = pa.table(cols)
    # competitor (linear / mono / loop) best explanation of the same spectrum
    n = len(sp)
    lin_score = np.zeros(n)
    lin_kind = np.full(n, -1, np.int8)
    lin_row = np.full(n, -1, np.int64)
    if lin is not None and lin.num_rows:
        ls = lin.column("spec_idx").to_numpy()
        lr = lin.column("rank").to_numpy()
        m1 = lr == 1
        lin_score[ls[m1]] = lin.column("hyperscore").to_numpy()[m1]
        mode = lin.column("mode").to_numpy()[m1]
        off = lin.column("hyp_offset").to_numpy()[m1]
        kind = np.where(mode == 0, KIND_LINEAR, np.where(np.isclose(off, linker.loop_mass, atol=1e-3), KIND_LOOP,
                                                         KIND_MONO))
        lin_kind[ls[m1]] = kind
        lin_row[ls[m1]] = lin.column("pf_row").to_numpy()[m1]
    tbl = tbl.append_column("competitor_score", pa.array(lin_score[s]))
    tbl = tbl.append_column("competitor_kind", pa.array(lin_kind[s]))
    tbl = tbl.append_column("competitor_pf_row", pa.array(lin_row[s]))
    tbl = tbl.append_column("xl_wins", pa.array(B["score"] > lin_score[s]))
    return tbl


# ------------------------------------------------------------------------------------------------ exhaustive
def exhaustive_pair_count(sp, gi, peaks, xp: XLParams, linker: Linker, specs=None):
    """Size of the exhaustive non-cleavable search space: per spectrum, #unordered (a, b) row pairs with
    m(a) + m(b) + linker within tolerance of some M_k, and #distinct rows taking part. Counted on the GPU."""
    dev = gi.dev
    s_ok, z_ok, pm_ok = _precursor(sp, peaks["valid"], xp)
    if specs is not None:
        m = np.isin(s_ok, specs)
        s_ok, z_ok, pm_ok = s_ok[m], z_ok[m], pm_ok[m]
    mass = torch.from_numpy(gi.mass32_np.astype(np.float64)).to(dev)
    n_pairs = np.zeros(len(s_ok), np.int64)
    n_rows = np.zeros(len(s_ok), np.int64)
    for i in range(len(s_ok)):
        part = torch.zeros(len(mass), dtype=torch.bool, device=dev)
        tot = 0
        for k in range(xp.isotope_errors[0], xp.isotope_errors[1] + 1):
            M = pm_ok[i] - k * float(NEUTRON)
            need = M - linker.mass - mass
            tol = M * xp.precursor_ppm * 1e-6
            lo = torch.searchsorted(mass, need - tol, side="left")
            hi = torch.searchsorted(mass, need + tol, side="right")
            c = (hi - lo).clamp(min=0)
            valid = (mass >= xp.min_pep_mass) & (need >= xp.min_pep_mass - 1.0)
            c = torch.where(valid, c, torch.zeros_like(c))
            tot += int(c.sum())            # ordered pairs (a, b) incl. a == b once
            part |= c > 0
        n_pairs[i] = (tot + 0) // 2 + (tot % 2)
        n_rows[i] = int(part.sum())
    return s_ok, n_pairs, n_rows


def search_exhaustive_noncleavable(sp, gi, l1, p, xp: XLParams, linker: Linker, specs: np.ndarray, lin=None):
    """Reference: score EVERY row that takes part in a mass-matching pair (both roles), no K1 shortlist, for the
    spectra `specs`; then the same pairing/ranking as `search`. Returns (CSM table, Stats)."""
    st = E.Stats()
    peaks = E.k0_preprocess(sp, p, gi.dev, st)
    s_ok, z_ok, pm_ok = _precursor(sp, peaks["valid"], xp)
    m = np.isin(s_ok, specs)
    s_ok, z_ok, pm_ok = s_ok[m], z_ok[m], pm_ok[m]
    pm_of = np.full(len(sp), np.nan)
    z_of = np.zeros(len(sp), np.int64)
    pm_of[s_ok], z_of[s_ok] = pm_ok, z_ok
    mass = gi.mass32_np.astype(np.float64)
    hs_l, hk_l, lo_l, hi_l, sh_l = [], [], [], [], []
    for k in range(xp.isotope_errors[0], xp.isotope_errors[1] + 1):
        # every row scored once per (spectrum, k) with its own exact shift: hypotheses = (spectrum, k, row)
        Mk = pm_ok - k * float(NEUTRON)
        for i in range(len(s_ok)):
            need = Mk[i] - linker.mass - mass
            tol = Mk[i] * xp.precursor_ppm * 1e-6
            lo = np.searchsorted(mass, need - tol, side="left")
            hi = np.searchsorted(mass, need + tol, side="right")
            rows = np.flatnonzero((hi > lo) & (mass >= xp.min_pep_mass) & (need >= xp.min_pep_mass - 1.0))
            hs_l.append(np.full(len(rows), s_ok[i])); hk_l.append(np.full(len(rows), k))
            lo_l.append(gi.mass32_np[rows]); hi_l.append(gi.mass32_np[rows])
            sh_l.append(np.full(len(rows), Mk[i]))
    h_spec = np.concatenate(hs_l).astype(np.int64)
    h_k = np.concatenate(hk_l)
    h_lo, h_hi = np.concatenate(lo_l), np.concatenate(hi_l)
    h_shift = np.concatenate(sh_l)
    hh, rr, hs_, mt = score_windows(gi, peaks, p, h_spec, _maxfz(z_of[h_spec], p), h_lo, h_hi, None, st,
                                    h_obs=h_shift)
    dev = gi.dev
    ev = torch.where(mt >= xp.min_matched_per_peptide, hs_, torch.zeros_like(hs_))
    spec = torch.from_numpy(h_spec).to(dev)[hh]
    hk = torch.from_numpy(h_k.astype(np.int64)).to(dev)[hh]
    w = _best_per_key((spec * (len(mass) + 1) + rr) * 8 + hk, ev)
    spec, rr, ev, mt, hk = (x[w].cpu().numpy() for x in (spec, rr, ev, mt, hk))
    w = slice(None)
    C = {"spec": spec[w], "row": rr[w], "evid": ev[w], "extra": {"matched": mt[w], "iso": hk[w]},
         "pm_of": pm_of, "z_of": z_of}
    P = _pairs_noncleavable(C, gi, linker, xp, l1.pf_id)
    B = best_pairs(P, l1.pf_id, xp.report)
    st.counts["k2_scored"] = int(st.counts.get("k2_window_candidates", 0))
    return csm_table(sp, l1, P, B, C, lin, linker), st
