"""Family top-2 retention (DESIGN §4 L1 "family tags and family-max entrapment", §6.3 replicates).

One search over a multi-replicate index (vulcan.index.build.build_multi_decoy: T, E_s, dec_j(T), dec_j(E_s) for every
entrapment seed s and decoy set j) keeps, per spectrum, the top-2 candidates of every *retention group*, so the
(decoy seed j x entrapment seed s) competition of each replicate is reconstructed exactly from one search:

  * K1 shortlist: top-K per spectrum **per K1 family** (forward family f in {T, E_s} x side in {forward, dec_j}),
    instead of one global top-K. A separate search of replicate (j, s) with the same rule shortlists exactly the
    same rows of its four families (K1 counts are per row; nothing else competes inside a family).
  * Pair closure over every decoy set. A candidate x is *admissible* for decoy context j iff x or its set-j partner
    is in its family's top-K, i.e. exactly the rows the separate replicate search would score.
  * Canonical hypothesis per (spectrum, row): the lowest-index precursor hypothesis whose mass window contains the
    row (falls back to the closing source's hypothesis); depends only on the spectrum and the row's mass, not on the DB.
  * Retention group g = ((f * J) + j) * 2 + is_decoy_side; top-2 per (spectrum, g) by (hyperscore desc, row asc),
    the same order as engine.psm_table.

Reconstruction of replicate (j, s): merge groups {T, E_s} x {fwd, dec} in context j and take the top-2 by
(hyperscore desc, row asc). Validated against separate searches by scripts/m2b/family_validate.py.

Label use: this module reads labels.arrow ONLY to derive the opaque per-row family/side codes used for retention
(the sanctioned family-tag exception of DESIGN §4 L1). Retention is symmetric in the target/decoy roles. Features
that depend on the whole DB (scored_candidates / poisson) are NOT reconstructable this way and are not emitted.
"""
from __future__ import annotations

import numpy as np
import pyarrow as pa
import torch

from . import engine as E


class FamilyIndex:
    """Per-row family codes on the device. fwd_fam: 0 = T, 1 + i = E of the i-th entrapment seed (sorted);
    side: -1 forward, j = decoy set j. partners: (J, n) pair partner per decoy set (self where none)."""

    def __init__(self, fwd_fam: np.ndarray, side: np.ndarray, partners: list[np.ndarray], n_fwd_fam: int, dev,
                 seeds=None):
        self.J = len(partners)
        self.F = int(n_fwd_fam)
        self.seeds = seeds or {}
        self.dev = dev
        self.fwd_fam = torch.from_numpy(fwd_fam.astype(np.int64)).to(dev)
        self.side = torch.from_numpy(side.astype(np.int64)).to(dev)
        self.partners = torch.stack([torch.from_numpy(np.array(p, np.int64)) for p in partners]).to(dev)
        self.kfam = self.fwd_fam * (self.J + 1) + (self.side + 1)       # K1 family
        assert self.F * (self.J + 1) < 256 and self.F * self.J * 2 < 256

    @staticmethod
    def codes_from_labels(lab: pa.Table):
        fam = lab.column("family").to_numpy()
        seed = lab.column("entrap_seed").to_numpy()
        is_e = (fam % 2) == 1
        eseeds = sorted(set(seed[is_e].tolist()))
        fwd_fam = np.zeros(len(fam), np.int64)
        for i, s in enumerate(eseeds):
            fwd_fam[is_e & (seed == s)] = 1 + i
        if "decoy_set" in lab.column_names:
            side = lab.column("decoy_set").to_numpy().astype(np.int64)
        else:
            side = np.where(lab.column("is_decoy").to_numpy(zero_copy_only=False), 0, -1)
        return fwd_fam, side, eseeds

    @classmethod
    def from_l1(cls, l1, gi, groups=None):
        from vulcan.index.build import read_labels
        lab = read_labels(l1.path)
        assert np.array_equal(lab.column("pf_id").to_numpy(), l1.pf_id), "labels/peptidoforms row order differ"
        fwd_fam, side, eseeds = cls.codes_from_labels(lab)
        decoy_seeds = l1.meta.get("params", {}).get("decoy_seeds") or [l1.meta.get("params", {}).get("decoy_seed")]
        return cls(fwd_fam, side, l1.partners, 1 + len(eseeds), gi.dev,
                   {"entrap_seeds": eseeds, "decoy_seeds": decoy_seeds})

    def group(self, f: int, j: int, decoy: bool) -> int:
        return (f * self.J + j) * 2 + int(decoy)

    def replicate_groups(self, j: int, s_index: int) -> list[int]:
        """Groups of replicate (decoy set j, entrapment family 1 + s_index)."""
        return [self.group(f, j, d) for f in (0, 1 + s_index) for d in (False, True)]


def _stable_order(*keys):
    """argsort by keys[0] (major) ... keys[-1] (minor), all ascending, stable (torch)."""
    o = torch.argsort(keys[-1], stable=True)
    for k in reversed(keys[:-1]):
        o = o[torch.argsort(k[o], stable=True)]
    return o


def select_family(o_h, o_r, o_c, d_spec, hyp, fi: FamilyIndex, p: E.SearchParams, dev):
    """Per-family top-K, closure over all decoy sets, canonical hypothesis. Returns sorted-by-(spec,row) candidates:
    (h, r, k1count, in_topk, spec)."""
    spec = d_spec[o_h.long()].long()
    r = o_r.long()
    c = o_c.long()
    kf = fi.kfam[r]
    low = ((65535 - c) << 27) | ((2 ** 27 - 1) - r)       # count desc, row desc (engine.select_topk order)
    o = _stable_order(spec * 256 + kf, low)
    g = (spec * 256 + kf)[o]
    first = torch.searchsorted(g, g, side="left")
    rank = torch.arange(len(g), device=dev) - first
    sel = o[rank < p.topk]
    h, r, c = o_h[sel].long(), r[sel], c[sel]
    hs, rs, cs, ts = [h], [r], [c], [torch.ones_like(r)]
    for j in range(fi.J):
        pr = fi.partners[j][r]
        hs.append(h)
        rs.append(pr)
        cs.append(torch.full_like(c, -1))
        ts.append(torch.zeros_like(r))
    h, r, c, t = torch.cat(hs), torch.cat(rs), torch.cat(cs), torch.cat(ts)
    sp = d_spec[h].long()
    key = (sp << 27) | r
    uk, inv = torch.unique(key, sorted=True, return_inverse=True)
    n = len(uk)
    hsrc = torch.full((n,), 2 ** 62, dtype=torch.int64, device=dev).scatter_reduce_(0, inv, h, reduce="amin")
    cc = torch.full((n,), -1, dtype=torch.int64, device=dev).scatter_reduce_(0, inv, c, reduce="amax")
    tk = torch.zeros(n, dtype=torch.int64, device=dev).scatter_reduce_(0, inv, t, reduce="amax")
    us = uk >> 27
    ur = uk & (2 ** 27 - 1)
    # canonical hypothesis: lowest-index hypothesis of the spectrum whose row window contains the row
    hspec = torch.from_numpy(hyp["spec"].astype(np.int64)).to(dev)
    hlo = torch.from_numpy(hyp["lo"]).to(dev)
    hhi = torch.from_numpy(hyp["hi"]).to(dev)
    h0 = torch.searchsorted(hspec, us, side="left")
    h1 = torch.searchsorted(hspec, us, side="right")
    hcan = torch.full((n,), -1, dtype=torch.int64, device=dev)
    kmax = int((h1 - h0).max()) if n else 0
    for k in range(kmax - 1, -1, -1):          # descending so the lowest valid k wins
        hk = h0 + k
        ok = hk < h1
        hk = torch.where(ok, hk, torch.zeros_like(hk))
        ok &= (hlo[hk] <= ur) & (ur < hhi[hk])
        hcan = torch.where(ok, hk, hcan)
    hcan = torch.where(hcan >= 0, hcan, hsrc)
    return hcan, ur, cc, tk.bool(), us, uk


def search_family(sp: E.Spectra, gi: E.GpuIndex, l1, fi: FamilyIndex, p: E.SearchParams):
    """Closed search with family top-2 retention. Returns (table, stats). Table: one row per
    (spectrum, retention group, grank in {1, 2})."""
    st = E.Stats()
    dev = gi.dev
    with E.Timer(st, "k0"):
        peaks = E.k0_preprocess(sp, p, dev, st)
    with E.Timer(st, "hyp"):
        hyp = E.make_hypotheses(sp, peaks["valid"], gi.mass32_np, p)
    with E.Timer(st, "k1"):
        o_h, o_r, o_c, scored, matched, d_spec, d_maxfz = E.run_k1(gi, hyp, peaks, p, st)
    with E.Timer(st, "topk"):
        h, r, c, intop, spec, ukey = select_family(o_h, o_r, o_c, d_spec, hyp, fi, p, dev)
    with E.Timer(st, "k2"):
        k2 = E.run_k2(gi, h, r, d_spec, d_maxfz, peaks, p)
    with E.Timer(st, "features"):
        tbl = family_table(sp, l1, fi, hyp, h, r, intop, spec, ukey, k2, p)
    st.counts.update(n_spectra=len(sp), n_searched=int(peaks["valid"].sum()), n_hypotheses=len(hyp["spec"]),
                     n_rescored=int(len(r)), n_rows=tbl.num_rows)
    return tbl, st


def family_table(sp, l1, fi: FamilyIndex, hyp, h, r, intop, spec, ukey, k2, p) -> pa.Table:
    dev = r.device
    ok = (k2["matched_b"].long() + k2["matched_y"].long()) >= p.min_matched_peaks
    side = fi.side[r]
    ff = fi.fwd_fam[r]
    # admissibility per decoy context j: x or its set-j partner was shortlisted in its own family
    idx_l, g_l = [], []
    for j in range(fi.J):
        pj = fi.partners[j][r]
        pk = (spec << 27) | pj
        pos = torch.searchsorted(ukey, pk).clamp(max=len(ukey) - 1)
        p_in = (ukey[pos] == pk) & intop[pos]
        adm = ok & (intop | p_in)
        fwd = adm & (side < 0)
        dec = adm & (side == j)
        for m, d in ((fwd, 0), (dec, 1)):
            ii = torch.nonzero(m).squeeze(1)
            idx_l.append(ii)
            g_l.append((ff[ii] * fi.J + j) * 2 + d)
    idx = torch.cat(idx_l)
    g = torch.cat(g_l)
    hs = k2["hyperscore"][idx]
    s = spec[idx]
    rr = r[idx]
    o = _stable_order(s * 256 + g, -hs, rr)
    idx, g, hs, s, rr = idx[o], g[o], hs[o], s[o], rr[o]
    grp = s * 256 + g
    first = torch.searchsorted(grp, grp, side="left")
    rank = torch.arange(len(grp), device=dev) - first
    keep = rank < 2
    idx, g, s, rr, rank = idx[keep], g[keep], s[keep], rr[keep], rank[keep]
    hh = h[idx].cpu().numpy()
    rows = rr.cpu().numpy()
    sn = s.cpu().numpy()
    return pa.table({
        "run": pa.array([sp.run] * len(sn)),
        "spec_idx": sn.astype(np.int64), "scan": sp.scan[sn],
        "fgroup": g.cpu().numpy().astype(np.int16), "grank": (rank.cpu().numpy() + 1).astype(np.int8),
        "pf_row": rows.astype(np.int64), "pf_id": l1.pf_id[rows],
        "hyperscore": k2["hyperscore"][idx].cpu().numpy(),
        "charge": hyp["z"][hh].astype(np.int8), "isotope_error": hyp["iso"][hh].astype(np.int8),
        "matched_b": k2["matched_b"][idx].cpu().numpy(), "matched_y": k2["matched_y"][idx].cpu().numpy(),
    })


def reconstruct(tbl: pa.Table, groups: list[int], n: int = 2):
    """Replicate competition from a family table: per spectrum the top-n over `groups` by (hyperscore desc, row asc).
    Returns a dict of numpy arrays sorted by (spec_idx, rank)."""
    g = tbl.column("fgroup").to_numpy()
    m = np.isin(g, groups)
    s = tbl.column("spec_idx").to_numpy()[m]
    hs = tbl.column("hyperscore").to_numpy()[m]
    row = tbl.column("pf_row").to_numpy()[m]
    pid = tbl.column("pf_id").to_numpy()[m]
    o = np.lexsort((row, -hs, s))
    s, hs, row, pid = s[o], hs[o], row[o], pid[o]
    first = np.searchsorted(s, s, side="left")
    rank = np.arange(len(s)) - first
    k = rank < n
    return {"spec_idx": s[k], "rank": rank[k] + 1, "hyperscore": hs[k], "pf_id": pid[k], "pf_row": row[k]}


class SubIndex:
    """Row subset of an L1 index (order preserved, so mass order and row tie-breaks are preserved) with the pair
    partner of one decoy set remapped: the 'separate search' index of one replicate."""

    def __init__(self, l1, keep: np.ndarray, j: int):
        rows = np.nonzero(keep)[0]
        remap = np.full(len(l1), -1, np.int64)
        remap[rows] = np.arange(len(rows))
        self.rows = rows
        self.path = l1.path
        self.meta = l1.meta
        for k in ("pf_id", "mass32", "seq_off", "length", "mask", "n_var", "missed"):
            setattr(self, k, getattr(l1, k)[rows])
        self.residues = l1.residues
        pj = remap[l1.partners[j][rows]]
        assert (pj >= 0).all(), "partner outside the replicate subset"
        self.partner = pj
        self.partners = [pj]

    def __len__(self):
        return len(self.rows)
