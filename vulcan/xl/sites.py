"""Link-site localization for reported CSMs (L3, label-blind; DESIGN L3 "localization operator" applied to XL).

K2's shifted scoring places the partner/stub mass on whichever residue maximises the hyperscore (any residue).
For the reported CSMs we re-localize on the CPU over the linker's linkable residues only (linker.sites; '^' =
peptide N-terminus, which gives the same b/y ions as residue 0): for each candidate site, b_k carries the shift iff
k >= site and y_k (residues k+1..L-1) iff site > k; fragment charges 1..min(z-1, max_fz). The site with the most
matched ions wins, then the higher matched intensity, then K over S/T/Y/N-term, then the lower position.
For cleavable linkers the shift set is {partner + linker, light stub, heavy stub} and a site is credited with the
best of them per ion. Only peaks and sequences are used.
"""
from __future__ import annotations

import numpy as np

from vulcan.ptm.localize import RESIDUE

PROTON = 1.007276467
H2O = 18.010565


def residue_masses(seq: str, mask: int, static: dict, var: dict) -> np.ndarray:
    return np.array([RESIDUE[c] + static.get(c, 0.0) + (var.get(c, 0.0) if (mask >> k) & 1 else 0.0)
                     for k, c in enumerate(seq)])


def candidate_sites(seq: str, linker_sites: str) -> list[int]:
    s = [k for k, c in enumerate(seq) if c in linker_sites]
    if "^" in linker_sites and 0 not in s:
        s = [0] + s
    return s


def localize(seq, mask, shifts, mz, inten, z, linker_sites, static, var, tol_ppm=20.0, max_fz=3):
    """Return (site, n_matched, intensity) for one peptide in one spectrum (mz sorted ascending)."""
    rm = residue_masses(seq, mask, static, var)
    L = len(seq)
    cb = np.cumsum(rm)[:-1]                       # b_k neutral (k = 0..L-2)
    tot = rm.sum() + H2O
    cy = tot - cb                                 # y ion complementary to b_k: residues k+1..L-1
    zs = np.arange(1, max(2, min(z - 1, max_fz) + 1))
    sites = candidate_sites(seq, linker_sites)
    shifts = np.atleast_1d(np.asarray(shifts, np.float64))

    def hit(neutral):                              # neutral: (..., ) masses -> matched intensity per ion, 0 if none
        out = np.zeros(neutral.shape)
        for c in zs:
            q = (neutral + c * PROTON) / c
            lo = np.searchsorted(mz, q * (1 - tol_ppm * 1e-6), side="left")
            hi = np.searchsorted(mz, q * (1 + tol_ppm * 1e-6), side="right")
            has = hi > lo
            v = np.where(has, inten[np.minimum(lo, len(inten) - 1)], 0.0)
            out = np.maximum(out, v)
        return out
    bu, yu = hit(cb), hit(cy)
    bs = np.max([hit(cb + s) for s in shifts], axis=0)
    ys = np.max([hit(cy + s) for s in shifts], axis=0)
    best = None
    k = np.arange(L - 1)
    for st in sites:
        b = np.where(k >= st, bs, bu)
        y = np.where(st > k, ys, yu)
        n = int((b > 0).sum() + (y > 0).sum())
        it = float(b.sum() + y.sum())
        key = (n, it, 1 if seq[st] == "K" else 0, -st)
        if best is None or key > best[0]:
            best = (key, st, n, it)
    return best[1], best[2], best[3]


def localize_csms(tbl, sp, l1, linker, static: dict, var: dict, tol_ppm=20.0, max_fz=3):
    """Sites (0-based within peptide) for every row of a CSM table. Returns (site_a, site_b, n_a, n_b)."""
    n = tbl.num_rows
    sa, sb = np.zeros(n, np.int16), np.zeros(n, np.int16)
    na, nb = np.zeros(n, np.int16), np.zeros(n, np.int16)
    spec = tbl.column("spec_idx").to_numpy()
    ra, rb = tbl.column("pf_row_a").to_numpy(), tbl.column("pf_row_b").to_numpy()
    z = tbl.column("charge").to_numpy()
    pm = tbl.column("expmass").to_numpy()
    iso = tbl.column("isotope_error").to_numpy()
    for i in range(n):
        s = spec[i]
        o0, o1 = sp.offsets[s], sp.offsets[s + 1]
        mz = sp.mz[o0:o1].astype(np.float64)
        it = sp.intensity[o0:o1].astype(np.float64)
        o = np.argsort(mz)
        mz, it = mz[o], it[o]
        M = pm[i] - iso[i] * 1.00335
        for which, r, sarr, narr in ((0, ra[i], sa, na), (1, rb[i], sb, nb)):
            seq = l1.sequence(int(r))
            shifts = [M - float(l1.mass32[r])]
            if linker.cleavable:
                shifts += list(linker.doublet)
            st, nm, _ = localize(seq, int(l1.mask[r]), shifts, mz, it, int(z[i]), linker.sites, static, var,
                                 tol_ppm, max_fz)
            sarr[i], narr[i] = st, nm
    return sa, sb, na, nb
