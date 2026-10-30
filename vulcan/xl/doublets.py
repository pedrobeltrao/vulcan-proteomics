"""MS2 doublet detection for MS-cleavable linkers (DSSO: alkene/thiol stubs, delta 31.97208 Da).

For every MS2 spectrum and fragment charge c < precursor charge, find peak pairs (i, j) with
mz_j - mz_i = delta / c (+- tol ppm of mz_i). For c >= 2 the lighter peak must have its +1 isotope
(mz_i + 1.00335 / c) present. Each doublet gives the unmodified peptide mass
    P = c * (mz_i - proton) - light_stub.
Vectorised over all spectra (one global sorted f64 key, torch on the search device). Label-free (L2): only peaks
and precursor masses are used.
"""
from __future__ import annotations

import numpy as np
import torch

PROTON = 1.0072764
NEUTRON = 1.00335


def detect(sp, linker, tol_ppm: float = 15.0, top_peaks: int = 300, max_per_spec: int = 10,
           min_pep_mass: float = 400.0, z_min: int = 3, require_iso: bool = True, device="cuda",
           precursor_ppm: float = 10.0, iso_range: tuple = (0, 3)) -> dict:
    """Return dict of numpy arrays (spec, mass, c, inten, paired) for up to `max_per_spec` doublets per spectrum,
    ranked by (has a complementary doublet, min(intensity_i, intensity_j)) desc (ties: lower mass). Spectra with
    precursor charge < z_min are skipped."""
    light, heavy = linker.doublet
    delta = heavy - light
    dev = torch.device(device)
    n = len(sp)
    lens = torch.from_numpy(np.diff(sp.offsets)).to(dev)
    spec_of = torch.repeat_interleave(torch.arange(n, device=dev), lens)
    mz = torch.from_numpy(sp.mz).to(dev).double()
    it = torch.from_numpy(sp.intensity).to(dev).double()
    z = torch.from_numpy(sp.charge.astype(np.int64)).to(dev)
    M = (torch.from_numpy(sp.precursor_mz).to(dev).double() - PROTON) * z
    ok = torch.isfinite(mz) & (z[spec_of] >= z_min)
    # top `top_peaks` per spectrum by intensity (stable: ties keep m/z order)
    o = torch.argsort(-it, stable=True)
    o = o[torch.argsort(spec_of[o], stable=True)]
    first = torch.searchsorted(spec_of[o], spec_of[o], side="left")
    rank = torch.arange(len(o), device=dev) - first
    keep = torch.zeros(len(mz), dtype=torch.bool, device=dev)
    keep[o[rank < top_peaks]] = True
    keep &= ok
    s, m, i = spec_of[keep], mz[keep], it[keep]
    key = s.double() * 1.0e5 + m              # store order is m/z ascending within spectrum -> key sorted
    oo = torch.argsort(key, stable=True)
    s, m, i, key = s[oo], m[oo], i[oo], key[oo]
    out = {"spec": [], "mass": [], "c": [], "inten": []}
    zmax = int(z.max()) if n else 0
    for c in range(1, max(zmax, 1)):
        cm = z[s] > c
        d = delta / c
        tol = m * tol_ppm * 1e-6
        lo = torch.searchsorted(key, key + d - tol, side="left")
        hi = torch.searchsorted(key, key + d + tol, side="right")
        has = (hi > lo) & cm
        if require_iso and c >= 2:
            ilo = torch.searchsorted(key, key + NEUTRON / c - tol, side="left")
            ihi = torch.searchsorted(key, key + NEUTRON / c + tol, side="right")
            has &= ihi > ilo
        idx = torch.nonzero(has).squeeze(1)
        if len(idx) == 0:
            continue
        cnt = hi[idx] - lo[idx]
        starts = torch.cumsum(cnt, 0) - cnt
        jj = torch.repeat_interleave(lo[idx], cnt) + (torch.arange(int(cnt.sum()), device=dev)
                                                      - torch.repeat_interleave(starts, cnt))
        seg = torch.repeat_interleave(torch.arange(len(idx), device=dev), cnt)
        best = torch.zeros(len(idx), dtype=torch.float64, device=dev).scatter_reduce_(0, seg, i[jj], "amax")
        P = c * (m[idx] - PROTON) - light
        sc = torch.minimum(i[idx], best)
        good = (P >= min_pep_mass) & (P <= M[s[idx]] - linker.mass - min_pep_mass + 1.0)
        out["spec"].append(s[idx][good]); out["mass"].append(P[good])
        out["c"].append(torch.full((int(good.sum()),), c, device=dev)); out["inten"].append(sc[good])
    if not out["spec"]:
        return {k: np.zeros(0) for k in list(out) + ["paired"]}
    r = {k: torch.cat(v) for k, v in out.items()}
    # complementary doublets: P_i + P_j + linker = M - k*NEUTRON (k in iso_range, tol precursor_ppm of M; j may be
    # i itself = homodimer). Doublets of fragment ions that carry a stub also show the 31.97 Da spacing, so the
    # precursor-level doublets are ranked first by having a complement, then by intensity.
    o = torch.argsort(r["mass"], stable=True)
    o = o[torch.argsort(r["spec"][o], stable=True)]
    r = {k: v[o] for k, v in r.items()}
    key = r["spec"].double() * 1.0e5 + r["mass"]
    paired = torch.zeros(len(key), dtype=torch.bool, device=dev)
    for k in range(iso_range[0], iso_range[1] + 1):
        Mk = M[r["spec"]] - k * NEUTRON
        need = Mk - linker.mass - r["mass"]
        tol = Mk * precursor_ppm * 1e-6
        lo = torch.searchsorted(key, r["spec"].double() * 1.0e5 + need - tol, side="left")
        hi = torch.searchsorted(key, r["spec"].double() * 1.0e5 + need + tol, side="right")
        paired |= hi > lo
    r["paired"] = paired
    o = torch.argsort(r["mass"], stable=True)
    o = o[torch.argsort(-r["inten"][o], stable=True)]
    o = o[torch.argsort((~r["paired"][o]).to(torch.int8), stable=True)]
    o = o[torch.argsort(r["spec"][o], stable=True)]
    r = {k: v[o] for k, v in r.items()}
    first = torch.searchsorted(r["spec"], r["spec"], side="left")
    sel = (torch.arange(len(o), device=dev) - first) < max_per_spec
    return {"spec": r["spec"][sel].cpu().numpy().astype(np.int64), "mass": r["mass"][sel].cpu().numpy(),
            "c": r["c"][sel].cpu().numpy().astype(np.int8), "inten": r["inten"][sel].cpu().numpy(),
            "paired": r["paired"][sel].cpu().numpy()}
