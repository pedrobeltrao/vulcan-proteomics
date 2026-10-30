"""Decoy-amino-acid false localization rate (DESIGN.md §6.5).

Reference: Ramsbottom, Prakash, Perez-Riverol, ..., Jones. "Method for
independent estimation of the false localization rate for phosphoproteomics."
J Proteome Res 21:1603 (2022), doi:10.1021/acs.jproteome.1c00827. Alanine and
leucine were the most reliable decoy residues. Aggregation: Camacho et al. JPR
2023, doi:10.1021/acs.jproteome.2c00823.

Estimator (pre-registered):

    FLR = n_pA * (N_STY / N_A) / n_pSTY

* n_pA   = number of reported sites localized to the decoy residue (A primary;
           L secondary, with I counted as L),
* n_pSTY = number of reported sites localized to S, T or Y,
* N_STY, N_A = number of candidate acceptor residues (S+T+Y, resp. the decoy
           residue) summed over the localized peptides in the reported set.

Rationale: a random (false) localization picks an acceptor roughly uniformly among
the candidates, so false STY sites ~= n_pA * N_STY / N_A. Requires localization
scoring blind to residue identity given the masses (DESIGN §4 L3).

Reporting strata (separately): mono- vs multi-phospho peptides, and sites that are
<= 2 residues from another candidate acceptor of the same peptide vs not.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = ["acceptor_counts", "decoy_aa_flr", "flr_curve", "flr_strata", "STY"]

STY = frozenset("STY")


def _norm(seq: str, decoy_residue: str) -> str:
    s = seq.upper()
    if decoy_residue == "L":
        s = s.replace("I", "L")
    return s


def acceptor_counts(peptides: Sequence[str], decoy_residue: str = "A") -> tuple[int, int]:
    """(N_STY, N_decoy) summed over the peptide sequences (I counted as L for L)."""
    n_sty = 0
    n_dec = 0
    for p in peptides:
        s = _norm(p, decoy_residue)
        n_sty += sum(s.count(c) for c in "STY")
        n_dec += s.count(decoy_residue)
    return n_sty, n_dec


def decoy_aa_flr(site_residue: Sequence[str], site_peptide_idx: Sequence[int], peptides: Sequence[str],
                 decoy_residue: str = "A") -> dict:
    """FLR estimate for a reported set of sites.

    Parameters
    ----------
    site_residue : residue each reported site was localized to ('S','T','Y' or the
        decoy residue; 'I' is read as 'L' when decoy_residue == 'L').
    site_peptide_idx : index into ``peptides`` of the peptide carrying the site.
    peptides : sequences of the localized peptides (each counted once).
    """
    res = np.array([_norm(r, decoy_residue) for r in site_residue]) if len(site_residue) else np.array([], dtype=str)
    used = sorted(set(int(i) for i in site_peptide_idx))
    n_sty, n_dec = acceptor_counts([peptides[i] for i in used], decoy_residue)
    n_pA = int(np.count_nonzero(res == decoy_residue))
    n_pSTY = int(np.count_nonzero(np.isin(res, list(STY))))
    if n_pSTY == 0:
        flr = 0.0 if n_pA == 0 else float("inf")
    elif n_dec == 0:
        flr = float("nan")
    else:
        flr = n_pA * (n_sty / n_dec) / n_pSTY
    return {"flr": flr, "n_pA": n_pA, "n_pSTY": n_pSTY, "N_STY": n_sty, "N_A": n_dec,
            "decoy_residue": decoy_residue}


def flr_curve(loc_score, site_residue, site_peptide_idx, peptides, decoy_residue: str = "A",
              monotone: bool = True) -> np.ndarray:
    """FLR estimate at each site's own localization-score threshold (sites with
    score >= s); optionally made monotone like a q-value. Peptide acceptor counts
    are accumulated over the distinct peptides present at each threshold.
    Returns values in the input order (decoy-residue sites get values too)."""
    sc = np.asarray(loc_score, float)
    n = sc.shape[0]
    res = np.array([_norm(r, decoy_residue) for r in site_residue])
    pidx = np.asarray(site_peptide_idx, dtype=np.int64)
    order = np.argsort(-sc, kind="stable")
    seen = set()
    n_sty = n_dec = n_pA = n_pSTY = 0
    raw = np.empty(n)
    vals = []
    for k, i in enumerate(order):
        if pidx[i] not in seen:
            seen.add(pidx[i])
            a, b = acceptor_counts([peptides[pidx[i]]], decoy_residue)
            n_sty += a
            n_dec += b
        n_pA += res[i] == decoy_residue
        n_pSTY += res[i] in STY
        vals.append((n_pA, n_pSTY, n_sty, n_dec))
    # evaluate at the end of tie runs
    ss = sc[order]
    f = np.empty(n)
    k = n - 1
    while k >= 0:
        j = k
        while j > 0 and ss[j - 1] == ss[k]:
            j -= 1
        a, b, c, d = vals[k]
        v = (a * (c / d) / b) if (b > 0 and d > 0) else (0.0 if a == 0 else 1.0)
        f[j:k + 1] = min(v, 1.0)
        k = j - 1
    if monotone:
        f = np.minimum.accumulate(f[::-1])[::-1]
    raw[order] = f
    return raw


def flr_strata(n_phospho: Sequence[int], site_pos: Sequence[int], peptides: Sequence[str],
               site_peptide_idx: Sequence[int], decoy_residue: str = "A", max_gap: int = 2) -> dict:
    """Stratum label per site: 'mono'/'multi' and whether another candidate
    acceptor (S/T/Y or the decoy residue) lies within ``max_gap`` residues."""
    out_mm = np.where(np.asarray(n_phospho) > 1, "multi", "mono")
    close = np.zeros(len(site_pos), dtype=bool)
    acc = set("STY") | {decoy_residue}
    for k, (pos, pi) in enumerate(zip(site_pos, site_peptide_idx)):
        s = _norm(peptides[pi], decoy_residue)
        lo, hi = max(0, pos - max_gap), min(len(s), pos + max_gap + 1)
        close[k] = any(s[j] in acc for j in range(lo, hi) if j != pos)
    return {"phospho_class": out_mm, "close_site": close}
