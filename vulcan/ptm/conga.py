"""CONGA-style FDR control for open / mass-offset search (DESIGN.md §4 L4; Freestone, Noble, Keich et al.,
"Group-walk, a rigorous approach to group-wise FDR analysis by target-decoy competition", Bioinformatics 2022, and
the CONGA open-search procedure, Freestone et al. 2023 J Proteome Res).

Pipeline (all stages before the final competition are label-blind):
1. Candidate PSMs: the top `top_n` PSMs per spectrum of a combined narrow + mass-offset/open search.
2. Neighbour filter (label-blind): within a spectrum, in score order, drop a PSM whose peptide shares > `frac` of its
   b/y ions (unshifted, singly charged, binned at `bin_da`) with any higher-scoring kept PSM of the same spectrum.
3. Stem-level competition: stem = stripped sequence (I->L); best PSM per (stem, side), then each target stem competes
   with its own paired decoy stem (pair key = forward stem). This is `vulcan.stats.tdc.double_competition` with
   unit = stem and every PSM its own "spectrum" (stage 1 of CONGA already chose which PSMs of a spectrum survive).
4. Groups = mass-shift bin of the pair winner (narrow bin 0 = |delta| < 0.5 Da, i.e. closed + isotope errors; else
   round(delta) Da). The partition is frozen from label-blind pair counts (label_blind_pair_counts) and merged with
   the pre-registered rule `merge_adjacent` (adjacent bins by mass, >= 200 pair competitions; narrow kept separate
   from shifted bins by ordering narrow first). q-values by Group-walk (K = 40); pooled TDC reported alongside.
"""
from __future__ import annotations

import numpy as np

from vulcan.stats import groupwalk as GW
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.tdc import double_competition, tdc_qvalues

PROTON = 1.007276467
H2O = 18.010564684


def by_ion_bins(seq: str, res_mass: np.ndarray, nterm: float, bin_da: float) -> set:
    pre = np.cumsum(res_mass)[:-1] + nterm
    tot = res_mass.sum() + nterm + H2O
    ions = np.concatenate([pre, tot - pre]) + PROTON
    return set(np.round(ions / bin_da).astype(np.int64).tolist())


def neighbour_filter(spec_key, score, ion_sets, frac: float = 0.05) -> np.ndarray:
    """Boolean keep-mask. ion_sets[i] = set of binned b/y ions of PSM i's peptide."""
    spec_key = np.asarray(spec_key)
    score = np.asarray(score, np.float64)
    keep = np.zeros(len(score), bool)
    order = np.lexsort((-score, spec_key))
    i = 0
    n = len(order)
    while i < n:
        j = i
        while j < n and spec_key[order[j]] == spec_key[order[i]]:
            j += 1
        kept = []
        for t in order[i:j]:
            a = ion_sets[t]
            ok = True
            for u in kept:
                b = ion_sets[u]
                if len(a & b) > frac * min(len(a), len(b)):
                    ok = False
                    break
            if ok:
                kept.append(t)
                keep[t] = True
        i = j
    return keep


def shift_bin(delta: np.ndarray) -> np.ndarray:
    d = np.asarray(delta, np.float64)
    return np.where(np.abs(d) < 0.5, 0, np.round(d)).astype(np.int64)


def conga(psm_id, stem, pair_stem, score, is_decoy, delta, alpha_list=(0.01, 0.05), seed: int = 11,
          min_size: int = GW.MIN_STRATUM_SIZE, K: int = 40):
    """Stem-level competition + Group-walk over mass-shift bins. Inputs are per surviving PSM (after the neighbour
    filter). Returns dict with winners, their groups, q (Group-walk) and q_pooled (plain TDC), counts."""
    stem = np.asarray(stem)
    unit = np.array([f"{s}|{int(d)}" for s, d in zip(stem, is_decoy)])
    tb = hash_tiebreak(stem, seed)
    r = double_competition(np.asarray(psm_id), unit, np.asarray(pair_stem), score, is_decoy, seed=seed,
                           tiebreak=tb)
    w = r.row
    g = shift_bin(np.asarray(delta)[w])
    counts = GW.label_blind_pair_counts(np.asarray(pair_stem)[w], g)
    shifted = sorted(k for k in counts if k != 0)
    mp = GW.merge_adjacent(shifted, counts, min_size) if shifted else {}
    grp = np.array(["narrow" if x == 0 else "s" + "_".join(map(str, mp[x])) for x in g])
    q = GW.group_walk(r.score, r.is_decoy, grp, K=K, seed=seed)
    qp = tdc_qvalues(r.score, r.is_decoy)
    out = {"row": w, "group": grp, "q": q, "q_pooled": qp, "is_decoy": r.is_decoy, "score": r.score,
           "n_groups": int(len(set(grp.tolist()))), "group_pair_counts": {str(k): int(v) for k, v in counts.items()}}
    for a in alpha_list:
        acc = (~r.is_decoy) & (q <= a)
        out[f"stems_gw_{a}"] = int(acc.sum())
        out[f"stems_pooled_{a}"] = int(((~r.is_decoy) & (qp <= a)).sum())
        out[f"offset_bins_gw_{a}"] = int(len(set(g[acc].tolist())))
    return out
