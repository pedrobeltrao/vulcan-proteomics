"""Cross-link FDR (DESIGN.md §4 L4, XL-MS).

References
----------
* Walzthoeni et al. Nat Methods 9:901 (2012): FDR = (TD - DD) / TT.
* Fischer & Rappsilber, Anal Chem 89:3829 (2017), xiFDR: FDR at the level of the
  claim (CSM -> peptide pair -> residue pair -> PPI), merging before estimating.
* Lenz et al. Nat Commun 12:3564 (2021): self and heteromeric links separately.

Implemented
-----------
At a score threshold s over the entities of one level and one link class:

    FDR_TDDD(s) = max(TD(s) - DD(s), 0) / max(TT(s), 1)     (clipped at 0)
    FDR_TD(s)   = TD(s) / max(TT(s), 1)                      (conservative)

TT/TD/DD = number of target-target, target-decoy, decoy-decoy entities with score
>= s. q-values are monotone (min over lower thresholds), capped at 1.
Optional ``plus_one`` adds 1 to the numerator (not used by xiFDR; no finite-sample
guarantee is claimed for either form). ``min_tt``: a threshold is only eligible
if it accepts at least ``min_tt`` TT entities (placeholder for the xiFDR-style
minimum-count rule, which DESIGN does not specify numerically; default 0).

Merging: entities are keyed per level, the entity score is the best member score.
Self vs heteromeric: a link is *self* when both sides map to the same base
protein (a target protein and its decoy share a base, so decoy-A x A is self).
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .tdc import _clean_scores

__all__ = ["xl_class", "is_self_link", "xl_qvalues", "merge_level", "xl_level_fdr", "LEVELS"]

LEVELS = ("csm", "peptide_pair", "residue_pair", "ppi")


def xl_class(decoy_a, decoy_b) -> np.ndarray:
    """Number of decoy sides: 0 = TT, 1 = TD, 2 = DD."""
    return np.asarray(decoy_a, bool).astype(np.int8) + np.asarray(decoy_b, bool).astype(np.int8)


def is_self_link(prot_a, prot_b, base_of: Mapping | None = None) -> np.ndarray:
    pa, pb = np.asarray(prot_a), np.asarray(prot_b)
    if base_of is not None:
        pa = np.array([base_of.get(x, x) for x in pa])
        pb = np.array([base_of.get(x, x) for x in pb])
    return pa == pb


def xl_qvalues(score, n_decoy, plus_one: bool = False, min_tt: int = 0):
    """Return (q_tddd, q_td) for every entity (decoy-containing entities included)."""
    s = _clean_scores(score)
    c = np.asarray(n_decoy, dtype=np.int8)
    n = s.shape[0]
    if n == 0:
        return np.zeros(0), np.zeros(0)
    order = np.argsort(-s, kind="stable")
    ss, cc = s[order], c[order]
    TT = np.cumsum(cc == 0)
    TD = np.cumsum(cc == 1)
    DD = np.cumsum(cc == 2)
    # counts at the end of tie runs
    tstart = np.r_[True, ss[1:] != ss[:-1]]
    tid = np.cumsum(tstart) - 1
    ends = np.r_[np.flatnonzero(tstart)[1:] - 1, n - 1]
    TT, TD, DD = TT[ends[tid]], TD[ends[tid]], DD[ends[tid]]
    add = 1.0 if plus_one else 0.0
    f1 = np.maximum(TD - DD + add, 0) / np.maximum(TT, 1)
    f2 = (TD + add) / np.maximum(TT, 1)
    bad = TT < min_tt
    f1 = np.where(bad, 1.0, np.minimum(f1, 1.0))
    f2 = np.where(bad, 1.0, np.minimum(f2, 1.0))
    q1 = np.minimum.accumulate(f1[::-1])[::-1]
    q2 = np.minimum.accumulate(f2[::-1])[::-1]
    o1, o2 = np.empty(n), np.empty(n)
    o1[order], o2[order] = q1, q2
    return o1, o2


def _unordered(a, b):
    a = np.asarray(a).astype(str)
    b = np.asarray(b).astype(str)
    lo = np.where(a <= b, a, b)
    hi = np.where(a <= b, b, a)
    return lo, hi


def merge_level(level: str, score, prot_a, prot_b, decoy_a, decoy_b, pep_a=None, pep_b=None,
                site_a=None, site_b=None, csm_id=None):
    """Merge CSMs to entities of ``level``; entity score = best member score.

    Returns dict with arrays: key, score, n_decoy, prot_a, prot_b (one row per entity).
    """
    s = _clean_scores(score)
    pa, pb = np.asarray(prot_a).astype(str), np.asarray(prot_b).astype(str)
    da, db = np.asarray(decoy_a, bool), np.asarray(decoy_b, bool)
    n = s.shape[0]
    if level == "csm":
        key = np.asarray(csm_id if csm_id is not None else np.arange(n)).astype(str)
    elif level == "peptide_pair":
        x, y = _unordered(np.char.add(np.char.add(np.asarray(pep_a).astype(str), "|"), da.astype(str)),
                          np.char.add(np.char.add(np.asarray(pep_b).astype(str), "|"), db.astype(str)))
        key = np.char.add(np.char.add(x, "~"), y)
    elif level == "residue_pair":
        x, y = _unordered(np.char.add(np.char.add(pa, ":"), np.asarray(site_a).astype(str)),
                          np.char.add(np.char.add(pb, ":"), np.asarray(site_b).astype(str)))
        key = np.char.add(np.char.add(x, "~"), y)
    elif level == "ppi":
        x, y = _unordered(pa, pb)
        key = np.char.add(np.char.add(x, "~"), y)
    else:
        raise ValueError(level)
    order = np.lexsort((-s, key))
    ks = key[order]
    first = np.r_[True, ks[1:] != ks[:-1]]
    w = order[first]
    return {"key": key[w], "score": s[w], "n_decoy": xl_class(da[w], db[w]), "prot_a": pa[w], "prot_b": pb[w]}


def xl_level_fdr(level: str, score, prot_a, prot_b, decoy_a, decoy_b, base_of: Mapping | None = None,
                 plus_one: bool = False, min_tt: int = 0, **kw) -> dict:
    """Merge to ``level`` then compute q-values separately for self and
    heteromeric links. Returns {"self": {...}, "hetero": {...}}, each with the
    merged entities and q_tddd / q_td."""
    ent = merge_level(level, score, prot_a, prot_b, decoy_a, decoy_b, **kw)
    selfm = is_self_link(ent["prot_a"], ent["prot_b"], base_of)
    out = {}
    for name, m in (("self", selfm), ("hetero", ~selfm)):
        sub = {k: v[m] for k, v in ent.items()}
        q1, q2 = xl_qvalues(sub["score"], sub["n_decoy"], plus_one=plus_one, min_tt=min_tt)
        sub["q_tddd"], sub["q_td"] = q1, q2
        out[name] = sub
    return out
