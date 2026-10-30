"""Protein-level FDR: Picked Protein Group FDR and classic protein TDC.

References
----------
* Savitski, Wilhelm, Hahne, Kuster, Bantscheff. "A scalable approach for protein
  false discovery rate estimation in large proteomic data sets." MCP 14:2394
  (2015). Picked target-decoy protein FDR.
* The, Samaras, Kuster, Wilhelm. "Well-calibrated protein-group FDR for large
  proteomic data sets: the picked protein group FDR." MCP 21:100437 (2022);
  reference implementation github.com/kusterlab/picked_group_fdr (competition.py,
  grouping.py, observed_peptides.py; read 2026-10-26).

Procedure implemented (DESIGN.md §4 L4)
---------------------------------------
1. Subset grouping over observed peptides: proteins with identical peptide sets
   form one group; a protein whose peptide set is a subset of another protein's
   set is merged into the superset protein with the most peptides.
   Razor/Occam assignment is never used.
2. Shared peptides (mapping to >1 group) are discarded for scoring. The group
   score is the **best peptide score** (higher is better), not a product of PEPs.
   Groups without a unique peptide are dropped. (The reference implementation's
   connected-component regrouping of such groups is not reproduced.)
3. Picked group competition: groups are shuffled (seeded) and sorted by score,
   descending (regular groups before obsolete ones at equal score). Walking down,
   a group is dropped if any member's target/decoy pair key was already "seen";
   otherwise it is kept and the pair keys of its leading proteins (most unique
   peptides) are marked seen.
4. q-values: (D+1)/max(T,1), monotone, decoys counted first at ties (vulcan TDC
   convention; the reference uses (D+1)/(T+1)).
5. Rescued subsets (optional, default on): the protein score at the
   ``rescue_alpha`` (1%) protein-FDR cutoff is used to filter peptides; grouping
   is redone on those peptides only, which merges subset groups that were split
   only by low-scoring peptides. Groups of the first pass whose proteins do not
   appear in the rescued grouping are kept as "obsolete" groups: they compete and
   count in the FDR but are not reported. Scores are then recomputed from all
   peptides and steps 3-4 rerun.
   Rescued subsets have only empirical support (The 2022); DESIGN requires
   protein-level entrapment validation.

Classic protein FDR (comparison only; anticonservative at scale): the same
grouping and scores, all target and decoy groups ranked together, no picking.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from ._util import make_rng
from .tdc import tdc_qvalues

__all__ = ["subset_grouping", "score_groups", "picked_group_fdr", "classic_protein_fdr",
           "ProteinGroupResult", "picked_protein_fdr_simple"]


@dataclass
class ProteinGroupResult:
    groups: list            # list of tuples of protein ids
    score: np.ndarray
    is_decoy: np.ndarray
    q: np.ndarray
    is_obsolete: np.ndarray

    def accepted(self, alpha: float) -> list:
        return [g for g, q, d, o in zip(self.groups, self.q, self.is_decoy, self.is_obsolete)
                if q <= alpha and not d and not o]


def subset_grouping(pep2prot: Mapping) -> list:
    """Subset grouping of proteins from a peptide -> proteins map."""
    prot2pep = defaultdict(set)
    for pep, prots in pep2prot.items():
        for p in prots:
            prot2pep[p].add(pep)
    parent = {p: p for p in prot2pep}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in sorted(prot2pep, key=lambda x: (-len(prot2pep[x]), str(x))):
        peps = prot2pep[p]
        cand = None
        for pep in peps:
            s = set(pep2prot[pep])
            cand = s if cand is None else cand & s
            if len(cand) == 1:
                break
        sup = sorted((c for c in cand if c != p), key=lambda x: (-len(prot2pep[x]), str(x)))
        if sup:
            a, b = find(sup[0]), find(p)
            if a != b:
                parent[b] = a
    groups = defaultdict(list)
    for p in prot2pep:
        groups[find(p)].append(p)
    return [tuple(sorted(g, key=lambda x: (-len(prot2pep[x]), str(x)))) for g in
            sorted(groups.values(), key=lambda g: str(sorted(map(str, g))))]


def score_groups(groups: Sequence, pep2prot: Mapping, pep_score: Mapping):
    """Best unique-peptide score per group and leading proteins.

    Returns (scores (float array, -inf if no unique peptide), leading list of tuples).
    Proteins not in any group are ignored when deciding uniqueness.
    """
    prot2g = {}
    for i, g in enumerate(groups):
        for p in g:
            prot2g[p] = i
    best = np.full(len(groups), -np.inf)
    counts = [defaultdict(int) for _ in groups]
    for pep, prots in pep2prot.items():
        gs = {prot2g[p] for p in prots if p in prot2g}
        if len(gs) != 1:
            continue
        (i,) = gs
        s = pep_score[pep]
        if s > best[i]:
            best[i] = s
        for p in prots:
            if p in prot2g:
                counts[i][p] += 1
    leading = []
    for i, g in enumerate(groups):
        if counts[i]:
            m = max(counts[i].values())
            leading.append(tuple(p for p in g if counts[i].get(p, 0) == m))
        else:
            leading.append(tuple(g))
    return best, leading


def _compete_groups(groups, scores, leading, obsolete, pair_of, is_decoy_prot, rng, picked=True):
    n = len(groups)
    perm = rng.permutation(n)
    # sort by score desc, regular before obsolete, then random
    keys = np.lexsort((perm, np.asarray(obsolete, int), -np.asarray(scores)))
    seen = set()
    keep = []
    for i in keys:
        if not np.isfinite(scores[i]):
            continue
        if picked:
            if any(pair_of[p] in seen for p in groups[i]):
                continue
            seen.update(pair_of[p] for p in leading[i])
        keep.append(i)
    keep = np.asarray(keep, dtype=np.int64)
    gl = [groups[i] for i in keep]
    sc = np.asarray(scores)[keep] if len(keep) else np.zeros(0)
    dec = np.array([all(is_decoy_prot[p] for p in g) for g in gl], dtype=bool)
    obs = np.asarray(obsolete, bool)[keep] if len(keep) else np.zeros(0, bool)
    q = tdc_qvalues(sc, dec) if len(keep) else np.zeros(0)
    return ProteinGroupResult(groups=gl, score=sc, is_decoy=dec, q=q, is_obsolete=obs)


def picked_group_fdr(pep2prot: Mapping, pep_score: Mapping, pair_of: Mapping, is_decoy_prot: Mapping,
                     rescue: bool = True, rescue_alpha: float = 0.01, seed=None, grouping: str = "subset",
                     picked: bool = True) -> ProteinGroupResult:
    """Picked Protein Group FDR.

    Parameters
    ----------
    pep2prot : peptide -> list of proteins (target peptides map to target
        proteins, decoy peptides to decoy proteins).
    pep_score : peptide -> best score (higher better). Peptides absent get no score.
    pair_of : protein -> pair key, identical for a target protein and its decoy.
    is_decoy_prot : protein -> bool.
    grouping : "subset" (default) or "none" (each protein its own group =
        picked protein FDR of Savitski 2015).
    picked : False gives the classic protein FDR on the same groups.
    """
    rng = make_rng(seed)
    obs_map = {k: v for k, v in pep2prot.items() if k in pep_score}
    if grouping == "subset":
        groups = subset_grouping(obs_map)
    elif grouping == "none":
        groups = sorted({(p,) for prots in obs_map.values() for p in prots}, key=str)
    else:
        raise ValueError(grouping)
    scores, leading = score_groups(groups, obs_map, pep_score)
    res = _compete_groups(groups, scores, leading, [False] * len(groups), pair_of, is_decoy_prot, rng, picked)
    if not rescue or grouping == "none":
        return res
    acc = (res.q <= rescue_alpha) & ~res.is_decoy
    cutoff = res.score[acc].min() if acc.any() else (res.score.min() if len(res.score) else -np.inf)
    filt = {k: v for k, v in obs_map.items() if pep_score[k] >= cutoff}
    new_groups = subset_grouping(filt)
    in_new = {p for g in new_groups for p in g}
    obsolete_groups = [g for g in groups if not any(p in in_new for p in g)]
    s_new, l_new = score_groups(new_groups, obs_map, pep_score)
    old_idx = {g: i for i, g in enumerate(groups)}
    s_obs = np.array([scores[old_idx[g]] for g in obsolete_groups]) if obsolete_groups else np.zeros(0)
    l_obs = [leading[old_idx[g]] for g in obsolete_groups]
    all_groups = list(new_groups) + list(obsolete_groups)
    all_scores = np.concatenate([s_new, s_obs])
    all_leading = list(l_new) + l_obs
    obsolete = [False] * len(new_groups) + [True] * len(obsolete_groups)
    return _compete_groups(all_groups, all_scores, all_leading, obsolete, pair_of, is_decoy_prot, rng, picked)


def classic_protein_fdr(pep2prot, pep_score, pair_of, is_decoy_prot, seed=None, grouping="subset"):
    """Classic (non-picked) protein-group TDC on the same groups (comparison only)."""
    return picked_group_fdr(pep2prot, pep_score, pair_of, is_decoy_prot, rescue=False, seed=seed,
                            grouping=grouping, picked=False)


def picked_protein_fdr_simple(target_score, decoy_score, block=None, picked: bool = True, seed=None):
    """Vectorised picked (or classic) protein FDR for one-protein-per-group data.

    ``target_score[i]``/``decoy_score[i]`` are the best-peptide scores of target
    protein i and its decoy (-inf if no peptide). Picked: each pair keeps its
    higher-scoring member (random ties). Classic: both enter. ``block`` keeps
    replicates separate. Returns (score, is_decoy, q, pair_index, block).
    """
    t = np.asarray(target_score, float)
    d = np.asarray(decoy_score, float)
    n = t.shape[0]
    b = np.zeros(n, dtype=np.int64) if block is None else np.asarray(block)
    idx = np.arange(n)
    if picked:
        rng = make_rng(seed)
        coin = rng.random(n) < 0.5
        dwin = (d > t) | ((d == t) & coin)
        s = np.where(dwin, d, t)
        isd = dwin
        keep = np.isfinite(s)
        s, isd, pi, bb = s[keep], isd[keep], idx[keep], b[keep]
    else:
        s = np.concatenate([t, d])
        isd = np.r_[np.zeros(n, bool), np.ones(n, bool)]
        pi = np.r_[idx, idx]
        bb = np.r_[b, b]
        keep = np.isfinite(s)
        s, isd, pi, bb = s[keep], isd[keep], pi[keep], bb[keep]
    q = tdc_qvalues(s, isd, block=bb)
    return s, isd, q, pi, bb
