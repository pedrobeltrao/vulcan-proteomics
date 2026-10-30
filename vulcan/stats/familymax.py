"""Family-max single-pass entrapment (DESIGN.md v0.3 §4 L1, §6.3).

Index entries carry a family tag: ``T`` (targets), ``E_s`` (entrapment, seed s),
``dec_d(T)`` and ``dec_d(E_s)`` (decoys of seed d), ``foreign`` (natural
foreign-proteome entrapment, with its own decoy families). The search keeps the
top-2 candidates per family per precursor hypothesis. Any replicate database,
i.e. a union of families such as {T, E_2001, dec_1001(T), dec_1001(E_2001)}, can
then be competed post hoc: the replicate's best and second-best candidate for a
hypothesis are always among the top-2 of the included families, so the winner and
the runner-up score (needed for the label-free Δscore feature) are exact.

Ties are broken by a storage-order-free key (e.g. hash(sequence, seed)): the
smallest ``tiebreak`` wins among equal scores, so reconstruction and a separate
search give identical winners.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np

from ._util import as_codes

__all__ = ["top_k_per_family", "reconstruct_competition", "ReplicateSpec", "Reconstructed"]


def top_k_per_family(spectrum_id, hypothesis_id, family, score, tiebreak, k: int = 2) -> np.ndarray:
    """Row indices of the top-k candidates per (spectrum hypothesis, family).

    This is what the search kernels keep (K1/K2 output). Order: by hypothesis,
    family, descending score, ascending tiebreak.
    """
    sp = as_codes(spectrum_id)
    hy = as_codes(hypothesis_id)
    fa = as_codes(family)
    s = np.asarray(score, float)
    tb = np.asarray(tiebreak)
    order = np.lexsort((tb, -s, fa, hy, sp))
    key_sp, key_hy, key_fa = sp[order], hy[order], fa[order]
    newgrp = np.r_[True, (key_sp[1:] != key_sp[:-1]) | (key_hy[1:] != key_hy[:-1]) | (key_fa[1:] != key_fa[:-1])]
    gid = np.cumsum(newgrp) - 1
    starts = np.flatnonzero(newgrp)
    rank = np.arange(order.size) - starts[gid]
    return order[rank < k]


@dataclass
class ReplicateSpec:
    """Which families form one replicate database and their roles."""

    targets: tuple          # e.g. ("T",)
    entrapment: tuple       # e.g. ("E_2001",)  (target-side, counted as entrapment)
    decoys: tuple           # e.g. ("dec_1001(T)", "dec_1001(E_2001)")
    entrapment_decoys: tuple = ()  # subset of decoys that are decoys of entrapment

    @property
    def families(self) -> tuple:
        return tuple(self.targets) + tuple(self.entrapment) + tuple(self.decoys)


@dataclass
class Reconstructed:
    spectrum: np.ndarray     # spectrum key per winner
    row: np.ndarray          # index of the winning candidate row (into the input arrays)
    score: np.ndarray
    runner_up: np.ndarray    # second-best score in the replicate DB (-inf if none)
    family: np.ndarray
    is_decoy: np.ndarray
    is_entrapment: np.ndarray


def reconstruct_competition(spectrum_id, family, score, tiebreak, spec: ReplicateSpec) -> Reconstructed:
    """Spectrum-level competition (all hypotheses of a spectrum compete) for one
    replicate database, from candidate rows (e.g. the top-2-per-family table)."""
    fam = np.asarray(family).astype(str)
    keep = np.isin(fam, np.array(spec.families, dtype=str))
    idx = np.flatnonzero(keep)
    sp = as_codes(np.asarray(spectrum_id))[idx]
    s = np.asarray(score, float)[idx]
    tb = np.asarray(tiebreak)[idx]
    order = np.lexsort((tb, -s, sp))
    sps = sp[order]
    first = np.r_[True, sps[1:] != sps[:-1]] if sps.size else np.zeros(0, bool)
    second = np.r_[False, ~first[1:]] & np.r_[False, first[:-1]] if sps.size else np.zeros(0, bool)
    win = order[first]
    ru = np.full(win.size, -np.inf)
    # second element of each segment
    seg = np.cumsum(first) - 1
    sec_pos = np.flatnonzero(second)
    ru[seg[sec_pos]] = s[order[sec_pos]]
    rows = idx[win]
    f = fam[rows]
    is_dec = np.isin(f, np.array(spec.decoys, dtype=str))
    is_ent = np.isin(f, np.array(spec.entrapment, dtype=str))
    return Reconstructed(spectrum=np.asarray(spectrum_id)[rows], row=rows, score=s[win], runner_up=ru,
                         family=f, is_decoy=is_dec, is_entrapment=is_ent)


def replicate_specs(decoy_seeds: Iterable[int], entrap_seeds: Iterable[int]) -> Mapping:
    """Paired-arm replicate specs for every (decoy seed, entrapment seed)."""
    out = {}
    for d in decoy_seeds:
        for s in entrap_seeds:
            out[(d, s)] = ReplicateSpec(targets=("T",), entrapment=(f"E_{s}",),
                                        decoys=(f"dec_{d}(T)", f"dec_{d}(E_{s})"),
                                        entrapment_decoys=(f"dec_{d}(E_{s})",))
    return out
