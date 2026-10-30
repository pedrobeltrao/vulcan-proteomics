"""E7 chimeric competition (DESIGN §4 L4; pre-registered in results/2026-10-28_E7_chimeric_v2/config.yaml).

1. Stage 1: per query unit (run, spectrum, precursor group) the best PSM by the score (hash tie-break).  Each
   co-isolated precursor keeps its own winner.
2. CONGA neighbour filter (label-blind, vulcan.ptm.conga.neighbour_filter): within a spectrum, in descending score
   order over the group winners, drop a winner whose peptide shares > 5% of its b/y ions (of the smaller ion set) with a
   higher-scoring kept winner.  Ion set: singly charged b_i and y_i with i >= 3 (the K1 fragment-index convention
   min_ion_index = 2: b1, b2, y1, y2 are excluded because y1/y2 coincide by construction between tryptic peptides
   ending in the same residue), binned at 0.02 Da.
3. Paired double competition over the survivors (peptidoform units, pair = tdc_pair) -> pair winners.
4. Strata = chimeric rank of the pair winner pooled {1, 2, >=3} (groupwalk.pool_chimeric_rank), merged by
   groupwalk.merge_adjacent over label-blind pair counts (>= 200); q by Group-walk (K = 40); pooled TDC alongside;
   q_sep = plain TDC separately within each merged stratum (amendment E7-A1: every stratum list controlled at alpha).
Labels enter only at steps 3-4 (as in every L4 competition).
"""
from __future__ import annotations

import numpy as np

from vulcan.ptm.conga import neighbour_filter
from vulcan.stats import groupwalk as GW
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.tdc import compete, double_competition, tdc_qvalues

PROTON = 1.007276466812
H2O = 18.010564684
MONO = {"G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414, "T": 101.047679,
        "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578,
        "K": 128.094963, "E": 129.042593, "M": 131.040485, "H": 137.058912, "F": 147.068414, "R": 156.101111,
        "Y": 163.06332, "W": 186.079313, "U": 150.95363, "O": 237.147727}
BIN_DA = 0.02
FRAC = 0.05
MIN_ION_INDEX = 3


class IonSets:
    """Cached b/y ion-bin sets per index row (label-blind: residues + mod mask only)."""

    def __init__(self, l1, static_mods=None, var_mods=None, bin_da=BIN_DA, min_ion=MIN_ION_INDEX):
        static_mods = {"C": 57.021464} if static_mods is None else static_mods
        var_mods = {"M": 15.9949} if var_mods is None else var_mods
        self.l1, self.bin, self.min_ion = l1, bin_da, min_ion
        self.aa = np.zeros(128)
        self.var = np.zeros(128)
        for a, m in MONO.items():
            self.aa[ord(a)] = m + static_mods.get(a, 0.0)
        for a, m in var_mods.items():
            if a != "^":
                self.var[ord(a)] = m
        self.cache = {}

    def get(self, row: int) -> frozenset:
        s = self.cache.get(row)
        if s is not None:
            return s
        l1 = self.l1
        o, n = int(l1.seq_off[row]), int(l1.length[row])
        r = np.asarray(l1.residues[o:o + n], np.int64)
        mask = int(np.int64(l1.mask[row]))
        bits = np.array([(mask >> k) & 1 for k in range(n)], bool)
        m = self.aa[r] + np.where(bits, self.var[r], 0.0)
        pre = np.cumsum(m)[:-1]
        tot = m.sum() + H2O
        k = self.min_ion
        b = pre[k - 1:] + PROTON                      # b_i, i >= k  (pre[i-1] = b_i neutral)
        y = (tot - pre)[: max(0, n - k)] + PROTON     # y_j = tot - pre[n-j-1]; j >= k  <=> index <= n-1-k
        ions = np.concatenate([b, y])
        s = frozenset(np.round(ions / self.bin).astype(np.int64).tolist())
        self.cache[row] = s
        return s

    def sets(self, rows):
        return [self.get(int(r)) for r in rows]


def group_winners(qkey, score, tiebreak):
    """Stage 1: best PSM per query unit (run, spectrum, group)."""
    return compete(qkey, score, tiebreak=tiebreak)


def filter_winners(spec_key, score, rows, ion_sets: IonSets, frac=FRAC):
    """CONGA neighbour filter over group winners of each spectrum -> keep mask."""
    return neighbour_filter(spec_key, score, ion_sets.sets(rows), frac)


def chimeric_competition(spec_key, group, chim_rank, unit, pair, score, is_decoy, rows, ion_sets: IonSets,
                         seed: int = 11, K: int = 40, pool_from: int = 3, min_size: int = GW.MIN_STRATUM_SIZE,
                         neighbour: bool = True) -> dict:
    """All arrays per PSM.  Returns dict with the double-competition result `res` (rows index the INPUT arrays),
    stratum per pair winner, q_gw (Group-walk), q_pooled, merge map and diagnostics."""
    spec_key = np.asarray(spec_key)
    unit = np.asarray(unit)
    score = np.asarray(score, np.float64)
    tb = hash_tiebreak(unit, seed)
    qkey = np.char.add(np.char.add(spec_key.astype(str), "#"), np.asarray(group).astype(str)) if spec_key.dtype.kind in "UO" \
        else spec_key.astype(np.int64) * 64 + np.asarray(group, np.int64)
    w1 = group_winners(qkey, score, tb)
    if neighbour:
        keep = filter_winners(spec_key[w1], score[w1], np.asarray(rows)[w1], ion_sets)
        w = w1[keep]
    else:
        keep = np.ones(len(w1), bool)
        w = w1
    res = double_competition(qkey[w], unit[w], np.asarray(pair)[w], score[w], np.asarray(is_decoy)[w], seed=seed,
                             tiebreak=tb[w])
    win = w[res.row]                                   # input-array index of each pair winner
    stratum = GW.pool_chimeric_rank(np.asarray(chim_rank)[win].astype(np.int64), pool_from)
    counts = GW.label_blind_pair_counts(res.pair, stratum)
    mp = GW.merge_adjacent(list(range(1, pool_from + 1)), counts, min_size)
    merged = np.array(["_".join(map(str, mp[int(s)])) for s in stratum])
    q_gw = GW.group_walk(res.score, res.is_decoy, merged, K=K, seed=seed)
    q_sep = tdc_qvalues(res.score, res.is_decoy, block=merged)   # amendment E7-A1: TDC within each stratum
    return {"res": res, "win": win, "survivors": w, "stratum": stratum, "merged": merged, "q_gw": q_gw, "q_pooled": res.q, "q_sep": q_sep,
            "pair_counts": {int(k): int(v) for k, v in counts.items()},
            "merge": {int(k): list(map(int, v)) for k, v in mp.items()},
            "n_group_winners": int(len(w1)), "n_after_neighbour_filter": int(len(w)),
            "n_dropped_by_neighbour_filter": int((~keep).sum())}
