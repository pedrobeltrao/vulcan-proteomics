"""L6 repository-scale FDR (DESIGN.md §4 L6, §7 E9): global paired TDC over the best evidence across runs.

Procedure (one frozen, label-blind score shared by every run; one pooling group, fixed before estimation):

1. **Spectrum level, per run.** All candidates of a spectrum compete (``psm_winners``); spectra are never pooled
   across runs, so this step is identical to the per-run analysis.
2. **Best evidence per unit across runs.** Among the PSM winners of *all* runs, the best PSM per unit (peptidoform
   for DDA, (peptide, charge) for precursor level) is kept, separately for the target and decoy side.
3. **Global paired competition.** Each target unit competes with its own paired decoy unit (shared ``pair_id``).
4. **Global q.** (D+1)/T, monotone, over the pair winners (the vulcan TDC convention, ``tdc.tdc_qvalues``).

Validity: for an absent (null) target unit, the best score over runs of the target and of its paired decoy are
exchangeable when each run's own target/decoy scores are (the max of exchangeable pairs is exchangeable), so the
pair-competition coin stays fair however many runs are pooled. Calibration across runs therefore affects *power*
(which run's evidence wins, how decoys from bad runs crowd the tail), not validity; it is checked with
:func:`score_comparability` (G4-style KS on decoys, decided on dev/training data only, §6.4).

Naive merge (the E9 comparator): per-run double competition at q_run <= alpha, then the union of the accepted units.
Its FDP grows with N because every run contributes its own ~alpha false units and false units differ across runs,
while true units are shared.

Outputs (:func:`global_outputs`): a master list (unit, best run, score, global q), per-run lists (unit, run q,
global q of the unit, globally accepted flag) and, per alpha, the expected false count q x N (= the TDC estimate
D+1 at the threshold) next to entrapment FDP bounds supplied by the caller.

Picked protein groups at repository scale are :func:`vulcan.stats.protein.picked_group_fdr` applied to the global
pair winners (best peptide score per group), exactly as per run.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._util import as_codes, segment_starts
from .tdc import compete, tdc_qvalues

__all__ = ["psm_winners", "runner_up_stats", "calibrate", "GlobalResult", "global_competition", "run_competitions", "naive_union",
           "score_comparability", "global_outputs"]


def psm_winners(spectrum_id, score, tiebreak) -> np.ndarray:
    """Step 1 for one run: winner row per spectrum (max score; ties by the storage-order-free key)."""
    return compete(spectrum_id, score, tiebreak=tiebreak)


def runner_up_stats(spectrum_id, score, q_lo: float = 0.5, q_hi: float = 0.9):
    """LABEL-BLIND per-run null proxy for cross-run calibration: the runner-up (2nd best) candidate score of every
    spectrum that has >= 2 candidates.  Runner-ups are overwhelmingly random matches of either label, and the
    statistic never looks at labels, so a calibration built from it keeps TDC validity (a label swap of any null
    pair leaves it unchanged; the per-run map is monotone and applied to targets and decoys alike).

    Returns (location, scale) = (quantile q_lo, quantile q_hi - quantile q_lo) of the runner-up scores.
    """
    g = as_codes(spectrum_id)
    s = np.asarray(score, np.float64)
    o = np.lexsort((-s, g))
    gs = g[o]
    first = segment_starts(gs)
    pos = np.arange(len(o)) - np.repeat(np.flatnonzero(first), np.diff(np.r_[np.flatnonzero(first), len(o)]))
    r2 = s[o][pos == 1]
    if len(r2) < 50:
        return 0.0, 1.0
    lo, hi = np.quantile(r2, [q_lo, q_hi])
    return float(lo), float(max(hi - lo, 1e-6))


def calibrate(score, run_index, loc, scale):
    """Per-run monotone affine map (score - loc[run]) / scale[run]."""
    loc, scale = np.asarray(loc, float), np.asarray(scale, float)
    return (np.asarray(score, float) - loc[run_index]) / scale[run_index]


@dataclass
class GlobalResult:
    """One row per competed pair unit (pair-competition winners)."""

    row: np.ndarray        # index into the (PSM-winner) input table
    unit: np.ndarray
    pair: np.ndarray
    run: np.ndarray        # run of the best evidence
    score: np.ndarray
    is_decoy: np.ndarray
    q: np.ndarray
    extra: dict = field(default_factory=dict)

    def accepted(self, alpha: float) -> np.ndarray:
        return (self.q <= alpha) & ~self.is_decoy


def _dense(x) -> np.ndarray:
    """Dense int64 codes 0..k-1 (hash ids such as u64 pf_id are not safe to combine arithmetically)."""
    return np.unique(np.asarray(x), return_inverse=True)[1].astype(np.int64).ravel()


def _unit_pair(rows_unit, rows_dec, rows_pair, rows_score, rows_tb):
    """Steps 2-3 on PSM-winner rows: best per (unit, side), then pair competition. Returns winning row indices."""
    key = _dense(rows_unit) * 2 + rows_dec.astype(np.int64)
    w2 = compete(key, rows_score, tiebreak=rows_tb)
    w3 = w2[compete(rows_pair[w2], rows_score[w2], tiebreak=rows_tb[w2])]
    return w3


def global_competition(run, unit, pair, score, is_decoy, tiebreak) -> GlobalResult:
    """Steps 2-4 over PSM winners pooled from all runs (step 1 must already have been applied per run).

    All arrays are row-aligned over the pooled PSM winners. ``tiebreak`` = hash(unit, seed) (storage-order free).
    """
    run, unit, pair = np.asarray(run), np.asarray(unit), np.asarray(pair)
    score = np.asarray(score, np.float64)
    is_decoy = np.asarray(is_decoy, bool)
    tb = np.asarray(tiebreak)
    w = _unit_pair(unit, is_decoy, pair, score, tb)
    q = tdc_qvalues(score[w], is_decoy[w])
    return GlobalResult(row=w, unit=unit[w], pair=pair[w], run=run[w], score=score[w], is_decoy=is_decoy[w], q=q)


def run_competitions(run, unit, pair, score, is_decoy, tiebreak) -> GlobalResult:
    """Per-run double competition (steps 2-4 inside each run) for every run at once.

    Returns one row per (run, competed pair unit) with the *run-level* q-value (vectorised over runs).
    """
    run = np.asarray(run)
    rc = _dense(run)
    unit, pair = np.asarray(unit), np.asarray(pair)
    score = np.asarray(score, np.float64)
    is_decoy = np.asarray(is_decoy, bool)
    tb = np.asarray(tiebreak)
    uc = _dense(unit)
    nu = int(uc.max()) + 1 if len(unit) else 1
    ukey = rc * (2 * nu) + uc * 2 + is_decoy.astype(np.int64)
    w2 = compete(ukey, score, tiebreak=tb)
    pc = _dense(pair)
    npc = int(pc.max()) + 1 if len(pc) else 1
    w3 = w2[compete(rc[w2] * npc + pc[w2], score[w2], tiebreak=tb[w2])]
    q = tdc_qvalues(score[w3], is_decoy[w3], block=rc[w3])
    return GlobalResult(row=w3, unit=unit[w3], pair=pair[w3], run=run[w3], score=score[w3], is_decoy=is_decoy[w3], q=q)


def naive_union(per_run: GlobalResult, alpha: float):
    """Union of the per-run lists at q_run <= alpha: (unit, max score over the runs that accepted it, n_runs)."""
    acc = per_run.accepted(alpha)
    u, s = per_run.unit[acc], per_run.score[acc]
    if not len(u):
        return u, s, np.zeros(0, np.int64)
    o = np.lexsort((-s, u))
    first = segment_starts(u[o])
    counts = np.diff(np.r_[np.flatnonzero(first), len(o)])
    return u[o][first], s[o][first], counts


def _ks(a, b):
    a, b = np.sort(a), np.sort(b)
    x = np.concatenate([a, b])
    return float(np.max(np.abs(np.searchsorted(a, x, side="right") / len(a) - np.searchsorted(b, x, side="right") / len(b))))


def score_comparability(run, score, is_decoy, strata=None, min_n: int = 200, seed: int = 11, max_ref: int = 200_000):
    """G4-style comparability of a frozen score across runs, on DECOYS only (DESIGN §6.4 G4; dev/training data only).

    For every run (and stratum) with >= min_n decoys: KS distance between the run's decoy scores and the pooled decoy
    scores of all other runs (subsampled to max_ref).  Reference null: KS between two random halves of the same
    run's decoys (what pure sampling noise gives at that n); ``excess = KS - null_KS``.  Returns a dict of per-run
    rows and summary quantiles.
    """
    rng = np.random.default_rng(seed)
    run = np.asarray(run)
    score = np.asarray(score, float)
    d = np.asarray(is_decoy, bool)
    st = np.zeros(len(run), np.int64) if strata is None else as_codes(strata)
    strata_vals = np.asarray(strata)[np.unique(st, return_index=True)[1]] if strata is not None else np.array(["all"])
    rows = []
    for k, sv in enumerate(np.unique(st)):
        m = d & (st == sv)
        r_m, s_m = run[m], score[m]
        for r in np.unique(r_m):
            x = s_m[r_m == r]
            if len(x) < min_n:
                continue
            ref = s_m[r_m != r]
            if len(ref) < min_n:
                continue
            if len(ref) > max_ref:
                ref = rng.choice(ref, max_ref, replace=False)
            perm = rng.permutation(len(x))
            null = _ks(x[perm[: len(x) // 2]], x[perm[len(x) // 2:]])
            ks = _ks(x, ref)
            rows.append({"run": str(r), "stratum": str(strata_vals[k]), "n_decoys": int(len(x)), "ks": ks,
                         "null_ks_halves": null, "median": float(np.median(x)), "q90": float(np.quantile(x, 0.9)),
                         "q99": float(np.quantile(x, 0.99))})
    out = {"rows": rows}
    if rows:
        ks = np.array([r["ks"] for r in rows])
        ex = ks - np.array([r["null_ks_halves"] for r in rows])
        q90 = np.array([r["q90"] for r in rows])
        out["summary"] = {"n_run_strata": len(rows), "ks_median": float(np.median(ks)),
                          "ks_q90": float(np.quantile(ks, 0.9)), "ks_max": float(ks.max()),
                          "excess_median": float(np.median(ex)), "excess_q90": float(np.quantile(ex, 0.9)),
                          "decoy_q90_score_iqr_across_runs": float(np.quantile(q90, 0.75) - np.quantile(q90, 0.25))}
    return out


def global_outputs(glob: GlobalResult, per_run: GlobalResult, alphas=(0.001, 0.01, 0.05)):
    """Master list, per-run lists and expected false counts.

    master: dict of arrays over globally competed TARGET units (unit, best_run, score, q_global).
    per_run_lists: dict of arrays over per-run competed TARGET units (run, unit, q_run, q_global of the unit;
        NaN if the unit's global pair competition was won by its decoy).
    expected_false: per alpha, N accepted and the TDC expected false count (D+1 at the threshold = q x N bound).
    """
    t = ~glob.is_decoy
    master = {"unit": glob.unit[t], "best_run": glob.run[t], "score": glob.score[t], "q_global": glob.q[t]}
    o = np.argsort(master["unit"], kind="stable")
    mu, mq = master["unit"][o], master["q_global"][o]
    pt = ~per_run.is_decoy
    pu = per_run.unit[pt]
    pos = np.searchsorted(mu, pu)
    pos_c = np.minimum(pos, max(len(mu) - 1, 0))
    hit = (len(mu) > 0) & (mu[pos_c] == pu) if len(mu) else np.zeros(len(pu), bool)
    qg = np.where(hit, mq[pos_c] if len(mu) else np.nan, np.nan)
    per = {"run": per_run.run[pt], "unit": pu, "q_run": per_run.q[pt], "q_global": qg}
    exp = {}
    for a in alphas:
        acc = glob.accepted(a)
        n = int(acc.sum())
        if n:
            s_cut = glob.score[acc].min()
            d = int(((glob.score >= s_cut) & glob.is_decoy).sum())
        else:
            d = 0
        exp[str(a)] = {"n_accepted": n, "decoys_at_threshold": d, "expected_false_q_times_n": a * n,
                       "tdc_estimate_false": float(min(d + 1, n)) if n else 0.0,
                       "per_run_mean_globally_accepted": float(np.bincount(as_codes(per["run"]),
                                                                           weights=(per["q_global"] <= a)).mean())
                       if len(per["run"]) else 0.0}
    return master, per, exp
