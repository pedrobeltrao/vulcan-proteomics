"""Amendment A3 (DESIGN §6.0c = reviews/amendment_A3_draft.md rev 1) gate statistics for v2_sep and ce_s.

Docstring = specification of the NEW statistics (A3.7); every A1/A2 rule is reused unchanged from
vulcan.gates.a2_stats, vulcan.gates.stats and vulcan.stats.g2 (imported, never redefined).  Where the A3 text is open,
the stricter reading is taken (A2.9.5); each choice is also listed in scripts/a3/IMPLEMENTATION_NOTES.md.

Verdict vocabulary (as A2): "pass", "fail", "not-evaluable", "not-evaluable (invalid pool)", plus
"not-evaluable (implementation)" (V3 mismatch, A3.2.4).

V6  G3-NF, neighbour-filter symmetry (A3.2.4 row V6).
    Unit: every stage-1 group winner of an A ∪ B null spectrum with >= 2 precursor groups (after (c-final)).
    Class 1 = winners of rank-1 units, class 2 = winners of rank >= 2 units (rank = the unit's chimeric rank).
    kept = the winner survives the CONGA neighbour filter.  Per class, two bins: the top decile of the winner's
    v2 score within the class (score >= the class's 90% quantile, numpy linear, both labels pooled; label-free
    edges), and the rest.  Per bin: Delta = P(kept | target winner) - P(kept | decoy winner); 90% two-sided Newcombe
    hybrid-score interval (A1 newcombe_diff); the bin passes iff lo >= -0.01 and hi <= +0.01 (closed, as A2.8).
    Evaluable iff BOTH bins of the class have n_T >= n_min AND n_D >= n_min (n_min per label per bin from
    config/a3.yaml = A3.1.9); otherwise the class is not-evaluable (blocks its tier; not final).  Class verdict =
    pass iff both bins pass, else fail.  Two-sided (A3.5.4).
V12 top-5 truncation symmetry (A3.2.4 row V12).  Per A spectrum (rank 1): the exhaustive v2-best target and v2-best
    decoy of the single-hypothesis window W (every candidate with >= min_matched_peaks); kept = that candidate's
    peptidoform row is among the production top-5 rows of the spectrum (K1 top-K + hyperscore top-5, i.e. the whole
    pruning chain a top-K rescorer sees).  Bins = A2.8 bins on the exhaustive v2-best score max(best T, best D)
    (deciles; top five each need n_T, n_D >= n_min; bottom half merged upward; leftover -> not-evaluable), the TOST
    of A2.8 with recall = kept, n_min from config (A3.1.9 by the V6 rule).  Implemented by calling
    vulcan.gates.a2_stats.g3_gate_a2 with best_*_k = best_* where kept else -inf (identical arithmetic).  The same
    statistic with the ce_s logit is ce_s's V12.
V8  G2-original on secondary units: vulcan.stats.g2.g2_differential unchanged; not-evaluable iff
    n < g2_sample_size() (the §6.4 budget, A3.2.4).  Stricter reading: n counts null competitions (units) in the
    row's pool, as G2-original's own unit.
V9  G5 (point rule, AUC <= 0.55, vulcan.gates.stats.g5_gate unchanged).  (c) is evaluable iff >= 200 target and
    >= 200 decoy winners; (a), (b) as A2 (always evaluated).
V13 tail-G2 on A-syn: vulcan.gates.a2_stats.tail_g2_a2 unchanged (n_tail >= 200; UCB_abs <= 0.60; UCB_diff <= 0.10;
    cluster bootstrap B = 10,000 seed 83; hyperscore pool-validity control).  Tail membership is assembled by the
    data layer (unit winner >= tau of the merged stratum containing the unit's rank 2 or 3, per mixture and pipeline).
V14 G1 on A-syn: vulcan.gates.a2_stats.g1_gate_a2 unchanged (n_top >= 7,200, hyperscore control).
Tiers (A3.2.4-A3.2.6).  A row verdict is pass / fail / not-evaluable*.  FS rows (V10): an evaluable fail blocks
    (counted as a fail of tier 1), not-evaluable never blocks (A2.4 rule 3).  A tier's outcome:
        fail            if any of its rows fails (final under P4);
        not-evaluable   else if any of its rows is not evaluable (blocks; not final);
        pass            else.
    Tier 2 requires tier 1 to pass (tier 2 ⊃ tier 1).  Licensed list: tier 2 pass -> full v2_sep list; tier 1 pass
    only -> stratum-1 list only; else none.
"""
from __future__ import annotations

import numpy as np

from vulcan.gates import a2_stats as A2S
from vulcan.gates import stats as GS

PASS, FAIL = "pass", "fail"
NE = A2S.NOT_EVALUABLE
NE_INVALID = A2S.INVALID
NE_IMPL = "not-evaluable (implementation)"
TOL = 0.01
CONF = 0.90


def is_ne(v: str) -> bool:
    return isinstance(v, str) and v.startswith(NE)


# ------------------------------------------------------------------------------------------------ V6
def nf_symmetry_class(kept, is_target, score, n_min: int, tol: float = TOL, conf: float = CONF,
                      top_frac: float = 0.10) -> dict:
    """V6 for one class.  Arrays per stage-1 group winner of the class."""
    kept = np.asarray(kept, bool)
    t = np.asarray(is_target, bool)
    s = np.asarray(score, float)
    if s.size == 0:
        return {"verdict": NE, "reason": "no winners", "n": 0}
    edge = float(np.quantile(s, 1 - top_frac))
    top = s >= edge
    bins = {}
    ok, evaluable = True, True
    for name, m in (("top_decile", top), ("rest", ~top)):
        nT, nD = int((m & t).sum()), int((m & ~t).sum())
        r = A2S.tost_bin(int((m & t & kept).sum()), nT, int((m & ~t & kept).sum()), nD, tol, conf)
        r["n_min"] = n_min
        r["evaluable"] = bool(nT >= n_min and nD >= n_min)
        evaluable &= r["evaluable"]
        ok &= r["pass"]
        bins[name] = r
    verdict = NE if not evaluable else (PASS if ok else FAIL)
    return {"verdict": verdict, "edge_top_decile": edge, "bins": bins, "n": int(s.size), "n_min": n_min,
            "drop_rate": float(1 - kept.mean())}


def nf_symmetry(kept, is_target, score, unit_rank, n_min: int, **kw) -> dict:
    r = np.asarray(unit_rank)
    return {"class1": nf_symmetry_class(np.asarray(kept)[r == 1], np.asarray(is_target)[r == 1],
                                        np.asarray(score)[r == 1], n_min, **kw),
            "class2": nf_symmetry_class(np.asarray(kept)[r >= 2], np.asarray(is_target)[r >= 2],
                                        np.asarray(score)[r >= 2], n_min, **kw)}


# ------------------------------------------------------------------------------------------------ V12
def truncation_symmetry(best_t, best_d, kept_t, kept_d, n_min: int, tol: float = TOL, conf: float = CONF) -> dict:
    """V12: arrays per A spectrum (best_* = exhaustive v2-best (or ce_s-best) score of that label, -inf if none)."""
    bt, bd = np.asarray(best_t, float), np.asarray(best_d, float)
    kt = np.where(np.asarray(kept_t, bool), bt, -np.inf)
    kd = np.where(np.asarray(kept_d, bool), bd, -np.inf)
    r = A2S.g3_gate_a2(bt, bd, kt, kd, n_min=n_min, tol=tol, conf=conf)
    r["statistic"] = "recall_T - recall_D of the exhaustive best candidate within the production top-5"
    return r


# ------------------------------------------------------------------------------------------------ V8 / V9
def g2_original(win_fixed, win_learned, score_fixed, n_min: int | None = None) -> dict:
    from vulcan.stats import g2 as G2
    n = int(np.asarray(win_fixed).size)
    if n_min is not None and n < n_min:
        return {"verdict": NE, "n": n, "n_min": n_min}
    if n == 0:
        return {"verdict": NE, "n": 0}
    r = G2.g2_differential(np.asarray(win_fixed, bool), np.asarray(win_learned, bool), np.asarray(score_fixed, float))
    return {"verdict": PASS if r["pass"] else FAIL, "n": n, "n_min": n_min, "detail": r}


def g5(score, win_t, alarm: float = 0.55, min_winners: int | None = None) -> dict:
    w = np.asarray(win_t, bool)
    nT, nD = int(w.sum()), int((~w).sum())
    if min_winners is not None and (nT < min_winners or nD < min_winners):
        return {"verdict": NE, "nT": nT, "nD": nD, "min_winners": min_winners}
    if nT == 0 or nD == 0:
        return {"verdict": NE, "nT": nT, "nD": nD}
    r = GS.g5_gate(score, w, alarm=alarm)
    return {"verdict": PASS if r["pass"] else FAIL, "nT": nT, "nD": nD, "detail": r}


# ------------------------------------------------------------------------------------------------ tiers
def tier_outcome(verdicts: dict, fs_keys=()) -> dict:
    """verdicts: row -> verdict string.  FS rows (fs_keys): evaluable fail -> fail; not-evaluable ignored."""
    fails, nes = [], []
    for k, v in verdicts.items():
        if k in fs_keys:
            if v == FAIL:
                fails.append(k)
            continue
        if v == FAIL:
            fails.append(k)
        elif v != PASS:
            nes.append(k)
    out = FAIL if fails else (NE if nes else PASS)
    return {"outcome": out, "failed_rows": sorted(fails), "not_evaluable_rows": sorted(nes)}


def licensed(t1: dict, t2: dict) -> dict:
    if t1["outcome"] != PASS:
        lst = "none"
    elif t2["outcome"] == PASS:
        lst = "full v2_sep list (all strata)"
    else:
        lst = "v2_sep stratum-1 list only (strata >= 2 not gate-validated)"
    return {"tier1": t1["outcome"], "tier2": t2["outcome"] if t1["outcome"] == PASS else "requires tier 1",
            "licensed_list": lst}


# ------------------------------------------------------------------------------------------------ power helpers
def tost_bin_pass_prob(n: int, p: float, n_rep: int, rng, tol: float = TOL, conf: float = CONF) -> float:
    """MC pass probability of one TOST bin for a symmetric process (both labels at rate p, n per label)."""
    kT = rng.binomial(n, p, size=n_rep)
    kD = rng.binomial(n, p, size=n_rep)
    ok = 0
    for a, b in zip(kT, kD):
        ok += A2S.tost_bin(int(a), n, int(b), n, tol, conf)["pass"]
    return ok / n_rep


def n_min_symmetric(p: float, target: float = 0.90, n_rep: int = 4000, seed: int = 131, lo: int = 500,
                    hi: int = 200_000, step: int = 100) -> dict:
    """Smallest n per label per bin (multiple of `step`) at which a symmetric process at rate p passes one TOST bin
    with probability >= target (bisection on a monotone MC curve with common random numbers per n)."""
    def prob(n):
        return tost_bin_pass_prob(n, p, n_rep, np.random.default_rng(seed))
    a, b = lo, hi
    if prob(b) < target:
        return {"n_min": None, "p": p, "note": f"not reached at n = {hi}"}
    while b - a > step:
        m = (a + b) // 2
        if prob(m) >= target:
            b = m
        else:
            a = m
    b = int(np.ceil(b / step) * step)
    while prob(b) < target:          # MC curve is not exactly monotone: step up until the target is met at n_min
        b += step
    return {"n_min": b, "p": p, "pass_prob_at_n_min": prob(b), "target": target, "n_rep": n_rep}
