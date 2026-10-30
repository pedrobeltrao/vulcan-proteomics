"""T1 driver, Amendment A3 (DESIGN §6.0c, reviews/amendment_A3_draft.md A3.3): v2_sep as the pre-registered
SECONDARY T1 arm.  Additive only: the v1 primary (A1) path, its comparators and every number the driver computed
before A3 are untouched (A3.1.6 dry run: v1 columns bitwise identical with and without these additions).

What A3 adds to a driver invocation (`--a3`; forced for the sealed T1 look, A3.3.8 "inclusion is not outcome-dependent"):
  methods      + vulcan_v2_sep (numerator of family S; chimeric search + frozen v2 model, per-stratum TDC q_sep)
               + every A3.3.3 comparator not already in the A1 set:
                   C2 = sage_chim2+pred+mokapot (E7b fair comparator, R = 2, frozen by results/2026-10-28_E7b_fair_speed/
                        config.yaml, 7a2d6ee), and sage_chim2 (Sage's own LDA, all ranks).
                 R = 5 is NOT included: E7b's pre-registered rule ran R = 5 beyond D0/D1 only if it beat R = 2 at matched
                 FDP 1% on BOTH D0 and D1; it did not on D1 (6,238.3 vs 6,238.9; E7b NOTES §1), so R = 5 did not qualify.
                 No wide-window Sage chimera variant was pre-registered (E7b config: wide_window false), so v2_sep claim
                 wording carries "Sage chimera at narrow precursor tolerance".
  numerators   "vulcan,vulcan_v2_sep": `noninferiority` stays v1's (primary; it now also lists v1 vs v2_sep, C2 and
               sage_chim2, which A3 Q2 reports in the primary section); `noninferiority_by_numerator.vulcan_v2_sep` is
               family S's input (v2_sep vs every comparator, A3.3.3).
  `a3_secondary` (this module, computed after evaluate from the same outputs; nothing is re-searched):
    loo        A3.3.4 few-cluster backstop: the six (= R) leave-one-run-out POINT ratios (same estimator as
               evaluate.noninferiority with run r removed: run-level IDs summed over the other runs, mean over decoy
               shuffles; matched FDP from the arm-P curve of the other runs) and their minimum, per numerator x
               comparator x level x threshold.  With one run it is undefined (None).
    strata     A3.3.6 entrapment verdicts for v2_sep per MERGED stratum (the merged block that contains chimeric rank
               k in {1, 2, >=3}, per replicate, label-blind merge_adjacent >= 200 pairs, as vulcan.chimeric.compete)
               at alpha 1% / 5%, arms P (15) and C (5), CI conventions (a) (deciding) and (b); counts restricted to the
               stratum with the paired target's acceptance and score from the FULL list (scripts/e7/strata.py
               paired_stratum / combined_stratum, imported unchanged).  Raw strata are reported only.  Production IDs by
               stratum.  The merge maps of the production replicates are recorded.
    void       A3.3.6 void flags (inputs to the claim table, not a decision): any "violated" verdict in arm P or C, in
               the overall v2_sep peptidoform list or any merged stratum, at 1% or 5%, under the deciding convention
               (a) (and, reported, (b)); and "violated" of the v2_sep protein-group overall list (stops S at S3).
    sage_rank  A3.3.6 (reported): arm-P FDP of sage_chim2 and C2 split by the Sage rank of the peptidoform's best PSM
               (scripts/e7b/sage_rank_strata.py, imported unchanged).
    no_ms1     A3.3.2: the number of runs without MS1 hills (their chimeric search equals the closed search).
Family S itself (fixed sequence S1..S5, IUT over every comparator, by the LB AND the LOO minimum) is decided from these
numbers by the analysis agent, per A3.3.5; this module does not decide it.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from . import evaluate as EV

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

SECONDARY = "vulcan_v2_sep"
PRIMARY = "vulcan"
A1_METHODS = ["vulcan", "sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot"]
A3_ADDED_METHODS = [SECONDARY, "sage_chim2", "sage_chim2+pred+mokapot"]
A3_METHODS = A1_METHODS + A3_ADDED_METHODS
A3_COMPARATORS = [m for m in A3_METHODS if m not in (PRIMARY, SECONDARY)]
A3_NUMERATORS = f"{PRIMARY},{SECONDARY}"
SAGE_CHIM_R = (2,)                      # E7b R rule: R = 5 did not qualify (D1)
STRATA = (1, 2, 3)
STRATUM_ALPHAS = (0.01, 0.05)


# ------------------------------------------------------------------------------------------------ LOO backstop
def _own_q(curve, lv, a, runs_idx):
    k = int(np.flatnonzero(np.isclose(EV.GRID, a))[0])
    return float(curve["ids"][lv][runs_idx][:, :, k].sum(0).mean())


def _matched(curve, lv, f, runs_idx):
    c = curve["cnt"][lv][runs_idx]
    x = curve["ids"][lv][runs_idx]
    fdp = EV.fdp_paired(c.sum(0)).mean(axis=(0, 1))
    return EV.matched_ids(fdp, x.sum(0).mean(0), f)[0]


def loo_ratios(curves: dict, numerator: str, comparators, levels=EV.LEVELS, own_alphas=(0.01, 0.05, 0.001),
               fmatch=EV.FMATCH) -> dict:
    """Leave-one-run-out point ratios numerator / comparator (A3.3.4).  -> {comp: {level: {own_q|matched_fdp:
    {thr: {"full": point ratio on all runs, "loo": [ratio without run r], "loo_min": min}}}}}."""
    if numerator not in curves:
        return {}
    R = curves[numerator]["ids"]["sequence"].shape[0]
    allr = np.arange(R)
    subsets = [np.delete(allr, r) for r in range(R)] if R > 1 else []
    out = {}
    for m in comparators:
        if m == numerator or m not in curves:
            continue
        out[m] = {}
        for lv in levels:
            d = {"own_q": {}, "matched_fdp": {}}
            for kind, thrs, fn in (("own_q", own_alphas, _own_q), ("matched_fdp", fmatch, _matched)):
                for t in thrs:
                    def ratio(idx):
                        o, c = fn(curves[numerator], lv, t, idx), fn(curves[m], lv, t, idx)
                        return (o / c) if c else None
                    loo = [ratio(s) for s in subsets]
                    vals = [x for x in loo if x is not None]
                    d[kind][str(t)] = {"full": ratio(allr), "loo": loo,
                                       "loo_min": (min(vals) if vals and len(vals) == len(loo) else None)}
            out[m][lv] = d
    return out


def check_loo_full_equals_noninferiority(loo: dict, ni: dict) -> int:
    """The all-runs point ratio of loo_ratios must equal evaluate.noninferiority's point ratio (same estimator)."""
    n = 0
    for m, byl in loo.items():
        for lv, d in byl.items():
            for kind in ("own_q", "matched_fdp"):
                for t, v in d[kind].items():
                    ref = ni.get(m, {}).get(lv, {}).get(kind, {}).get(t, {}).get("ratio")
                    if ref is None and v["full"] is None:
                        continue
                    if ref is None or v["full"] is None or not np.isclose(ref, v["full"], rtol=1e-12, atol=0):
                        raise AssertionError(f"LOO full ratio != noninferiority ratio ({m}, {lv}, {kind}, {t}): "
                                             f"{v['full']} vs {ref}")
                    n += 1
    return n


# ------------------------------------------------------------------------------------------------ strata (A3.3.6)
def _e7_strata():
    p = os.path.join(ROOT, "scripts", "e7")
    if p not in sys.path:
        sys.path.insert(0, p)
    import strata as S     # scripts/e7/strata.py (paired_stratum, combined_stratum), unchanged
    return S


def _arm(reps_, al, est, r):
    from vulcan.stats import entrapment as ENT
    out = {}
    for cv in ENT.CI_CONVENTIONS:
        e = ENT.evaluate_arm(reps_, al, est, r=r, n_boot=EV.N_BOOT, seed=EV.SEED, ci_convention=cv)
        out[cv] = {x: e[x] for x in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")}
    return out


def strata_v2_sep(work, runs, decoy_seeds, entrap_seeds, log=print) -> dict:
    from vulcan.stats import entrapment as ENT

    from . import ours_chimeric as OC
    from .dbs import parse_rep, replicate_names
    S = _e7_strata()
    reps = replicate_names(decoy_seeds, entrap_seeds)
    kinds = ("merged", "raw")
    P = {kd: {k: {al: [] for al in STRATUM_ALPHAS} for k in STRATA} for kd in kinds}
    C = {kd: {k: {al: [] for al in STRATUM_ALPHAS} for k in STRATA} for kd in kinds}
    ids = {kd: {k: {str(al): [] for al in STRATUM_ALPHAS} for k in STRATA} for kd in kinds}
    merges, rC = {}, None
    for rep in reps:
        kind, j, s = parse_rep(rep)
        rl = EV.RepLabels(os.path.join(work, "index", rep))
        tab, cc = OC.chim_units(os.path.join(work, "ours_chim", rep), rl, runs, "score_v2", return_detail=True,
                                qcol="q_sep")
        merged = np.asarray(cc["merged"])
        assert len(merged) == len(tab)
        raw = tab.stratum.values.copy()
        views = {}
        for k in STRATA:
            lab = "_".join(map(str, cc["merge"][k]))
            tm = tab.copy()
            tm["stratum"] = np.where(merged == lab, k, 0)
            views[("merged", k)] = tm
            tr = tab.copy()
            tr["stratum"] = np.where(raw == k, k, 0)
            views[("raw", k)] = tr
        if kind == "prod":
            merges[rep] = {"merge": cc["merge"], "pair_counts": cc["pair_counts"],
                           "n_group_winners": cc["n_group_winners"],
                           "n_dropped_by_neighbour_filter": cc["n_dropped_by_neighbour_filter"]}
            for (kd, k), tv in views.items():
                fw = ~tv.is_decoy & (tv.fam == 0) & (tv.stratum == k)
                for al in STRATUM_ALPHAS:
                    ids[kd][k][str(al)].append(int((fw & (tv.q <= al)).sum()))
        elif kind == "P":
            for (kd, k), tv in views.items():
                for al in STRATUM_ALPHAS:
                    P[kd][k][al].append(ENT.Replicate(S.paired_stratum(tv, al, k), s, j))
        else:
            rC = rl.n_fwd_pf[-1] / rl.n_fwd_pf[0]
            for (kd, k), tv in views.items():
                for al in STRATUM_ALPHAS:
                    C[kd][k][al].append(ENT.Replicate(S.combined_stratum(tv, al, k), -1, j))
        log(f"a3 strata {rep}")
    out = {"merge_maps_prod": merges, "arm_C_r": rC}
    for kd in kinds:
        out[kd] = {}
        for k in STRATA:
            rec = {"ids_prod_dataset_level": {al: {"per_decoy_seed": v, "mean": float(np.mean(v))}
                                              for al, v in ids[kd][k].items()}, "verdicts": {}}
            for al in STRATUM_ALPHAS:
                v = {"P": _arm(P[kd][k][al], al, "paired", 1.0), "C": _arm(C[kd][k][al], al, "combined", rC)}
                v["overall"] = {cv: ENT.combine_arms({"P": v["P"][cv]["verdict"], "C": v["C"][cv]["verdict"]},
                                                     "static_low_capacity") for cv in ENT.CI_CONVENTIONS}
                rec["verdicts"][str(al)] = v
            out[kd][str(k) if k < 3 else ">=3"] = rec
    return out


def void_flags(metrics: dict, strata: dict) -> dict:
    """A3.3.6 inputs: 'violated' in any arm of any v2_sep peptidoform list (overall + merged strata) at 1% / 5%."""
    from vulcan.stats import entrapment as ENT
    rec = metrics["methods"].get(SECONDARY, {})
    out = {}
    for cv in ENT.CI_CONVENTIONS:
        hits = []
        for al in ("0.01", "0.05"):
            v = rec.get("verdicts", {}).get("peptidoform", {}).get(al, {})
            for arm in ("P", "C"):
                if v.get(arm, {}).get(cv, {}).get("verdict") == "violated":
                    hits.append(f"overall/{arm}/{al}")
            for k, sv in strata.get("merged", {}).items():
                for arm in ("P", "C"):
                    if sv["verdicts"][al][arm][cv]["verdict"] == "violated":
                        hits.append(f"merged_stratum_{k}/{arm}/{al}")
        prot = [al for al in ("0.01", "0.05")
                if rec.get("verdicts", {}).get("protein", {}).get(al, {}).get("P", {}).get(cv, {}).get("verdict")
                == "violated"]
        out[cv] = {"peptidoform_violated": hits, "any_peptidoform_violated": bool(hits),
                   "protein_group_overall_violated": prot}
    out["deciding_convention"] = EV.CI_CONVENTION_DECIDING
    return out


def sage_rank_split(work, decoy_seeds, entrap_seeds, methods, log=print) -> dict:
    """A3.3.6 (reported): Sage-rank split of arm-P FDP for the Sage chimera comparators (E7b read-out, unchanged)."""
    p = os.path.join(ROOT, "scripts", "e7b")
    if p not in sys.path:
        sys.path.insert(0, p)
    import sage_rank_strata as SR        # noqa: E402  scripts/e7b/sage_rank_strata.py, unchanged
    from vulcan.stats import entrapment as ENT

    from .dbs import parse_rep, replicate_names
    S = _e7_strata()
    reps = replicate_names(decoy_seeds, entrap_seeds)
    out = {}
    for m in methods:
        P = {k: {al: [] for al in STRATUM_ALPHAS} for k in (1, 2)}
        ids = {k: {al: [] for al in STRATUM_ALPHAS} for k in (1, 2)}
        for rep in reps:
            kind, j, s = parse_rep(rep)
            if kind == "C":
                continue
            rl = EV.RepLabels(os.path.join(work, "index", rep))
            tab = EV._finish(SR.units(SR.load(m, os.path.join(work, "cmp"), rep)), rl, [])
            for k in (1, 2):
                for al in STRATUM_ALPHAS:
                    if kind == "prod":
                        ids[k][al].append(int((~tab.is_decoy & (tab.fam == 0) & (tab.stratum == k) &
                                               (tab.q <= al)).sum()))
                    else:
                        P[k][al].append(ENT.Replicate(S.paired_stratum(tab, al, k), s, j))
        rec = {}
        for k in (1, 2):
            rk = {}
            for al in STRATUM_ALPHAS:
                e = ENT.evaluate_arm(P[k][al], al, "paired", r=1.0, n_boot=EV.N_BOOT, seed=EV.SEED, ci_convention="a")
                rk[str(al)] = {"prod_ids_mean": float(np.mean(ids[k][al])),
                               **{x: e[x] for x in ("mean_fdp", "ci95", "verdict", "mean_lb", "mean_n_t")}}
            rec["rank_1" if k == 1 else "rank_ge2"] = rk
        out[m] = rec
        log(f"a3 sage rank split {m}")
    return out


def no_ms1_runs(chim_search_meta: dict | None) -> dict:
    """A3.3.2: runs whose chimeric search had no secondary group (no MS1 hills): chimeric search == closed search."""
    if not chim_search_meta:
        return {"n": None}
    first = sorted(chim_search_meta)[0]
    runs = chim_search_meta[first].get("runs", {})
    none = sorted(r for r, v in runs.items() if int(v.get("n_secondary_groups", 0)) == 0)
    return {"n": len(none), "runs": none, "n_runs": len(runs), "from_replicate": first}


def secondary(work, runs, decoy_seeds, entrap_seeds, metrics, curves, methods, meta, log=print) -> dict:
    out = {"amendment": "A3 (DESIGN §6.0c) A3.3; secondary arm, separate error family S; v1 primary unchanged",
           "secondary": SECONDARY, "comparators": [m for m in A3_COMPARATORS if m in methods],
           "sage_chim_R": list(SAGE_CHIM_R),
           "claim_wording_comparator": "Sage chimera at narrow precursor tolerance (no wide-window variant pre-registered)"}
    comps = [m for m in methods if m != SECONDARY]
    out["loo"] = {SECONDARY: loo_ratios(curves, SECONDARY, comps), PRIMARY: loo_ratios(curves, PRIMARY,
                                                                                         [m for m in methods if m != PRIMARY])}
    nb = metrics.get("noninferiority_by_numerator", {})
    out["loo_full_ratio_checks"] = {o: check_loo_full_equals_noninferiority(out["loo"][o], nb.get(o, {}))
                                    for o in out["loo"] if o in nb}
    if SECONDARY in methods and SECONDARY in curves:
        out["strata"] = strata_v2_sep(work, runs, decoy_seeds, entrap_seeds, log)
        out["void_flags"] = void_flags(metrics, out["strata"])
    sr = [m for m in methods if m.startswith("sage_chim")]
    if sr:
        out["sage_rank_split"] = sage_rank_split(work, decoy_seeds, entrap_seeds, sr, log)
    out["no_ms1_runs"] = no_ms1_runs(meta.get("chim_search"))
    return out
