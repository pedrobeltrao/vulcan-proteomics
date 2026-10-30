#!/usr/bin/env python3
"""Amendment A3 step 5 (part a): gates for v2_sep (tiers 1 and 2) and ce_s (licenses nothing) on DN2.  Runs ONCE,
in the single step-5 session, AFTER gates_a2.py has run for v1 and v1's results are committed (A3.6 step 4).

Normative text: reviews/amendment_A3_draft.md rev 1 (DESIGN §6.0c).  Statistics: vulcan.gates.a3_stats (new rows; its
docstring is the specification), vulcan.gates.a2_stats / vulcan.gates.stats / vulcan.stats.g2 (unchanged A1/A2 rules).
Config: config/a3.yaml.  Stricter-reading choices: scripts/a3/IMPLEMENTATION_NOTES.md.

Phases (each refuses to run unless the previous phase's record is committed and unmodified, the tree is the A3
freeze commit and clean, and every model / offsets hash matches):
  --finalize   (label-free; after step 4b)  v2_sep (c-final) membership: start from A2's committed final membership
               (v1's); additionally drop a null spectrum if, for ANY of its units (real or A-syn), the stage-1 winner or
               its pair partner (l1.partner) under the learned (v2), reference (hyperscore) or noRT pipeline is A1.3-near
               the project's reference (same HomologFilter reference set as A2).  ce_s (c-final): v1's minus spectra
               whose ce_s PSM-level winner or pair partner is A1.3-near.  Training-overlap guard (A3.2.3): a DN2
               project whose accession is in train.tsv / train_m7.tsv / train_v2_folds.tsv is dropped from the v2_sep
               and ce_s pools (label-free; the change in the number of projects is recorded).  Writes
               <work>/a3_final_membership.npz + --record json (sha256, label-free counts only).
  --gates      v2_sep rows V1, V2, V6-V10, V12-V14 and ce_s rows (labels read here, step 5 only).
  --g3         V3 (G3(A), recomputed and asserted byte-identical to v1's committed G3(A)), V4 and V5 (exhaustive vs
               production top-K over the CHIMERIC query units: primary group = A2 hypotheses; secondary groups = the
               production virtual spectra of the step-4b search (source-run MS1 groups, isotope errors 0..1); the
               exhaustive search uses exactly these hypotheses (gates_a2.search_exhaustive_hyp on each virtual
               spectrum set); per unit: best target / best decoy hyperscore, exhaustive vs production; pools A∪B all
               units (V4) and A∪B secondary units (V5); A2.8 rule unchanged).
Verdicts: vulcan.gates.a3_stats.tier_outcome / licensed.  Everything is reported whatever the outcome (A3.2.8).

The PURE layer (mixture_units, evaluate_v2sep, evaluate_ce_s, overall) is tested on synthetic pools and in an
end-to-end synthetic-label run (tests/a3/test_a3_gates_e2e.py).  The data layer is marked `pragma: no cover`; A3.6
step 3 requires an end-to-end run on the real DN2 label-free inputs (finalize, counts) and a synthetic-label dry run on
the real pool structure before any per-label statistic.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "a2"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3e"))

import numpy as np  # noqa: E402

from vulcan.gates import a2_stats as A2S  # noqa: E402
from vulcan.gates import a3_stats as A3S  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import compete  # noqa: E402

A2W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")
A3W = os.environ.get("VULCAN_A3_WORK", "/tmp/vulcan_a3")
MIX = os.path.join(A2W, "mix")
S4B = os.path.join(A3W, "step4b", "psms")
PIPES = ("v2", "hyperscore", "v2_noRT")          # learned, reference, noRT (reported)
NULL_PROV = {20: "A", 21: "B", 22: "FS", 1: "OFF"}
SEED = 11
ALPHA = 0.01


# ------------------------------------------------------------------------------------------------ pure: one mixture
def unit_key(spec, group):
    return np.asarray(spec, np.int64) * 64 + np.asarray(group, np.int64)


def stage1(spec, group, score, pf_id):
    """Stage-1 winner (row index) per query unit (spectrum, group), production hash tie-break (pf_id, seed 11)."""
    return compete(unit_key(spec, group), np.asarray(score, np.float64), tiebreak=hash_tiebreak(pf_id, SEED))


def stratum_thresholds(cc, alpha: float = ALPHA) -> dict:
    """tau per merged stratum label: min accepted-target pair-winner score with q_sep <= alpha (+inf if none), and
    the map raw rank stratum {1, 2, 3} -> merged label (this pipeline's own merge map)."""
    res = cc["res"]
    merged = np.asarray(cc["merged"])
    acc = (np.asarray(cc["q_sep"]) <= alpha) & ~np.asarray(res.is_decoy)
    tau = {}
    for m in sorted(set(merged.tolist()) | {"_".join(map(str, v)) for v in cc["merge"].values()}):
        sel = acc & (merged == m)
        tau[m] = float(np.min(res.score[sel])) if sel.any() else float("inf")
    rank_to = {int(k): "_".join(map(str, v)) for k, v in cc["merge"].items()}
    return {"tau": tau, "rank_to_merged": rank_to}


def mixture_units(chim: dict, asyn: dict | None, is_decoy, pair, prov, project, keep, ions, logit_fn=None,
                  pipes=PIPES, alpha: float = ALPHA) -> dict:
    """Per-unit gate inputs of ONE mixture (labels used: step 5 only).

    chim: arrays per PSM of the mixture's chimeric table: spec_idx, group, chim_rank, pf_id, pf_row, and one score
          column per pipeline (pipes); asyn: the same for the A-syn table (or None).
    is_decoy / pair: label arrays indexed by pf_row.  prov / project / keep: per mixture spectrum ((c-final) keep).
    ions: IonSets-like (.sets(rows)) for the neighbour filter.  logit_fn(psm_rows) -> sequence-classifier logit of
          the given chim rows (None -> NaN).
    Returns {"units": {...}, "asyn": {...}, "v6": {...}, "thresholds": {...}} with, for null units of kept spectra:
          pool, project, spec, rank, and per pipeline win_t / score / tail1 / logit (stage-1 winners)."""
    from vulcan.chimeric.compete import chimeric_competition
    spec = np.asarray(chim["spec_idx"], np.int64)
    grp = np.asarray(chim["group"], np.int64)
    rank = np.asarray(chim["chim_rank"], np.int64)
    rows = np.asarray(chim["pf_row"], np.int64)
    pfid = np.asarray(chim["pf_id"])
    dec = np.asarray(is_decoy, bool)[rows]
    pr = np.asarray(pair)[rows]
    prov = np.asarray(prov)
    keep = np.asarray(keep, bool)
    nullspec = np.isin(prov, list(NULL_PROV)) & keep
    out = {"units": {}, "asyn": {}, "v6": {}, "thresholds": {}}
    ref = None
    for p in pipes:
        sc = np.asarray(chim[p], np.float64)
        cc = chimeric_competition(spec, grp, rank, pfid, pr, sc, dec, rows, ions, seed=SEED)
        th = stratum_thresholds(cc, alpha)
        out["thresholds"][p] = {**th, "merge": cc["merge"], "pair_counts": cc["pair_counts"]}
        w = stage1(spec, grp, sc, pfid)
        o = np.argsort(unit_key(spec[w], grp[w]), kind="stable")
        w = w[o]
        m = nullspec[spec[w]]
        w = w[m]
        if ref is None:
            ref = unit_key(spec[w], grp[w])
            out["units"] = {"pool": np.array([NULL_PROV[int(x)] for x in prov[spec[w]]], dtype=object),
                            "project": np.asarray(project, dtype=object)[spec[w]], "spec": spec[w], "group": grp[w],
                            "rank": np.minimum(rank[w], 3)}
        elif not np.array_equal(ref, unit_key(spec[w], grp[w])):
            raise RuntimeError("unit sets differ between pipelines")
        tau_unit = np.array([th["tau"][th["rank_to_merged"][int(r)]] for r in np.minimum(rank[w], 3)])
        out["units"][p] = {"win_t": ~dec[w], "score": sc[w], "tail1": sc[w] >= tau_unit,
                           "logit": logit_fn(w) if logit_fn is not None else np.full(w.size, np.nan),
                           "row": rows[w], "partner_row": None}
        if p != "hyperscore":                 # V6 per learned pipeline (v2 decides; v2_noRT reported)
            w1 = stage1(spec, grp, sc, pfid)
            kept = np.isin(w1, cc["survivors"])
            nwin = np.bincount(spec[w1], minlength=int(spec.max()) + 1 if spec.size else 0)
            sel = nullspec[spec[w1]] & np.isin(prov[spec[w1]], [20, 21]) & (nwin[spec[w1]] >= 2)
            out["v6"][p] = {"kept": kept[sel], "is_target": ~dec[w1][sel], "score": sc[w1][sel],
                            "rank": np.minimum(rank[w1][sel], 3), "chim_rank": rank[w1][sel],
                            "project": np.asarray(project, dtype=object)[spec[w1][sel]]}
        if asyn is not None and len(asyn["spec_idx"]):
            aspec = np.asarray(asyn["spec_idx"], np.int64)
            agrp = np.asarray(asyn["group"], np.int64)
            arank = np.asarray(asyn["chim_rank"], np.int64)
            arow = np.asarray(asyn["pf_row"], np.int64)
            asc = np.asarray(asyn[p], np.float64)
            wa = stage1(aspec, agrp, asc, np.asarray(asyn["pf_id"]))
            wa = wa[np.argsort(unit_key(aspec[wa], agrp[wa]), kind="stable")]
            wa = wa[keep[aspec[wa]] & (prov[aspec[wa]] == 20)]
            ta = np.array([th["tau"][th["rank_to_merged"][int(min(r, 3))]] for r in arank[wa]])
            if "spec" not in out["asyn"]:
                out["asyn"].update({"spec": aspec[wa], "rank": arank[wa],
                                    "project": np.asarray(project, dtype=object)[aspec[wa]]})
            out["asyn"][p] = {"win_t": ~np.asarray(is_decoy, bool)[arow[wa]], "score": asc[wa], "tail1": asc[wa] >= ta,
                              "logit": asyn["logit_fn"](wa) if asyn.get("logit_fn") else np.full(wa.size, np.nan)}
    return out


def concat_units(parts: list) -> dict:
    """Concatenate mixture_units outputs over mixtures."""
    def cat(get):
        xs = [get(q) for q in parts]
        xs = [x for x in xs if x is not None and len(x)]
        return np.concatenate(xs) if xs else np.zeros(0)
    U = {k: cat(lambda q, k=k: q["units"][k]) for k in ("pool", "project", "spec", "rank")}
    for p in PIPES:
        if p in parts[0]["units"]:
            U[p] = {k: cat(lambda q, k=k, p=p: q["units"][p][k]) for k in ("win_t", "score", "tail1", "logit")}
    A = {}
    if any(q["asyn"] for q in parts):
        A = {k: cat(lambda q, k=k: q["asyn"].get(k)) for k in ("rank", "project")}
        for p in PIPES:
            if any(p in q["asyn"] for q in parts):
                A[p] = {k: cat(lambda q, k=k, p=p: q["asyn"][p][k] if p in q["asyn"] else None)
                        for k in ("win_t", "score", "tail1", "logit")}
    V6 = {p: {k: cat(lambda q, k=k, p=p: q["v6"].get(p, {}).get(k)) for k in ("kept", "is_target", "score", "rank",
                                                                             "chim_rank", "project")}
          for p in PIPES if any(p in q["v6"] for q in parts)}
    return {"units": U, "asyn": A, "v6": V6}


# ------------------------------------------------------------------------------------------------ pure: verdicts
def _tail(cfg_row, cluster, L, F, B):
    return A2S.tail_g2_a2(cluster, L["tail1"], L["win_t"], F["tail1"], F["win_t"], B=B,
                          seed=cfg_row.get("bootstrap_seed", 83), n_min=cfg_row["n_min"],
                          margin_frac=cfg_row["margin_abs"], margin_diff=cfg_row["margin_diff"])


def _g1(n_min, L, F):
    return A2S.g1_gate_a2(L["score"], L["logit"], F["score"], F["logit"], n_min=n_min)


def _sub(d, m):
    return {k: np.asarray(v)[m] for k, v in d.items()}


def evaluate_v2sep(D: dict, cfg: dict, g3: dict | None, v1_g3_A: dict | None, v12: dict | None,
                   learned: str = "v2", B: int | None = None) -> dict:
    """Every V row for one learned pipeline ("v2" decides; "v2_noRT" is reported with the same function).

    D: concat_units output.  g3: {"A", "AB_all", "AB_secondary"} g3_gate_a2 dicts of the chimeric G3 (None -> rows
    not evaluable).  v1_g3_A: v1's committed G3(A) dict (for V3's byte-identity).  v12: {"best_t", "best_d",
    "kept_t", "kept_d"} per A spectrum for the learned scorer (None -> not evaluable)."""
    rows_cfg = cfg["rows"]
    U, A = D["units"], D["asyn"]
    V6 = D["v6"][learned]
    B = B or rows_cfg["V1"]["bootstrap_B"]
    pool, proj, rank = np.asarray(U["pool"]), np.asarray(U["project"]), np.asarray(U["rank"])
    L, F = U[learned], U["hyperscore"]
    mA, mB, mFS = pool == "A", pool == "B", pool == "FS"
    mAB = mA | mB
    sec = mAB & (rank >= 2)
    r = {}
    r["V1"] = _tail(rows_cfg["V1"], proj[mA], _sub(L, mA), _sub(F, mA), B)
    r["V2"] = _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mA), _sub(F, mA))
    # V3: G3(A) recomputed; byte-identical to v1's committed G3(A)
    if g3 is None or v1_g3_A is None:
        r["V3"] = {"verdict": A3S.NE, "reason": "G3 not run"}
    else:
        same = json.dumps(g3["A"], sort_keys=True, default=float) == json.dumps(v1_g3_A, sort_keys=True, default=float)
        r["V3"] = dict(g3["A"]) if same else {"verdict": A3S.NE_IMPL, "reason": "G3(A) differs from v1's committed G3(A)"}
        r["V3"]["byte_identical_to_v1"] = same
    r["V4"] = g3["AB_all"] if g3 is not None else {"verdict": A3S.NE, "reason": "G3 not run"}
    r["V5"] = g3["AB_secondary"] if g3 is not None else {"verdict": A3S.NE, "reason": "G3 not run"}
    v6 = A3S.nf_symmetry(V6["kept"], V6["is_target"], V6["score"], V6["rank"],
                         rows_cfg["V6"]["n_min_per_label_per_bin"])
    r["V6.class1"], r["V6.class2"] = v6["class1"], v6["class2"]
    r["V7"] = A3S.g2_original(F["win_t"][mAB], L["win_t"][mAB], F["score"][mAB])
    from vulcan.stats.g2 import g2_sample_size
    n8 = rows_cfg["V8"]["n_min"]
    n8 = int(g2_sample_size()) if n8 == "g2_sample_size" else int(n8)
    r["V8"] = A3S.g2_original(F["win_t"][sec], L["win_t"][sec], F["score"][sec], n_min=n8)
    al = rows_cfg["V9"]
    r["V9.a"] = A3S.g5(L["score"][mAB], L["win_t"][mAB], alarm=al["alarm"])
    r["V9.b"] = A3S.g5(L["score"][mA], L["win_t"][mA], alarm=al["alarm"])
    r["V9.c"] = A3S.g5(L["score"][sec], L["win_t"][sec], alarm=al["alarm"],
                       min_winners=al["parts"]["c"]["min_winners_per_label"])
    fs_t = _tail(rows_cfg["V1"], proj[mFS], _sub(L, mFS), _sub(F, mFS), B)
    fs_g1 = _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mFS), _sub(F, mFS))
    fs_g5 = A3S.g5(L["score"][mFS], L["win_t"][mFS], alarm=al["alarm"], min_winners=rows_cfg["V10"]["fs_g5_min_winners"])
    fs_fail = any(x["verdict"] == A3S.FAIL for x in (fs_t, fs_g1, fs_g5))
    fs_all_ne = all(A3S.is_ne(x["verdict"]) for x in (fs_t, fs_g1, fs_g5))
    r["V10"] = {"verdict": A3S.FAIL if fs_fail else (A3S.NE if fs_all_ne else A3S.PASS), "tail_G2": fs_t, "G1": fs_g1,
                "G5": fs_g5}
    if v12 is None:
        r["V12"] = {"verdict": A3S.NE, "reason": "V12 inputs missing"}
    else:
        r["V12"] = A3S.truncation_symmetry(v12["best_t"], v12["best_d"], v12["kept_t"], v12["kept_d"],
                                           rows_cfg["V12"]["n_min_per_label_per_bin"])
    if A and learned in A:
        AL, AF = A[learned], A["hyperscore"]
        aproj = np.asarray(A["project"])
        r["V13"] = _tail(rows_cfg["V13"], aproj, AL, AF, B)
        r["V14"] = _g1(rows_cfg["V14"]["n_min_top"], AL, AF)
        rep = {}
        for k in (2, 3):
            mk = np.asarray(A["rank"]) == k
            rep[f"V13_rank{k}"] = _tail(rows_cfg["V13"], aproj[mk], _sub(AL, mk), _sub(AF, mk), B)["verdict"]
            rep[f"V14_rank{k}"] = _g1(rows_cfg["V14"]["n_min_top"], _sub(AL, mk), _sub(AF, mk))["verdict"]
    else:
        r["V13"] = {"verdict": A3S.NE, "reason": "no A-syn units"}
        r["V14"] = {"verdict": A3S.NE, "reason": "no A-syn units"}
        rep = {}
    # reported only (never used to pass)
    rep["B_tail_G2_by_stratum"] = {str(k): _tail(rows_cfg["V1"], proj[mB & (rank == k)], _sub(L, mB & (rank == k)),
                                                 _sub(F, mB & (rank == k)), B)["verdict"] for k in (1, 2, 3)}
    rep["B_tail_G2_secondary_pooled"] = _tail(rows_cfg["V1"], proj[mB & (rank >= 2)], _sub(L, mB & (rank >= 2)),
                                              _sub(F, mB & (rank >= 2)), B)["verdict"]
    rep["B_G1_by_stratum"] = {str(k): _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mB & (rank == k)),
                                          _sub(F, mB & (rank == k)))["verdict"] for k in (1, 2, 3)}
    rep["V1_per_project"] = {str(u): _tail(rows_cfg["V1"], proj[mA & (proj == u)], _sub(L, mA & (proj == u)),
                                           _sub(F, mA & (proj == u)), B)["a1_statistics"] for u in np.unique(proj[mA])}
    rep["V1_leave_one_project_out"] = {str(u): _tail(rows_cfg["V1"], proj[mA & (proj != u)], _sub(L, mA & (proj != u)),
                                                     _sub(F, mA & (proj != u)), B)["verdict"] for u in np.unique(proj[mA])}
    if A and learned in A:
        aproj = np.asarray(A["project"])
        rep["V13_leave_one_project_out"] = {str(u): _tail(rows_cfg["V13"], aproj[aproj != u],
                                                          _sub(A[learned], aproj != u), _sub(A["hyperscore"], aproj != u),
                                                          B)["verdict"] for u in np.unique(aproj)}
    rep["V13_per_project"] = {} if not (A and learned in A) else {
        str(u): _tail(rows_cfg["V13"], np.asarray(A["project"])[np.asarray(A["project"]) == u],
                      _sub(A[learned], np.asarray(A["project"]) == u), _sub(A["hyperscore"], np.asarray(A["project"]) == u),
                      B)["a1_statistics"] for u in np.unique(np.asarray(A["project"]))}
    c2 = np.asarray(V6["rank"]) >= 2
    rep["V6_class2_by_rank"] = {str(k): A3S.nf_symmetry_class(
        np.asarray(V6["kept"])[c2 & (np.asarray(V6["chim_rank"]) == k if k == 2 else np.asarray(V6["chim_rank"]) >= 3)],
        np.asarray(V6["is_target"])[c2 & (np.asarray(V6["chim_rank"]) == k if k == 2 else np.asarray(V6["chim_rank"]) >= 3)],
        np.asarray(V6["score"])[c2 & (np.asarray(V6["chim_rank"]) == k if k == 2 else np.asarray(V6["chim_rank"]) >= 3)],
        rows_cfg["V6"]["n_min_per_label_per_bin"])["verdict"] for k in (2, 3)}
    rep["V8_by_stratum"] = {str(k): A3S.g2_original(F["win_t"][mAB & (rank == k)], L["win_t"][mAB & (rank == k)],
                                                    F["score"][mAB & (rank == k)])["verdict"] for k in (2, 3)}
    rep["B_G1_secondary_pooled"] = _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mB & (rank >= 2)),
                                       _sub(F, mB & (rank >= 2)))["verdict"]
    mO = pool == "OFF"
    rep["A1_offsets"] = {"tail_G2": _tail(rows_cfg["V1"], proj[mO], _sub(L, mO), _sub(F, mO), B)["verdict"],
                         "G1": _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mO), _sub(F, mO))["verdict"],
                         "G5": A3S.g5(L["score"][mO], L["win_t"][mO], alarm=al["alarm"])["verdict"]}
    rep["V7_by_stratum"] = {str(k): A3S.g2_original(F["win_t"][mAB & (rank == k)], L["win_t"][mAB & (rank == k)],
                                                    F["score"][mAB & (rank == k)])["verdict"] for k in (1, 2, 3)}
    r["reported"] = rep
    tiers = cfg["tiers"]
    v = {k: r[k]["verdict"] for k in tiers["tier1"]}
    t1 = A3S.tier_outcome(v, fs_keys=tiers["fs_rows"])
    t2 = A3S.tier_outcome({k: r[k]["verdict"] for k in tiers["tier2"]})
    r["tier1"], r["tier2"] = t1, t2
    r["licensed"] = A3S.licensed(t1, t2)
    return r


def evaluate_ce_s(R: dict, P: dict, cfg: dict, g3_v1: dict | None, v12: dict | None, B: int | None = None) -> dict:
    """A3.2.7: ce_s on the closed pipeline (A2's v1 analysis with the ce_s logit).  R: pool / project per null
    spectrum; P["ce_s"], P["hyperscore"]: win_t / score / tail1 / logit aligned with R.  Licenses nothing."""
    rows_cfg = cfg["rows"]
    B = B or rows_cfg["V1"]["bootstrap_B"]
    pool, proj = np.asarray(R["pool"]), np.asarray(R["project"])
    L, F = P["ce_s"], P["hyperscore"]
    mA, mB, mFS = pool == "A", pool == "B", pool == "FS"
    mAB = mA | mB
    al = rows_cfg["V9"]["alarm"]
    r = {"tail_G2_A": _tail(rows_cfg["V1"], proj[mA], _sub(L, mA), _sub(F, mA), B),
         "G1_A": _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mA), _sub(F, mA)),
         "G2_original_AB": A3S.g2_original(F["win_t"][mAB], L["win_t"][mAB], F["score"][mAB]),
         "G2_original_A": A3S.g2_original(F["win_t"][mA], L["win_t"][mA], F["score"][mA]),
         "G5_AB": A3S.g5(L["score"][mAB], L["win_t"][mAB], alarm=al),
         "G5_A": A3S.g5(L["score"][mA], L["win_t"][mA], alarm=al)}
    fs_t = _tail(rows_cfg["V1"], proj[mFS], _sub(L, mFS), _sub(F, mFS), B)
    fs_g1 = _g1(rows_cfg["V2"]["n_min_top"], _sub(L, mFS), _sub(F, mFS))
    fs_g5 = A3S.g5(L["score"][mFS], L["win_t"][mFS], alarm=al, min_winners=rows_cfg["V10"]["fs_g5_min_winners"])
    fs_fail = any(x["verdict"] == A3S.FAIL for x in (fs_t, fs_g1, fs_g5))
    r["FS"] = {"verdict": A3S.FAIL if fs_fail else (A3S.NE if all(A3S.is_ne(x["verdict"]) for x in (fs_t, fs_g1, fs_g5))
                                                   else A3S.PASS), "tail_G2": fs_t, "G1": fs_g1, "G5": fs_g5}
    r["V12_ce_s"] = (A3S.truncation_symmetry(v12["best_t"], v12["best_d"], v12["kept_t"], v12["kept_d"],
                                             rows_cfg["V12"]["n_min_per_label_per_bin"]) if v12 is not None
                     else {"verdict": A3S.NE, "reason": "V12 inputs missing"})
    if g3_v1 is not None:
        r["G3_A_v1"], r["G3_AB_v1"] = g3_v1["A"], g3_v1["AB"]
    else:
        r["G3_A_v1"] = r["G3_AB_v1"] = {"verdict": A3S.NE, "reason": "v1 G3 missing"}
    ver = {k: r[k]["verdict"] for k in ("tail_G2_A", "G1_A", "G2_original_AB", "G2_original_A", "G5_AB", "G5_A",
                                        "FS", "V12_ce_s", "G3_A_v1", "G3_AB_v1")}
    r["verdict"] = A3S.tier_outcome(ver, fs_keys=("FS",))
    r["licenses"] = "nothing (A3.2.7); final for ce_s under P4"
    return r


def v12_inputs(spec_a, cand_spec, cand_score, cand_is_decoy, cand_in_top5) -> dict:
    """Per A spectrum (spec_a): exhaustive best target / decoy score of the given scorer, and kept = SOME production
    top-5 candidate of that label has a score equal to the exhaustive best (A2.8 / g3_gate_a2 counts recall by score
    equality, so ties cannot depend on candidate order; independent review minor 6)."""
    spec_a = np.asarray(spec_a, np.int64)
    pos = {int(s): i for i, s in enumerate(spec_a)}
    out = {"best_t": np.full(spec_a.size, -np.inf), "best_d": np.full(spec_a.size, -np.inf),
           "kept_t": np.zeros(spec_a.size, bool), "kept_d": np.zeros(spec_a.size, bool)}
    cs = np.asarray(cand_spec, np.int64)
    sc = np.asarray(cand_score, np.float64)
    dec = np.asarray(cand_is_decoy, bool)
    k5 = np.asarray(cand_in_top5, bool)
    inA = np.isin(cs, spec_a)
    for lab, d in (("t", ~dec), ("d", dec)):
        m = d & inA
        if not m.any():
            continue
        ii = np.array([pos[int(x)] for x in cs[m]])
        best = np.full(spec_a.size, -np.inf)
        np.maximum.at(best, ii, sc[m])
        top = np.full(spec_a.size, -np.inf)
        if (m & k5).any():
            np.maximum.at(top, np.array([pos[int(x)] for x in cs[m & k5]]), sc[m & k5])
        out["best_" + lab] = best
        out["kept_" + lab] = np.isfinite(best) & (top == best)
    return out


def overall(v2: dict, ce: dict | None) -> dict:
    return {"v2_sep": {"tier1": v2["tier1"], "tier2": v2["tier2"], "licensed": v2["licensed"]},
            "ce_s": (ce or {}).get("verdict"), "production_rule": "A3.2.6 (needs v1's A2 verdict, A3.1.8 and T1)"}


# ------------------------------------------------------------------------------------------------ guards
def committed_unmodified(path: str) -> bool:
    rel = os.path.relpath(path, ROOT)
    ls = subprocess.run(["git", "-C", ROOT, "ls-files", "--error-unmatch", rel], capture_output=True)
    st = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", rel], capture_output=True, text=True)
    return ls.returncode == 0 and not st.stdout.strip()


def require_committed(path: str, what: str) -> dict:
    if not committed_unmodified(path):
        raise SystemExit(f"gates_a3: {what} ({path}) is not committed/unmodified")
    return json.load(open(path))


# ------------------------------------------------------------------------------------------------ data layer
def _labels():  # pragma: no cover
    import m3_common as C
    import m3_data as MD
    return MD.labels(C.index_dir("human_contam", "prod"))


def _src_proj():  # pragma: no cover
    info = json.load(open(os.path.join(MIX, "mix_info.json")))
    return {s["id"]: s["accession"] for s in info["sources"]}, {s["id"]: s["run"] for s in info["sources"]}


def _read(path, cols=None):  # pragma: no cover
    import pyarrow.parquet as pq
    t = pq.read_table(path, columns=cols)
    return {c: t.column(c).to_numpy(zero_copy_only=False) for c in t.column_names}


def training_overlap(cfg) -> set:  # pragma: no cover
    import csv
    acc = set()
    for f in cfg["pipeline_v2_sep"]["training_overlap_files"]:
        rd = csv.DictReader(open(os.path.join(ROOT, f)), delimiter="\t")
        if "accession" not in (rd.fieldnames or []):
            raise SystemExit(f"gates_a3: training-overlap file {f} has no 'accession' column")
        for r in rd:
            acc.add(r["accession"])
    return acc


def finalize(cfg, record, a2_record):  # pragma: no cover (data)
    """Label-free v2_sep and ce_s (c-final) membership on top of A2's committed final membership (its npz must match
    the committed A2 record's sha256; independent review minor 4)."""
    a2rec = require_committed(a2_record, "A2 final membership record")
    if __import__("hashlib").sha256(open(os.path.join(A2W, "final_membership.npz"), "rb").read()).hexdigest() \
            != a2rec["sha256"]:
        raise SystemExit("gates_a3: A2 final membership on disk differs from its committed record")
    import own_search as OS
    from vulcan.gates import a2_config
    from vulcan.gates.homolog import HomologFilter
    from vulcan.index.build import L1Index
    import m3_common as C
    cfg2 = a2_config.load()
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    fm = np.load(os.path.join(A2W, "final_membership.npz"))
    src_proj, _ = _src_proj()
    overlap = training_overlap(cfg)
    refs = {}
    for acc in sorted(set(src_proj.values())):
        upid = json.load(open(os.path.join(A2W, "admission", acc, "counts.json")))["upid"]
        _, src = OS.build_own_index(upid, cfg2)
        refs[acc] = HomologFilter(src["own"] + src["contam"])
    keep_v2, keep_ce, info = {}, {}, {"dropped_v2": 0, "dropped_ce_s": 0, "overlap_projects": sorted(overlap &
                                                                                                     set(src_proj.values()))}
    ce = _ce_model(cfg)
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        srun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
        k2 = fm[name].copy()
        chim = _read(os.path.join(S4B, name + ".chim.scored.parquet"))
        rows_by_spec = []
        for p in PIPES:
            w = stage1(chim["spec_idx"], chim["group"], chim[p], chim["pf_id"])
            rows_by_spec.append((chim["spec_idx"][w], chim["pf_row"][w]))
        fa = os.path.join(S4B, name + ".asyn.scored.parquet")
        if os.path.exists(fa):
            asyn = _read(fa)
            for p in PIPES:
                w = stage1(asyn["spec_idx"], asyn["group"], asyn[p], asyn["pf_id"])
                rows_by_spec.append((asyn["spec_idx"][w], asyn["pf_row"][w]))
        sp_all = np.concatenate([s for s, _ in rows_by_spec])
        rw_all = np.concatenate([r for _, r in rows_by_spec])
        sp_all = np.concatenate([sp_all, sp_all])
        rw_all = np.concatenate([rw_all, np.asarray(l1.partner)[rw_all]])
        near = np.zeros(prov.size, bool)
        for acc, hf in refs.items():
            ks = [k for k, v in src_proj.items() if v == acc]
            m = np.isin(srun[sp_all], ks) & np.isin(prov[sp_all], list(NULL_PROV))
            if m.any():
                rr = np.unique(rw_all[m])
                bad = rr[hf.near([l1.sequence(int(r)) for r in rr])]
                hit = m & np.isin(rw_all, bad)
                near[sp_all[hit]] = True
        drop_overlap = np.isin(srun, [k for k, v in src_proj.items() if v in overlap])
        kv = k2 & ~near & ~(drop_overlap & np.isin(prov, list(NULL_PROV)))
        keep_v2[name] = kv
        info["dropped_v2"] += int((k2 & ~kv).sum())
        kc = k2 & ~(drop_overlap & np.isin(prov, list(NULL_PROV)))
        kc &= ~_ce_near(name, prov, srun, src_proj, refs, l1, ce)
        keep_ce[name] = kc
        info["dropped_ce_s"] += int((k2 & ~kc).sum())
    np.savez(os.path.join(A3W, "a3_final_membership.npz"), **{"v2_" + k: v for k, v in keep_v2.items()},
             **{"ce_" + k: v for k, v in keep_ce.items()})
    rec = {"sha256": __import__("hashlib").sha256(open(os.path.join(A3W, "a3_final_membership.npz"), "rb").read()).hexdigest(),
           **info}
    json.dump(rec, open(record, "w"), indent=1)
    return rec


def _ce_model(cfg):  # pragma: no cover
    from vulcan.encoders.score import load_model
    return load_model(os.path.join(ROOT, cfg["scorers"]["ce_s"]["model"]))


_CE_CACHE = {}


def ce_scores(name, ce):  # pragma: no cover
    """ce_s logits of the step-4 closed top-5 table of mixture `name` (as M7; deterministic, label-free)."""
    if name not in _CE_CACHE:
        import m3_common as C
        from vulcan.encoders.score import score_table
        from vulcan.index.build import L1Index
        l1 = L1Index(C.index_dir("human_contam", "prod"))
        psm = _read(os.path.join(A2W, "psms", "A2mix", name + ".psms.parquet"))
        mix = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
        lg, _ = score_table(ce, l1, {"pf_row": psm["pf_row"], "charge": psm["charge"], "spec_idx": psm["spec_idx"]}, mix)
        _CE_CACHE[name] = (psm, np.asarray(lg, np.float64))
    return _CE_CACHE[name]


def _ce_near(name, prov, srun, src_proj, refs, l1, ce):  # pragma: no cover
    from vulcan.stats.tdc import psm_competition
    psm, lg = ce_scores(name, ce)
    w = psm_competition(psm["spec_idx"], lg, tiebreak=hash_tiebreak(psm["pf_id"], SEED))
    sp_w, rw = psm["spec_idx"][w], psm["pf_row"][w]
    near = np.zeros(prov.size, bool)
    for acc, hf in refs.items():
        ks = [k for k, v in src_proj.items() if v == acc]
        m = np.isin(srun[sp_w], ks) & np.isin(prov[sp_w], list(NULL_PROV))
        if m.any():
            rr = np.unique(np.r_[rw[m], np.asarray(l1.partner)[rw[m]]])
            bad = rr[hf.near([l1.sequence(int(r)) for r in rr])]
            near[sp_w[m & (np.isin(rw, bad) | np.isin(np.asarray(l1.partner)[rw], bad))]] = True
    return near


def collect(cfg):  # pragma: no cover (data; labels read here, step 5 only)
    import m3_common as C
    import seqclf
    from vulcan.chimeric.compete import IonSets
    from vulcan.index.build import L1Index
    clf = seqclf.load()
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    ions = IonSets(l1)
    lab = _labels()
    fm = np.load(os.path.join(A3W, "a3_final_membership.npz"))
    src_proj, src_run = _src_proj()
    parts, v12p = [], {"spec": [], "cs": [], "v2": [], "dec": [], "k5": []}
    v12c = {"spec": [], "cs": [], "ce": [], "dec": [], "k5": []}
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        srun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
        project = np.array([src_proj.get(int(s), "host") for s in srun], dtype=object)
        chim = _read(os.path.join(S4B, name + ".chim.scored.parquet"))
        fa = os.path.join(S4B, name + ".asyn.scored.parquet")
        asyn = _read(fa) if os.path.exists(fa) else None
        if asyn is not None:
            asyn["logit_fn"] = lambda w, a=asyn: clf.logit_rows(l1, a["pf_row"][w], a["rt_pred"][w])
        parts.append(mixture_units(chim, asyn, lab["is_decoy"], lab["tdc_pair"], prov, project, fm["v2_" + name], ions,
                                   logit_fn=lambda w, c=chim: clf.logit_rows(l1, c["pf_row"][w], c["rt_pred"][w])))
        mi = len(parts) - 1                          # mixture index -> unique spectrum keys across mixtures
        v = _read(os.path.join(S4B, name + ".v12.parquet"))
        ka = fm["v2_" + name]
        m = ka[v["spec_idx"]]
        v12p["cs"].append((mi << 32) + v["spec_idx"][m].astype(np.int64))
        v12p["v2"].append(v["v2"][m])
        v12p["dec"].append(lab["is_decoy"][v["pf_row"][m]])
        v12p["k5"].append(v["in_top5"][m])
        v12p["spec"].append((mi << 32) + np.flatnonzero((prov == 20) & ka).astype(np.int64))
        kc = fm["ce_" + name]                        # ce_s's own (c-final) membership (independent review M2)
        mc = kc[v["spec_idx"]]
        v12c["cs"].append((mi << 32) + v["spec_idx"][mc].astype(np.int64))
        v12c["ce"].append(v["ce_s"][mc])
        v12c["dec"].append(lab["is_decoy"][v["pf_row"][mc]])
        v12c["k5"].append(v["in_top5"][mc])
        v12c["spec"].append((mi << 32) + np.flatnonzero((prov == 20) & kc).astype(np.int64))
    D = concat_units(parts)
    cat = {k: np.concatenate(x) for k, x in v12p.items() if x}
    v12_v2 = v12_inputs(cat["spec"], cat["cs"], cat["v2"], cat["dec"], cat["k5"])
    cc_ = {k: np.concatenate(x) for k, x in v12c.items() if x}
    v12_ce = v12_inputs(cc_["spec"], cc_["cs"], cc_["ce"], cc_["dec"], cc_["k5"])
    return D, v12_v2, v12_ce


def collect_ce_s(cfg):  # pragma: no cover (data)
    import m3_common as C
    import seqclf
    from vulcan.gates import stats as GS
    from vulcan.index.build import L1Index
    from vulcan.stats.tdc import double_competition, psm_competition
    clf = seqclf.load()
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    lab = _labels()
    fm = np.load(os.path.join(A3W, "a3_final_membership.npz"))
    src_proj, _ = _src_proj()
    ce = _ce_model(cfg)
    R = {"pool": [], "project": []}
    P = {p: {k: [] for k in ("win_t", "score", "tail1", "logit")} for p in ("ce_s", "hyperscore")}
    feat_dir = os.path.join(A2W, "psms", "A2mix")
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        srun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
        psm, lg = ce_scores(name, ce)
        rtp = _read(os.path.join(feat_dir, name + ".feat.parquet"), ["rt_pred"])["rt_pred"]
        rows = psm["pf_row"]
        dec = lab["is_decoy"][rows]
        tb = hash_tiebreak(psm["pf_id"], SEED)
        scores = {"ce_s": lg, "hyperscore": psm["hyperscore"].astype(np.float64)}
        w0 = psm_competition(psm["spec_idx"], scores["hyperscore"], tiebreak=tb)
        w0 = w0[np.argsort(psm["spec_idx"][w0], kind="stable")]
        sp_w = psm["spec_idx"][w0]
        m = fm["ce_" + name][sp_w] & np.isin(prov[sp_w], list(NULL_PROV))
        R["pool"].append(np.array([NULL_PROV[int(x)] for x in prov[sp_w][m]], dtype=object))
        R["project"].append(np.array([src_proj.get(int(s), "host") for s in srun[sp_w][m]], dtype=object))
        for p, s in scores.items():
            w = psm_competition(psm["spec_idx"], s, tiebreak=tb)
            w = w[np.argsort(psm["spec_idx"][w], kind="stable")]
            r = double_competition(psm["spec_idx"].astype(np.int64), psm["pf_id"], lab["tdc_pair"][rows], s, dec,
                                   seed=SEED, tiebreak=tb)
            thr = GS.tail_threshold(r.score, r.is_decoy, r.q, ALPHA)
            P[p]["win_t"].append(~dec[w][m])
            P[p]["score"].append(s[w][m])
            P[p]["tail1"].append((s[w] >= thr)[m])
            P[p]["logit"].append(clf.logit_rows(l1, rows[w][m], rtp[w][m]))
    return ({k: np.concatenate(v) for k, v in R.items()},
            {p: {k: np.concatenate(v) for k, v in d.items()} for p, d in P.items()})


def g3_chimeric(cfg):  # pragma: no cover (GPU/data)
    """V3/V4/V5: exhaustive vs production top-K per chimeric query unit (best target / decoy hyperscore)."""
    from dataclasses import replace

    import m3_common as C
    import m3_data as MD
    import torch
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m2"))
    from e2_correctness import _best_per_label, _rescored
    from gates_a2 import search_exhaustive_hyp
    import step4b_chimeric as S4
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import make_hypotheses, pre_hypotheses
    from vulcan.index.build import L1Index
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    isdec = MD.labels(idx)["is_decoy"]
    p0 = E.SearchParams()
    p2 = replace(p0, isotope_errors=(0, 1))
    gi = E.GpuIndex(l1, p0, C.STATIC, C.VAR)
    fm2 = np.load(os.path.join(A2W, "final_membership.npz"))
    fm3 = np.load(os.path.join(A3W, "a3_final_membership.npz"))
    info, src = S4.admitted_sources()
    hosts = {"store": [C.E2.store_path(r) for r in C.DEV["D0"]["raw"]]}
    hosts.update({h: os.path.basename(C.E2.store_path(r))[:-len(".vspec.arrow")] for h, r in enumerate(C.DEV["D0"]["raw"])})
    accA = {k: [] for k in ("bt_x", "bd_x", "bt_k", "bd_k", "pool")}       # v1-identical G3 (A2 membership) for V3
    acc = {k: [] for k in ("bt_x", "bd_x", "bt_k", "bd_k", "pool", "sec")}  # chimeric units, v2 membership
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        sp = pickle.load(open(f, "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        srun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
        single = np.load(os.path.join(MIX, name + ".single.npy"))
        valid, _ = pre_hypotheses(sp, gi, p0)
        hyp = a2_hypotheses(sp, valid, gi.mass32_np, p0, single)
        n = len(sp)
        sel2 = np.flatnonzero(np.isin(prov, [20, 21]) & fm2[name])
        _, _, dx = search_exhaustive_hyp(sp, gi, l1, p0, hyp, sel2)
        sx, hx, rx = _rescored(dx, l1, p0)
        btx, bdx = _best_per_label(sx, hx, isdec[rx], n)
        _, _, dk = E.search(sp, gi, l1, p0, pre_hyp=(valid, hyp))
        sk, hk, rk = _rescored(dk, l1, p0)
        btk, bdk = _best_per_label(sk, hk, isdec[rk], n)
        for k, v in (("bt_x", btx), ("bd_x", bdx), ("bt_k", btk), ("bd_k", bdk)):
            accA[k].append(v[sel2])
        accA["pool"].append(np.where(prov[sel2] == 20, "A", "B"))
        sel3 = np.flatnonzero(np.isin(prov, [20, 21]) & fm3["v2_" + name])
        for k, v in (("bt_x", btx), ("bd_x", bdx), ("bt_k", btk), ("bd_k", bdk)):
            acc[k].append(v[sel3])
        acc["pool"].append(np.where(prov[sel3] == 20, "A", "B"))
        acc["sec"].append(np.zeros(sel3.size, bool))
        groups, _ = S4.mixture_groups(sp, prov, srun, src, hosts)
        sec = groups["sec"]
        ms = fm3["v2_" + name][sec["spec"]] & np.isin(prov[sec["spec"]], [20, 21])
        if ms.any():
            secm = {k: v[ms] for k, v in sec.items()}
            vsp = S4.subset_spectra(sp, secm["spec"], pmz=secm["mz"], charge=secm["charge"])
            vv, _ = pre_hypotheses(vsp, gi, p2)
            vh = make_hypotheses(vsp, vv, gi.mass32_np, p2)
            nv = len(vsp)
            _, _, vdx = search_exhaustive_hyp(vsp, gi, l1, p2, vh, np.arange(nv))
            a1, a2_, a3_ = _rescored(vdx, l1, p2)
            vbtx, vbdx = _best_per_label(a1, a2_, isdec[a3_], nv)
            _, _, vdk = E.search(vsp, gi, l1, p2, pre_hyp=(vv, vh))
            b1, b2, b3 = _rescored(vdk, l1, p2)
            vbtk, vbdk = _best_per_label(b1, b2, isdec[b3], nv)
            for k, v in (("bt_x", vbtx), ("bd_x", vbdx), ("bt_k", vbtk), ("bd_k", vbdk)):
                acc[k].append(v)
            acc["pool"].append(np.where(prov[secm["spec"]] == 20, "A", "B"))
            acc["sec"].append(np.ones(nv, bool))
        torch.cuda.empty_cache()
    gg = cfg["rows"]["V4"]
    AA = {k: np.concatenate(v) for k, v in accA.items()}
    mA = AA["pool"] == "A"
    a2cfg = __import__("vulcan.gates.a2_config", fromlist=["load"]).load()["gates"]["g3"]
    g3_A = A2S.g3_gate_a2(*(np.where(mA, AA[k], -np.inf) for k in ("bt_x", "bd_x", "bt_k", "bd_k")),
                          n_min=a2cfg["n_min_per_label_per_bin"], tol=a2cfg["tol"])
    X = {k: np.concatenate(v) for k, v in acc.items()}
    s = X["sec"]
    return {"A": g3_A,
            "AB_all": A2S.g3_gate_a2(X["bt_x"], X["bd_x"], X["bt_k"], X["bd_k"], n_min=gg["n_min"], tol=gg["tol"]),
            "AB_secondary": A2S.g3_gate_a2(*(np.where(s, X[k], -np.inf) for k in ("bt_x", "bd_x", "bt_k", "bd_k")),
                                           n_min=cfg["rows"]["V5"]["n_min"], tol=cfg["rows"]["V5"]["tol"])}


def main():  # pragma: no cover
    from vulcan.gates import a3_config
    ap = argparse.ArgumentParser()
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--g3", action="store_true")
    ap.add_argument("--record", required=True, help="results/<date>_A2_gates/a3_final_membership.json")
    ap.add_argument("--step4b-record", required=True, help="committed step-4b record (A3.2.2)")
    ap.add_argument("--a2-record", required=True, help="committed A2 final_membership.json (results/<date>_A2_admission)")
    ap.add_argument("--v1-results", help="committed gates_a2.py output of v1 (required for --gates/--g3, A3.6 4b)")
    ap.add_argument("--out")
    a = ap.parse_args()
    cfg = a3_config.load()
    a3_config.require_freeze_commit(cfg)
    a3_config.check_models(cfg)
    a3_config.check_offsets(cfg)
    require_committed(a.step4b_record, "step-4b record")
    if a.finalize:
        rec = finalize(cfg, a.record, a.a2_record)
        print(json.dumps(rec))
        return
    rec = require_committed(a.record, "A3 final membership record")
    if __import__("hashlib").sha256(open(os.path.join(A3W, "a3_final_membership.npz"), "rb").read()).hexdigest() \
            != rec["sha256"]:
        raise SystemExit("gates_a3: final membership on disk differs from the committed record")
    if a.gates and not a.g3:
        raise SystemExit("gates_a3: --gates requires --g3 in the same run (V3-V5 belong to tier 1/2)")
    if not a.v1_results:
        raise SystemExit("gates_a3: --v1-results (committed) is required: v1's results are committed first (A3.6)")
    v1 = require_committed(a.v1_results, "v1 gates_a2 results")
    out = {}
    g3 = g3_chimeric(cfg) if a.g3 else None
    if g3 is not None:
        out["g3_chimeric"] = g3
    if a.gates:
        D, v12_v2, v12_ce = collect(cfg)
        v1g3 = v1.get("g3")
        out["v2_sep"] = evaluate_v2sep(D, cfg, g3, (v1g3 or {}).get("A"), v12_v2, learned="v2")
        out["v2_sep_noRT_reported"] = evaluate_v2sep(D, cfg, g3, (v1g3 or {}).get("A"), None, learned="v2_noRT")
        out["v2_sep_noRT_reported"]["note"] = "reported only; V12 is not computed for noRT (exhaustive scoring is v2)"
        R, P = collect_ce_s(cfg)
        out["ce_s"] = evaluate_ce_s(R, P, cfg, v1g3, v12_ce)
        out["overall"] = overall(out["v2_sep"], out["ce_s"])
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
