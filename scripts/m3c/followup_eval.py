#!/usr/bin/env python3
"""3c follow-up evaluation (dev only; estimators, units and verdict code unchanged).

    python scripts/m3c/followup_eval.py ours
        3c eval (scripts/m3c/eval_itcid.py, unchanged) for hyperscore + the frozen itcid_linear_pred on the dev
        replicates whose features were computed from the pinned prediction store ($VULCAN_3C_FU/psms/<ds>_<kind>)
        -> followup/eval_pinned.json
    python scripts/m3c/followup_eval.py ratios
        ours / comparator in SEQUENCE units (I/L-collapsed; the unit of 3c's Sage rows), as bench/t1/evaluate.py:
        R = run-level IDs summed over runs (a sequence counts in run r if accepted at the dataset-level q and it has
        >= 1 top-1 PSM in r); (i) at each method's own q <= alpha; (ii) at matched arm-P FDP f (each method's own-q
        threshold moved until its arm-P paired FDP estimate, mean over the 3 entrapment seeds of run-summed paired
        counts at sequence level (pair key = target sequence), reaches f; first crossing, linear interpolation);
        95% CI = percentile bootstrap (2000, seed 11) over runs (same draw for both methods; the FDP curve is
        re-estimated in every draw).  3c has one decoy set (pseudo-reverse for ours, M1e rev_ for Sage), so J = 1.
        Comparators: sage+pred+mokapot (strongest, A1.9), sage+mokapot(seeded), sage (LDA).
        -> followup/ratios.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
import common as C  # noqa: E402
import pinned_cache as PC  # noqa: E402

while C.ROOT in sys.path:
    sys.path.remove(C.ROOT)
sys.path.insert(0, C.ROOT)                        # before scripts/m3c (its bench.py shadows the bench package)
import types as _types  # noqa: E402
_bench = _types.ModuleType("bench")                # scripts/m3c/bench.py would shadow the namespace package
_bench.__path__ = [os.path.join(C.ROOT, "bench")]
sys.modules.setdefault("bench", _bench)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCORER = "itcid_linear_pred"
PKINDS = ("E2001", "E2002", "E2003")
SAGE_P = {"E2001": "P2001", "E2002": "P2002", "E2003": "P2003"}


def ours_eval(datasets=("D4", "D10cf")):
    import eval_itcid as EV
    C.psm_dir = PC.fu_psm_dir                       # features from the pinned store; PSM tables are 3c's (symlinks)
    sys.argv = ["eval_itcid.py", "--datasets", *datasets, "--scorers", "hyperscore", SCORER, "--out",
                os.path.join(PC.FU_RES, "eval_pinned.json")]
    EV.main()


# ------------------------------------------------------------------------------------------------ unit tables
def ours_table(ds, kind, runs):
    """Sequence-unit table (bench.t1.evaluate columns) of our frozen scorer for one replicate search."""
    import evaluate as M3E                           # scripts/m3/evaluate.py (Search, double competition wiring)
    from vulcan.rescore.static import StaticModel
    M3E.C = C
    C.psm_dir = PC.fu_psm_dir
    M3E.FAM = {k: v["family"] for k, v in C.DEV.items()}
    s = M3E.Search(ds, kind)
    model = StaticModel.load(os.path.join(C.RES, "models", SCORER + ".pkl"))
    res = s.compete(s.static_score(model))
    rows = s.psm["pf_row"][res.row]
    from vulcan.index.build import L1Index
    l1 = L1Index(s.index)
    seq = np.array([l1.sequence(int(r)).replace("I", "L") for r in rows], dtype=object)
    fam = s.lab["family"][rows] % 2                   # 0 T, 1 E (decoys carry their partner's family + 2)
    tseq = np.array([t.replace("I", "L") for t in s.lab["target_seq"][rows]], dtype=object)
    rix = {r: i for i, r in enumerate(runs)}
    w = res.psm_winners
    bit = np.array([1 << rix[r] for r in s.psm["run"][w]], np.int64)
    ub = pd.DataFrame({"u": s.unit[w], "b": bit}).groupby("u").b.agg(lambda x: int(np.bitwise_or.reduce(x.values)))
    df = pd.DataFrame({"seq": seq, "is_decoy": res.is_decoy, "q": res.q, "score": res.score, "fam": fam,
                       "target_seq": tseq,
                       "runbits": ub.reindex(s.unit[res.row]).fillna(0).astype(np.int64).values})
    return seq_units(df)


def seq_units(df):
    g = df.groupby(["seq", "is_decoy"], sort=False)
    return g.agg(q=("q", "min"), score=("score", "max"),
                 runbits=("runbits", lambda x: int(np.bitwise_or.reduce(x.values))),
                 fam=("fam", "first"), target_seq=("target_seq", "first")).reset_index()


_maps = {}


def m1e_map(db):
    """seq (I/L) -> (fam 0/1/-1, target pair key) for forward peptides; decoy seq -> forward seq (M1e peptide DB)."""
    if db not in _maps:
        p = os.path.join(C.ROOT, "data", "cache", "M1e", "entrapment", "db", f"human_contam.{db}.peptides.tsv")
        t = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        fwd = t[t.kind.isin(["T", "E"])]
        fmap = dict(zip(fwd.seq, zip(np.where(fwd.kind == "T", 0, 1 if db != "C" else -1),
                                     np.where(fwd.pair_key != "", fwd.pair_key, fwd.seq))))
        dec = t[~t.kind.isin(["T", "E"])]
        _maps[db] = (fmap, dict(zip(dec.seq, dec.decoy_of)))
    return _maps[db]


def comparator_table(method, ds, db, runs):
    import analyze as A                                # scripts/baselines/analyze.py (strip = I/L sequence)
    from bench.t1 import comparators as CMP
    rix = {r: i for i, r in enumerate(runs)}
    rows = []
    if method == "sage":
        t = pd.read_csv(os.path.join(C.WORK, "comp", "sage", ds, db, "results.sage.tsv"), sep="\t",
                        usecols=["peptide", "label", "peptide_q", "sage_discriminant_score", "filename", "rank"])
        t = t[t["rank"] == 1]
        for p, lab, q, s, f in zip(t.peptide, t.label, t.peptide_q, t.sage_discriminant_score, t.filename):
            rows.append((A.strip(p), lab == -1, float(q), float(s), f[:-4]))
    else:
        d = os.path.join(PC.FU_WORK, "comp", {"sage+pred+mokapot": "sage_pred_mokapot",
                                              "sage+mokapot(seeded)": "sage_mokapot_seeded"}[method], ds, db)
        pep = {}
        for fn, dec in (("mokapot.peptides.txt", False), ("mokapot.decoy.peptides.txt", True)):
            for x in CMP.read_ragged(os.path.join(d, fn)):
                pep[(x["Peptide"], dec)] = (float(x["mokapot q-value"]), float(x["mokapot score"]))
        for fn, dec in (("mokapot.psms.txt", False), ("mokapot.decoy.psms.txt", True)):
            for x in CMP.read_ragged(os.path.join(d, fn)):
                k = (x["Peptide"], dec)
                if k in pep:
                    rows.append((A.strip(x["Peptide"]), dec, pep[k][0], pep[k][1], x["FileName"][:-4]))
    df = pd.DataFrame(rows, columns=["seq", "is_decoy", "q", "score", "run"])
    df["runbits"] = [1 << rix[r] for r in df.run]
    df = df.groupby(["seq", "is_decoy"], sort=False).agg(
        q=("q", "min"), score=("score", "max"),
        runbits=("runbits", lambda x: int(np.bitwise_or.reduce(x.values)))).reset_index()
    # orientation (label-blind, from the method's own q), as bench.t1.evaluate._orient
    if len(df) > 10 and np.corrcoef(df.q.rank(), df.score.rank())[0, 1] > 0:
        df["score"] = -df["score"]
    if db == "prod":
        df["fam"], df["target_seq"] = 0, df.seq
        return df
    fmap, dmap = m1e_map(db)
    fam, ts = [], []
    for s, d in zip(df.seq, df.is_decoy):
        f = dmap.get(s, s) if d else s
        v = fmap.get(f)
        if v is None:
            raise RuntimeError(f"{method} {ds} {db}: peptide {s} not in the M1e DB")
        fam.append(int(v[0]))
        ts.append(v[1])
    df["fam"], df["target_seq"] = fam, ts
    return df


def ratios():
    from bench.t1 import evaluate as EV
    EV.LEVELS = ("sequence",)
    methods = ["sage+pred+mokapot", "sage+mokapot(seeded)", "sage"]
    out = {"unit": "sequence (I/L-collapsed); run-level IDs summed over runs", "scorer": SCORER,
           "bootstrap": "2000 percentile, seed 11, runs resampled (J = 1)", "datasets": {}}
    for ds in ("D4", "D10cf"):
        runs = sorted(os.path.basename(s)[:-len(".vspec.arrow")] for s in C.DEV[ds]["stores"]())
        R = len(runs)
        curves = {}
        dsl = {"ids_dataset_level_seq": {}}
        for m in ["vulcan"] + methods:
            ids = np.zeros((R, 1, len(EV.GRID)))
            cnt = np.zeros((R, 1, 3, len(EV.GRID), 4))
            for si, k in enumerate(("prod",) + PKINDS):
                tab = ours_table(ds, k, runs) if m == "vulcan" else comparator_table(m, ds, "prod" if k == "prod"
                                                                                     else SAGE_P[k], runs)
                if k == "prod":
                    ids[:, 0] = EV.ids_per_run(tab, R, "sequence")
                    dsl["ids_dataset_level_seq"][m] = {str(a): EV.ids_dataset(tab, "sequence", a)
                                                      for a in (0.001, 0.01, 0.05)}
                else:
                    cnt[:, 0, si - 1] = EV.paired_grid(tab, "sequence", R)
            curves[m] = {"ids": {"sequence": ids}, "cnt": {"sequence": cnt}}
            print(ds, m, dsl["ids_dataset_level_seq"][m], flush=True)
        fake = {"methods": {m: {"verdicts": {"sequence": {}}} for m in methods}}
        ni = EV.noninferiority(fake, curves, ["vulcan"] + methods)
        for m in methods:
            for k, v in ni[m]["sequence"]["own_q"].items():
                for drop in ("comparator_verdict_deciding", "counts_for_noninferiority",
                             "comparator_verdict_conservative"):
                    v.pop(drop, None)
        fdp = {m: {"fdp_P_mean_at_1pct_own_q": float(EV.fdp_paired(curves[m]["cnt"]["sequence"].sum(0)).mean(
            axis=(0, 1))[int(np.flatnonzero(np.isclose(EV.GRID, 0.01))[0])])} for m in curves}
        out["datasets"][ds] = {"runs": runs, **dsl, "ratios": ni, "arm_P_seq_fdp": fdp}
        json.dump(out, open(os.path.join(PC.FU_RES, "ratios.json"), "w"), indent=1, default=float)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["ours", "ratios"])
    ap.add_argument("--datasets", nargs="+", default=["D4", "D10cf"])
    a = ap.parse_args()
    sys.path.insert(0, os.path.join(C.ROOT, "scripts", "baselines"))
    if a.what == "ours":
        ours_eval(a.datasets)
    else:
        ratios()
