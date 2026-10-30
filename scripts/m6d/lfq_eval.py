#!/usr/bin/env python3
"""M6d: full LFQ (link -> align -> opt-in MBR -> XIC re-quant [-> charge siblings]) on a dev set / run subset, version
m6 (M6 settings) or m6d (fix), scored with the unchanged §6.8 metrics: precursor epsilon / CV / completeness,
peptidoform epsilon (charge states summed), protein Top3 / MaxLFQ / directLFQ with DE pAUC, and - on D9 - epsilon
on the peptidoforms matched with Sage LFQ.  MBR uses the M6 defaults unchanged (MBRParams(), donors at ID q <= 0.1%).
    scripts/m2b/vpy.sh scripts/m6d/lfq_eval.py DATASET SUBSET VERSION OUT.json [--mbr] [--no-directlfq]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m6"))
sys.path.append(os.path.join(ROOT, "data", "cache", "pylibs_m6"))
import compare_common as CC  # noqa: E402
import datasets as DS  # noqa: E402
import evalq as E  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.quant import lfq as Q  # noqa: E402
from vulcan.quant import metrics as MET  # noqa: E402


def species_of_prot(p):
    import evaluate_lfq as EV
    return EV.species_from_proteins(p)


def protein_eval(Mn, M_raw, pk, A, B, expected, workdir, directlfq=True):
    out, T = {}, {}
    psp = pd.Series({p: species_of_prot(p) for p in set(pk)})
    P = {}
    t0 = time.perf_counter(); P["top3"] = Q.topn(Mn, pk, 3); T["top3_s"] = time.perf_counter() - t0
    t0 = time.perf_counter(); P["maxlfq"] = Q.maxlfq(Mn, pk); T["maxlfq_s"] = time.perf_counter() - t0
    if directlfq:
        t0 = time.perf_counter()
        try:
            P["directlfq"] = Q.directlfq_rollup(M_raw, pk, workdir)
        except Exception as e:   # pragma: no cover
            print("directLFQ failed", repr(e))
        T["directlfq_s"] = time.perf_counter() - t0
    for k, X in P.items():
        out[k] = {"n": int(len(X)), "eps_min2": MET.epsilon(X, psp, A, B, expected, min(2, len(A), len(B)))["median_abs_eps"],
                  "cv_min3": MET.cv_stats(X, A, B, 3)["median_cv"] if min(len(A), len(B)) >= 3 else None,
                  "de": MET.de_pauc(X, psp, A, B) if min(len(A), len(B)) >= 2 else None}
    return out, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("subset")
    ap.add_argument("version", choices=["m6", "m6d"])
    ap.add_argument("out")
    ap.add_argument("--mbr", action="store_true")
    ap.add_argument("--no-directlfq", action="store_true")
    ap.add_argument("--tq", type=float, default=0.01)
    ap.add_argument("--save-long", default=None, help="write the long (run, precursor, intensity, provenance) table")
    a = ap.parse_args()
    cfg = DS.get(a.dataset)
    runs, A, B = cfg.subset(a.subset)
    T = {}
    t0 = time.perf_counter()
    ids = pd.read_parquet(os.path.join(cfg.ids, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(cfg.ids, "peptides.parquet")).set_index("pf_id")
    acc = ids[ids["pass_0.01"] & ids["run"].isin(runs)]
    ms1, tabs = {}, {}
    for run in runs:
        r = Q.load_run_ms1(cfg.ms1, cfg.store(run), run)
        ms1[run] = r
        tabs[run] = Q.precursor_table(Q.link_psms(acc[acc["run"] == run], r), r)
    T["link_s"] = time.perf_counter() - t0
    dec = pep.index[pep["is_decoy"]].to_numpy()
    t0 = time.perf_counter()
    al = Q.Alignment({r: t[~t["pf_id"].isin(dec)] for r, t in tabs.items()})
    T["align_s"] = time.perf_counter() - t0
    tr = None
    if a.mbr:
        t0 = time.perf_counter()
        dd = ids[ids["pass_0.001"] & ids["run"].isin(runs)]
        dk = set(zip(dd["pf_id"], dd["charge"].astype(np.int8)))
        tr, _ = Q.match_between_runs(tabs, ms1, al, dec, donor_keys=dk, params=Q.MBRParams())
        T["mbr_s"] = time.perf_counter() - t0
    stores = {r: cfg.store(r) for r in runs}
    xp = Q.XICParams() if a.version == "m6d" else None
    t0 = time.perf_counter()
    xa = Q.xic_requant(tabs, ms1, stores, al, transfers=tr, params=xp)
    T["xic_s"] = time.perf_counter() - t0
    parts = []
    for run in runs:
        t = tabs[run]
        parts.append(pd.DataFrame({"run": run, "pf_id": t["pf_id"].to_numpy(), "charge": t["charge"].to_numpy(),
                                   "intensity": np.where(t["feature_id"].to_numpy() >= 0, xa[run], np.nan),
                                   "provenance": "ms2"}))
    n_sib = 0
    if a.version == "m6d" and xp.siblings:
        t0 = time.perf_counter()
        for run in runs:
            tabs[run]["intensity"] = xa[run]          # anchor = most intense quantified charge (XIC value)
        sib = Q.charge_siblings(tabs, ms1, stores, xp, exclude_pf=dec)
        T["siblings_s"] = time.perf_counter() - t0
        n_sib = int(len(sib))
        parts.append(sib[["run", "pf_id", "charge", "intensity", "provenance"]])
    if tr is not None:
        x = tr[(~tr["is_decoy_donor"]) & (tr["transfer_q"] <= a.tq)]
        parts.append(pd.DataFrame({"run": x["run"], "pf_id": x["pf_id"], "charge": x["charge"], "intensity": x["xic"],
                                   "provenance": "mbr"}))
    L = pd.concat(parts, ignore_index=True)
    L = L[~L["pf_id"].isin(dec) & np.isfinite(L["intensity"]) & (L["intensity"] > 0)]
    # an MS2 or sibling value wins over a transfer for the same (run, precursor)
    L["_p"] = L["provenance"].map({"ms2": 0, "ms2_sibling": 1, "mbr": 2})
    L = L.sort_values("_p").drop_duplicates(["run", "pf_id", "charge"])
    L["pf"] = pep["modseq"].reindex(L["pf_id"].to_numpy()).map(CC.norm).to_numpy()
    L["key"] = L["pf"] + "/" + L["charge"].astype(str)
    L["prot"] = pep["proteins"].reindex(L["pf_id"].to_numpy()).to_numpy()
    if a.save_long:
        L.drop(columns=["_p"]).to_parquet(a.save_long)
    res = {"dataset": a.dataset, "subset": a.subset, "runs": runs, "version": a.version, "mbr": a.mbr,
           "n_sibling_values": n_sib, "provenance_counts": L["provenance"].value_counts().to_dict()}
    Mp, Pp = E.matrix(L, "precursor")
    op, Mn = E.score(Mp, Pp, A, B, cfg.expected)
    res["precursor"] = op
    Mf, Pf = E.matrix(L, "peptidoform", "sum")
    res["peptidoform"], _ = E.score(Mf, Pf, A, B, cfg.expected)
    # protein level from the precursor matrix (same post-processing as M6 evaluate_lfq)
    Mraw = Mp.reindex(columns=A + B).replace(0, np.nan)
    Mraw = Mraw[Mraw.notna().any(axis=1)]
    sp = E.species(Pp.reindex(Mraw.index))
    Mraw = Mraw[sp.isin(list(cfg.expected)).to_numpy()]
    Mn2, _ = MET.median_normalize(Mraw)
    pk = Pp.reindex(Mraw.index).map(lambda s: ";".join(sorted(p for p in str(s).split(";") if p)))
    t0 = time.perf_counter()
    res["protein"], Tp = protein_eval(Mn2, Mraw, pk, A, B, cfg.expected,
                                      f"/tmp/vulcan_m6d/dlfq_{a.dataset}_{a.subset}_{a.version}_{int(a.mbr)}",
                                      directlfq=not a.no_directlfq and min(len(A), len(B)) >= 2)
    T.update(Tp)
    if a.dataset == "d9":
        import comparators as CMP
        Ms, Ps = CMP.sage("/tmp/vulcan_m6/sage/D9")
        Ms.index = [CC.norm(i) for i in Ms.index]
        Ms = Ms.groupby(level=0).sum(min_count=1).reindex(columns=A + B).replace(0, np.nan)
        Ps.index = [CC.norm(i) for i in Ps.index]
        Ps = Ps[~Ps.index.duplicated()]
        ge1 = lambda M: set(M.index[(M[A].notna().sum(axis=1) >= 1) & (M[B].notna().sum(axis=1) >= 1)])
        Mf2 = Mf.reindex(columns=A + B)
        keep = sorted(ge1(Mf2) & ge1(Ms))
        fixed = json.load(open("/tmp/vulcan_m6d/d9_keep_m6.json")) if os.path.exists("/tmp/vulcan_m6d/d9_keep_m6.json") else None
        if a.version == "m6" and not a.mbr and a.subset == "all":
            json.dump(keep, open("/tmp/vulcan_m6d/d9_keep_m6.json", "w"))
            fixed = keep
        mt = {}
        for name, rows in (("matched_this_version", keep), ("matched_m6_set", fixed)):
            if rows is None:
                continue
            so, _ = E.score(Ms, Ps, A, B, cfg.expected, rows=rows, norm_rows=rows)
            oo, _ = E.score(Mf2, Pf, A, B, cfg.expected, rows=rows, norm_rows=rows)
            mt[name] = {"n": oo["n"], "ours": oo, "sage": so}
        res["vs_sage"] = mt
    T["total_s"] = sum(v for k, v in T.items() if k.endswith("_s"))
    res["timing_s"] = T
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    short = {"prec_eps": op["eps"], "prec_cv": op.get("cv", op.get("cv2")), "compl": op["completeness"], "n": op["n"],
             "pf_eps": res["peptidoform"]["eps"]}
    if "vs_sage" in res:
        short.update({k: (v["ours"]["eps"], v["sage"]["eps"]) for k, v in res["vs_sage"].items()})
    short["pauc_dlfq"] = (res["protein"].get("directlfq") or {}).get("de", {}) and res["protein"]["directlfq"]["de"]["pauc_fpr5pct_std"]
    print(json.dumps(short, default=float))


if __name__ == "__main__":
    main()
