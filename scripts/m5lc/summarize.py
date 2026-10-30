#!/usr/bin/env python3
"""Assemble results/2026-10-27_M5_lip_cf/metrics.json from the M5 LiP/CF work dir.
    python scripts/m5lc/summarize.py"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.hw import fingerprint  # noqa: E402
from vulcan.lip import halftryptic as H  # noqa: E402
from vulcan.lip.digest import read_fasta  # noqa: E402

W = C.WORK


def J(p):
    return json.load(open(p))


def search_stats(d):
    s = J(os.path.join(d, "stats.json"))
    warm = s["passes"][-1]
    return {"setup": s["setup"], "gpu_max_mem_bytes": s["gpu_max_mem_bytes"], "params": s["params"],
            "n_runs": len(warm),
            "warm_per_run_wall_s_median": float(np.median([r["wall_s_total"] for r in warm])),
            "warm_per_run_search_s_median": float(np.median([r["search_wall_s"]["total"] for r in warm])),
            "warm_per_run_core_s_median": float(np.median([r["cpu_core_s"] for r in warm])),
            "warm_total_wall_s": float(sum(r["wall_s_total"] for r in warm)),
            "n_spectra_median": float(np.median([r["counts"]["n_spectra"] for r in warm])),
            "process_wall_s": s["process_wall_s"], "hw": s["hw"]}


def sage_summary(group, db):
    d = os.path.join(W, "sage", f"semi_{group}_t32")
    t = pd.read_csv(os.path.join(d, "results.sage.tsv"), sep="\t")
    tg = t[t["label"] == 1]
    psm = tg[tg["spectrum_q"] <= 0.01]
    pep = tg[tg["peptide_q"] <= 0.01].drop_duplicates("peptide")
    stripped = pep["peptide"].str.replace(r"\[[^\]]*\]", "", regex=True).str.replace("I", "L")
    calls = []
    for s, prots in zip(stripped, pep["proteins"]):
        acc = [p.split("|")[1] if p.count("|") >= 2 else p for p in str(prots).split(";")]
        idx = [db.idx[a] for a in acc if a in db.idx]
        calls.append(H.call(s, idx, db)["call"])
    calls = pd.Series(calls).value_counts().to_dict()
    tm = [json.loads(x) for x in open(os.path.join(C.RES, "timings.jsonl")) if f"M5_sage_semi_D7_{group}_t32" in x]
    return {"psms_spectrum_q01": int(len(psm)), "peptides_peptide_q01": int(len(pep)),
            "peptides_semi_enzymatic_flag": int(pep["semi_enzymatic"].astype(bool).sum()) if "semi_enzymatic" in pep else None,
            "peptides_by_our_call": calls, "peptide_set": set(stripped),
            "timing": {k: tm[-1][k] for k in ("wall_s", "core_s", "max_rss_mb", "lock_wait_s")} if tm else None}


def main():
    out = {"hw": fingerprint(), "dev_only": True, "sealed_touched": False}
    # ------------------------------------------------------------------ LiP
    meta = J(os.path.join(C.LIP_INDEX, "meta.json"))
    lip = {"index": {k: meta[k] for k in ("n_targets", "n_peptidoforms", "n_peptidoforms_by_family", "timing_s",
                                         "semi", "n_forward_without_decoy", "n_decoy_shuffled", "content_hash")},
           "search": search_stats(os.path.join(W, "psms", "d7_semi"))}
    si = J(os.path.join(W, "psms", "d7_semi", "score_info.json"))
    lip["scoring"] = {"model": si["model"], "per_run_features_s_median": float(np.median(
        [r["timing_s"]["features_total"] for r in si["runs"].values()])), "python_import_s": si["python_import_s"],
        "model_load_s": si["model_load_s"]}
    for g in ("lip", "trp"):
        lip[f"ids_{g}"] = J(os.path.join(W, "ids", g, "ids_summary.json"))["groups"]["all"]
    lq = J(os.path.join(W, "lip", "lip_quant.json"))
    ic = os.path.join(W, "lip", "ids_calls.json")      # half-tryptic calls on ALL identified peptidoforms (1%)
    if os.path.exists(ic):
        lip["calls_all_identified_1pct"] = J(ic)
    lip["quant"] = lq
    acc, seqs = read_fasta([C.LIP_FASTA])
    db = H.ProteinDB(acc, seqs, H.parse_uniprot_features(os.path.join(W, "truth", "yeast_uniprot_features.tsv")))
    ours_lip = pd.read_parquet(os.path.join(W, "lip", "lip_long.parquet"))
    sage = {}
    for g in ("LiP", "TrP"):
        s = sage_summary(g, db)
        ourset = set(ours_lip.loc[(ours_lip["digest"] == g) & (ours_lip["provenance"] == "ms2"), "target_seq"])
        s["overlap_stripped_with_ours_quantified"] = len(s.pop("peptide_set") & ourset)
        s["ours_quantified_stripped"] = len(ourset)
        sage[g] = s
    lip["sage_semi"] = sage
    rj = os.path.join(W, "lip", "msstatslip_run.json")
    if os.path.exists(rj):
        lip["msstatslip"] = J(rj)
        m = pd.read_csv(os.path.join(W, "lip", "msstatslip_Adjusted_LiP_model.csv"))
        u = pd.read_csv(os.path.join(W, "lip", "msstatslip_LiP_model.csv"))
        t = pd.read_csv(os.path.join(W, "lip", "msstatslip_TrP_model.csv"))
        lip["msstatslip"]["n_peptides_tested"] = int(len(m))
        lip["msstatslip"]["adjusted_sig_adjp05"] = int((m["adj.pvalue"] < 0.05).sum())
        lip["msstatslip"]["unadjusted_sig_adjp05"] = int((u["adj.pvalue"] < 0.05).sum())
        lip["msstatslip"]["trp_proteins_tested"] = int(len(t))
        cols = ["PeptideSequence", "log2FC", "pvalue", "adj.pvalue", "fully_TRI", "NSEMI_TRI", "CSEMI_TRI",
                "StartPos", "EndPos"]
        sn = m[m["ProteinName"].astype(str) == C.SNCA].sort_values("StartPos")
        snu = u[u["ProteinName"].astype(str) == C.SNCA].sort_values("StartPos")
        lip["snca"] = {"adjusted": sn[cols].to_dict("records"), "unadjusted": snu[cols].to_dict("records"),
                       "trp": t[t["Protein"].astype(str) == C.SNCA].drop(columns=["issue"], errors="ignore").to_dict("records")}
        # expected structural change: fibril core (NAC region, residues ~37-99) protected in F1 vs M1
        core = (sn["StartPos"] < 99) & (sn["EndPos"] > 37)
        sig = sn["adj.pvalue"] < 0.05
        lip["snca"]["summary"] = {
            "peptides_tested": int(len(sn)), "in_core_37_99": int(core.sum()), "sig_in_core": int((core & sig).sum()),
            "sig_outside_core": int((~core & sig).sum()), "sig_in_core_down": int((core & sig & (sn["log2FC"] < 0)).sum()),
            "sig_in_core_up": int((core & sig & (sn["log2FC"] > 0)).sum())}
    out["lip"] = lip
    # ------------------------------------------------------------------ CF
    cf = {"search": search_stats(os.path.join(W, "psms", "cf")),
          "ids": {r: g["0.01"]["peptidoforms"] for r, g in J(os.path.join(W, "ids", "cf", "ids_summary.json"))["groups"].items()},
          "e10": J(os.path.join(W, "cf", "e10.json"))}
    cmeta = J(os.path.join(C.CF_INDEX, "meta.json"))
    cf["index"] = {k: cmeta[k] for k in ("n_targets", "n_peptidoforms", "n_peptidoforms_by_family", "timing_s")}
    out["cf"] = cf
    # ------------------------------------------------------------------ ingest / MS1 per run
    ing, m1 = [], []
    for f in glob.glob(os.path.join(W, "logs", "convert_*.json")):
        r = J(f)
        tv = r.get("time_v") or {}
        ing.append({"run": r["run"], "wall_s": r["wall_s"], "copy_s": r["copy_s"], "core_s": tv.get("core_s")})
    for f in glob.glob(os.path.join(W, "logs", "ms1_*.time.json")):
        r = J(f)
        tv = r.get("time_v") or {}
        m1.append({"run": r["run"], "wall_s": tv.get("wall_s_timev"), "core_s": tv.get("core_s")})
    out["ingest_runs"] = ing
    out["ms1_runs"] = m1
    out["ingest_median"] = {k: float(np.nanmedian([x[k] or np.nan for x in ing])) for k in ("wall_s", "core_s", "copy_s")} if ing else None
    out["ms1_median"] = {k: float(np.nanmedian([x[k] or np.nan for x in m1])) for k in ("wall_s", "core_s")} if m1 else None
    os.makedirs(C.RES, exist_ok=True)
    json.dump(out, open(os.path.join(C.RES, "metrics.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps({"lip_sage": lip["sage_semi"], "snca": lip.get("snca", {}).get("summary"),
                      "ingest": out["ingest_median"], "ms1": out["ms1_median"]}, indent=1, default=str))


if __name__ == "__main__":
    main()
