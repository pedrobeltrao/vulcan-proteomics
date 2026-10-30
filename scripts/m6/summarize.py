#!/usr/bin/env python3
"""Collect M6 results into results/2026-10-27_M6_quant/metrics.json (+ hwinfo.json)."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402

W = C.WORK
J = lambda p: json.load(open(p))  # noqa: E731


def slim(rep):
    """Headline numbers of an evaluate_lfq report."""
    if rep is None:
        return None
    p = rep["precursor"]
    out = {"n_any": p["completeness"]["n_any"], "n_ge1_per_condition": p["completeness"]["n_ge1_per_condition"],
           "n_3_per_condition": p["completeness"]["n_3_per_condition"],
           "completeness": p["completeness"]["completeness"],
           "median_abs_eps_min1": p["epsilon_min1"]["median_abs_eps"], "mean_abs_eps_min1": p["epsilon_min1"]["mean_abs_eps"],
           "median_abs_eps_min3": p["epsilon_min3"]["median_abs_eps"],
           "median_abs_eps_min1_unnormalized": p["epsilon_min1_unnormalized"]["median_abs_eps"],
           "eps_per_species_min1": {s: {k: v[k] for k in ("n", "median_abs_eps", "median_log2fc")}
                                    for s, v in p["epsilon_min1"]["per_species"].items()},
           "median_cv_min3": p["cv_min3"]["median_cv"], "n_cv": p["cv_min3"]["A"]["n"] + p["cv_min3"]["B"]["n"]}
    if rep.get("protein"):
        out["protein"] = {k: {"n_any": v["completeness"]["n_any"], "completeness": v["completeness"]["completeness"],
                              "median_cv_min3": v["cv_min3"]["median_cv"],
                              "median_abs_eps_min2": v["epsilon_min2"]["median_abs_eps"],
                              "de": v["de"]} for k, v in rep["protein"].items()}
    return out


def species_space(index_dir):
    """Species-unique target peptide sequences in the index (for the H-run natural-entrapment ratio r)."""
    with ipc.open_file(os.path.join(index_dir, "proteins.arrow")) as rd:
        t = rd.read_all()
    pep = t.column("peptide").to_pylist()
    off = t.column("prot_off").to_numpy()
    ln = t.column("prot_len").to_numpy()
    pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
    acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
    sp_of = np.array([a.rsplit("_", 1)[-1] if "Cont_" not in a else "CONT" for a in acc])
    cnt = {}
    for i in range(len(pep)):
        s = set(sp_of[pidx[off[i]:off[i] + ln[i]]])
        k = s.pop() if len(s) == 1 else "MULTI"
        cnt[k] = cnt.get(k, 0) + 1
    return cnt


def h_run_ids(tag, r):
    d = os.path.join(W, "ids", tag)
    ids = pd.read_parquet(os.path.join(d, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(d, "peptides.parquet")).set_index("pf_id")
    out = {}
    for al in ("0.001", "0.01", "0.05"):
        a = ids[ids["pass_" + al] & ~ids["is_decoy"] & ids["run"].str.contains("_Human_")].drop_duplicates("pf_id")
        sp = pep["species"].reindex(a["pf_id"]).to_numpy()
        ct = pep["contaminant"].reindex(a["pf_id"]).to_numpy().astype(bool)
        ye = int((np.isin(sp, ["YEAST", "ECOLI"]) & ~ct).sum())
        n = int(len(a))
        out[al] = {"peptidoforms_in_H_runs": n, "YE_unique": ye,
                   "combined_est_FDP (N_YE (1+1/r)/N)": ye * (1 + 1 / r) / max(n, 1)}
    return out


def main():
    M = {"experiment": "M6 quant (dev only: D9, D4)", "date": "2026-10-27"}
    M["ids"] = {"D9_hye_static_gbm_pred": J(os.path.join(W, "ids", "d9_hye", "ids_summary.json")),
                "D9_hye_hyperscore": J(os.path.join(W, "ids", "d9_hye_hs", "ids_summary.json")),
                "D9_all9_static_gbm_pred": J(os.path.join(W, "ids", "d9_all", "ids_summary.json")),
                "D4_hyperscore": J(os.path.join(W, "d4", "ids", "ids_summary.json"))}
    t = pd.read_csv(os.path.join(W, "sage", "D9", "results.sage.tsv"), sep="\t",
                    usecols=["peptide", "charge", "label", "peptide_q", "spectrum_q"])
    a = t[(t["label"] == 1) & (t["peptide_q"] <= 0.01)]
    M["ids"]["D9_hye_sage_own_q"] = {"peptides_q1pct": int(a["peptide"].nunique()),
                                     "precursors_q1pct": int(a[["peptide", "charge"]].drop_duplicates().shape[0])}
    sp = species_space(C.D9_INDEX)
    r = (sp.get("YEAST", 0) + sp.get("ECOLI", 0)) / sp["HUMAN"]
    M["ids"]["D9_H_run_natural_entrapment"] = {
        "note": "dev sanity check, not a §6.3 verdict: yeast/E. coli-unique peptidoforms accepted in the 3 human-only "
                "runs (q pooled over 9 runs); r = YE/HUMAN species-unique target peptide sequences in the index",
        "species_unique_peptides": sp, "r": r,
        "static_gbm_pred": h_run_ids("d9_all", r), "hyperscore": h_run_ids("d9_all_hs", r)}
    # ---------------- LFQ
    L = {}
    ours = J(os.path.join(W, "lfq", "d9_hye_qtot_d01", "lfq_metrics.json"))
    L["ours_mbr_off"] = slim(ours["hye"]["mbr_off"])
    L["ours_mbr_off_peptidoform"] = slim(ours["hye"]["mbr_off_peptidoform"])
    L["ours_linking"] = ours["linking"]
    L["ours_alignment"] = ours["alignment"]
    L["ours_ids"] = ours["ids"]
    for v, lab in (("qtot_d01", "donors q<=0.1%, transfer_q (total)"), ("qtot_d1", "donors q<=1%, transfer_q (total)"),
                   ("qpeak_d1", "donors q<=1%, transfer_q_peak (peak decoys only)")):
        h = J(os.path.join(W, "lfq", f"d9_hye_{v}", "lfq_metrics.json"))
        for k, rr in h["hye"].items():
            if k.startswith("mbr_on"):
                L[f"ours_mbr_{v}_{k[7:]}"] = {"setting": lab, "n_transfers": rr["n_transfers_accepted"],
                                             "precursor": slim(rr["metrics"]),
                                             "peptidoform": slim(rr["metrics_peptidoform"])}
        L[f"mbr_competition_{v}"] = h.get("mbr_competition")
    ar = J(os.path.join(W, "lfq", "d9_hye_area_qtot_d01", "lfq_metrics.json"))
    L["ours_arm_feature_area_mbr_off"] = slim(ar["hye"]["mbr_off"])
    L["ours_arm_feature_area_mbr_on_q0.01"] = slim(ar["hye"]["mbr_on_q0.01"]["metrics"])
    if os.path.exists(os.path.join(C.RES, "eval_common.json")):
        L["matched_analytes"] = J(os.path.join(C.RES, "eval_common.json"))
    sg = J(os.path.join(C.RES, "eval_sage_d9.json")) if os.path.exists(os.path.join(C.RES, "eval_sage_d9.json")) else None
    L["sage_lfq"] = slim(sg["peptidoform"]) if sg else None
    for k in ("ours_hye", "sage_hye"):
        p = os.path.join(C.RES, f"eval_flashlfq_{k}.json")
        if os.path.exists(p):
            f = J(p)
            for mb in ("mbr_off", "mbr_on"):
                L[f"flashlfq_{k}_{mb}"] = {"precursor": slim(f[mb]["precursor"]), "peptidoform": slim(f[mb]["peptidoform"]),
                                           "detection_types": f[mb]["detection_types"]}
    M["lfq_d9"] = L
    # ---------------- MBR truth (human-only acceptor runs)
    T = {}
    for v in ("qtot_d01", "qtot_d1", "qpeak_d1"):
        d = J(os.path.join(W, "lfq", f"d9_all_{v}", "lfq_metrics.json"))
        T[v] = {"truth": d["mbr_truth_human_only"], "competition": d["mbr_competition"], "mbr_info": d["mbr"]}
    for v, sub in (("v1_allow_claimed", "d9_all_qtot_d01"), ("v1_allow_claimed", "d9_all_qtot_d1"),
                   ("v2_claimed_pooledq", "d9_all_qtot_d01")):
        d = J(os.path.join(W, "lfq", v, sub, "lfq_metrics.json"))
        T[f"{v}_{sub[7:]}"] = {"truth": d["mbr_truth_human_only"], "competition": d["mbr_competition"]}
    T["default"] = "qtot_d01 (donors at ID q <= 0.1%, claimed features excluded, per-acceptor-run transfer q, threshold 1%)"
    p = os.path.join(C.RES, "eval_flashlfq_all9_truth.json")
    if os.path.exists(p):
        T["flashlfq_ours_all9"] = J(p)
    M["mbr_truth_d9"] = T
    # ---------------- TMT
    M["tmt_d4"] = {"metrics": J(os.path.join(W, "d4", "tmt_metrics.json")),
                   "agreement": J(os.path.join(W, "d4", "tmt_agreement.json")),
                   "search": J(os.path.join(W, "d4", "psms", "stats.json"))}
    # ---------------- timings
    tim = [json.loads(x) for x in open(os.path.join(C.RES, "timings.jsonl"))]
    M["timings"] = {x["label"]: {k: x.get(k) for k in ("wall_s", "core_s", "max_rss_mb", "threads_env", "lock_wait_s",
                                                       "loadavg_before", "rc")} for x in tim}
    try:
        from vulcan.hw import fingerprint
        hw = fingerprint()
        json.dump(hw, open(os.path.join(C.RES, "hwinfo.json"), "w"), indent=1, default=str)
        M["hw"] = {k: hw.get(k) for k in ("hostname", "cpu_model", "n_affinity", "gpu", "gpus") if k in hw}
    except Exception as e:  # pragma: no cover
        M["hw"] = {"error": repr(e)}
    json.dump(M, open(os.path.join(C.RES, "metrics.json"), "w"), indent=1, default=float)
    print("wrote metrics.json")


if __name__ == "__main__":
    main()
