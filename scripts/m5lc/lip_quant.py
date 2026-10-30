#!/usr/bin/env python3
"""LiP-MS on D7 (dev): accepted peptidoforms (scripts/m5lc/ids.py, one ids dir per digest group) -> MS1 link ->
RT alignment -> MBR **within the digest group only** (never LiP<->TrP; DESIGN L5) -> fixed-window XIC re-quant ->
half-tryptic calls with biological exclusions + in-source flags (vulcan.lip.halftryptic) -> peptide x run matrix
with cut sites + MSstatsLiP input (LiP / TrP tables in MSstats format, FULL_PEPTIDE = protein_peptide).

    python scripts/m5lc/lip_quant.py --ids-lip DIR --ids-trp DIR --out DIR [--alpha 0.01] [--transfer-q 0.01]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

t_start = time.perf_counter()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.lip import halftryptic as H  # noqa: E402
from vulcan.lip.digest import read_fasta  # noqa: E402
from vulcan.quant import lfq as Q  # noqa: E402


def group_quant(ids_dir, alpha, donor_alpha, tq, T, tag):
    ids = pd.read_parquet(os.path.join(ids_dir, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(ids_dir, "peptides.parquet")).set_index("pf_id")
    acc = ids[ids[f"pass_{alpha}"]].copy()
    runs = sorted(acc["run"].unique())
    t0 = time.perf_counter()
    ms1, tabs, link = {}, {}, {}
    for run in runs:
        r = Q.load_run_ms1(os.path.join(C.WORK, "ms1"), C.store(run), run)
        ms1[run] = r
        L = Q.link_psms(acc[acc["run"] == run], r)
        tabs[run] = Q.precursor_table(L, r)
        link[run] = dict(r.stats)
    T[f"{tag}_link_s"] = time.perf_counter() - t0
    dec_pf = pep.index[pep["is_decoy"]].to_numpy()
    t0 = time.perf_counter()
    al = Q.Alignment({r: t[~t["pf_id"].isin(dec_pf)] for r, t in tabs.items()})
    T[f"{tag}_align_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    dd = ids[ids[f"pass_{donor_alpha}"]]
    dk = set(zip(dd["pf_id"], dd["charge"].astype(np.int8)))
    tr, info = Q.match_between_runs(tabs, ms1, al, dec_pf, donor_keys=dk, params=Q.MBRParams())
    T[f"{tag}_mbr_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    xa = Q.xic_requant(tabs, ms1, {r: C.store(r) for r in runs}, al, transfers=tr, params=Q.XICParams())
    for run in runs:
        tabs[run]["intensity_area"] = tabs[run]["intensity"]
        tabs[run]["intensity"] = xa[run]
    tr["intensity"] = tr["xic"]
    T[f"{tag}_xic_s"] = time.perf_counter() - t0
    allp = pd.concat(tabs.values(), ignore_index=True)
    allp = allp[~allp["pf_id"].isin(dec_pf) & (allp["feature_id"] >= 0)]
    x = tr[(~tr["is_decoy_donor"]) & (tr["transfer_q"] <= tq)].copy()
    cols = ["run", "pf_id", "charge", "intensity", "rt_apex", "provenance", "transfer_q"]
    long = pd.concat([allp[cols], x[cols]], ignore_index=True)
    long = long[np.isfinite(long["intensity"]) & (long["intensity"] > 0)]
    # MS2-only identification counts per run and class
    pid = pep.loc[~pep["is_decoy"]]
    tgt = acc[~acc["is_decoy"]].merge(pid[["cls", "target_seq"]], left_on="pf_id", right_index=True)
    stats = {"runs": runs, "linking": link, "alignment": al.summary(), "mbr": info,
             "n_transfers_accepted": int(len(x)), "transfers_per_run": x.groupby("run").size().to_dict(),
             "peptidoforms_identified": int(tgt["pf_id"].nunique()),
             "peptidoforms_by_class": tgt.drop_duplicates("pf_id")["cls"].value_counts().sort_index().to_dict(),
             "per_run_peptidoforms": tgt.groupby("run")["pf_id"].nunique().to_dict(),
             "per_run_by_class": {r: g.drop_duplicates("pf_id")["cls"].value_counts().sort_index().to_dict()
                                  for r, g in tgt.groupby("run")}}
    return long, pep, stats, acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids-lip", required=True)
    ap.add_argument("--ids-trp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--donor-alpha", type=float, default=0.001)
    ap.add_argument("--transfer-q", type=float, default=0.01)
    ap.add_argument("--features", default=os.path.join(C.WORK, "truth", "yeast_uniprot_features.tsv"))
    ap.add_argument("--insource-tol", type=float, default=0.1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = {}
    res = {}
    longs, peps = {}, {}
    for tag, d in (("LiP", a.ids_lip), ("TrP", a.ids_trp)):
        longs[tag], peps[tag], res[tag], _ = group_quant(d, a.alpha, a.donor_alpha, a.transfer_q, T, tag)
        longs[tag]["digest"] = tag
    # ------------------------------------------------ peptide annotation (half-tryptic calls, cut sites)
    t0 = time.perf_counter()
    acc_, seqs = read_fasta([C.LIP_FASTA])
    feats = H.parse_uniprot_features(a.features)
    db = H.ProteinDB(acc_, seqs, feats)
    allpep = pd.concat([peps["LiP"], peps["TrP"]])
    allpep = allpep[~allpep.index.duplicated()]
    allpep = allpep[~allpep["is_decoy"] & ~allpep["is_entrapment"]]
    used = set(longs["LiP"]["pf_id"]) | set(longs["TrP"]["pf_id"])
    allpep = allpep[allpep.index.isin(used)].copy()
    pidx = [[db.idx[p] for p in s.split(";") if p in db.idx] if s else [] for s in allpep["proteins"].map(
        lambda s: ";".join(x.split("|")[1] if x.count("|") >= 2 else x for x in s.split(";")) if s else "")]
    allpep["prot_idx"] = pidx
    allpep["n_proteins"] = [len(x) for x in pidx]
    ann = H.annotate(allpep.reset_index()[["pf_id", "target_seq", "prot_idx"]], db).set_index("pf_id")
    allpep = allpep.join(ann.drop(columns=["target_seq", "prot_idx"]))
    T["annotate_s"] = time.perf_counter() - t0
    # ------------------------------------------------ long table + in-source flags + TrP presence
    long = pd.concat([longs["LiP"], longs["TrP"]], ignore_index=True)
    long = long.join(allpep[["modseq", "target_seq", "orig_seq", "protein", "n_proteins", "start", "end", "call",
                             "cut_site", "cut_side", "bio_reason", "cls"]], on="pf_id", how="inner")
    long["condition"] = long["run"].map(lambda r: C.D7_DESIGN[r][1])
    long["replicate"] = long["run"].map(lambda r: C.D7_DESIGN[r][2])
    t0 = time.perf_counter()
    long["insource_flag"] = H.insource_flags(long[long["provenance"] == "ms2"].reset_index(drop=True)
                                             .set_index(long.index[long["provenance"] == "ms2"]),
                                             a.insource_tol).reindex(long.index).fillna(False).astype(bool)
    T["insource_s"] = time.perf_counter() - t0
    trp_seqs = set(long.loc[(long["digest"] == "TrP") & (long["provenance"] == "ms2"), "target_seq"])
    long["half_tryptic_seen_in_TrP"] = (long["call"] == "half_tryptic") & long["target_seq"].isin(trp_seqs)
    long.to_parquet(os.path.join(a.out, "lip_long.parquet"))
    # peptide x run matrix with cut sites (LiP and TrP separately; MS2 + transfers with provenance)
    for tag in ("LiP", "TrP"):
        x = long[long["digest"] == tag]
        M = x.pivot_table(index=["protein", "orig_seq", "modseq", "charge", "call", "cut_site", "cut_side", "start",
                                 "end"], columns="run", values="intensity", aggfunc="sum")
        M.reset_index().to_csv(os.path.join(a.out, f"matrix_{tag}.tsv.gz"), sep="\t", index=False)
        Pm = x.pivot_table(index=["protein", "orig_seq", "modseq", "charge"], columns="run", values="provenance",
                           aggfunc="first")
        Pm.reset_index().to_csv(os.path.join(a.out, f"provenance_{tag}.tsv.gz"), sep="\t", index=False)
    # ------------------------------------------------ MSstatsLiP input
    def msstats(x, tag):
        x = x[x["n_proteins"] == 1].copy()       # proteotypic only (protein-unique peptides)
        nmod = x["modseq"].str.count(r"\[")
        df = pd.DataFrame({
            "ProteinName": x["protein"], "PeptideSequence": x["orig_seq"], "PrecursorCharge": x["charge"],
            "FragmentIon": np.where(nmod > 0, "MS1_ox" + nmod.astype(str), "MS1"), "ProductCharge": 0,
            "IsotopeLabelType": "L", "Condition": x["condition"],
            "BioReplicate": x["condition"] + ("_LiP" if tag == "LiP" else "_TrP") + "_" + x["replicate"].astype(str),
            "Run": x["run"], "Fraction": 1, "Intensity": x["intensity"],
            "FULL_PEPTIDE": x["protein"] + "_" + x["orig_seq"]})
        # one feature per (peptide, charge, mod count) and run: positional Met-ox isomers are summed
        keys = [c for c in df.columns if c != "Intensity"]
        return df.groupby(keys, as_index=False)["Intensity"].sum()
    for tag in ("LiP", "TrP"):
        ms = msstats(long[long["digest"] == tag], tag)
        if tag == "TrP":
            ms = ms.drop(columns=["FULL_PEPTIDE"])
        ms.to_csv(os.path.join(a.out, f"msstatslip_{tag}.csv"), index=False)
    # ------------------------------------------------ stats
    out = {"groups": res, "alpha": a.alpha, "donor_alpha": a.donor_alpha, "transfer_q": a.transfer_q,
           "no_LiP_TrP_transfer": "MBR runs separately inside each digest group (acceptors and donors of one group)"}
    pl = long[(long["digest"] == "LiP")]
    ms2 = pl[pl["provenance"] == "ms2"]
    up = allpep[allpep.index.isin(set(ms2["pf_id"]))]
    calls = up["call"].value_counts().to_dict()
    ht = up[up["call"] == "half_tryptic"]
    out["LiP_calls_peptidoforms_ms2"] = calls
    out["LiP_half_tryptic"] = {
        "peptidoforms": int(len(ht)), "cut_sites_unique": int(ht[["protein", "cut_site", "cut_side"]].drop_duplicates().shape[0]),
        "proteins_with_cut": int(ht["protein"].nunique()),
        "by_side": ht["cut_side"].value_counts().to_dict(),
        "search_class_vs_call": pd.crosstab(up["cls"], up["call"]).to_dict(),
        "bio_excluded_reasons": up.loc[up["call"] == "bio_terminal", "bio_reason"].value_counts().to_dict(),
        "frac_semi_peptidoforms_rescued_as_bio": float((up["call"] == "bio_terminal").sum() / max((up["cls"] > 0).sum(), 1)),
        "insource_flagged_peptidoforms": int(ms2.loc[ms2["insource_flag"], "pf_id"].nunique()),
        "insource_flagged_frac": float(ms2.loc[ms2["call"] == "half_tryptic"].groupby("pf_id")["insource_flag"].any().mean()),
        "half_tryptic_also_identified_in_TrP": int(ms2.loc[ms2["half_tryptic_seen_in_TrP"], "pf_id"].nunique()),
    }
    trp = long[(long["digest"] == "TrP") & (long["provenance"] == "ms2")]
    tp = allpep[allpep.index.isin(set(trp["pf_id"]))]
    out["TrP_calls_peptidoforms_ms2"] = tp["call"].value_counts().to_dict()
    out["matrix_rows"] = {tag: int(long[long["digest"] == tag][["pf_id", "charge"]].drop_duplicates().shape[0])
                          for tag in ("LiP", "TrP")}
    T["total_s"] = time.perf_counter() - t_start
    out["timing_s"] = T
    json.dump(out, open(os.path.join(a.out, "lip_quant.json"), "w"), indent=1, default=lambda o: o if not hasattr(o, "item") else o.item())
    print(json.dumps({k: out[k] for k in ("LiP_calls_peptidoforms_ms2", "LiP_half_tryptic", "TrP_calls_peptidoforms_ms2",
                                          "matrix_rows", "timing_s")}, indent=1, default=str))


if __name__ == "__main__":
    main()
