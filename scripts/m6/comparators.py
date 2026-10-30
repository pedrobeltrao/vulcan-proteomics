#!/usr/bin/env python3
"""D9 comparators through the same evaluation as ours (scripts/m6/evaluate_lfq.py):
  sage      : Sage 0.14.6 built-in LFQ (lfq.tsv; peptide level, charge states combined; q_value <= 0.01; Sage LFQ
              integrates every 1%-peptide in every file, i.e. it always borrows identity across runs)
  flashlfq  : FlashLFQ 2.1.5 QuantifiedPeaks.tsv (precursor level; MSMS-only and MSMS+MBR variants)
    python scripts/m6/comparators.py sage DIR OUT.json
    python scripts/m6/comparators.py flashlfq DIR OUT.json
    python scripts/m6/comparators.py flashlfq_truth DIR OUT.json   (9-run FlashLFQ: MBR into human-only runs)
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import evaluate_lfq as EV  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def sage(d, q=0.01):
    t = pd.read_csv(os.path.join(d, "lfq.tsv"), sep="\t")
    t = t[(t["q_value"] <= q)]
    t = t[~t["proteins"].astype(str).str.split(";").map(lambda ps: all(p.startswith("rev_") for p in ps))]
    files = [c for c in t.columns if c.endswith(".mzML")]
    M = t.set_index("peptide")[files]
    M.columns = [c[:-5] for c in files]
    M = M.groupby(level=0).sum(min_count=1)
    prots = t.drop_duplicates("peptide").set_index("peptide")["proteins"]
    return M, prots


def flashlfq(d, mbr):
    t = pd.read_csv(os.path.join(d, "QuantifiedPeaks.tsv"), sep="\t")
    t = t[~t["Decoy Peptide"].astype(str).eq("True") & ~t["Random RT"].astype(str).eq("True")]   # PIP-ECHO decoys
    det = t["Peak Detection Type"].astype(str).str.upper()
    if not mbr:
        t = t[det.str.startswith("MSMS")]
    else:   # MSMS peaks + MBR peaks at PIP-ECHO q <= 1% (FlashLFQ's own transfer FDR)
        t = t[det.str.startswith("MSMS") | (det.str.contains("MBR") & (pd.to_numeric(t["PIP Q-Value"], errors="coerce") <= 0.01))]
    t = t[t["Peak intensity"] > 0]
    t["key"] = t["Full Sequence"].astype(str) + "/" + t["Precursor Charge"].astype(str)
    M = t.pivot_table(index="key", columns="File Name", values="Peak intensity", aggfunc="sum")
    prots = t.drop_duplicates("key").set_index("key")["Protein Group"]
    pfM = M.copy()
    pfM.index = [k.rsplit("/", 1)[0] for k in pfM.index]
    pfM = pfM.groupby(level=0).sum(min_count=1)
    pp = prots.copy()
    pp.index = [k.rsplit("/", 1)[0] for k in pp.index]
    pp = pp[~pp.index.duplicated()]
    return M, prots, pfM, pp, t


def main():
    kind, d, out = sys.argv[1], sys.argv[2], sys.argv[3]
    res = {}
    wd = os.path.join(C.WORK, "eval", f"work_{kind}")
    if kind == "sage":
        M, prots = sage(d)
        res["peptidoform"], _, _ = EV.evaluate(M, prots, wd, "sage")
    else:
        for mbr in (False, True):
            M, prots, pfM, pp, t = flashlfq(d, mbr)
            hye = [c for c in M.columns if C.condition(c) in "AB"]
            tag = "mbr_on" if mbr else "mbr_off"
            res[tag] = {}
            res[tag]["precursor"], _, _ = EV.evaluate(M[hye], prots, wd, f"flash_{tag}")
            res[tag]["peptidoform"], _, _ = EV.evaluate(pfM[hye], pp, wd, f"flash_{tag}_pf", rollups=())
            res[tag]["detection_types"] = t["Peak Detection Type"].value_counts().to_dict()
    json.dump(res, open(out, "w"), indent=1, default=float)
    print(out)




def flashlfq_truth(d, out):
    """FlashLFQ MBR into the 3 human-only runs (ours_all9): yeast/E. coli MBR peaks there are false transfers."""
    t = pd.read_csv(os.path.join(d, "QuantifiedPeaks.tsv"), sep="\t")
    t = t[~t["Decoy Peptide"].astype(str).eq("True") & ~t["Random RT"].astype(str).eq("True")]
    t["sp"] = t["Protein Group"].map(EV.species_from_proteins)
    H = t["File Name"].astype(str).str.contains("_Human_")
    det = t["Peak Detection Type"].astype(str)
    mbr = det.str.upper().str.contains("MBR") & (t["Peak intensity"] > 0)
    res = {"columns": list(t.columns), "detection_types": det.value_counts().to_dict()}
    qcols = ["PIP Q-Value"]
    sets = {"all_mbr": mbr}
    for c in qcols:
        for a in (0.01, 0.05):
            sets[f"{c}<={a}"] = mbr & (pd.to_numeric(t[c], errors="coerce") <= a)
    for k, m in sets.items():
        x = t[H & m]
        n_ye = int(x["sp"].isin(["YEAST", "ECOLI"]).sum())
        n_hu = int((x["sp"] == "HUMAN").sum())
        res[k] = {"mbr_peaks_in_H": int(len(x)), "YE_false": n_ye, "HUMAN": n_hu,
                  "false_fraction_lower_bound": n_ye / max(n_ye + n_hu, 1)}
    json.dump(res, open(out, "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float)[:2000])


if __name__ == "__main__":
    if sys.argv[1] == "flashlfq_truth":
        flashlfq_truth(sys.argv[2], sys.argv[3])
    else:
        main()
