#!/usr/bin/env python3
"""D4 agreement checks: (1) our reporter intensities vs Sage 0.14.6 TMT11 level-3 reporters on the same MS2 scans;
(2) our accepted PSMs vs Sage's and vs the submitter's SEQUEST/LDA identifications (mzIdentML passThreshold=true,
data/truth/PXD016766, fetched by scripts/data/fetch_tier_s.py --manifest d4_ref).  The submitter files carry no
reporter intensities, so quant agreement is only checked against Sage.
    python scripts/m6/tmt_agree.py
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from vulcan.quant.tmt import CHANNELS  # noqa: E402

D4 = os.path.join(C.WORK, "d4")


def strip(s):
    return re.sub(r"\[.*?\]", "", str(s)).replace("-", "").replace("I", "L")   # I/L isobaric


def parse_mzid(path):
    txt = open(path).read()
    peps = dict(re.findall(r'<Peptide id="([^"]+)">\s*<PeptideSequence>([A-Z]+)</PeptideSequence>', txt))
    out = []
    for m in re.finditer(r'<SpectrumIdentificationResult [^>]*spectrumID="index=\d+ scan=(\d+)"[^>]*>(.*?)'
                         r'</SpectrumIdentificationResult>', txt, re.S):
        scan, body = int(m.group(1)), m.group(2)
        s = re.search(r'<SpectrumIdentificationItem [^>]*rank="1"[^>]*peptide_ref="([^"]+)"[^>]*passThreshold="(\w+)"', body)
        if s and s.group(2) == "true":
            out.append((scan, peps.get(s.group(1), "")))
    return pd.DataFrame(out, columns=["scan", "seq"])


def main():
    Q = pd.read_parquet(os.path.join(D4, "tmt_psms.parquet"))
    Q["seq"] = Q["modseq"].map(strip)
    res = {}
    # ---- Sage reporters
    S = pd.read_csv(os.path.join(C.WORK, "sage", "D4", "tmt.tsv"), sep="\t")
    S["run"] = S["filename"].str.replace(".mzML", "", regex=False)
    S["scan"] = S["scannr"].str.extract(r"scan=(\d+)")[0].astype(int)
    m = Q.merge(S, on=["run", "scan"], how="inner")
    ours = m[[f"raw_{c}" for c in CHANNELS]].to_numpy()
    sage = m[[f"tmt_{i}" for i in range(1, 12)]].to_numpy()
    ok = (ours > 0) & (sage > 0)
    lo, ls = np.log2(ours[ok]), np.log2(sage[ok])
    rel = np.abs(ours[ok] / sage[ok] - 1)
    sh_o = ours / np.where(ours.sum(1, keepdims=True) > 0, ours.sum(1, keepdims=True), np.nan)
    sh_s = sage / np.where(sage.sum(1, keepdims=True) > 0, sage.sum(1, keepdims=True), np.nan)
    both = np.isfinite(sh_o).all(1) & np.isfinite(sh_s).all(1)
    res["reporters_vs_sage"] = {
        "psms_matched_by_scan": int(len(m)), "channel_values_both_nonzero": int(ok.sum()),
        "pearson_log2": float(stats.pearsonr(lo, ls)[0]), "spearman": float(stats.spearmanr(lo, ls)[0]),
        "frac_within_1pct": float((rel < 0.01).mean()), "frac_within_5pct": float((rel < 0.05).mean()),
        "zero_in_ours_nonzero_in_sage": int(((ours == 0) & (sage > 0)).sum()),
        "zero_in_sage_nonzero_in_ours": int(((sage == 0) & (ours > 0)).sum()),
        "channel_share_pearson": float(stats.pearsonr(sh_o[both].ravel(), sh_s[both].ravel())[0])}
    # ---- IDs vs Sage
    T = pd.read_csv(os.path.join(C.WORK, "sage", "D4", "results.sage.tsv"), sep="\t",
                    usecols=["peptide", "label", "peptide_q", "spectrum_q", "filename", "scannr"])
    T = T[(T["label"] == 1) & (T["peptide_q"] <= 0.01) & (T["spectrum_q"] <= 0.01)].copy()
    T["run"] = T["filename"].str.replace(".mzML", "", regex=False)
    T["scan"] = T["scannr"].str.extract(r"scan=(\d+)")[0].astype(int)
    T["seq"] = T["peptide"].map(strip)
    mm = Q[["run", "scan", "seq"]].merge(T[["run", "scan", "seq"]], on=["run", "scan"], how="outer",
                                         suffixes=("_o", "_s"), indicator=True)
    both = mm["_merge"] == "both"
    res["ids_vs_sage"] = {"ours_psms": int(len(Q)), "sage_psms": int(len(T)), "same_scan": int(both.sum()),
                          "same_scan_same_seq": int((both & (mm["seq_o"] == mm["seq_s"])).sum()),
                          "ours_only": int((mm["_merge"] == "left_only").sum()),
                          "sage_only": int((mm["_merge"] == "right_only").sum()),
                          "sage_ids_note": "Sage LDA-rescored q; ours = plain hyperscore TDC (no IT-CID model)"}
    # ---- IDs vs submitter
    sub = []
    for f in sorted(glob.glob(os.path.join(C.ROOT, "data", "truth", "PXD016766", "*.mzIdentML"))):
        run = re.sub(r"^\d+_", "", os.path.basename(f))[:-len(".mzIdentML")]
        d = parse_mzid(f)
        d["run"] = run
        sub.append(d)
    U = pd.concat(sub, ignore_index=True)
    U["seq"] = U["seq"].str.replace("I", "L")
    Qs = Q[["run", "scan", "seq"]].assign(seq=lambda x: x["seq"].str.replace("I", "L"))
    mu = Qs.merge(U, on=["run", "scan"], how="outer", suffixes=("_o", "_u"), indicator=True)
    b = mu["_merge"] == "both"
    res["ids_vs_submitter"] = {"submitter_psms_passThreshold": int(len(U)), "ours_psms": int(len(Q)),
                               "same_scan": int(b.sum()), "same_scan_same_seq": int((b & (mu["seq_o"] == mu["seq_u"])).sum()),
                               "ours_only": int((mu["_merge"] == "left_only").sum()),
                               "submitter_only": int((mu["_merge"] == "right_only").sum()),
                               "submitter_peptides": int(U["seq"].nunique()), "ours_peptides": int(Qs["seq"].nunique())}
    json.dump(res, open(os.path.join(D4, "tmt_agreement.json"), "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
