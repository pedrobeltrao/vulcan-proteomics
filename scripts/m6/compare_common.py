#!/usr/bin/env python3
"""Accuracy at matched analytes: epsilon / CV of ours (MBR off and on) vs Sage LFQ vs FlashLFQ on the peptidoforms
that every tool quantifies (I/L and the static CAM mark normalised away).  ProteoBench's own plot trades depth
against epsilon; this removes the depth difference.
    python scripts/m6/compare_common.py OUT.json
"""
from __future__ import annotations

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import comparators as CMP  # noqa: E402
import evaluate_lfq as EV  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.quant import metrics as MET  # noqa: E402


def norm(s):
    s = re.sub(r"\[\+57\.02\d*\]", "", str(s)).replace("I", "L")
    return re.sub(r"\[\+15\.99\d*\]", "[ox]", s)


def ours_matrix(lfq_dir, ids_dir, tq=None):
    pr = pd.read_parquet(os.path.join(lfq_dir, "precursors.parquet"))
    pep = pd.read_parquet(os.path.join(ids_dir, "peptides.parquet")).set_index("pf_id")
    pr["pf"] = pep["modseq"].reindex(pr["pf_id"].to_numpy()).map(norm).to_numpy()
    parts = [pr[["run", "pf", "intensity", "pf_id"]]]
    if tq is not None:
        tr = pd.read_parquet(os.path.join(lfq_dir, "transfers.parquet"))
        tr = tr[~tr["is_decoy_donor"] & (tr["transfer_q"] <= tq)].copy()
        tr["pf"] = pep["modseq"].reindex(tr["pf_id"].to_numpy()).map(norm).to_numpy()
        parts.append(tr[["run", "pf", "intensity", "pf_id"]])
    x = pd.concat(parts)
    M = x.pivot_table(index="pf", columns="run", values="intensity", aggfunc="sum")
    prots = pd.Series(pep["proteins"].reindex(x.drop_duplicates("pf")["pf_id"].to_numpy()).to_numpy(),
                      index=x.drop_duplicates("pf")["pf"].to_numpy())
    return M, prots


def main():
    out = sys.argv[1]
    tabs = {}
    lfq = os.path.join(C.WORK, "lfq", "d9_hye_qtot_d01")
    ids = os.path.join(C.WORK, "ids", "d9_hye")
    tabs["ours_mbr_off"] = ours_matrix(lfq, ids)
    tabs["ours_mbr_q1pct"] = ours_matrix(lfq, ids, 0.01)
    M, p = CMP.sage(os.path.join(C.WORK, "sage", "D9"))
    M.index = [norm(i) for i in M.index]
    M = M.groupby(level=0).sum(min_count=1)
    p.index = [norm(i) for i in p.index]
    tabs["sage_lfq"] = (M, p[~p.index.duplicated()])
    for k in ("ours_hye", "sage_hye"):
        d = os.path.join(C.WORK, "flashlfq_untimed", k)
        if os.path.exists(os.path.join(d, "QuantifiedPeaks.tsv")):
            for mbr in (False, True):
                _, _, pfM, pp, _ = CMP.flashlfq(d, mbr)
                pfM.index = [norm(i) for i in pfM.index]
                pfM = pfM.groupby(level=0).sum(min_count=1)
                pp.index = [norm(i) for i in pp.index]
                tabs[f"flashlfq_{k}_{'mbr_on' if mbr else 'mbr_off'}"] = (pfM, pp[~pp.index.duplicated()])
    hye = [f"LFQ_Astral_DDA_15min_50ng_Condition_{c}_REP{r}" for c in "AB" for r in (1, 2, 3)]
    A, B = hye[:3], hye[3:]
    # analytes quantified (>= 1 value per condition) by every tool
    keep = None
    tabs = {k: (M.reindex(columns=hye).replace(0, np.nan), p) for k, (M, p) in tabs.items()}   # 0 = missing
    for k, (M, p) in tabs.items():
        M = M.reindex(columns=hye)
        ok = set(M.index[(M[A].notna().sum(1) >= 1) & (M[B].notna().sum(1) >= 1)])
        keep = ok if keep is None else keep & ok
    res = {"n_common": len(keep), "tools": list(tabs)}
    for k, (M, p) in tabs.items():
        M = M.reindex(columns=hye).loc[sorted(keep)]
        sp = p.reindex(M.index).map(EV.species_from_proteins)
        Mn, _ = MET.median_normalize(M)
        e = MET.epsilon(Mn, sp, A, B, C.EXPECTED_LOG2FC, 1)
        res[k] = {"median_abs_eps": e["median_abs_eps"], "per_species": {s: v["median_abs_eps"] for s, v in e["per_species"].items()},
                  "n": e["n"], "median_cv_min3": MET.cv_stats(Mn, A, B, 3)["median_cv"],
                  "n_complete": int(M.notna().all(axis=1).sum())}
    json.dump(res, open(out, "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
