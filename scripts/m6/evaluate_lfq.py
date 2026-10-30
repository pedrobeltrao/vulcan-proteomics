"""Common LFQ evaluation (§6.8) for ours / Sage LFQ / FlashLFQ on D9 HYE: identical post-processing for every tool
(per-run label-free median normalisation, Top3 / MaxLFQ / directLFQ protein roll-up on the same precursor matrix,
ProteoBench epsilon, CV, completeness, DE pAUC).  Also reports the un-normalised precursor epsilon (ProteoBench
uses tool output as-is)."""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.quant import lfq as Q  # noqa: E402
from vulcan.quant import metrics as MET  # noqa: E402

SP_OK = set(C.EXPECTED_LOG2FC)


def species_from_proteins(prots: str) -> str:
    ps = [p for p in str(prots).split(";") if p and "Cont_" not in p and not p.startswith("rev_")]
    if any("Cont_" in p for p in str(prots).split(";")):
        return "CONT"
    sp = {p.rsplit("_", 1)[-1] for p in ps}
    return sp.pop() if len(sp) == 1 else ("MULTI" if sp else "NONE")


def protein_key(prots: str) -> str:
    return ";".join(sorted(p for p in str(prots).split(";") if p))


def evaluate(M: pd.DataFrame, prots: pd.Series, workdir: str, label: str, rollups=("top3", "maxlfq", "directlfq")):
    """M: analyte x run linear intensity (NaN = missing), columns = run names; prots: analyte -> protein string."""
    runs = list(M.columns)
    A = [r for r in runs if C.condition(r) == "A"]
    B = [r for r in runs if C.condition(r) == "B"]
    M = M[A + B].replace(0, np.nan)
    M = M[M.notna().any(axis=1)]
    sp = prots.reindex(M.index).map(species_from_proteins)
    keep = sp.isin(SP_OK)
    M, sp = M[keep], sp[keep]
    pk = prots.reindex(M.index).map(protein_key)
    t0 = time.perf_counter()
    Mn, shifts = MET.median_normalize(M)
    timing = {"normalize_s": time.perf_counter() - t0}
    P = {}
    for r in rollups:
        t0 = time.perf_counter()
        if r == "top3":
            P[r] = Q.topn(Mn, pk, 3)
        elif r == "maxlfq":
            P[r] = Q.maxlfq(Mn, pk)
        elif r == "directlfq":
            try:
                P[r] = Q.directlfq_rollup(M, pk, os.path.join(workdir, f"directlfq_{label}"))
            except Exception as e:  # pragma: no cover
                print("directLFQ failed:", repr(e), flush=True)
        timing[f"rollup_{r}_s"] = time.perf_counter() - t0
    psp = pd.Series({p: species_from_proteins(p) for p in set(pk)})
    rep = MET.full_report(Mn, pk, sp, P, psp, A, B, C.EXPECTED_LOG2FC)
    rep["precursor"]["epsilon_min1_unnormalized"] = MET.epsilon(M, sp, A, B, C.EXPECTED_LOG2FC, 1)
    rep["normalization_log2_shifts"] = shifts
    rep["timing_s"] = timing
    rep["n_by_species"] = sp.value_counts().to_dict()
    return rep, Mn, P
