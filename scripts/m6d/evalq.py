"""M6d evaluation helpers: build precursor / peptidoform matrices from harness intensities and score them with the
unchanged §6.8 metrics (vulcan.quant.metrics; ProteoBench epsilon, CV, completeness)."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m6"))
import compare_common as CC  # noqa: E402
import evaluate_lfq as EV  # noqa: E402

from vulcan.quant import metrics as MET  # noqa: E402


def long_table(obj, vals: dict, extra: pd.DataFrame | None = None):
    """vals: run -> intensity array aligned with obj['tabs'][run] rows.  extra: optional long rows
    (run, pf_id, charge, intensity) appended (e.g. accepted MBR transfers or sibling-charge extractions)."""
    parts = []
    for run, v in vals.items():
        t = obj["tabs"][run]
        parts.append(pd.DataFrame({"run": run, "pf_id": t["pf_id"].to_numpy(), "charge": t["charge"].to_numpy(),
                                   "intensity": np.asarray(v, float)}))
    if extra is not None and len(extra):
        parts.append(extra[["run", "pf_id", "charge", "intensity"]])
    L = pd.concat(parts, ignore_index=True)
    L = L[~L["pf_id"].isin(obj["dec"]) & np.isfinite(L["intensity"]) & (L["intensity"] > 0)]
    pep = obj["pep"]
    L["pf"] = pep["modseq"].reindex(L["pf_id"].to_numpy()).map(CC.norm).to_numpy()
    L["key"] = L["pf"] + "/" + L["charge"].astype(str)
    L["prot"] = pep["proteins"].reindex(L["pf_id"].to_numpy()).to_numpy()
    return L


def species(prots: pd.Series) -> pd.Series:
    return prots.map(EV.species_from_proteins)


def matrix(L, level="precursor", charge_rule="sum"):
    """level precursor: key = pf/charge.  level peptidoform: charge_rule 'sum' (all charges, as M6), 'dominant'
    (only the charge state quantified in the most runs; ties: higher summed intensity), 'consistent' (sum over the
    charge states present in every run where the peptidoform is quantified)."""
    if level == "precursor":
        M = L.pivot_table(index="key", columns="run", values="intensity", aggfunc="sum")
        P = L.drop_duplicates("key").set_index("key")["prot"]
        return M, P
    if charge_rule == "dominant":
        c = L.groupby(["pf", "charge"]).agg(n=("run", "nunique"), s=("intensity", "sum")).reset_index()
        c = c.sort_values(["pf", "n", "s"], ascending=[True, False, False]).drop_duplicates("pf")
        L = L.merge(c[["pf", "charge"]], on=["pf", "charge"])
    elif charge_rule == "consistent":
        nr = L.groupby("pf")["run"].nunique().rename("nr")
        nc = L.groupby(["pf", "charge"])["run"].nunique().rename("nc").reset_index().merge(nr, on="pf")
        good = nc[nc["nc"] == nc["nr"]][["pf", "charge"]]
        L = L.merge(good, on=["pf", "charge"])
    M = L.pivot_table(index="pf", columns="run", values="intensity", aggfunc="sum")
    P = L.drop_duplicates("pf").set_index("pf")["prot"]
    return M, P


def score(M, P, A, B, expected, rows=None, norm_rows=None):
    """Median-normalise (on norm_rows if given, else all rows), restrict to rows, return the §6.8 numbers."""
    cols = A + B
    M = M.reindex(columns=cols).replace(0, np.nan)
    M = M[M.notna().any(axis=1)]
    sp = species(P.reindex(M.index))
    M = M[sp.isin(list(expected)).to_numpy()]
    _, sh = MET.median_normalize(M if norm_rows is None else M.loc[M.index.intersection(norm_rows)])
    Mn = M / pd.Series({c: 2 ** s for c, s in sh.items()})
    if rows is not None:
        Mn = Mn.reindex(sorted(set(rows) & set(Mn.index)))
    sp = species(P.reindex(Mn.index))
    e1 = MET.epsilon(Mn, sp, A, B, expected, 1)
    out = {"n": e1["n"], "eps": e1["median_abs_eps"],
           "eps_sp": {s: v["median_abs_eps"] for s, v in e1["per_species"].items()},
           "lfc_sp": {s: v["median_log2fc"] for s, v in e1["per_species"].items()}}
    if len(A) >= 3 and len(B) >= 3:
        out["eps_min3"] = MET.epsilon(Mn, sp, A, B, expected, 3)["median_abs_eps"]
        out["cv"] = MET.cv_stats(Mn, A, B, 3)["median_cv"]
    else:
        out["cv2"] = MET.cv_stats(Mn, A, B, 2)["median_cv"]
    out["completeness"] = MET.completeness(Mn, A, B)["completeness"]
    return out, Mn
