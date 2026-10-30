#!/usr/bin/env python3
"""M6d part 1: diagnose the D9 matched-peptide accuracy gap to Sage LFQ (dev).  Sage's per-peptide intensities are the
paired reference.  Writes OUT/diag_d9.json and OUT/diag_tables/*.tsv (+ /tmp/vulcan_m6d/diag_long.parquet).
    scripts/m2b/vpy.sh scripts/m6d/diag_d9.py OUT
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m6"))
import common as C  # noqa: E402
import comparators as CMP  # noqa: E402
import compare_common as CC  # noqa: E402
import evalq as E  # noqa: E402
import harness as H  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.ms1 import averagine as AV  # noqa: E402
from vulcan.quant import metrics as MET  # noqa: E402

OUT = sys.argv[1]
os.makedirs(os.path.join(OUT, "diag_tables"), exist_ok=True)
W = "/tmp/vulcan_m6"
RUNS = [f"LFQ_Astral_DDA_15min_50ng_Condition_{c}_REP{r}" for c in "AB" for r in (1, 2, 3)]
A, B = RUNS[:3], RUNS[3:]
EXP = C.EXPECTED_LOG2FC
res = {}

obj = H.prepare(f"{W}/ids/d9_hye", RUNS, f"{W}/ms1", C.store, cache="/tmp/vulcan_m6d/d9_hye_obj.pkl")
tabs = obj["tabs"]
# ---------------------------------------------------------------- reference: Sage peptide matrix
Ms, Ps = CMP.sage(f"{W}/sage/D9")
Ms.index = [CC.norm(i) for i in Ms.index]
Ms = Ms.groupby(level=0).sum(min_count=1).reindex(columns=RUNS).replace(0, np.nan)
Ps.index = [CC.norm(i) for i in Ps.index]
Ps = Ps[~Ps.index.duplicated()]

base = H.requant(obj, C.store, tabs)
def iso_sum(res_, cols):   # cols: isotope offsets to sum
    j = [list(H.ISO_OFFS).index(c) for c in cols]
    return {r: v[0][:, j].sum(1) for r, v in res_.items()}

L0 = E.long_table(obj, iso_sum(base, [0, 1, 2]))
Mo, Po = E.matrix(L0, "peptidoform", "sum")
def ge1(M):
    return set(M.index[(M[A].notna().sum(axis=1) >= 1) & (M[B].notna().sum(axis=1) >= 1)])
KEEP = sorted(ge1(Mo) & ge1(Ms))
res["n_matched"] = len(KEEP)
SP = E.species(Ps.reindex(KEEP))

def eps_rows(M, rows, norm_rows=KEEP):
    o, Mn = E.score(M, Ps.reindex(M.index).fillna(Po.reindex(M.index)), A, B, EXP, rows=rows, norm_rows=norm_rows)
    return o, Mn

# ---------------------------------------------------------------- 1. decomposition of the matched gap
dec = {}
o_s, Msn = eps_rows(Ms, KEEP)
o_o, Mon = eps_rows(Mo, KEEP)
dec["sage"], dec["ours_M6"] = o_s, o_o
Msm = Ms.reindex(KEEP).where(Mo.reindex(KEEP).notna())
dec["sage_masked_to_our_missingness"] = eps_rows(Msm, KEEP)[0]
full = [k for k in KEEP if Mo.loc[k].notna().all()]
dec["complete_6of6_rows"] = {"n": len(full), "ours": eps_rows(Mo, full)[0]["eps"], "sage": eps_rows(Ms, full)[0]["eps"]}
# charge-set consistency
cs = L0.groupby(["pf", "run"])["charge"].apply(lambda s: tuple(sorted(s)))
incons = set((cs.groupby(level=0).nunique() > 1).pipe(lambda s: s[s].index))
ki = [k for k in KEEP if k in incons]
kc = [k for k in KEEP if k not in incons]
dec["charge_set_varies_across_runs"] = {"n": len(ki), "ours": eps_rows(Mo, ki)[0], "sage": eps_rows(Ms, ki)[0]}
dec["charge_set_constant"] = {"n": len(kc), "ours": eps_rows(Mo, kc)[0], "sage": eps_rows(Ms, kc)[0]}
for rule in ("dominant", "consistent"):
    M2, _ = E.matrix(L0, "peptidoform", rule)
    dec[f"ours_charge_rule_{rule}"] = eps_rows(M2, KEEP)[0]
res["decomposition"] = dec

# ---------------------------------------------------------------- 2. per-species and intensity dependence (paired)
lg_o = np.log2(Mon.reindex(KEEP))
lg_s = np.log2(Msn.reindex(KEEP))
fc = lambda L: L[A].mean(axis=1) - L[B].mean(axis=1)
eo = fc(lg_o) - SP.map(EXP)
es = fc(lg_s) - SP.map(EXP)
I = lg_s.mean(axis=1)
q = pd.qcut(I, 5, labels=False)
rows = []
for sp in ("HUMAN", "YEAST", "ECOLI", "ALL"):
    for b in range(5):
        m = (q == b) & ((SP == sp) if sp != "ALL" else True)
        m &= eo.notna() & es.notna()
        rows.append({"species": sp, "sage_intensity_quintile": b + 1, "n": int(m.sum()),
                     "ours_median_abs_eps": float(eo[m].abs().median()), "sage_median_abs_eps": float(es[m].abs().median()),
                     "ours_median_eps": float(eo[m].median()), "sage_median_eps": float(es[m].median()),
                     "ours_complete_frac": float(Mo.reindex(KEEP)[m.to_numpy()].notna().all(axis=1).mean())})
T2 = pd.DataFrame(rows)
T2.to_csv(os.path.join(OUT, "diag_tables", "species_by_intensity.tsv"), sep="\t", index=False)
res["species_by_intensity"] = rows

# ---------------------------------------------------------------- 3-5. integration variants (precursor level and matched)
def var_eval(vals, label, extra=None):
    L = E.long_table(obj, vals, extra)
    Mp, Pp = E.matrix(L, "precursor")
    op, _ = E.score(Mp, Pp, A, B, EXP)
    Mpf, _ = E.matrix(L, "peptidoform", "sum")
    om = eps_rows(Mpf, KEEP)[0]
    Mpc, _ = E.matrix(L, "peptidoform", "dominant")
    omc = eps_rows(Mpc, KEEP)[0]
    return {"variant": label, "prec_n": op["n"], "prec_eps": op["eps"], "prec_eps_H": op["eps_sp"]["HUMAN"],
            "prec_eps_Y": op["eps_sp"]["YEAST"], "prec_eps_E": op["eps_sp"]["ECOLI"], "prec_lfc_Y": op["lfc_sp"]["YEAST"],
            "prec_lfc_E": op["lfc_sp"]["ECOLI"], "prec_cv": op["cv"], "matched_eps_sum": om["eps"],
            "matched_eps_dominant": omc["eps"], "matched_eps_dominant_E": omc["eps_sp"]["ECOLI"]}

var = []
isos = {"M0": [0], "M1": [1], "M0-1": [0, 1], "M0-2 (M6)": [0, 1, 2], "M0-3": [0, 1, 2, 3], "M0-4": [0, 1, 2, 3, 4]}
for k, cols in isos.items():
    var.append(dict(var_eval(iso_sum(base, cols), f"iso {k}"), group="isotopes"))
# averagine-weighted: estimate the mono-equivalent from M0..M2 by least squares against the envelope
def avg_weighted(res_):
    out = {}
    for r, (a, _, _) in res_.items():
        t = tabs[r]
        env = np.stack([AV.envelope(m, 3) for m in t["calcmass"].to_numpy()])
        X = a[:, 1:4]
        out[r] = (X * env).sum(1) / (env * env).sum(1) * env.sum(1)   # projection onto the envelope, scaled to M0..2
    return out
var.append(dict(var_eval(avg_weighted(base), "iso M0-2 averagine projection"), group="isotopes"))
for wm in (0.5, 0.75, 1.5, 2.0, 3.0):
    rr = H.requant(obj, C.store, tabs, width_mult=wm)
    var.append(dict(var_eval(iso_sum(rr, [0, 1, 2]), f"window x{wm}"), group="window"))
for hw in (0.05, 0.1, 0.2):
    rr = H.requant(obj, C.store, tabs, fixed_half_width=hw)
    var.append(dict(var_eval(iso_sum(rr, [0, 1, 2]), f"fixed half-width {hw} min"), group="window"))
# centring on the aligned consensus RT instead of the run's feature apex
al = obj["align"]
cons = {}
for r in RUNS:
    t = tabs[r]
    ok = np.isfinite(t["rt_apex"].to_numpy())
    rc = np.full(len(t), np.nan)
    rc[ok] = al.to_consensus(r, t["rt_apex"].to_numpy()[ok])
    cons[r] = pd.Series(rc, index=pd.MultiIndex.from_arrays([t["pf_id"], t["charge"]]))
consm = pd.concat(cons, axis=1).median(axis=1)
tabs_pred = {}
for r in RUNS:
    t = tabs[r].copy()
    c = consm.reindex(pd.MultiIndex.from_arrays([t["pf_id"], t["charge"]])).to_numpy()
    t["rt_pred"] = al.to_run(r, c)
    tabs_pred[r] = t
rr = H.requant(obj, C.store, tabs_pred, centre="pred")
var.append(dict(var_eval(iso_sum(rr, [0, 1, 2]), "centre = aligned consensus RT"), group="window"))
# feature hill area (M6 arm)
var.append(dict(var_eval({r: tabs[r]["intensity"].to_numpy() for r in RUNS}, "feature hill area (variable bounds)"), group="window"))
for tk, tm, sm in ((3.0, 3.0, True), (5.0, 5.0, False), (10.0 / 0.45, 10.0, False), (1.5, 2.0, False)):
    rr = H.requant(obj, C.store, tabs, tol_k=tk, tol_min=tm, mode_sum=sm)
    var.append(dict(var_eval(iso_sum(rr, [0, 1, 2]), f"tol max({tk:.1f} SD,{tm} ppm) {'sum' if sm else 'max'}"), group="peak"))
rr = H.requant(obj, C.store, tabs, mz_source="theory")
var.append(dict(var_eval(iso_sum(rr, [0, 1, 2]), "m/z from peptide (not feature)"), group="peak"))
var.insert(0, dict(var_eval(iso_sum(base, [0, 1, 2]), "M6 baseline (XIC M0-2, x1 window, apex)"), group="baseline"))
TV = pd.DataFrame(var)
TV.to_csv(os.path.join(OUT, "diag_tables", "integration_variants.tsv"), sep="\t", index=False)
res["integration_variants"] = var

# ---------------------------------------------------------------- 6. per-value paired residuals vs Sage (single-charge peptides)
a0 = base
recs = []
for r in RUNS:
    t = tabs[r]
    a, nsc, apx = a0[r]
    f = obj["ms1"][r].feat
    fid = t["feature_id"].to_numpy()
    ok = fid >= 0
    env = np.stack([AV.envelope(m, 4) for m in t["calcmass"].to_numpy()])
    obs = a[:, 1:5]
    cos = (obs * env).sum(1) / np.sqrt((obs ** 2).sum(1) * (env ** 2).sum(1) + 1e-300)
    rc = np.full(len(t), np.nan)
    rc[ok] = al.to_consensus(r, t["rt_apex"].to_numpy()[ok])
    recs.append(pd.DataFrame({
        "run": r, "pf_id": t["pf_id"].to_numpy(), "charge": t["charge"].to_numpy(), "n_psms": t["n_psms"].to_numpy(),
        "feature_id": fid, "x3": a[:, 1:4].sum(1), "m_minus1_ratio": a[:, 0] / np.maximum(a[:, 1], 1e-9),
        "xic_iso_cos": cos, "feat_iso_cos": t["iso_cos"].to_numpy(), "feat_trace_cos": t["trace_cos"].to_numpy(),
        "rt_apex": t["rt_apex"].to_numpy(), "rt_cons": rc, "xic_apex": apx, "n_scans": nsc,
        "feat_n_iso": np.where(ok, f["n_isotopes"].to_numpy()[np.where(ok, fid, 0)], 0),
        "dt_psm_apex": t["rt_psm"].to_numpy() - t["rt_apex"].to_numpy()}))
D = pd.concat(recs, ignore_index=True)
D = D[~D["pf_id"].isin(obj["dec"]) & np.isfinite(D["x3"]) & (D["x3"] > 0)]
D["pf"] = obj["pep"]["modseq"].reindex(D["pf_id"].to_numpy()).map(CC.norm).to_numpy()
D["rt_cons_med"] = D.groupby(["pf_id", "charge"])["rt_cons"].transform("median")
D["n_runs_prec"] = D.groupby(["pf_id", "charge"])["run"].transform("nunique")
D["d_rt_cons"] = D["rt_cons"] - D["rt_cons_med"]
sd_rt = {r: obj["align"].maps[r][3] for r in RUNS}
D["d_rt_cons_sd"] = D["d_rt_cons"] / D["run"].map(sd_rt)
# paired residual: our normalised log2 minus Sage's, centred per peptide (removes the per-peptide scale difference)
single = D.groupby("pf")["charge"].nunique()
single = set(single.index[single == 1]) & set(KEEP)
X = D[D["pf"].isin(single)].copy()
sh_o = {c: float(np.log2(Mo[c] / Mon[c]).dropna().median()) for c in RUNS}
sh_s = {c: float(np.log2(Ms[c] / Msn[c]).dropna().median()) for c in RUNS}
X["lo"] = np.log2(X["x3"]) - X["run"].map(sh_o)
X["ls"] = [np.log2(Ms.at[p, r]) - sh_s[r] if np.isfinite(Ms.at[p, r]) else np.nan for p, r in zip(X["pf"], X["run"])]
X["d"] = X["lo"] - X["ls"]
X["resid"] = X["d"] - X.groupby("pf")["d"].transform("median")
X = X[np.isfinite(X["resid"])]
X["n_runs_pf"] = X.groupby("pf")["run"].transform("nunique")
X["species"] = SP.reindex(X["pf"].to_numpy()).to_numpy()
X["log2_sage"] = X["ls"]
X.to_parquet("/tmp/vulcan_m6d/diag_long.parquet")   # 8.6 MB, not committed
big = X["n_runs_pf"] >= 3
summ = {"n_values": int(big.sum()), "median_abs_resid": float(X.loc[big, "resid"].abs().median()),
        "frac_abs_resid_gt_0.5": float((X.loc[big, "resid"].abs() > 0.5).mean()),
        "frac_abs_resid_gt_1": float((X.loc[big, "resid"].abs() > 1).mean())}
def strat(col, bins, labels=None):
    out = []
    cut = pd.cut(X.loc[big, col], bins, labels=labels)
    for b, g in X[big].groupby(cut, observed=True):
        out.append({"factor": col, "bin": str(b), "n": int(len(g)), "median_abs_resid": float(g["resid"].abs().median()),
                    "frac_gt_0.5": float((g["resid"].abs() > 0.5).mean()), "median_resid": float(g["resid"].median())})
    return out
st = []
st += strat("d_rt_cons_sd", [-np.inf, -5, -3, -1, 1, 3, 5, np.inf])
st += strat("m_minus1_ratio", [-np.inf, 0.05, 0.15, 0.3, 0.6, np.inf])
st += strat("xic_iso_cos", [-np.inf, 0.8, 0.9, 0.95, 0.98, 0.99, np.inf])
st += strat("feat_iso_cos", [-np.inf, 0.8, 0.9, 0.95, 0.99, np.inf])
st += strat("feat_n_iso", [0, 2, 3, 4, 20])
st += strat("n_scans", [0, 3, 5, 8, 12, 1e9])
st += strat("log2_sage", [-np.inf, 18, 20, 22, 24, np.inf])
st += strat("dt_psm_apex", [-np.inf, -0.1, -0.03, 0.03, 0.1, np.inf])
TS = pd.DataFrame(st)
TS.to_csv(os.path.join(OUT, "diag_tables", "paired_residual_strata.tsv"), sep="\t", index=False)
res["paired_residuals"] = {"summary": summ, "strata": st}
# linking: how often is the run's feature apex off the aligned consensus by > 3 SD (wrong feature)?
D3 = D[D["n_runs_prec"] >= 3]
res["linking"] = {"precursor_values": int(len(D3)),
                  "frac_apex_off_consensus_gt3sd": float((D3["d_rt_cons_sd"].abs() > 3).mean()),
                  "frac_apex_off_consensus_gt0.1min": float((D3["d_rt_cons"].abs() > 0.1).mean()),
                  "frac_xic_apex_ne_feature_apex_gt0.05min": float(((D3["xic_apex"] - D3["rt_apex"]).abs() > 0.05).mean()),
                  "alignment_resid_sd_min": sd_rt}
json.dump(res, open(os.path.join(OUT, "diag_d9.json"), "w"), indent=1, default=float)
pd.set_option("display.width", 250)
print(json.dumps(res["decomposition"], indent=1, default=float)[:3000])
print(TV.to_string())
print(T2[T2.species == "ALL"].to_string())
print(TS.to_string())
print(json.dumps(summ), json.dumps(res["linking"], default=float))
