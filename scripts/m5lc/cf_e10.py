#!/usr/bin/env python3
"""E10 on D10cf (dev): CF-MS transfer modes.  Input: per-fraction ids (scripts/m5lc/ids.py --per-run) over an index
with paired shuffled entrapment (seed 2001, r = 1).

Modes (identifications are identical in all modes: run-specific per-fraction q <= alpha, MS2 only):
  none           MS2-identified values only
  local          + fraction-local transfers (donors in adjacent fractions f+-1), transfer q <= tq (DESIGN default)
  global         + unconstrained MBR (donors in any fraction), transfer q <= tq
  global_naive   + unconstrained MBR without a transfer-FDR filter (any target peak in the window; MaxQuant-like)
All transfers use the same profile-agnostic score and windows (vulcan.cf.transfer); only the donor set differs.
Intensities: fixed-window XIC (M6d XICParams), no normalization of any kind.

Metrics per mode:
  * matrix size, values by provenance;
  * elution peak width per precursor and per protein (fractions >= 50% of the profile max; fractions with any
    value), and the inflation vs 'none' on the same precursors/proteins;
  * CORUM co-complex PR-AUC: protein pairs with both proteins in CORUM (OmniPath CORUM subset), positive = share a
    complex; score = Pearson r of the raw protein elution profiles (protein = sum of its proteotypic precursors;
    missing = 0); reported on each mode's own eligible proteins (>= 3 non-zero fractions) and on the pairs eligible
    in every mode;
  * entrapment null ("false co-elution"): entrapment precursors (family E) present in the matrix; combined
    entrapment FDP estimate of matrix values 2 N_E / (N_T + N_E) (r = 1) split by provenance; entrapment profiles
    vs all protein profiles: fraction of entrapment profiles with max r >= 0.8 to some real protein (false co-elution),
    vs the same rate for real proteotypic precursors of proteins not co-complexed... (reported as counts).
    python scripts/m5lc/cf_e10.py --ids DIR --out DIR
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

from vulcan.cf.transfer import CFTransferParams, cf_transfers  # noqa: E402
from vulcan.quant import lfq as Q  # noqa: E402

MODES = ("none", "local", "global", "global_naive")


def corum_pairs(path, proteins):
    c = pd.read_csv(path, sep="\t")
    c = c[c["sources"].str.contains("CORUM")]
    comp = [set(x.split("_")) for x in c["components"]]
    pset = set(proteins)
    members = set().union(*comp) & pset
    pos = set()
    for s in comp:
        s = sorted(s & pset)
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                pos.add((s[i], s[j]))
    return members, pos, len(c)


def width_stats(M):
    """M: entity x fraction (ordered) matrix, NaN = missing."""
    X = M.fillna(0).to_numpy()
    mx = X.max(1, keepdims=True)
    half = (X >= 0.5 * np.where(mx > 0, mx, np.inf)).sum(1)
    nz = (X > 0).sum(1)
    return pd.DataFrame({"fwhm_fracs": half, "n_fracs": nz}, index=M.index)


def pr_auc(y, s):
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(y, s)) if y.sum() and (~y.astype(bool)).sum() else float("nan")


def profile_corr(A):
    X = A.fillna(0).to_numpy(float)
    X = X - X.mean(1, keepdims=True)
    n = np.linalg.norm(X, axis=1, keepdims=True)
    X = X / np.where(n > 0, n, 1)
    return X @ X.T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--donor-alpha", type=float, default=0.001)
    ap.add_argument("--transfer-q", type=float, default=0.01)
    ap.add_argument("--corum", default=os.path.join(C.WORK, "truth", "omnipath_corum.tsv"))
    ap.add_argument("--min-fracs", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = {}
    ids = pd.read_parquet(os.path.join(a.ids, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(a.ids, "peptides.parquet")).set_index("pf_id")
    acc = ids[ids[f"pass_{a.alpha}"]].copy()
    runs = sorted(acc["run"].unique(), key=C.cf_fraction)
    frac = {r: C.cf_fraction(r) for r in runs}
    # ---------------- link + precursor tables (per fraction)
    t0 = time.perf_counter()
    ms1, tabs = {}, {}
    for run in runs:
        r = Q.load_run_ms1(os.path.join(C.WORK, "ms1"), C.store(run), run)
        ms1[run] = r
        tabs[run] = Q.precursor_table(Q.link_psms(acc[acc["run"] == run], r), r)
    T["link_s"] = time.perf_counter() - t0
    dec_pf = pep.index[pep["is_decoy"]].to_numpy()
    t0 = time.perf_counter()
    al = Q.Alignment({r: t[~t["pf_id"].isin(dec_pf)] for r, t in tabs.items()})
    T["align_s"] = time.perf_counter() - t0
    dd = ids[ids[f"pass_{a.donor_alpha}"]]
    dk = set(zip(dd["pf_id"], dd["charge"].astype(np.int8)))
    trs, infos = {}, {}
    for m in ("local", "global"):
        t0 = time.perf_counter()
        trs[m], infos[m] = cf_transfers(tabs, ms1, al, dec_pf, frac, CFTransferParams(mode=m), donor_keys=dk)
        T[f"transfer_{m}_s"] = time.perf_counter() - t0
    # ---------------- XIC re-quant (one pass over MS2 tables + all candidate transfers)
    t0 = time.perf_counter()
    allc = pd.concat([trs["local"].assign(mode="local"), trs["global"].assign(mode="global")], ignore_index=True)
    allc = allc[allc["feature_id"] >= 0].copy()
    xa = Q.xic_requant(tabs, ms1, {r: C.store(r) for r in runs}, al, transfers=allc, params=Q.XICParams())
    for run in runs:
        tabs[run]["intensity"] = xa[run]
    allc["intensity"] = allc["xic"]
    T["xic_s"] = time.perf_counter() - t0
    # ---------------- long tables per mode
    ms2 = pd.concat(tabs.values(), ignore_index=True)
    ms2 = ms2[~ms2["pf_id"].isin(dec_pf) & (ms2["feature_id"] >= 0)]
    ms2 = ms2.assign(provenance="ms2", mode="none")
    tg = allc[~allc["is_decoy_donor"]]
    sel = {"none": ms2.iloc[:0],
           "local": tg[(tg["mode"] == "local") & (tg["transfer_q"] <= a.transfer_q)],
           "global": tg[(tg["mode"] == "global") & (tg["transfer_q"] <= a.transfer_q)],
           "global_naive": tg[(tg["mode"] == "global")]}
    cols = ["run", "pf_id", "charge", "intensity", "provenance"]
    pep_t = pep[~pep["is_decoy"]]
    ent_pf = set(pep_t.index[pep_t["is_entrapment"]])
    one = pep_t[~pep_t["is_entrapment"] & (pep_t["proteins"].str.count(";") == 0) & (pep_t["proteins"] != "")]
    prot_of = one["proteins"].map(lambda s: s.split("|")[1] if s.count("|") >= 2 else s)
    corum_members, corum_pos, n_cplx = corum_pairs(a.corum, set(prot_of))
    res = {"runs": runs, "fractions": [frac[r] for r in runs], "alpha": a.alpha, "donor_alpha": a.donor_alpha,
           "transfer_q": a.transfer_q, "transfer_info": infos,
           "corum": {"n_complexes_corum_subset": n_cplx, "n_member_proteins_in_proteome": len(corum_members),
                     "n_positive_pairs_possible": len(corum_pos)},
           "ids_per_fraction": acc[~acc["is_decoy"]].groupby("run")["pf_id"].nunique().to_dict(),
           "entrapment_ids_per_fraction": acc[~acc["is_decoy"] & acc["pf_id"].isin(ent_pf)].groupby("run")["pf_id"].nunique().to_dict(),
           "modes": {}}
    order = runs
    prof_by_mode, pw_by_mode, ppw_by_mode = {}, {}, {}
    for m in MODES:
        L = pd.concat([ms2[cols], sel[m][cols]], ignore_index=True)
        L = L[np.isfinite(L["intensity"]) & (L["intensity"] > 0)]
        L["key"] = L["pf_id"].astype(str) + "/" + L["charge"].astype(str)
        M = L.pivot_table(index="key", columns="run", values="intensity", aggfunc="sum").reindex(columns=order)
        kpf = M.index.map(lambda k: np.uint64(k.split("/")[0]))
        is_e = np.isin(np.array(kpf, dtype=np.uint64), np.array(sorted(ent_pf), dtype=np.uint64))
        Lv = L.assign(is_e=L["pf_id"].isin(ent_pf))
        nT = Lv[~Lv["is_e"]].groupby("provenance").size().to_dict()
        nE = Lv[Lv["is_e"]].groupby("provenance").size().to_dict()
        fdp = {p: 2 * nE.get(p, 0) / max(nT.get(p, 0) + nE.get(p, 0), 1) for p in set(nT) | set(nE)}
        # precursor widths (target, non-entrapment)
        Mt = M[~is_e]
        pw = width_stats(Mt)
        pw_by_mode[m] = pw
        # protein profiles: sum of proteotypic precursors
        kp = pd.Series(np.array(Mt.index.map(lambda k: np.uint64(k.split("/")[0])), dtype=np.uint64), index=Mt.index)
        prot = kp.map(lambda u: prot_of.get(u, None))
        Pm = Mt[prot.notna().to_numpy()].groupby(prot.dropna().to_numpy()).sum(min_count=1)
        prof_by_mode[m] = Pm
        ppw = width_stats(Pm)
        ppw_by_mode[m] = ppw
        # entrapment null profiles
        Me = M[is_e]
        rec = {"matrix_precursors": int(len(M)), "values_by_provenance_target": nT, "values_by_provenance_entrapment": nE,
               "entrapment_combined_FDP_by_provenance": fdp,
               "entrapment_combined_FDP_all_values": 2 * sum(nE.values()) / max(sum(nT.values()) + sum(nE.values()), 1),
               "entrapment_precursors": int(len(Me)),
               "median_precursor_fwhm_fracs": float(pw["fwhm_fracs"].median()),
               "median_precursor_n_fracs": float(pw["n_fracs"].median()),
               "proteins": int(len(Pm)), "median_protein_fwhm_fracs": float(ppw["fwhm_fracs"].median()),
               "median_protein_n_fracs": float(ppw["n_fracs"].median())}
        # CORUM PR-AUC on this mode's own eligible proteins
        elig = [p for p in Pm.index if p in corum_members and (Pm.loc[p].fillna(0) > 0).sum() >= a.min_fracs]
        if len(elig) > 2:
            R = profile_corr(Pm.loc[elig])
            iu = np.triu_indices(len(elig), 1)
            pa_ = np.array(elig)[iu[0]]
            pb_ = np.array(elig)[iu[1]]
            y = np.array([(x, z) in corum_pos or (z, x) in corum_pos for x, z in zip(pa_, pb_)])
            rec["corum_own"] = {"proteins": len(elig), "pairs": int(len(y)), "positives": int(y.sum()),
                                "pr_auc": pr_auc(y, R[iu]), "baseline_prevalence": float(y.mean())}
        # false co-elution from entrapment profiles: max |r| to any real protein profile (>= min_fracs values)
        Pe = Pm[(Pm.fillna(0) > 0).sum(axis=1) >= a.min_fracs]
        Ee = Me[(Me.fillna(0) > 0).sum(axis=1) >= a.min_fracs]
        if len(Ee) and len(Pe):
            X = profile_corr(pd.concat([Ee, Pe]))
            cross = X[:len(Ee), len(Ee):]
            mx = cross.max(1)
            rec["entrapment_null"] = {"entrapment_profiles_ge_min_fracs": int(len(Ee)),
                                      "frac_with_max_r_ge_0.8": float((mx >= 0.8).mean()),
                                      "n_with_max_r_ge_0.8": int((mx >= 0.8).sum()),
                                      "n_pairs_r_ge_0.8": int((cross >= 0.8).sum()),
                                      "median_max_r": float(np.median(mx))}
        else:
            rec["entrapment_null"] = {"entrapment_profiles_ge_min_fracs": int(len(Ee))}
        # transfer placement: distance (fractions) of each transferred value to the nearest MS2 detection of the
        # same precursor and of the same protein; >= 3 fractions away = an isolated, off-peak value (a secondary
        # "peak" that MS2 never saw), the co-elution artifact the DESIGN hypothesis is about
        if len(sel[m]):
            x = sel[m][~sel[m]["pf_id"].isin(ent_pf)]
            d_prec, d_prot = [], []
            ms2f = ms2.assign(f=ms2["run"].map(frac).astype(int))
            pf_fr = ms2f.groupby(["pf_id", "charge"])["f"].apply(np.array).to_dict()
            ms2f["prot"] = ms2f["pf_id"].map(prot_of)
            pr_fr = ms2f.dropna(subset=["prot"]).groupby("prot")["f"].apply(np.unique).to_dict()
            for pfid, z, run in zip(x["pf_id"], x["charge"], x["run"]):
                f = frac[run]
                d_prec.append(np.abs(pf_fr[(pfid, z)] - f).min() if (pfid, z) in pf_fr else 99)
                pr = prot_of.get(pfid)
                d_prot.append(np.abs(pr_fr[pr] - f).min() if pr in pr_fr else 99)
            d_prec, d_prot = np.array(d_prec), np.array(d_prot)
            rec["transfer_placement"] = {
                "n_target_transfers": int(len(x)),
                "frac_ge3_from_same_precursor_ms2": float((d_prec >= 3).mean()),
                "frac_ge3_from_same_protein_ms2": float((d_prot[d_prot < 99] >= 3).mean()) if (d_prot < 99).any() else None,
                "median_dist_precursor": float(np.median(d_prec))}
            # PIP-ECHO estimate of the false fraction among the accepted set: decoy-window wins at the same acceptance
            allm = tg[tg["mode"] == ("global" if m.startswith("global") else m)]
            if m == "global_naive":
                rec["transfer_placement"]["est_false_frac_rt_decoy_window_has_peak"] = float(
                    np.isfinite(allm["Sd_rt"]).sum() / max(np.isfinite(allm["St"]).sum(), 1))
        res["modes"][m] = rec
        M.to_parquet(os.path.join(a.out, f"precursor_matrix_{m}.parquet"))
        Pm.to_parquet(os.path.join(a.out, f"protein_matrix_{m}.parquet"))
    # ---------------- width inflation on common entities, PR-AUC on common pairs
    base = pw_by_mode["none"].index
    for m in MODES:
        w = pw_by_mode[m].reindex(base)
        b = pw_by_mode["none"]
        res["modes"][m]["precursor_width_inflation_common"] = {
            "n": int(len(base)), "mean_fwhm_ratio": float((w["fwhm_fracs"] / b["fwhm_fracs"]).mean()),
            "mean_n_fracs_ratio": float((w["n_fracs"] / b["n_fracs"]).mean()),
            "median_n_fracs": float(w["n_fracs"].median()), "median_fwhm": float(w["fwhm_fracs"].median())}
        pb = ppw_by_mode["none"]
        pwm = ppw_by_mode[m].reindex(pb.index)
        res["modes"][m]["protein_width_inflation_common"] = {
            "n": int(len(pb)), "mean_fwhm_ratio": float((pwm["fwhm_fracs"] / pb["fwhm_fracs"]).mean()),
            "mean_n_fracs_ratio": float((pwm["n_fracs"] / pb["n_fracs"]).mean())}
    common = None
    for m in MODES:
        Pm = prof_by_mode[m]
        e = {p for p in Pm.index if p in corum_members and (Pm.loc[p].fillna(0) > 0).sum() >= a.min_fracs}
        common = e if common is None else common & e
    common = sorted(common)
    iu = np.triu_indices(len(common), 1)
    pa_, pb_ = np.array(common)[iu[0]], np.array(common)[iu[1]]
    y = np.array([(x, z) in corum_pos or (z, x) in corum_pos for x, z in zip(pa_, pb_)])
    for m in MODES:
        R = profile_corr(prof_by_mode[m].loc[common])
        res["modes"][m]["corum_common"] = {"proteins": len(common), "pairs": int(len(y)), "positives": int(y.sum()),
                                           "pr_auc": pr_auc(y, R[iu]), "baseline_prevalence": float(y.mean())}
    T["total_s"] = time.perf_counter() - t_start
    res["timing_s"] = T
    allc.to_parquet(os.path.join(a.out, "transfer_candidates.parquet"))
    json.dump(res, open(os.path.join(a.out, "e10.json"), "w"), indent=1, default=float)
    for m in MODES:
        r = res["modes"][m]
        print(m, json.dumps({k: r.get(k) for k in ("matrix_precursors", "entrapment_combined_FDP_by_provenance",
                                                   "median_protein_fwhm_fracs", "corum_own", "corum_common",
                                                   "entrapment_null", "precursor_width_inflation_common")}, default=float))


if __name__ == "__main__":
    main()
