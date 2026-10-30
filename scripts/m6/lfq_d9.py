#!/usr/bin/env python3
"""Vulcan LFQ on D9 (dev): accepted PSMs (scripts/m6/ids.py) -> MS1 features (vulcan.ms1) -> RT alignment ->
optional MBR (transfer q-values) -> normalisation / roll-up -> §6.8 metrics.

    python scripts/m6/lfq_d9.py --ids DIR --out DIR [--mbr] [--alpha 0.01] [--transfer-q 0.01]

Writes OUT/precursors.parquet (run x precursor, provenance ms2|mbr, transfer_q), OUT/transfers.parquet (all MBR
candidates incl. decoy columns), OUT/lfq_metrics.json (HYE metrics MBR off / on, alignment, linking, timings,
and - if human-only runs are in the experiment - the MBR false-transfer truth check).
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
import evaluate_lfq as EV  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.quant import lfq as Q  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ms1", default=os.path.join(C.WORK, "ms1"))
    ap.add_argument("--mbr", action="store_true")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--transfer-q", type=float, nargs="+", default=[0.01, 0.05])
    ap.add_argument("--donor-alpha", type=float, default=0.001,
                    help="MBR donor ID level (default 0.1%%: with 1%% donors the propagated ID error alone is ~1%%)")
    ap.add_argument("--qcol", default="transfer_q", help="transfer_q (total) | transfer_q_peak")
    ap.add_argument("--allow-claimed", action="store_true", help="MBR may transfer onto features claimed by MS2 IDs")
    ap.add_argument("--pooled-q", action="store_true", help="transfer q pooled over acceptor runs (default per run)")
    ap.add_argument("--intensity", default="xic", choices=["xic", "area"],
                    help="xic: fixed-window XIC re-integration (default); area: feature hill areas")
    ap.add_argument("--xic", default="m6d", choices=["m6", "m6d"],
                    help="m6d (default since M6d): averagine-projected M0..M2, 1.5x window, charge-sibling completion; "
                         "m6: plain M0..M2 sum, 1x window, no siblings (reproduces results/2026-10-27_M6_quant)")
    ap.add_argument("--no-eval", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = {}
    t0 = time.perf_counter()
    ids = pd.read_parquet(os.path.join(a.ids, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(a.ids, "peptides.parquet")).set_index("pf_id")
    col = f"pass_{a.alpha}"
    acc = ids[ids[col]].copy()
    runs = sorted(acc["run"].unique())
    T["load_ids_s"] = time.perf_counter() - t0
    # ---------------- link + precursor tables
    t0 = time.perf_counter()
    ms1, tabs, link_stats = {}, {}, {}
    for run in runs:
        r = Q.load_run_ms1(a.ms1, C.store(run), run)
        ms1[run] = r
        L = Q.link_psms(acc[acc["run"] == run], r)
        tabs[run] = Q.precursor_table(L, r)
        link_stats[run] = dict(r.stats)
        link_stats[run]["precursors"] = int(len(tabs[run]))
        link_stats[run]["precursors_with_feature"] = int((tabs[run]["feature_id"] >= 0).sum())
    T["link_s"] = time.perf_counter() - t0
    dec_pf = pep.index[pep["is_decoy"]].to_numpy()
    tgt_tabs = {r: t[~t["pf_id"].isin(dec_pf)] for r, t in tabs.items()}
    t0 = time.perf_counter()
    al = Q.Alignment(tgt_tabs)
    T["align_s"] = time.perf_counter() - t0
    allp = pd.concat(tabs.values(), ignore_index=True)
    out = {"runs": runs, "alpha": a.alpha, "linking": link_stats, "alignment": al.summary()}
    # ---------------- MBR
    tr = None
    if a.mbr:
        t0 = time.perf_counter()
        dk = None
        if a.donor_alpha is not None and a.donor_alpha != a.alpha:
            dd = ids[ids[f"pass_{a.donor_alpha}"]]
            dk = set(zip(dd["pf_id"], dd["charge"].astype(np.int8)))
        tr, info = Q.match_between_runs(tabs, ms1, al, dec_pf, donor_keys=dk,
                                        params=Q.MBRParams(exclude_claimed=not a.allow_claimed, per_run_q=not a.pooled_q))
        info["exclude_claimed"] = not a.allow_claimed
        info["per_run_q"] = not a.pooled_q
        info["donor_alpha"] = a.donor_alpha or a.alpha
        info["qcol"] = a.qcol
        T["mbr_s"] = time.perf_counter() - t0
        out["mbr"] = info
        tr.to_parquet(os.path.join(a.out, "transfers.parquet"))
    # ---------------- intensities
    if a.intensity == "xic":
        t0 = time.perf_counter()
        xp = Q.XICParams() if a.xic == "m6d" else None
        xa = Q.xic_requant(tabs, ms1, {r: C.store(r) for r in runs}, al, transfers=tr, params=xp)
        for run in runs:
            tabs[run]["intensity_area"] = tabs[run]["intensity"]
            tabs[run]["intensity"] = xa[run]
        if tr is not None:
            tr["intensity_area"] = tr["intensity"]
            tr["intensity"] = tr["xic"]
            tr.to_parquet(os.path.join(a.out, "transfers.parquet"))
        allp = pd.concat(tabs.values(), ignore_index=True)
        if xp is not None and xp.siblings:     # within-run charge-sibling values (provenance ms2_sibling, never IDs)
            sib = Q.charge_siblings(tabs, ms1, {r: C.store(r) for r in runs}, xp, exclude_pf=dec_pf)
            sib["feature_id"] = -2          # no feature: located from the run's own ID of another charge state
            allp = pd.concat([allp, sib], ignore_index=True)
            out["n_sibling_values"] = int(len(sib))
        T["xic_requant_s"] = time.perf_counter() - t0
    out["intensity"] = a.intensity
    # ---------------- precursor matrix
    def key(df):
        ms = pep["modseq"].reindex(df["pf_id"].to_numpy()).to_numpy()
        return pd.Series(ms, index=df.index).astype(str) + "/" + df["charge"].astype(str)
    tg = allp[~allp["pf_id"].isin(dec_pf) & ((allp["feature_id"] >= 0) | (allp["feature_id"] == -2))].copy()
    tg["key"] = key(tg)
    tg[["run", "pf_id", "charge", "key", "intensity", "rt_apex", "provenance", "transfer_q", "n_psms",
        "feature_id"]].to_parquet(os.path.join(a.out, "precursors.parquet"))
    prots = pd.Series(pep["proteins"].reindex(tg["pf_id"].to_numpy()).to_numpy(), index=tg["key"].to_numpy())
    prots = prots[~prots.index.duplicated()]
    M0 = tg.pivot_table(index="key", columns="run", values="intensity", aggfunc="sum")
    # identification counts (precursors / peptidoforms with an accepted PSM; MS2 only, never transfers)
    tall = allp[~allp["pf_id"].isin(dec_pf) & (allp["provenance"] == "ms2")]
    out["ids"] = {"precursors_identified_any_run": int(tall[["pf_id", "charge"]].drop_duplicates().shape[0]),
                  "precursors_with_ms1_feature_any_run": int(tg[["pf_id", "charge"]].drop_duplicates().shape[0]),
                  "peptidoforms_identified": int(tall["pf_id"].nunique()),
                  "per_run_precursors": tall.groupby("run").size().to_dict(),
                  "per_run_quantified": tg.groupby("run").size().to_dict()}
    hye = [r for r in runs if C.condition(r) in "AB"]
    res = {}
    if not a.no_eval and len(hye) == 6:
        t0 = time.perf_counter()
        res["mbr_off"], _, _ = EV.evaluate(M0[hye], prots, a.out, "mbr_off")
        # peptidoform level (sum over charge states; comparable with Sage lfq.tsv combine_charge_states)
        pf = M0[hye].copy()
        pf.index = [k.rsplit("/", 1)[0] for k in pf.index]
        pf = pf.groupby(level=0).sum(min_count=1)
        pprot = prots.copy()
        pprot.index = [k.rsplit("/", 1)[0] for k in pprot.index]
        pprot = pprot[~pprot.index.duplicated()]
        res["mbr_off_peptidoform"], _, _ = EV.evaluate(pf, pprot, a.out, "mbr_off_pf", rollups=())
        T["eval_s"] = time.perf_counter() - t0
    if tr is not None:
        for tq in a.transfer_q:
            ok = (~tr["is_decoy_donor"]) & (tr[a.qcol] <= tq)
            x = tr[ok].copy()
            x["key"] = key(x)
            Mt = x.pivot_table(index="key", columns="run", values="intensity", aggfunc="sum")
            M1 = M0.combine_first(Mt) if len(Mt) else M0
            p2 = pd.concat([prots, pd.Series(pep["proteins"].reindex(x["pf_id"].to_numpy()).to_numpy(),
                                              index=x["key"].to_numpy())])
            p2 = p2[~p2.index.duplicated()]
            rec = {"n_transfers_accepted": int(ok.sum()),
                   "per_run": x.groupby("run").size().to_dict()}
            if not a.no_eval and len(hye) == 6:
                rec["metrics"], _, _ = EV.evaluate(M1.reindex(columns=hye), p2, a.out, f"mbr_q{tq}")
                pf = M1.reindex(columns=hye).copy()
                pf.index = [k.rsplit("/", 1)[0] for k in pf.index]
                pf = pf.groupby(level=0).sum(min_count=1)
                pp = p2.copy()
                pp.index = [k.rsplit("/", 1)[0] for k in pp.index]
                pp = pp[~pp.index.duplicated()]
                rec["metrics_peptidoform"], _, _ = EV.evaluate(pf, pp, a.out, f"mbr_q{tq}_pf", rollups=())
            res[f"mbr_on_q{tq}"] = rec
        # decoy-model internals
        tgt = ~tr["is_decoy_donor"]
        St, Sr, Sm = tr["St"], tr["Sd_rt"], tr["Sd_mz"]
        out["mbr_competition"] = {
            "candidates": int(tgt.sum()), "decoy_donor_candidates": int((~tgt).sum()),
            "target_wins_both": int((tgt & np.isfinite(St) & (St > Sr) & (St > Sm)).sum()),
            "rt_decoy_wins": int((tgt & np.isfinite(Sr) & (Sr >= St)).sum()),
            "mz_decoy_wins": int((tgt & np.isfinite(Sm) & (Sm >= St)).sum())}
        # truth: human-only runs as acceptors (any YEAST/ECOLI-unique precursor transferred there is false)
        H = [r for r in runs if C.condition(r) == "H"]
        if H:
            sp = pep["species"].reindex(tr["pf_id"].to_numpy()).to_numpy()
            cont = pep["contaminant"].reindex(tr["pf_id"].to_numpy()).to_numpy().astype(bool)
            inH = tr["run"].isin(H).to_numpy() & tgt.to_numpy() & ~cont
            ye = np.isin(sp, ["YEAST", "ECOLI"])
            hu = sp == "HUMAN"
            truth = {"note": "acceptor = human-only runs; donors = all runs; YEAST/ECOLI-unique precursors are absent"}
            for tq in a.transfer_q:
                accq = (tr[a.qcol] <= tq).to_numpy()
                n_acc = int((inH & accq & (ye | hu)).sum())
                n_ye = int((inH & accq & ye).sum())
                n_ye_c = int((inH & ye).sum())
                n_hu_c = int((inH & hu).sum())
                # decoy-estimated false transfers inside the pure-null YE subset (exchangeability check)
                win = (np.isfinite(St) & (St > Sr) & (St > Sm)).to_numpy()
                thr = tr.loc[(tr[a.qcol] <= tq).to_numpy(), "St"].min() if accq.any() else np.inf
                d_rt = int((inH & ye & np.isfinite(Sr.to_numpy()) & (Sr.to_numpy() >= St.to_numpy())
                            & (Sr.to_numpy() >= thr)).sum())
                d_mz = int((inH & ye & np.isfinite(Sm.to_numpy()) & (Sm.to_numpy() >= St.to_numpy())
                            & (Sm.to_numpy() >= thr)).sum())
                truth[f"q{tq}"] = {
                    "accepted_in_H": n_acc, "accepted_YE_in_H (false)": n_ye,
                    "YE_candidates_in_H": n_ye_c, "HUMAN_candidates_in_H": n_hu_c,
                    "false_fraction_lower_bound": n_ye / max(n_acc, 1),
                    "YE_null_acceptance_rate": n_ye / max(n_ye_c, 1),
                    "HUMAN_acceptance_rate": int((inH & accq & hu).sum()) / max(n_hu_c, 1),
                    "decoy_wins_in_YE_null_subset_rt_mz": [d_rt, d_mz],
                    "YE_target_wins_above_thr": int((inH & ye & win & (St.to_numpy() >= thr)).sum()),
                    # entrapment-style estimate assuming false transfers spread over candidates in proportion
                    "combined_est_FDP_in_H": (n_ye * (1 + n_hu_c / max(n_ye_c, 1))) / max(n_acc, 1),
                }
            # Lim-style: YE precursors identified (MS2) in A/B runs that get transferred into >= 1 H run
            hye_ye = tall[tall["run"].isin(hye)]
            yek = set(zip(hye_ye["pf_id"], hye_ye["charge"]))
            yek = {k for k in yek if pep["species"].get(k[0]) in ("YEAST", "ECOLI") and not pep["contaminant"].get(k[0])}
            for tq in a.transfer_q:
                x = tr[inH & (tr[a.qcol] <= tq).to_numpy()]
                tk = set(zip(x["pf_id"], x["charge"]))
                truth[f"q{tq}"]["frac_YE_precursors_transferred_into_any_H"] = len(yek & tk) / max(len(yek), 1)
            # MS2-level false IDs in H runs (identification, not MBR)
            hi = tall[tall["run"].isin(H)]
            hsp = pep["species"].reindex(hi["pf_id"].to_numpy()).to_numpy()
            hct = pep["contaminant"].reindex(hi["pf_id"].to_numpy()).to_numpy().astype(bool)
            truth["ms2_ids_in_H"] = {"precursors": int(len(hi)), "YE_unique": int((np.isin(hsp, ["YEAST", "ECOLI"]) & ~hct).sum()),
                                     "HUMAN": int(((hsp == "HUMAN") & ~hct).sum())}
            out["mbr_truth_human_only"] = truth
    out["hye"] = res
    T["total_s"] = time.perf_counter() - t_start
    out["timing_s"] = T
    json.dump(out, open(os.path.join(a.out, "lfq_metrics.json"), "w"), indent=1, default=float)
    print(json.dumps({"ids": out["ids"], "timing": T}, default=float))


if __name__ == "__main__":
    main()
