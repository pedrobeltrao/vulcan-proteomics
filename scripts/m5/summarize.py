#!/usr/bin/env python3
"""Assemble results/2026-10-27_M5_xl/metrics.json from the M5 run outputs (run_m5.json, CSM tables -> L4 evaluate,
exhaustive_d3.json + its tables, comparators.json).

    python scripts/m5/summarize.py --out results/2026-10-27_M5_xl
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaluate as EV  # noqa: E402

WORK = "/tmp/vulcan_m5"
FASTA = {"D3b": f"{WORK}/fasta/D3b_ribo_human671.fasta", "D3": f"{WORK}/fasta/D3_cas9_contam.fasta"}

PUBLISHED = {
    "Beveridge2020_PXD014337_actual_FDR_at_5pct_nominal": {
        "pLink": 0.113, "StavroX": 0.052, "xi": 0.025, "Kojak(Percolator vs PeptideProphet)": [0.024, 0.32],
        "MeroX": 0.049, "XlinkX": 0.0,
        "source": "Beveridge et al. Nat Commun 11:742 (2020) as summarised in lit/04 A1; pooled over their "
                  "replicates/linkers, not exactly our single DSS R1 file; level as in the paper"},
    "Matzinger2022_PXD029252": {
        "statement": "all six engines (MS Annika, MaxLynx, MeroX, XlinkX, pLink 2, xiSEARCH) exceeded the 1% "
                     "estimated FDR; error higher for inter-protein links",
        "source": "Matzinger et al. Nat Commun 13:3975 (2022) via lit/04 A1 (per-file numbers not in our notes)"},
}


def headline(ev):
    out = {}
    for lvl, d in ev["levels"].items():
        for cl, c in d.items():
            for a in ("0.01", "0.05"):
                x = c[f"tddd_plus1@{a}"]
                y = c[f"td_plus1@{a}"]
                out[f"{lvl}.{cl}@{a}"] = {"tt": x["accepted_tt"], "true_fdp": x["true_fdp"],
                                          "reportable": x["reportable"],
                                          "conservative_td_plus1": {"tt": y["accepted_tt"], "true_fdp": y["true_fdp"]}}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run = json.load(open(f"{WORK}/out/run_m5.json"))
    M = {"experiment": "M5a XL (dev only: D3 = PXD014337 DSS R1, D3b = PXD029252 DSSO mainlib_rep2; T3 untouched)",
         "estimator": "FDR = (TD+1-DD)/TT clipped at 0 (vulcan.stats.xl plus_one), merge-then-estimate per level, "
                      "self vs heteromeric separately; conservative (TD+1)/TT alongside; a cell is reportable "
                      "only if accepted TT >= 50",
         "truth": "library groups (IMP-X-FDR tables); see scripts/m5/truth.py", "datasets": {},
         "hw": run["hw"], "process": run["process"], "published_comparators": PUBLISHED}
    for which, R in run["datasets"].items():
        csm = glob.glob(f"{WORK}/out/{which}/*.csm.parquet")[0]
        pm = EV.ProteinMap(f"{WORK}/index_timed/{which}", FASTA[which])
        inp, t = EV.ours_inputs(csm, f"{WORK}/index_timed/{which}", homo=False)
        ev = EV.evaluate(which, pm=pm, label="vulcan", **inp)
        ev["homo_family"] = EV.homo_summary(which, pm, *EV.ours_inputs(csm, f"{WORK}/index_timed/{which}", homo=True))
        json.dump(ev, open(os.path.join(a.out, f"eval_{which}.json"), "w"), indent=1)
        ps = R["search_passes"]
        M["datasets"][which] = {
            "linker": R["linker"], "headline": headline(ev), "homo_family": ev["homo_family"],
            "n_csm_input_hetero_family": ev["n_csm_input"], "n_tt_td_dd": [ev["n_tt"], ev["n_td"], ev["n_dd"]],
            "runtime": {"ingest": {k: v for k, v in R.get("ingest", {}).items() if k != "store_meta"},
                        "ingest_store_timing_s": R.get("ingest", {}).get("store_meta", {}).get("timing_s"),
                        "index_build": R["index_build"], "gpu_index_build": R["gpu_index_build"],
                        "load_store": R["load_store"],
                        "search_cold": {"wall_s": ps[0]["wall_s"], "core_s": ps[0]["core_s"]},
                        "search_warm": {"wall_s": ps[-1]["wall_s"], "core_s": ps[-1]["core_s"],
                                        "stage_wall_s": ps[-1]["stage_wall_s"], "stage_gpu_ms": ps[-1]["stage_gpu_ms"]},
                        "passes_identical": all(x.get("identical_to_pass1", True) for x in ps),
                        "sites_cpu": R["sites"], "gpu_max_mem_bytes": R["gpu_max_mem_bytes"]},
            "counts": ps[-1]["counts"], "n_rank1": R["n_rank1"], "n_rank1_xl_wins": R["n_rank1_xl_wins"],
            "competitor_kind_of_spectra_lost_to_linear(0)/mono(1)/loop(2)": R["competitor_kind_of_lost_spectra"]}
    ex = f"{WORK}/exh_final/exhaustive_d3.json"
    if os.path.exists(ex):
        E = json.load(open(ex))
        E.pop("hw", None)
        fdr = {}
        for f in sorted(glob.glob(f"{WORK}/exh_final/*.csm.parquet")):
            name = os.path.basename(f)[:-12]
            pm = EV.ProteinMap("/tmp/vulcan_m5/index/D3", FASTA["D3"])
            inp, _ = EV.ours_inputs(f, "/tmp/vulcan_m5/index/D3", homo=False)
            fdr[name] = headline(EV.evaluate("D3", pm=pm, **inp))
        E["fdr_by_search"] = {k: {kk: vv for kk, vv in v.items() if ".self@" in kk} for k, v in fdr.items()}
        M["D3_search_space_reduction"] = E
    cp = f"{WORK}/comparators.json"
    if os.path.exists(cp):
        M["comparators"] = json.load(open(cp))
    json.dump(M, open(os.path.join(a.out, "metrics.json"), "w"), indent=1, default=str)
    for w, d in M["datasets"].items():
        for k, v in d["headline"].items():
            if v["tt"]:
                print(w, k, v["tt"], v["true_fdp"], "R" if v["reportable"] else "-")


if __name__ == "__main__":
    main()
