#!/usr/bin/env python3
"""Assemble results/2026-10-27_M6d_quant_accuracy/metrics.json from the M6d work files (/tmp/vulcan_m6d)."""
import glob
import json
import os
import sys

OUT = sys.argv[1]
W = "/tmp/vulcan_m6d"
m = {"experiment": "M6d quant accuracy (dev only: D9, D11, D12)", "date": "2026-10-27",
     "diagnosis_d9": json.load(open(os.path.join(OUT, "diag_d9.json"))), "before_after": {}, "mbr_truth": {}}


def short(d):
    p, f = d["precursor"], d["peptidoform"]
    r = {"prec_n": p["n"], "prec_eps": p["eps"], "prec_eps_species": p["eps_sp"], "prec_lfc_species": p["lfc_sp"],
         "prec_cv": p.get("cv", p.get("cv2")), "prec_cv_kind": "min3" if "cv" in p else "min2",
         "prec_completeness": p["completeness"], "pf_n": f["n"], "pf_eps": f["eps"], "pf_cv": f.get("cv", f.get("cv2")),
         "provenance_counts": d["provenance_counts"], "protein": d["protein"], "timing_s": d["timing_s"]}
    if "vs_sage" in d:
        r["vs_sage"] = {k: {"n": v["n"], "ours_eps": v["ours"]["eps"], "sage_eps": v["sage"]["eps"],
                            "ours_eps_species": v["ours"]["eps_sp"], "sage_eps_species": v["sage"]["eps_sp"],
                            "ours_cv": v["ours"].get("cv", v["ours"].get("cv2")), "sage_cv": v["sage"].get("cv", v["sage"].get("cv2"))}
                        for k, v in d["vs_sage"].items()}
    return r


for f in sorted(glob.glob(os.path.join(W, "ev", "d*_*_m6*_o*.json"))):
    k = os.path.basename(f)[:-5]
    m["before_after"][k] = short(json.load(open(f)))
for f in sorted(glob.glob(os.path.join(W, "ev", "mbr_truth_*.json"))):
    m["mbr_truth"][os.path.basename(f)[10:-5]] = json.load(open(f))
for f in sorted(glob.glob(os.path.join(W, "ev", "same_set_*.json"))):
    m.setdefault("same_precursors_as_m6", {})[os.path.basename(f)[9:-5]] = json.load(open(f))
m["grid_d11_all"] = os.path.join("diag_tables", "grid_d11_all.tsv")
m["grid_d9_sel"] = os.path.join("diag_tables", "grid_d9_sel.tsv")
fe = {}
for n in ("d11", "d12"):
    p = os.path.join(W, n, "frontend_timing.json")
    if os.path.exists(p):
        fe[n] = json.load(open(p))
        s = os.path.join(W, n, "ids", "all", "ids_summary.json")
        fe[n]["ids_summary"] = json.load(open(s))
m["frontend_new_sets"] = fe
t = []
tp = os.path.join(OUT, "timings.jsonl")
if os.path.exists(tp):
    t = [json.loads(l) for l in open(tp)]
tm = {}
for k in ("m6", "m6d"):
    tm[k] = [json.load(open(f))["timing_s"] for f in sorted(glob.glob(os.path.join(W, "timing", f"{k}_*.json")))]
m["runtime"] = {"locked_records": t, "stage_timing_s": tm}
json.dump(m, open(os.path.join(OUT, "metrics.json"), "w"), indent=1, default=float)
print("ok", list(m["before_after"]))
