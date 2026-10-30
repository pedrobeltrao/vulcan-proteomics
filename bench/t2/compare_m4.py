#!/usr/bin/env python3
"""Dev reproduction check: a `t2_driver.py --dataset D2` run vs M4's committed E8 (results/2026-10-27_M4_ptm/e8.json).

    python bench/t2/compare_m4.py <dev out>/metrics.json [--work /tmp/vulcan_look_t2_<id>] [--m4-work /tmp/vulcan_m4]

Deterministic stages must be identical: search FASTA bytes, PSM tables (if both work dirs exist), peptidoform FDR counts,
truth (via eligible counts), pyAscore evaluations (A and L).  Our localizer: M4's run used the non-deterministic CUDA
index_add_ (IMPLEMENTATION_NOTES §3), so its numbers are compared with the driver's deterministic ones field by field and
the differences are printed; the acceptance band is the measured M4-mode run-to-run spread (--spread JSON, optional).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def flat(d, pre=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{pre}.{k}" if pre else str(k)))
    else:
        out[pre] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metrics")
    ap.add_argument("--work", default=None)
    ap.add_argument("--m4-work", default="/tmp/vulcan_m4")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    m = json.load(open(a.metrics))
    ref = json.load(open(os.path.join(ROOT, "results/2026-10-27_M4_ptm/e8.json")))
    rep = {"identical": {}, "ours_diffs": {}}
    fa = os.path.join(a.m4_work, "d2_phospho.fasta")
    if os.path.exists(fa):
        h = hashlib.sha256(open(fa, "rb").read()).hexdigest()
        rep["identical"]["search_fasta_sha256"] = h == m["meta"]["search_fasta_sha256"]
    if a.work and os.path.isdir(os.path.join(a.m4_work, "psm_d2_phospho")):
        import pyarrow.parquet as pq
        same = {}
        for acc, d in m["datasets"].items():
            for r in d["runs"]:
                x = pq.read_table(os.path.join(a.work, "psms", f"{r}.psms.parquet"))
                y = pq.read_table(os.path.join(a.m4_work, "psm_d2_phospho", f"{r}.psms.parquet"))
                same[r] = x.equals(y)
        rep["identical"]["psm_tables"] = same
    for acc, d in m["datasets"].items():
        r = ref["datasets"][acc]
        for k in ("n_rank1_psms", "n_peptidoforms_q01", "n_phospho_psms_localized", "n_mono", "n_multi",
                  "n_truth_eligible_psms", "search_rank1_site_correct_frac_on_eligible"):
            rep["identical"][f"{acc}.{k}"] = d.get(k) == r.get(k)
        rep["identical"][f"{acc}.search_mode_counts"] = d.get("search_mode_counts") == r.get("search_mode_counts")
        for k in ("pyascore_A", "pyascore_L"):
            rep["identical"][f"{acc}.{k}.eval"] = d[k]["eval"] == r[k]["eval"]
        for k in ("ours_A", "ours_L"):
            a_, b_ = flat(d[k]["eval"]), flat(r[k]["eval"])
            diffs = {x: [b_.get(x), a_.get(x)] for x in sorted(set(a_) | set(b_)) if a_.get(x) != b_.get(x)}
            rep["ours_diffs"][f"{acc}.{k}"] = {"n_fields": len(a_), "n_diff": len(diffs), "m4_vs_driver": diffs}
            # PI D35: the frozen localizer is non-deterministic; M4's committed value must lie in the driver's run band
            band = d.get("t2_band", {})
            rm = r[k]["eval"].get("mono", {})
            inb = {}
            for al in ("0.01", "0.05"):
                for parts, bk in (((f"all_thr_{al}", "true_flr"), "true_flr"),
                                  ((f"all_thr_{al}", "n_sty_eligible"), "n_accepted_eligible"),
                                  ((f"recall_at_true_flr_{al}",), "recall_at_true_flr")):
                    m4k = "/".join(parts)
                    v = rm
                    for part in parts:
                        v = v.get(part) if isinstance(v, dict) else None
                    b = band.get(f"{k}/{al}/{bk}")
                    if v is None or not b or "min" not in b:
                        inb[m4k] = None
                    else:
                        inb[m4k] = {"m4": v, "band": [b["min"], b["max"]], "n_runs": b["n"], "inside": b["min"] <= v <= b["max"]}
            rep.setdefault("m4_inside_band", {})[f"{acc}.{k}"] = inb
    ins = [x["inside"] for v in rep.get("m4_inside_band", {}).values() for x in v.values() if x]
    rep["m4_inside_band_all"] = all(ins) if ins else None
    rep["m4_inside_band_count"] = [sum(ins), len(ins)]
    rep["all_deterministic_identical"] = all(
        (all(v.values()) if isinstance(v, dict) else v) for v in rep["identical"].values())
    s = json.dumps(rep, indent=1, default=float)
    if a.out:
        open(a.out, "w").write(s)
    print(s[:20000])


if __name__ == "__main__":
    main()
