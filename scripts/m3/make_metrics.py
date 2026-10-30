#!/usr/bin/env python3
"""Assemble results/2026-10-27_M3_scoring/metrics.json from eval.json, gates.json, models/training_log.json,
predictor_library.json, bench_scoring.json, per-dir feat_info.json and the baselines (M1e metrics.json when present;
E2 ids.json for Sage's own LDA peptide_q)."""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

R = C.RES


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    ev = load(os.path.join(R, "eval.json"))
    out = {"experiment": "M3a / E4 rescoring (dev only: D0 PXD001468 x2, D1 PXD001077)", "date": "2026-10-27",
           "units": "peptidoforms (paired double competition units, as E2); ids_unique_sequences = distinct I/L-collapsed "
                    "sequences among them", "scorers": {}}
    for ds, d in ev.items():
        for sc, r in d.items():
            e = {a: {"overall": v["overall"], "P": {k: v["P"][k] for k in ("mean_fdp", "ci95", "mean_lb", "lb_ci95", "verdict",
                                                                         "consistent_with_control", "mean_n_t", "fdp_per_replicate")},
                     "C": {k: v["C"][k] for k in ("mean_fdp", "ci95", "mean_lb", "lb_ci95", "verdict", "r", "mean_n_t")}}
                 for a, v in r["entrapment"].items()}
            pe = {a: {k: v[k] for k in ("mean_fdp", "ci95", "mean_lb", "verdict", "mean_n_t")}
                  for a, v in r.get("protein_entrapment", {}).items()}
            out["scorers"].setdefault(sc, {})[ds] = {"model_class": r["model_class"], "ids": r["ids"],
                                                     "ids_unique_sequences": r["ids_unique_sequences"],
                                                     "proteins_picked_group": r["proteins"], "entrapment_peptide": e,
                                                     "entrapment_protein_P": pe, "reset": r.get("reset")}
    # baselines
    m1e = load(os.path.join(C.ROOT, "results", "2026-10-27_M1e_baselines", "metrics.json"))
    e2 = load(os.path.join(C.ROOT, "results", "2026-10-27_E2_gpu_search", "ids.json"))
    base = {"source": "results/2026-10-27_M1e_baselines/metrics.json (baselines agent; read-only)", "datasets": {}}
    for ds in ("D0", "D1"):
        b = {}
        if e2 and ds in e2:
            b["sage_own_lda_peptide_q_peptidoforms_at_1pct"] = e2[ds]["sage_own_q"]["peptides_at_1pct"]
        pipes = (m1e or {}).get("datasets", {}).get(ds, {}).get("pipelines", {})
        for p in ("sage", "sage+mokapot", "comet+mokapot"):
            if p in pipes:
                v = pipes[p]
                b[p] = {"ids_production": v["ids_production"], "proteins": v["protein_ids_production_picked_group"],
                        "verdicts_peptide": v["verdicts_peptide"], "verdicts_protein": v["verdicts_protein"]}
        b["status"] = "available" if pipes else "pending (baselines agent has not produced this dataset yet)"
        base["datasets"][ds] = b
    out["baselines"] = base
    out["gates"] = load(os.path.join(R, "gates.json"))
    out["training"] = load(os.path.join(R, "models", "training_log.json"))
    out["predictor_library"] = {k: v for k, v in (load(os.path.join(R, "predictor_library.json")) or {}).items() if k != "hw"}
    bench = load(os.path.join(R, "bench_scoring.json"))
    out["runtime"] = bench
    fi = {}
    for f in sorted(glob.glob(os.path.join(C.WORK, "psms", "*", "feat_info.json"))):
        j = json.load(open(f))
        fi[os.path.basename(os.path.dirname(f))] = {r: {k: v[k] for k in ("nce", "nce_grid_median_sa", "rt_mad_min",
                                                                           "rt_calib_n", "cache_stats", "timing_s", "n_psms")}
                                                    for r, v in j["runs"].items()}
    out["feature_runs"] = fi
    json.dump(out, open(os.path.join(R, "metrics.json"), "w"), indent=1, default=float)
    print("wrote metrics.json")


if __name__ == "__main__":
    main()
