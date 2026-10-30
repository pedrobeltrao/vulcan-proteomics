#!/usr/bin/env python3
"""Collate results/2026-10-28_M7_encoders/metrics.json from e2b.json, e6.json, gates_ce.json, e5.json, bench_ce.json,
dataset_meta.json and the model training log; apply the pre-stated go/no-go rules of NOTES §2 mechanically."""
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
R = os.path.join(ROOT, "results", "2026-10-28_M7_encoders")
J = lambda f: json.load(open(os.path.join(R, f)))  # noqa: E731


def main():
    e2b, e6, g, e5, b = J("e2b.json"), J("e6.json"), J("gates_ce.json"), J("e5.json"), J("bench_ce.json")
    ds, tr = J("dataset_meta.json"), J("models/ce_s.json")
    m = {"scope": "M7 encoder arm, dev only (training-role runs, D0, D1, D8 holdout first 10 runs, D7 arm M, M3e G2 nulls); "
                  "no sealed data; no DN2 data", "hw": J("hwinfo.json")}
    # E2b
    rows = [r for r in e2b["rows"] if r["dtype"] == "bf16"]
    m["E2b"] = {"bf16": {f"{r['size']}_K{r['K']}": {"params": r["params"], "spectra_per_s": r["full_spectra_per_s"],
                                                    "candidates_per_s": r["full_candidates_per_s"]} for r in rows},
                "decision": {"size": "s", "params": 1659164, "dtype": "bf16", "K_production": 5,
                             "rule": "largest size with <= 1.0 GPU-s per median run at the production top-5 table"}}
    # training
    m["training"] = {"n_runs_selected": 43, "n_runs_used_with_positives": sum(1 for v in ds["runs"].values() if v.get("n_pos")),
                     "n_spectra": ds["n_spectra"], "n_candidates": ds["n_candidates"],
                     "n_by_kind": dict(zip(["positive", "hard_topK", "shifted_precursor", "wrong_species"], ds["n_by_kind"])),
                     "n_spectra_by_fold": ds["n_spectra_by_fold"], "val_fold": 4, "raw_sha256_overlap": ds["raw_sha256"],
                     "spectrum_hash_dropped_in_eval": ds["dropped_in_eval"], "spectrum_hash_dup_within_train": ds["dropped_dup_within_train"],
                     "model_sha256": tr["sha256"], "best_epoch": tr["best_epoch"], "val": tr["epochs"][tr["best_epoch"]]["val"],
                     "train_wall_s": tr["train_wall_s"]}
    # E6
    E = {}
    for ds_ in ("D0", "D1"):
        E[ds_] = {sc: {"ids": r["ids"], "proteins": r["proteins"], "matched_ids": {k: v[0] for k, v in r["matched"].items()},
                       "verdict_1pct": {a: r["entrapment"]["0.01"]["a"][a]["verdict"] for a in ("P", "C")},
                       "fdp_1pct": {a: r["entrapment"]["0.01"]["a"][a]["mean_fdp"] for a in ("P", "C")},
                       "verdict_5pct": {a: r["entrapment"]["0.05"]["a"][a]["verdict"] for a in ("P", "C")}}
                  for sc, r in e6[ds_]["scorers"].items()}
        E[ds_]["ratio_ce_vs"] = e6[ds_]["ratio_ce_vs"]
    E["D8"] = {sc: {"ids": r["ids"], "matched_ids": {k: v[0] for k, v in r["matched"].items()},
                    "verdict_1pct": {a: r["entrapment"]["0.01"]["a"][a]["verdict"] for a in ("P", "C")},
                    "fdp_1pct": {a: r["entrapment"]["0.01"]["a"][a]["mean_fdp"] for a in ("P", "C")},
                    "verdict_5pct": {a: r["entrapment"]["0.05"]["a"][a]["verdict"] for a in ("P", "C")}}
               for sc, r in e6["D8"]["pooled"].items()}
    E["D8"]["ratio_ce_vs"] = e6["D8"]["ratio_ce_vs"]
    E["armM"] = {sc: {a: {k: r["entrapment"][a]["a"][k] for k in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95")}
                      for a in ("0.01", "0.05")} | {"memorization_1pct": r["memorization"]}
                 for sc, r in e6["D7_armM"]["scorers"].items()}
    m["E6"] = E
    gl = g["G2_G5"]["learned"]
    m["gates_ce"] = {"G0": g["G0"],
                     "G2_original": gl["ce"]["G2"], "G5": gl["ce"]["G5"],
                     "frac_null_winner_changed_vs_hyperscore": gl["ce"]["frac_winner_changed_vs_hyperscore"],
                     "pending": ["tail-G2", "G1", "G3"],
                     "DIAG_on_D12_nulls_not_verdicts": {"ce": {"tail_G2": gl["ce"]["DIAG_tail_G2_1pct_D12_nulls"],
                                                               "G1": gl["ce"]["DIAG_G1_D12_nulls"]},
                                                        "static_gbm_pred": {"tail_G2": gl["gbm"]["DIAG_tail_G2_1pct_D12_nulls"],
                                                                            "G1": gl["gbm"]["DIAG_G1_D12_nulls"]}},
                     "static_gbm_pred_reference": {"G2": gl["gbm"]["G2"]["pass"], "G5_auc": gl["gbm"]["G5"]["all"]["auc"]}}
    m["runtime_ce"] = {k: b[k] for k in ("median_gpu_s", "median_wall_s", "median_core_s", "candidates_per_gpu_s",
                                         "spectra_per_gpu_s", "one_time", "loadavg")}
    m["E5"] = {"open_recall": {k: e5["open_recall"][k] for k in ("recall_all", "recall_narrow", "recall_shifted", "median_rank",
                                                                 "n_reference", "n_shifted", "peptide_embedding_s_two_banks",
                                                                 "bruteforce_query_s_per_spectrum")},
               "k1_open_s_per_spectrum": sum(v["wall_s"]["total"] for v in e5["open_recall"]["k1_open_stats"].values())
               / sum(v["n_spectra"] for v in e5["open_recall"]["k1_open_stats"].values()),
               "g3_style": {K: {"pass": v["pass"], "all_nulls_pass": v["all_nulls"]["pass"],
                                "clean_species_pass": v["clean_species_nulls"]["pass"]} for K, v in e5["g3_symmetry"]["by_K"].items()},
               "DIAG_open_window_null_target_fraction": e5.get("DIAG_open_window_null_symmetry", {}).get("by_K")}
    # pre-stated rules
    r0, r1, r8 = E["D0"]["ratio_ce_vs"]["static_gbm_pred"]["P_0.01"], E["D1"]["ratio_ce_vs"]["static_gbm_pred"]["P_0.01"], \
        E["D8"]["ratio_ce_vs"]["static_gbm_pred"]["P_0.01"]
    rr = E["D8"]["ratio_ce_vs"]["reset_lda_pred"]["P_0.01"]["ratio"]
    c1 = r0 >= 1.0 and r1 >= 1.0 and r8["ratio"] >= 1.02 and r8["ci95"][0] > 1.0 and rr >= 1.0
    verdicts = [E[d]["ce_s"][f"verdict_{a}"][arm] for d in ("D0", "D1", "D8") for a in ("1pct", "5pct") for arm in ("P", "C")]
    verdicts += [E["armM"]["ce_s"][a]["verdict"] for a in ("0.01", "0.05")]
    c2 = "violated" not in verdicts and E["armM"]["ce_s"]["0.01"]["mean_fdp"] <= E["armM"]["static_gbm_pred"]["0.01"]["mean_fdp"] + 0.0025
    c3 = g["G0"]["pass"] and gl["ce"]["G2"]["pass"] and not gl["ce"]["G5"]["all"]["alarm"]
    c4 = b["median_gpu_s"] <= 1.0
    e5r = e5["open_recall"]
    g5a = max(e5r["recall_all"][k] for k in e5r["recall_all"]) >= 0.95 and max(e5r["recall_shifted"].values()) >= 0.95
    g5b = all(v["pass"] for v in e5["g3_symmetry"]["by_K"].values())
    m["go_no_go"] = {"E6": {"1_sensitivity": c1, "2_entrapment_incl_armM": c2, "3_gates_available": c3, "4_runtime": c4,
                            "verdict": "GO (conditional: tail-G2/G1/G3 pending)" if all((c1, c2, c3, c4)) else
                            ("NO-GO as default (validity OK; research arm)" if (c2 and c3) else "NO-GO (validity)")},
                     "E5": {"1_recall_ge_95pct": g5a, "2_g3_style_symmetry": g5b,
                            "verdict": "GO" if g5a and g5b else "NO-GO"}}
    json.dump(m, open(os.path.join(R, "metrics.json"), "w"), indent=1, default=float)
    print(json.dumps(m["go_no_go"], indent=1))


if __name__ == "__main__":
    main()
