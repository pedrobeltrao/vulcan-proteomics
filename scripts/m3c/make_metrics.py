#!/usr/bin/env python3
"""Assemble results/2026-10-27_3c_itcid/metrics.json from eval*.json, sage_eval.json, gates.json, bench.json,
check_da.json, koina_repeat.json, models/training_log.json, timings.jsonl and the feature/search stats."""
import glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

R = C.RES
CHOSEN = "itcid_linear_pred"


def J(n):
    p = os.path.join(R, n)
    return json.load(open(p)) if os.path.exists(p) else None


def summ(rec):
    out = {"ids_peptidoforms": rec["ids"], "ids_unique_seq_IL": rec["ids_unique_seq_IL"], "protein_groups": rec["proteins"],
           "verdicts": {}}
    for a in ("0.01", "0.05", "0.001"):
        e = rec["entrapment"][a]
        out["verdicts"][a] = {"P": {"mean_fdp": e["P"]["mean_fdp"], "ci95": e["P"]["ci95"], "verdict": e["P"]["verdict"],
                                    "consistent_with_control": e["P"].get("consistent_with_control")},
                              "C": {"mean_fdp": e["C"]["mean_fdp"], "ci95": e["C"]["ci95"], "verdict": e["C"]["verdict"]},
                              "overall": e["overall"],
                              "protein_P": {"mean_fdp": rec["protein_entrapment"][a]["mean_fdp"],
                                            "verdict": rec["protein_entrapment"][a]["verdict"]}}
    return out


def main():
    m = {"experiment": "3c IT-CID instrument-class static-first scorer (dev only)", "chosen_model": CHOSEN}
    tl = J("models/training_log.json")
    meta = json.load(open(os.path.join(R, "models", CHOSEN + ".pkl.json")))
    m["chosen_model_sha256"] = meta["sha256"]
    m["training"] = {"n_runs": tl["n_runs"], "runs": tl["runs"], "n_psms": tl["n_psms"],
                     "hyperscore_heldout": tl["hyperscore_psms_at_1pct"],
                     "variants": {k: {"hp": v["hyperparams"], "cv_heldout_best": max(v["cv"].values()),
                                      "sha256": v["sha256"]} for k, v in tl["variants"].items()}}
    ev = {}
    for f in ("eval.json", "eval_D4.json"):
        for ds, d in (J(f) or {}).items():
            ev[ds] = {sc: summ(rec) for sc, rec in d.items()}
    m["dev"] = ev
    se = J("sage_eval.json") or {}
    m["comparators"] = {ds: {p: {"peptides_IL_own_q": {a: r["ids"][a]["peptides"] for a in ("0.001", "0.01", "0.05")},
                                 "peptides_IL_vulcan_tdc": {a: r["ids"][a].get("peptides_vulcan_tdc") for a in
                                                            ("0.001", "0.01", "0.05")},
                                 "protein_groups": {a: r["proteins"][a] for a in ("0.001", "0.01", "0.05")},
                                 "verdicts": r["verdicts"]} for p, r in d.items()} for ds, d in se.items()}
    m["gates"] = J("gates.json")
    m["check_da"] = J("check_da.json")
    m["koina_repeat"] = J("koina_repeat.json")
    b = J("bench.json")
    m["runtime"] = {"bench_locked_warm_cache": b and {k: {"median": v["median"], "first_run_incl_cache_load":
                                                          v["runs"][0], "setup_s_unlocked": v["setup_s_unlocked"]}
                                                      for k, v in b["sets"].items()}}
    cold = {}
    for ds in ("D4", "D10cf"):
        fi = json.load(open(os.path.join(C.psm_dir(f"{ds}_prod"), "feat_info.json")))
        st = json.load(open(os.path.join(C.psm_dir(f"{ds}_prod"), "stats.json")))
        cold[ds] = {"n_runs": len(fi["runs"]),
                    "features_total_s_cold": sum(r["timing_s"]["total"] for r in fi["runs"].values()),
                    "of_which_koina_predict_s": sum(r["timing_s"]["of_which_predict"] for r in fi["runs"].values()),
                    "search_wall_s": sum(r["wall_s"] for r in st["runs"]), "search_core_s": sum(r["core_s"] for r in st["runs"])}
    m["runtime"]["first_contact_unlocked"] = cold
    tim = [json.loads(l) for l in open(os.path.join(R, "timings.jsonl"))] if os.path.exists(os.path.join(R, "timings.jsonl")) else []
    m["runtime"]["comparators_unlocked_16thr"] = {t["label"]: {"wall_s": t["wall_s"], "core_s": t.get("core_s"),
                                                               "rc": t["rc"]} for t in tim}
    json.dump(m, open(os.path.join(R, "metrics.json"), "w"), indent=1, default=float)
    print("wrote metrics.json")


if __name__ == "__main__":
    main()
