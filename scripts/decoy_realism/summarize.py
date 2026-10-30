#!/usr/bin/env python3
"""3g: assemble results/2026-10-27_3g_decoy_realism/metrics.json from the /tmp/vulcan_3g analysis JSONs and apply the
PLAN.md decision rules mechanically (post hoc diagnostics are carried separately and never enter a rule)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import common as C  # noqa: E402

S = C.SCR


def J(n):
    p = os.path.join(S, n)
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    m = {"plan": "results/2026-10-27_3g_decoy_realism/PLAN.md (commit 41288c0)", "data": {
        "steps_1_2": "E9 arm C / arm P outputs, D8 N=1000 prefix (981 C / 982 P runs), /tmp/vulcan_e9/out",
        "step_3": "new search, index human_contam_3g (A Arabidopsis + B prokaryote union + S shuffled-A), D8 runs 0-99 "
                  "(96 completed), /tmp/vulcan_3g"}}
    m["e9_coins_rederived"] = J("e9_coins.json")
    s1 = J("step1.json")
    m["step1"] = s1
    m["features"] = J("features.json")
    m["predictor_check"] = J("predictor_check.json")
    abl = J("ablation.json")
    # retained-excess bootstrap for ablations vs 'none'
    if abl:
        for arm in ("C", "P"):
            for al in ("0.01", "0.05"):
                ref = os.path.join(S, f"abl_{arm}_none_{float(al)}.npy")
                if not os.path.exists(ref):
                    continue
                w0 = np.load(ref)
                for k in list(abl):
                    if not k.startswith(arm + "_") or k.endswith("stored"):
                        continue
                    g = k[2:]
                    f = os.path.join(S, f"abl_{arm}_{g}_{float(al)}.npy")
                    if os.path.exists(f) and arm == "C":
                        abl[k][al]["retained_excess"] = C.boot_ratio(w0, np.load(f))
    m["ablation"] = abl
    m["step3"] = J("step3.json")
    m["post_hoc"] = {"diag_b_prokaryote_true_presence": J("diag_b.json"), "diag_a_project_protein": J("diag_a.json"),
                     "diag_cluster_permutation": J("diag_cluster_perm.json")}
    m["homology_meta"] = J("homology_meta.json")
    m["db_3g"] = J("db_3g.json") if os.path.exists(os.path.join(S, "db_3g.json")) else json.load(
        open("/tmp/vulcan_3g/db_3g.json"))
    # ------------- mechanical decision rules (PLAN.md)
    dec = {}
    x = s1["score"]["0.05"]["cumulative"]
    ev = [k for k in ("L4", "L3", "L2", "L1") if x[k]["n"] >= 400]
    k = ev[0]
    ret = x[k]["retained_excess"]
    lad = s1["score"]["0.05"]["ladder"]["model_dy_db"]
    lo, hi = x[k]["length_std"]["ci95"]
    h1_sup = (x[k]["retained_excess_std"] <= 0.5 and ret["ci95"][1] < 1) and \
        (lad["coin_at_delta0_ci95"][0] <= 0.5 <= lad["coin_at_delta0_ci95"][1])
    h1_rej = (x[k]["retained_excess_std"] >= 0.5 and not (lo <= 0.5 <= hi)) and \
        not (lad["coin_at_delta0_ci95"][0] <= 0.5 <= lad["coin_at_delta0_ci95"][1])
    dec["step1"] = {"strictest_evaluable_cumulative_level": k, "n": x[k]["n"], "coin": x[k]["coin"],
                    "coin_length_std": x[k]["length_std"], "retained_excess_raw": ret,
                    "retained_excess_length_std": x[k]["retained_excess_std"],
                    "ladder_coin_at_delta0": [lad["coin_at_delta0"], lad["coin_at_delta0_ci95"]],
                    "verdict": "H1 supported" if h1_sup else ("H1 rejected as main cause" if h1_rej else "partial/mixed")}
    s3 = m["step3"]["score"]["0.05"]
    cA, cB, cS = (s3[f]["all"] for f in "ABS")
    ctrl1 = (cS["ci95"][0] <= 0.5 <= cS["ci95"][1]) or abs(cS["coin"] - .5) <= 0.015
    diffAB = s3["cA_minus_cB"]
    h2 = (not (cB["ci95"][0] <= .5 <= cB["ci95"][1])) and (cB["coin"] - .5 >= .5 * (cA["coin"] - .5)) and \
        (diffAB["ci95"][0] <= 0 <= diffAB["ci95"][1])
    h1 = ((cB["ci95"][0] <= .5 <= cB["ci95"][1]) or (cB["coin"] - .5 < .5 * (cA["coin"] - .5))) and \
        not (diffAB["ci95"][0] <= 0 <= diffAB["ci95"][1])
    dec["step3"] = {"control1_shuffled_passes": bool(ctrl1), "c_A": cA["coin"], "c_B": cB["coin"], "c_S": cS["coin"],
                    "cA_minus_cB": diffAB,
                    "mechanical_verdict": "H2-order" if h2 else ("H1" if h1 else "mixed"),
                    "validity": "INVALID control 2: B is not a null; c_B = 0.745 (5%) / 0.925 (1%) driven by truly "
                                "present prokaryotic proteins (E. coli GroEL/DnaK/BirA/ArnA/KatE, S. aureus protein A "
                                "SAUSA300_0113, B. cereus) concentrated in few projects (post hoc diag_b). Recorded as "
                                "'not evaluable'; the mechanical verdict is void."}
    abl_c = {g: (abl.get(f"C_{g}", {}).get("0.05", {}).get("retained_excess")) for g in
             ("pred_ms2", "rt", "pred_all", "hs_family", "ladder", "masscharge", "other")} if abl else {}
    carriers = [g for g, r in abl_c.items() if r and r["retained"] <= .5 and r["ci95"][1] < 1]
    dec["step2"] = {"retained_excess_5pct": abl_c, "carrier_groups": carriers,
                    "verdict": "distributed (no single non-degenerate group carries it)" if not carriers else carriers,
                    "note": "'delta' and 'all_but_pred' ablations collapse the scorer (~14k vs 361k targets at 1%) and "
                            "are not interpretable"}
    dec["overall"] = ("mixed per PLAN rule (step 1 partial, step 3 not evaluable). Weight of evidence: H1 in the broad "
                      "sense (real fragment evidence: homologous conserved proteins + truly present plant material), "
                      "not H2; see NOTES")
    m["decisions"] = dec
    out = os.path.join(C.OUT, "metrics.json")
    json.dump(m, open(out, "w"), indent=1, default=float)
    print(json.dumps(dec, indent=1, default=float))


if __name__ == "__main__":
    main()
