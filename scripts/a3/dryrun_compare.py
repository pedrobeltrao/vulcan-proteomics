#!/usr/bin/env python3
"""A3.1.6 dry-run comparison (dev D0 / D1; nothing sealed).  Evidence that the A3 additions leave v1 untouched.

Two driver runs per dataset, both from clean git worktrees, sharing the db stage (the same replicate indexes / FASTAs,
built by the unchanged bench/t1/dbs.py):
  "without"  the driver at prereg-v0.4.2 (the A1/A2-tagged T1 driver, no E7/E7b/A3 code), A1 method set; stages
             ours (fresh cohort ingest, searches, features with a fresh prediction cache, frozen v1) + evaluate.
  "with"     the driver at the A3 code commit with --a3 (v2_sep + C2 + sage_chim2 added, both numerators,
             a3_secondary); every stage fresh (ours, chim, all comparators incl. a fresh TRFP mzML).
  The "without" evaluate reads the "with" run's comparator outputs (cmp/ symlinked), so every common number can be
  compared bit for bit; Sage's run-to-run non-reproducibility is measured separately against independent earlier runs.

Checks (exact float equality; NaN == NaN):
  (1) v1 artefacts: every <rep>/<run>.psms.parquet and .scored.parquet of `ours` byte-identical (sha256) between runs.
  (2) v1 columns of metrics.json identical: methods.vulcan, noninferiority (numerator vulcan) restricted to the A1
      comparators, matched_fdp_curves.vulcan, sage_predict_rt_control.  Also every A1 comparator record (same inputs,
      evaluate code with and without A3).
  (3) against independent earlier runs (other commits, other Sage invocations): v1 (methods.vulcan) and Comet vs the
      3f dry run (results/2026-10-27_T1driver_dryrun); Sage-based comparators: max |delta| of arm-P / arm-C mean FDP
      (documented tolerance 1e-4) and of the v1 / comparator ratios (+-0.0004), and the Sage-chimera comparators and
      v2_sep vs E7b (results/2026-10-28_E7b_fair_speed/t1).

    python scripts/a3/dryrun_compare.py --dataset D1 --with-work W --with-out O --without-work W2 --without-out O2 --out J
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
A1_CMP = ["sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot"]
SAGE_BASED = ["sage", "sage+mokapot", "sage+pred+mokapot"]
TOL_ARM_P = 1e-4
TOL_RATIO = 4e-4


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def diff(a, b, path="", out=None, numeric=None):
    """Exact recursive comparison -> list of differing paths; numeric collects |delta| of numeric leaves."""
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{path}/{k} (missing in {'with' if k not in a else 'without'})")
            else:
                diff(a[k], b[k], f"{path}/{k}", out, numeric)
    elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            out.append(f"{path} (len {len(a)} vs {len(b)})")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                diff(x, y, f"{path}[{i}]", out, numeric)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        same = (a == b) or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b))
        if not same:
            out.append(path)
        if numeric is not None and not (isinstance(a, float) and math.isnan(a)) and \
                not (isinstance(b, float) and math.isnan(b)):
            numeric.append((path, abs(float(a) - float(b))))
    elif a != b:
        out.append(path)
    return out


def ni_subset(m, comps, key="noninferiority"):
    return {c: m.get(key, {}).get(c) for c in comps}


def classify(paths):
    """diff paths -> value differences, keys missing in 'with', keys added in 'with' (A3 additions)."""
    added = [p for p in paths if p.endswith("(missing in without)")]
    missing = [p for p in paths if p.endswith("(missing in with)")]
    return {"value_diffs": [p for p in paths if p not in added and p not in missing], "missing_in_with": missing,
            "added_in_with": added}


def arm_means(rec, arms=("P", "C")):
    """{level/alpha/arm/conv: mean_fdp} of a method record."""
    out = {}
    for lv, byal in rec.get("verdicts", {}).items():
        for al, v in byal.items():
            for arm in arms:
                for cv, x in v.get(arm, {}).items():
                    out[f"{lv}/{al}/{arm}/{cv}"] = x["mean_fdp"]
    return out


def max_delta(d1, d2):
    ks = sorted(set(d1) & set(d2))
    if not ks:
        return None, None
    k = max(ks, key=lambda x: abs(d1[x] - d2[x]))
    return abs(d1[k] - d2[k]), k


def ratios(ni_m):
    out = {}
    for lv, d in (ni_m or {}).items():
        for kind in ("own_q", "matched_fdp"):
            for t, v in d.get(kind, {}).items():
                if v.get("ratio") is not None:
                    out[f"{lv}/{kind}/{t}"] = v["ratio"]
    return out


def ids_prod(rec):
    return {f"{lv}/{al}": v["mean"] for lv, d in rec.get("ids_dataset_level", {}).items() for al, v in d.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--with-work", required=True)
    ap.add_argument("--with-out", required=True)
    ap.add_argument("--without-work", required=True)
    ap.add_argument("--without-out", required=True)
    ap.add_argument("--a1-out", default=None, help="metrics of the SAME commit without --a3 (A1 method set)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    W = json.load(open(os.path.join(a.with_out, "metrics.json")))
    O = json.load(open(os.path.join(a.without_out, "metrics.json")))
    rep = {"dataset": a.dataset, "with": {"git_sha": W["meta"]["git_sha"], "dirty": W["meta"]["git_dirty_tracked"],
                                         "methods": list(W["methods"]), "a3": W["meta"].get("a3")},
           "without": {"git_sha": O["meta"]["git_sha"], "dirty": O["meta"]["git_dirty_tracked"],
                       "methods": list(O["methods"])}}
    # (1) v1 artefacts
    files = sorted(glob.glob(os.path.join(a.without_work, "ours", "*", "*.parquet")))
    bad, n = [], 0
    for f in files:
        rel = os.path.relpath(f, os.path.join(a.without_work, "ours"))
        g = os.path.join(a.with_work, "ours", rel)
        n += 1
        if not os.path.exists(g) or sha(f) != sha(g):
            bad.append(rel)
    rep["v1_artefacts"] = {"n_files": n, "n_not_identical": len(bad), "not_identical": bad[:50],
                           "identical": n > 0 and not bad}
    # (2) v1 columns of metrics.json
    chk = {
        "methods.vulcan": diff(W["methods"]["vulcan"], O["methods"]["vulcan"]),
        "noninferiority[vulcan vs A1 comparators]": diff(ni_subset(W, A1_CMP), ni_subset(O, A1_CMP)),
        "matched_fdp_curves.vulcan": diff(W["matched_fdp_curves"]["vulcan"], O["matched_fdp_curves"]["vulcan"]),
        "sage_predict_rt_control": diff(W.get("sage_predict_rt_control"), O.get("sage_predict_rt_control")),
        "A1 comparator records": diff({c: W["methods"][c] for c in A1_CMP}, {c: O["methods"][c] for c in A1_CMP}),
        "matched_fdp_curves[A1 comparators]": diff({c: W["matched_fdp_curves"][c] for c in A1_CMP},
                                                   {c: O["matched_fdp_curves"][c] for c in A1_CMP}),
        "noninferiority_numerator": [] if W.get("noninferiority_numerator") == "vulcan" else ["numerator"],
    }
    rep["v1_columns"] = {}
    for k, v in chk.items():
        c = classify(v)
        rep["v1_columns"][k] = {"identical_on_common_keys": not c["value_diffs"] and not c["missing_in_with"],
                                "n_value_diffs": len(c["value_diffs"]), "value_diffs": c["value_diffs"][:20],
                                "missing_in_with": c["missing_in_with"][:20], "n_keys_added_in_with": len(c["added_in_with"]),
                                "keys_added_in_with_examples": sorted({x.rsplit("/", 1)[-1] for x in c["added_in_with"]})[:10]}
    rep["v1_columns_bitwise_identical"] = all(x["identical_on_common_keys"] for x in rep["v1_columns"].values())
    # comparator artefacts recomputed by the "without" run (prediction features + every mokapot run) from the SAME
    # Sage / Comet searches: byte-identical => the A3 additions do not reach the A1 comparators' inputs
    cmpd = {}
    for d in ("sage_pred_mokapot", "sage_mokapot", "comet_mokapot"):
        root = os.path.join(a.without_work, "cmp", d)
        if not os.path.isdir(root) or os.path.realpath(root) == os.path.realpath(os.path.join(a.with_work, "cmp", d)):
            continue                     # not recomputed by the "without" run (shared directory): no evidence
        fs = sorted(f for f in glob.glob(os.path.join(root, "*", "*")) if os.path.isfile(f)
                    and not f.endswith((".log", ".log.time", ".json")))     # logs / timings / run info
        nb = [os.path.relpath(f, root) for f in fs
              if sha(f) != sha(os.path.join(a.with_work, "cmp", d, os.path.relpath(f, root)))]
        cmpd[d] = {"n_files": len(fs), "n_not_identical": len(nb), "not_identical": nb[:20]}
    rep["a1_comparator_artefacts_recomputed"] = cmpd
    rep["a1_comparator_artefacts_identical"] = bool(cmpd) and all(v["n_files"] > 0 and v["n_not_identical"] == 0
                                                                   for v in cmpd.values())
    if a.a1_out:
        A = json.load(open(os.path.join(a.a1_out, "metrics.json")))
        same = {"methods": diff({m: W["methods"][m] for m in A["methods"]}, A["methods"]),
                "noninferiority": diff({c: W["noninferiority"][c] for c in A["noninferiority"]}, A["noninferiority"]),
                "matched_fdp_curves": diff({m: W["matched_fdp_curves"][m] for m in A["matched_fdp_curves"]},
                                           A["matched_fdp_curves"])}
        rep["same_commit_without_a3"] = {"git_sha": A["meta"]["git_sha"], "methods": list(A["methods"]),
                                         **{k: {"identical": not v, "diff_paths": v[:20]} for k, v in same.items()}}
        rep["same_commit_without_a3_identical"] = all(not v for v in same.values())
    rep["added_by_a3"] = {"methods": [m for m in W["methods"] if m not in O["methods"]],
                          "noninferiority_comparators": [c for c in W["noninferiority"] if c not in O["noninferiority"]],
                          "top_level_keys": [k for k in W if k not in O]}
    # (3) independent earlier runs
    ind = {}
    p3f = os.path.join(ROOT, "results", "2026-10-27_T1driver_dryrun", a.dataset, "metrics.json")
    if os.path.exists(p3f):
        F = json.load(open(p3f))
        r = {"source": os.path.relpath(p3f, ROOT), "git_sha": F["meta"]["git_sha"],
             "v1_methods.vulcan_identical": not diff(W["methods"]["vulcan"], F["methods"]["vulcan"]),
             "v1_methods.vulcan_diff_paths": diff(W["methods"]["vulcan"], F["methods"]["vulcan"])[:20],
             "comet+mokapot_identical": not diff(W["methods"]["comet+mokapot"], F["methods"]["comet+mokapot"]),
             "sage_based": {}}
        for c in SAGE_BASED:
            dm, km = max_delta(arm_means(W["methods"][c], ("P",)), arm_means(F["methods"][c], ("P",)))
            dc, kc = max_delta(arm_means(W["methods"][c], ("C",)), arm_means(F["methods"][c], ("C",)))
            gated = {k: v for k, v in ratios(W["noninferiority"].get(c)).items() if k.endswith(("/0.01", "/0.05"))}
            dg, kg = max_delta(gated, ratios(F["noninferiority"].get(c)))
            dr, kr = max_delta(ratios(W["noninferiority"].get(c)), ratios(F["noninferiority"].get(c)))
            di, ki = max_delta(ids_prod(W["methods"][c]), ids_prod(F["methods"][c]))
            byl = {}
            for lv in ("peptidoform", "sequence", "protein"):
                pw = {k: v for k, v in arm_means(W["methods"][c], ("P",)).items() if k.startswith(lv + "/")}
                rw = {k: v for k, v in ratios(W["noninferiority"].get(c)).items() if k.startswith(lv + "/")}
                byl[lv] = {"max_abs_delta_arm_P_mean_fdp": max_delta(pw, arm_means(F["methods"][c], ("P",)))[0],
                           "max_abs_delta_v1_ratio": max_delta(rw, ratios(F["noninferiority"].get(c)))[0]}
            r["sage_based"][c] = {"by_level": byl,"max_abs_delta_arm_P_mean_fdp": dm, "at": km, "within_1e-4": dm is not None and dm <= TOL_ARM_P,
                                  "max_abs_delta_arm_C_mean_fdp (reported)": dc, "arm_C_at": kc,
                                  "max_abs_delta_v1_ratio_at_1%_5%": dg, "gated_ratio_at": kg,
                                  "max_abs_delta_v1_ratio": dr, "ratio_at": kr,
                                  "within_0.0004": dr is not None and dr <= TOL_RATIO,
                                  "max_abs_delta_prod_ids_mean": di, "ids_at": ki,
                                  "deciding_verdicts_identical": not diff(
                                      {lv: {al: v["deciding"] for al, v in d.items()} for lv, d in W["methods"][c]["verdicts"].items()},
                                      {lv: {al: v["deciding"] for al, v in d.items()} for lv, d in F["methods"][c]["verdicts"].items()})}
        ind["3f_dry_run"] = r
    for fn in ("metrics_nofm0001.json", "metrics.json"):
        pe = os.path.join(ROOT, "results", "2026-10-28_E7b_fair_speed", "t1", a.dataset, fn)
        if os.path.exists(pe):
            break
    if os.path.exists(pe):
        E = json.load(open(pe))
        r = {"source": os.path.relpath(pe, ROOT), "git_sha": E["meta"]["git_sha"], "methods": {}}
        for c in ("sage_chim2", "sage_chim2+pred+mokapot", "sage+pred+mokapot"):
            if c in E["methods"] and c in W["methods"]:
                dm, km = max_delta(arm_means(W["methods"][c], ("P",)), arm_means(E["methods"][c], ("P",)))
                di, ki = max_delta(ids_prod(W["methods"][c]), ids_prod(E["methods"][c]))
                r["methods"][c] = {"max_abs_delta_arm_P_mean_fdp": dm, "at": km,
                                   "within_1e-4": dm is not None and dm <= TOL_ARM_P,
                                   "max_abs_delta_prod_ids_mean": di, "ids_at": ki}
        if "vulcan_v2_sep" in E["methods"]:
            dv = diff(W["methods"]["vulcan_v2_sep"], E["methods"]["vulcan_v2_sep"])
            r["vulcan_v2_sep_record_identical"] = not dv
            r["vulcan_v2_sep_diff_paths"] = dv[:20]
            wn = W.get("noninferiority_by_numerator", {}).get("vulcan_v2_sep", {})
            en = E["noninferiority"] if E.get("noninferiority_numerator") == "vulcan_v2_sep" else {}
            r["v2_sep_ratios"] = {}
            for c in ("sage_chim2+pred+mokapot", "sage_chim2", "sage+pred+mokapot"):
                if c in wn and c in en:
                    dr, kr = max_delta(ratios(wn[c]), ratios(en[c]))
                    r["v2_sep_ratios"][c] = {"max_abs_delta": dr, "at": kr, "within_0.0004": dr is not None and dr <= TOL_RATIO,
                                             "with_pf_matched_0.01": wn[c]["peptidoform"]["matched_fdp"]["0.01"]["ratio"],
                                             "e7b_pf_matched_0.01": en[c]["peptidoform"]["matched_fdp"]["0.01"]["ratio"]}
        ind["e7b"] = r
    rep["independent_runs"] = ind
    sb = ind.get("3f_dry_run", {}).get("sage_based", {})
    rep["sage_based_within_documented_tolerance"] = bool(sb) and all(v["within_1e-4"] and v["within_0.0004"]
                                                                     for v in sb.values())
    json.dump(rep, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: rep.get(k) for k in ("v1_columns_bitwise_identical", "a1_comparator_artefacts_identical",
                                               "same_commit_without_a3_identical",
                                               "sage_based_within_documented_tolerance")}),
          rep["v1_artefacts"]["identical"], rep["v1_artefacts"]["n_files"])


if __name__ == "__main__":
    main()
