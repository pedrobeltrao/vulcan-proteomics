#!/usr/bin/env python3
"""A3.1.8(ii): arm M (DESIGN §6.3 (M), A1.5) for v2_sep on D7, chimeric search.  Dev only; no DN2, no sealed data.

Exactly the A1 arm-M design of scripts/m3e/armM.py (D7 yeast PXD031616 runs LM2486/2487/2490, one experiment; the
existing entrapment indexes yeast + contaminants + P37840 with human entrapment at r = 1, one per ENTRAP_SEEDS_DEV seed
2001-2003, $VULCAN_M3E_WORK/index/yeast_M<seed>, built by scripts/m3e/armM_prepare.py), with the v2_sep pipeline in
place of the closed search + v1:
  * MS1: vulcan.ms1 features + hypotheses of each run (E7 defaults) at the committed code -> $W/ms1/D7;
  * search: vulcan.chimeric.search.search_chimeric, E2 production SearchParams, report_psms 5, per seed index;
  * score: v2 features (vulcan.rescore.features_v2.FeatureBuilderV2, AlphaPeptDeep generic, PredCache per NCE, as
    bench/t1/ours_chimeric.ScorerV2) and the frozen v2 model (sha256 3f2c69cc...);
  * competition: vulcan.chimeric.compete.chimeric_competition over the 3 runs as one experiment (spectrum key =
    run index << 26 + spec_idx, as bench/t1/ours_chimeric.chim_units), q = q_sep (E7-A1, per merged stratum TDC);
  * arm M estimator: combined FDP_C = N_E (1 + 1/r) / (N_T + N_E) with r = 1.0 (peptide-level effective search space,
    as armM.py), units = peptidoform of the target-side pair winners; one Replicate per seed (decoy seed 1001);
    vulcan.stats.entrapment.evaluate_arm under CI conventions (a) (deciding, A1.6 / 3f) and (b) (reported);
  * lists: overall, and each raw chimeric stratum {1, 2, >=3} of the pair winner (counts restricted to the stratum,
    as scripts/e7/strata.py combined_stratum), alpha in {0.1%, 1%, 5%} (0.1% reported only).
Rule (A3.1.8(ii)): no "violated" verdict at 1% or 5%, overall or in any stratum.  Reported: the hyperscore pipeline
through the identical chimeric competition, the merge maps, stratum sizes.

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/a3/dev_armM.py --ms1 --search --score --eval --out <json>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_A3_WORK", "/tmp/vulcan_a3")
M3E = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
STORES = os.environ.get("VULCAN_M3_STORE_NULL", "/tmp/vulcan_m3/store_null")
RUNS = ["LM2486_DDA", "LM2487_DDA", "LM2490_DDA"]
SEEDS = (2001, 2002, 2003)
ALPHAS = (0.001, 0.01, 0.05)
GATED = (0.01, 0.05)
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
V2 = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2", "models", "static_v2.pkl")
V2_SHA = "3f2c69cc6cad002e42c51475dc1ad52b9bfbbc84f3dd51997b22cb24a70f0186"
MS1 = os.path.join(W, "ms1", "D7")


def idx(s):
    return os.path.join(M3E, "index", f"yeast_M{s}")


def pdir(s):
    return os.path.join(W, "armM", f"D7_M{s}")


def ms1():
    os.makedirs(MS1, exist_ok=True)
    for r in RUNS:
        if not os.path.exists(os.path.join(MS1, r + ".hypotheses.arrow")):
            subprocess.run([sys.executable, "-m", "vulcan.ms1", os.path.join(STORES, r + ".vspec.arrow"), MS1,
                            "--threads", "8"], check=True, cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT})
    return {r: hashlib.sha256(open(os.path.join(MS1, r + ".hypotheses.arrow"), "rb").read()).hexdigest() for r in RUNS}


def search():
    import pyarrow.parquet as pq
    import torch

    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(report_psms=5)
    sps = {r: E.load_store(os.path.join(STORES, r + ".vspec.arrow")) for r in RUNS}
    groups = {r: CS.precursor_groups(CS.read_hypotheses(os.path.join(MS1, r + ".hypotheses.arrow")), sps[r])
              for r in RUNS}
    for s in SEEDS:
        os.makedirs(pdir(s), exist_ok=True)
        if all(os.path.exists(os.path.join(pdir(s), r + ".chim.psms.parquet")) for r in RUNS):
            continue
        l1 = L1Index(idx(s))
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        for r in RUNS:
            out = os.path.join(pdir(s), r + ".chim.psms.parquet")
            if os.path.exists(out):
                continue
            tbl, st = CS.search_chimeric(sps[r], gi, l1, p, groups[r])
            pq.write_table(tbl, out)
            print("search", s, r, tbl.num_rows, st["n_secondary_groups"], flush=True)
        del gi
        torch.cuda.empty_cache()


def score():
    import pyarrow as pa
    import pyarrow.parquet as pq

    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD, check_allowed
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features_v2 import FeatureBuilderV2, model_psm, read_ms1_features
    from vulcan.rescore.static import StaticModel
    check_allowed("D7")
    assert hashlib.sha256(open(V2, "rb").read()).hexdigest() == V2_SHA
    v2 = StaticModel.load(V2)
    pred = APD()
    caches = {}

    def fac(nce):
        if nce not in caches:
            caches[nce] = PredCache(os.path.join(W, "predcache_D7"), pred, nce, "QE", STATIC, VAR)
        return caches[nce]
    fb = FeatureBuilderV2(fac, STATIC, VAR)
    sps = {r: load_store(os.path.join(STORES, r + ".vspec.arrow")) for r in RUNS}
    mf = {r: read_ms1_features(os.path.join(MS1, r + ".features.arrow")) for r in RUNS}
    for s in SEEDS:
        l1 = L1Index(idx(s))
        for r in RUNS:
            out = os.path.join(pdir(s), r + ".chim.scored.parquet")
            if os.path.exists(out):
                continue
            t = pq.read_table(os.path.join(pdir(s), r + ".chim.psms.parquet"))
            psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in t.column_names}
            feats, finfo = fb.run_v2(l1, psm, sps[r], mf[r])
            sv2 = v2.apply(model_psm(psm), feats).astype(np.float64)
            pq.write_table(pa.table({"spec_idx": psm["spec_idx"], "group": psm["group"], "chim_rank": psm["chim_rank"],
                                     "rank": psm["rank"], "pf_row": psm["pf_row"], "pf_id": psm["pf_id"],
                                     "charge": psm["charge"], "hyperscore": psm["hyperscore"], "score_v2": sv2}), out)
            json.dump({"nce": finfo["nce"]}, open(out + ".json", "w"), default=float)
            print("score", s, r, len(sv2), finfo["nce"], flush=True)


def evaluate(out_path, ms1_sha=None):
    import m3_data as MD
    import pyarrow.parquet as pq

    from vulcan.chimeric.compete import IonSets, chimeric_competition
    from vulcan.index.build import L1Index
    from vulcan.stats import entrapment as ENT
    dbi = json.load(open(os.path.join(M3E, "armM", "armM_db.json")))
    out = {"design": __doc__.split("Rule")[0].strip(), "runs": RUNS, "seeds": list(SEEDS), "r_eff": 1.0,
           "r_peptide_eff": {s: dbi["seeds"][str(s)]["r_peptide_eff"] for s in SEEDS}, "ms1_hypotheses_sha256": ms1_sha,
           "v2_sha256": V2_SHA, "scorers": {}, "competition": {}}
    per = {}
    for s in SEEDS:
        lab = MD.labels(idx(s))
        l1 = L1Index(idx(s))
        ions = IonSets(l1)
        parts = [pq.read_table(os.path.join(pdir(s), r + ".chim.scored.parquet")).to_pandas() for r in RUNS]
        for i, p in enumerate(parts):
            p["ri"] = i
        import pandas as pd
        psm = pd.concat(parts, ignore_index=True)
        row = psm.pf_row.values
        spec = psm.ri.values.astype(np.int64) * (1 << 26) + psm.spec_idx.values.astype(np.int64)
        for name, col in (("v2_sep", "score_v2"), ("hyperscore_sep", "hyperscore")):
            cc = chimeric_competition(spec, psm.group.values, psm.chim_rank.values, psm.pf_id.values,
                                      lab["tdc_pair"][row], psm[col].values.astype(np.float64), lab["is_decoy"][row],
                                      row, ions, seed=11)
            res = cc["res"]
            rows_w = row[cc["win"]]
            per.setdefault(name, {})[s] = {"is_decoy": res.is_decoy, "q": cc["q_sep"], "stratum": cc["stratum"],
                                           "ent": lab["is_entrapment"][rows_w], "unit": psm.pf_id.values[cc["win"]]}
            out["competition"].setdefault(name, {})[s] = {k: cc[k] for k in ("pair_counts", "merge", "n_group_winners",
                                                                             "n_after_neighbour_filter",
                                                                             "n_dropped_by_neighbour_filter")}
        print("eval loaded", s, flush=True)
    lists = {"overall": None, "1": 1, "2": 2, ">=3": 3}
    violated = []
    for name, d in per.items():
        rec = {}
        for ln, k in lists.items():
            rec[ln] = {}
            for a in ALPHAS:
                reps, n_tab = [], {}
                for s in SEEDS:
                    x = d[s]
                    ts = ~x["is_decoy"] if k is None else (~x["is_decoy"]) & (x["stratum"] == k)
                    acc = x["q"] <= a
                    reps.append(ENT.Replicate(ENT.combined_counts(x["unit"][ts], x["ent"][ts], acc[ts]), s, 1001))
                    n_tab[s] = {"N_T": int((acc & ts & ~x["ent"]).sum()), "N_E": int((acc & ts & x["ent"]).sum())}
                ev = {cv: ENT.evaluate_arm(reps, a, "combined", r=1.0, ci_convention=cv) for cv in ("a", "b")}
                rec[ln][str(a)] = {"counts": n_tab, **{cv: {kk: e[kk] for kk in ("mean_fdp", "ci95", "verdict", "mean_lb",
                                                                                 "lb_ci95", "mean_n_t")}
                                                        for cv, e in ev.items()}}
                if name == "v2_sep" and a in GATED and ev["a"]["verdict"] == "violated":
                    violated.append(f"{ln}@{a}")
                print(name, ln, a, round(ev["a"]["mean_fdp"], 5), ev["a"]["verdict"], "| b", ev["b"]["verdict"],
                      n_tab, flush=True)
        out["scorers"][name] = rec
    out["violated_v2_sep_gated"] = violated
    out["violated_v2_sep_gated_convention_b"] = [f"{ln}@{a}" for ln in lists for a in GATED
                                                 if out["scorers"]["v2_sep"][ln][str(a)]["b"]["verdict"] == "violated"]
    out["pass"] = not violated
    json.dump(out, open(out_path, "w"), indent=1, default=float)
    print("arm M v2_sep pass:", out["pass"], violated)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms1", action="store_true")
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--eval", action="store_true")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "2026-10-28_A3_preconditions", "armM.json"))
    a = ap.parse_args()
    t0 = time.perf_counter()
    sha = ms1() if (a.ms1 or a.eval) else None
    if a.search:
        search()
    if a.score:
        score()
    if a.eval:
        evaluate(a.out, sha)
    print("wall", round(time.perf_counter() - t0, 1), "s")


if __name__ == "__main__":
    main()
