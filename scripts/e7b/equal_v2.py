#!/usr/bin/env python3
"""E7b equivalence record for the fast v2_sep production path (dev only; results/2026-10-28_E7b_fair_speed).

For every replicate (5 prod, 15 arm P, 5 arm C) of D0, D1, D8h, on the STORED E7 chimeric PSM tables
(<work>/ours_chim/<rep>/<run>.chim.psms.parquet):
  features  vulcan.rescore.fast.FastFeatureBuilder.run(primary=group 0) on the consolidated copy of the dataset's own
            T1 prediction cache (data/cache/e7b_predlib/<ds>, GPU-resident), compared bitwise with the stored E7
            features (<run>.chim.feat.parquet, FeatureBuilderV2.run_v2 + PredCache) for every v1 feature column;
  scores    frozen v2 model (sha 3f2c69cc) on those features, bitwise vs the stored score_v2;
  q / IDs   dataset-level chimeric competition (spectrum key = run index * 2^26 + spec_idx, as bench/t1/ours_chimeric
            .chim_units): reference chimeric_competition on the stored scores vs chimeric_competition_fast on the new
            scores; winners, strata, q_sep, q_pooled compared exactly; target peptidoform (pf_id) sets at q_sep <= 0.1%,
            1%, 5% compared.
Labels (decoy / pair) are read only for the competition, as in every L4 step.  Writes --out JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

DATASETS = ("D0", "D1", "D8h")
ALPHAS = (0.001, 0.01, 0.05)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--datasets", nargs="+", default=list(DATASETS))
    ap.add_argument("--work-root", default=f"/tmp/vulcan_t1_{os.getuid()}")
    a = ap.parse_args()
    from vulcan import hw
    hw.pin_thread_pools(8)
    import pyarrow.parquet as pq
    import torch

    from bench.t1.ours_chimeric import MODEL_V2, PSM_COLS
    from vulcan.chimeric.compete import IonSets, chimeric_competition
    from vulcan.chimeric.fast_compete import chimeric_competition_fast
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.library import LibCache
    from vulcan.rescore import post as P
    from vulcan.rescore.fast import FastFeatureBuilder, FastStatic
    from vulcan.rescore.features_v2 import model_psm
    from vulcan.rescore.static import PRED_FEATURES as V1_FEATS, StaticModel
    from vulcan.stats._util import DECOY_SEEDS, ENTRAP_SEEDS_DEV
    from bench.t1.dbs import replicate_names
    st, vr = {"C": 57.021464}, {"M": 15.9949}
    v2 = StaticModel.load(MODEL_V2)
    v2_sha = hashlib.sha256(open(MODEL_V2, "rb").read()).hexdigest()
    assert v2_sha.startswith("3f2c69cc"), v2_sha
    fs = FastStatic(v2)
    pred = APD()
    out = {"model_v2_sha256": v2_sha, "v1_feature_columns": V1_FEATS, "datasets": {}, "hw": hw.fingerprint(),
           "loadavg_start": os.getloadavg()}
    reps = replicate_names(list(DECOY_SEEDS), list(ENTRAP_SEEDS_DEV))
    for ds in a.datasets:
        W = os.path.join(a.work_root, ds)
        lib = os.path.join(ROOT, "data", "cache", "e7b_predlib", ds)
        caches = {}

        def factory(nce):
            if nce not in caches:
                caches[nce] = LibCache(lib, pred, nce, "QE", st, vr, device="cuda")
            return caches[nce]
        fb = FastFeatureBuilder(factory, st, vr)
        runs = sorted(f[:-len(".chim.psms.parquet")] for f in os.listdir(os.path.join(W, "ours_chim", "prod_1001"))
                      if f.endswith(".chim.psms.parquet"))
        from vulcan.gpu.engine import load_store
        sps = {r: load_store(os.path.join(W, "stage", r + ".vspec.arrow")) for r in runs}
        rec = {"runs": runs, "replicates": {}}
        t_ds = time.time()
        for rep in reps:
            d = os.path.join(W, "ours_chim", rep)
            l1 = L1Index(os.path.join(W, "index", rep))
            rr = {"features_bitwise": True, "scores_bitwise": True, "n_psms": 0, "nce_equal": True,
                  "lib_misses": 0, "per_run": {}}
            new_scores = {}
            for run in runs:
                t = pq.read_table(os.path.join(d, f"{run}.chim.psms.parquet"))
                psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
                ref_f = pq.read_table(os.path.join(d, f"{run}.chim.feat.parquet"))
                ref_s = pq.read_table(os.path.join(d, f"{run}.chim.scored.parquet")).column("score_v2").to_numpy()
                ref_info = json.load(open(os.path.join(d, f"{run}.chim.feat.parquet.json")))
                feat, info = fb.run(l1, psm, sps[run], primary=psm["group"] == 0)
                s = fs.apply(model_psm(psm), feat)
                fe = {k: bool(np.asarray(feat[k], np.float32).tobytes() == ref_f.column(k).to_numpy().tobytes())
                      for k in V1_FEATS}
                se = bool(np.asarray(s, np.float64).tobytes() == np.asarray(ref_s, np.float64).tobytes())
                miss = int(info["cache_stats"].get("ms2_miss", 0) + info["cache_stats"].get("rt_miss", 0))
                rr["per_run"][run] = {"features": fe, "score": se, "nce": info["nce"], "ref_nce": ref_info["nce"],
                                      "misses": miss, "n": len(s)}
                rr["features_bitwise"] &= all(fe.values())
                rr["scores_bitwise"] &= se
                rr["nce_equal"] &= info["nce"] == ref_info["nce"]
                rr["lib_misses"] += miss
                rr["n_psms"] += len(s)
                new_scores[run] = s
            if rep.startswith("prod_"):
                ctx = P.IndexCtx(os.path.join(W, "index", rep))
                lab = ctx.lab
                parts = []
                for i, run in enumerate(runs):
                    t = pq.read_table(os.path.join(d, f"{run}.chim.scored.parquet")).to_pandas()
                    t["key"] = i * (1 << 26) + t.spec_idx.values.astype(np.int64)
                    t["new"] = new_scores[run]
                    parts.append(t)
                import pandas as pd
                T = pd.concat(parts, ignore_index=True)
                row = T.pf_row.values
                args = (T.key.values, T.group.values, T.chim_rank.values, T.pf_id.values, lab["tdc_pair"][row])
                t0 = time.time()
                ref = chimeric_competition(*args, T.score_v2.values, lab["is_decoy"][row], row, IonSets(l1), seed=11)
                t1 = time.time()
                new = chimeric_competition_fast(*args, T["new"].values, lab["is_decoy"][row], row, IonSets(l1),
                                                seed=11)
                t2 = time.time()
                eq = {k: bool(np.array_equal(ref[k], new[k])) for k in ("win", "survivors", "stratum", "merged",
                                                                          "q_sep", "q_pooled")}
                ids = {}
                for al in ALPHAS:
                    sr = set(T.pf_id.values[ref["win"]][(ref["q_sep"] <= al) & ~ref["res"].is_decoy].tolist())
                    sn = set(T.pf_id.values[new["win"]][(new["q_sep"] <= al) & ~new["res"].is_decoy].tolist())
                    ids[str(al)] = {"ref": len(sr), "new": len(sn), "identical": sr == sn,
                                    "jaccard": len(sr & sn) / max(1, len(sr | sn))}
                rr["competition_equal"] = eq
                rr["id_sets_q_sep"] = ids
                rr["competition_s"] = {"reference": t1 - t0, "fast": t2 - t1}
            rec["replicates"][rep] = rr
            print(ds, rep, {k: rr[k] for k in ("features_bitwise", "scores_bitwise", "nce_equal", "lib_misses",
                                               "n_psms")}, rr.get("id_sets_q_sep", {}).get("0.01"), flush=True)
            del l1
            torch.cuda.empty_cache()
        R = rec["replicates"].values()
        rec["all_features_bitwise"] = all(r["features_bitwise"] for r in R)
        rec["all_scores_bitwise"] = all(r["scores_bitwise"] for r in R)
        rec["all_nce_equal"] = all(r["nce_equal"] for r in R)
        rec["lib_misses_total"] = sum(r["lib_misses"] for r in R)
        rec["all_prod_id_sets_identical"] = all(v["identical"] for r in R if "id_sets_q_sep" in r
                                               for v in r["id_sets_q_sep"].values())
        rec["all_prod_competition_equal"] = all(all(r["competition_equal"].values()) for r in R
                                               if "competition_equal" in r)
        rec["wall_s"] = time.time() - t_ds
        out["datasets"][ds] = rec
        del sps, caches, fb
        torch.cuda.empty_cache()
    out["loadavg_end"] = os.getloadavg()
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({ds: {k: v for k, v in r.items() if k.startswith("all_") or k == "lib_misses_total"}
                      for ds, r in out["datasets"].items()}))


if __name__ == "__main__":
    main()
