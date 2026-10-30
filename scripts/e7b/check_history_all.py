#!/usr/bin/env python3
"""E7b: for EVERY D8h (replicate, run) whose fast-path features/scores are not bitwise equal to the stored E7 values
(equality_scoring.json), run the OLD slow path (FeatureBuilderV2.run_v2 + PredCache) on the same library copy and
compare it with the fast path (FastFeatureBuilder(primary) + LibCache), all v1 feature columns + v2 score bitwise.
slow == fast everywhere => the port is exact; the differences vs E7 come from the prediction cache's history (E7 scored
D8h in parallel processes that appended predictions of the same peptides with batch-dependent last bits)."""
import hashlib
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402


def main(out):
    import pyarrow.parquet as pq
    from bench.t1.ours_chimeric import MODEL_V2, PSM_COLS
    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.predict.library import LibCache
    from vulcan.rescore.fast import FastFeatureBuilder, FastStatic
    from vulcan.rescore.features_v2 import FeatureBuilderV2, model_psm, read_ms1_features
    from vulcan.rescore.static import PRED_FEATURES as V1, StaticModel
    st, vr = {"C": 57.021464}, {"M": 15.9949}
    W = f"/tmp/vulcan_t1_{os.getuid()}/D8h"
    lib = os.path.join(ROOT, "data", "cache", "e7b_predlib", "D8h")
    eq = json.load(open(os.path.join(ROOT, "results", "2026-10-28_E7b_fair_speed", "equality_scoring.json")))
    pairs = [(rep, run) for rep, x in eq["datasets"]["D8h"]["replicates"].items() for run, y in x["per_run"].items()
             if not y["score"] or not all(y["features"].values())]
    v2 = StaticModel.load(MODEL_V2)
    fs = FastStatic(v2)
    pred = APD()
    pc, lc = {}, {}
    slow = FeatureBuilderV2(lambda n: pc.setdefault(n, PredCache(lib, pred, n, "QE", st, vr)), st, vr)
    fast = FastFeatureBuilder(lambda n: lc.setdefault(n, LibCache(lib, pred, n, "QE", st, vr, device="cuda")), st, vr)
    sps, l1s, res = {}, {}, {}
    for rep, run in pairs:
        if run not in sps:
            sps[run] = load_store(os.path.join(W, "stage", run + ".vspec.arrow"))
        if rep not in l1s:
            l1s.clear()
            l1s[rep] = L1Index(os.path.join(W, "index", rep))
        l1 = l1s[rep]
        d = os.path.join(W, "ours_chim", rep)
        t = pq.read_table(os.path.join(d, f"{run}.chim.psms.parquet"))
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
        f_s, _ = slow.run_v2(l1, psm, sps[run], read_ms1_features(os.path.join(W, "ms1", run + ".features.arrow")))
        f_f, _ = fast.run(l1, psm, sps[run], primary=psm["group"] == 0)
        ref_s = pq.read_table(os.path.join(d, f"{run}.chim.scored.parquet")).column("score_v2").to_numpy()
        s_s = np.asarray(v2.apply(model_psm(psm), f_s), np.float64)
        s_f = np.asarray(fs.apply(model_psm(psm), f_f), np.float64)
        res[f"{rep}/{run}"] = {
            "slow_vs_fast_features_bitwise": all(np.asarray(f_s[k], np.float32).tobytes() ==
                                                 np.asarray(f_f[k], np.float32).tobytes() for k in V1),
            "slow_vs_fast_score_bitwise": s_s.tobytes() == s_f.tobytes(),
            "n_psms_score_differs_vs_E7": int((s_f != ref_s).sum()), "n_psms": int(len(s_f)),
            "max_abs_score_diff_vs_E7": float(np.abs(s_f - ref_s).max())}
        print(rep, run[:20], res[f"{rep}/{run}"], flush=True)
    o = {"n_pairs": len(pairs), "all_slow_equals_fast": all(v["slow_vs_fast_features_bitwise"] and
                                                            v["slow_vs_fast_score_bitwise"] for v in res.values()),
         "total_psms_differing_vs_E7": sum(v["n_psms_score_differs_vs_E7"] for v in res.values()),
         "total_psms_in_pairs": sum(v["n_psms"] for v in res.values()), "pairs": res}
    json.dump(o, open(out, "w"), indent=1, default=float)
    print({k: v for k, v in o.items() if k != "pairs"})


if __name__ == "__main__":
    main(sys.argv[1])
