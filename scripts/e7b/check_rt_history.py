#!/usr/bin/env python3
"""E7b: is the D8h run E9_0008 RT-feature difference (fast path vs stored E7 features) a property of the port or of the
prediction-cache history?  Same stored PSM table, same library copy (data/cache/e7b_predlib/D8h): the OLD slow path
(FeatureBuilderV2.run_v2 + PredCache) vs the fast path (FastFeatureBuilder(primary) + LibCache) vs the stored E7
features; plus the stored RT predictions' provenance (number of peptidoforms with >1 stored RT copy among the PSMs)."""
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402


def main(out):
    import pyarrow.parquet as pq
    from bench.t1.ours_chimeric import PSM_COLS
    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.predict.library import LibCache
    from vulcan.rescore.fast import FastFeatureBuilder
    from vulcan.rescore.features_v2 import FeatureBuilderV2, read_ms1_features
    st, vr = {"C": 57.021464}, {"M": 15.9949}
    W = f"/tmp/vulcan_t1_{os.getuid()}/D8h"
    lib = os.path.join(ROOT, "data", "cache", "e7b_predlib", "D8h")
    run, rep = "E9_0008__lgg1", "prod_1001"
    pred = APD()
    l1 = L1Index(os.path.join(W, "index", rep))
    sp = load_store(os.path.join(W, "stage", run + ".vspec.arrow"))
    d = os.path.join(W, "ours_chim", rep)
    t = pq.read_table(os.path.join(d, f"{run}.chim.psms.parquet"))
    psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
    ref = pq.read_table(os.path.join(d, f"{run}.chim.feat.parquet"))
    pc = {}
    slow = FeatureBuilderV2(lambda n: pc.setdefault(n, PredCache(lib, pred, n, "QE", st, vr)), st, vr)
    fs, _ = slow.run_v2(l1, psm, sp, read_ms1_features(os.path.join(W, "ms1", run + ".features.arrow")))
    lc = {}
    fast = FastFeatureBuilder(lambda n: lc.setdefault(n, LibCache(lib, pred, n, "QE", st, vr, device="cuda")), st, vr)
    ff, info = fast.run(l1, psm, sp, primary=psm["group"] == 0)
    rec = {}
    for k in ("rt_abs_resid", "rt_resid_z", "pred_sa"):
        a, b, r = np.asarray(fs[k], np.float32), np.asarray(ff[k], np.float32), ref.column(k).to_numpy()
        rec[k] = {"slow_vs_fast_bitwise": a.tobytes() == b.tobytes(), "slow_vs_E7_bitwise": a.tobytes() == r.tobytes(),
                  "fast_vs_E7_n_diff": int((b != r).sum()), "max_abs_diff_vs_E7": float(np.abs(b - r).max())}
    rtp_slow = pc[info["nce"]].lookup_rt(l1, psm["pf_row"])
    rtp_fast = lc[info["nce"]].lookup_rt(l1, psm["pf_row"])
    rec["rt_pred_slow_vs_fast_bitwise"] = np.asarray(rtp_slow).tobytes() == np.asarray(rtp_fast).tobytes()
    out_j = {"run": run, "replicate": rep, "nce": info["nce"], "columns": rec}
    json.dump(out_j, open(out, "w"), indent=1, default=float)
    print(json.dumps(out_j, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
