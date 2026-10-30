#!/usr/bin/env python3
"""3c runtime, production configuration, in ONE persistent process (setup excluded and reported separately):
per run = store load + GPU search (IT-CID, Da tolerance, top-5) + L3 features (Prosit cache warm: the dev
predictions are already cached; the cold first-contact Koina cost is reported from feat_info.json) + frozen model
apply + (pooled) paired double competition.  The timed loop holds `flock /tmp/vulcan_bench.lock` only for its own
duration (a few seconds per run); setup (imports, index upload) runs unlocked.

    python scripts/m3c/bench.py --model itcid_gbm_pred  -> results/2026-10-27_3c_itcid/bench.json
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402


def core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--n-runs", type=int, default=4)
    a = ap.parse_args()
    t_imp = time.perf_counter()
    import numpy as np
    import torch
    from vulcan import hw
    from vulcan.gpu import engine as E
    from vulcan.rescore.static import StaticModel
    import features as F
    t_imp = time.perf_counter() - t_imp
    model = StaticModel.load(os.path.join(C.RES, "models", a.model + ".pkl"))
    out = {"model": a.model, "model_sha256": model.meta.get("sha256"), "python_import_s": t_imp,
           "hw": hw.fingerprint() if hasattr(hw, "fingerprint") else None, "sets": {}}
    for ds in ("D4", "D10cf"):
        fam = C.DEV[ds]["family"]
        t0 = time.perf_counter()
        p = E.SearchParams(**C.SEARCH)
        l1 = C.load_l1(fam, "prod")
        st, var = C.mods(fam)
        gi = E.GpuIndex(l1, p, st, var)
        fb, _ = F.builder(fam)
        stores = C.DEV[ds]["stores"]()
        # warm-up run (CUDA kernels, cache shards load), excluded
        sp = E.load_store(stores[0])
        tbl, _, _ = E.search(sp, gi, l1, p)
        setup = time.perf_counter() - t0
        runs = []
        lock = open(os.environ.get("VULCAN_BENCH_LOCK", "/tmp/vulcan_bench.lock"), "a")
        tw = time.perf_counter()
        fcntl.flock(lock, fcntl.LOCK_EX)
        lock_wait = time.perf_counter() - tw
        try:
            for s in stores[1:1 + a.n_runs]:
                r = {}
                t, c = time.perf_counter(), core_s()
                sp = E.load_store(s)
                r["load_s"] = time.perf_counter() - t
                t1 = time.perf_counter()
                tbl, stt, _ = E.search(sp, gi, l1, p)
                torch.cuda.synchronize()
                r["search_s"] = time.perf_counter() - t1
                t1 = time.perf_counter()
                psm = {k: tbl.column(k).to_numpy(zero_copy_only=False) for k in
                       ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt", "delta_next", "delta_best",
                        "matched_peaks", "matched_b", "matched_y", "matched_intensity_pct", "average_ppm",
                        "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages",
                        "n_var", "isotope_error", "scored_candidates", "poisson")}
                feat, info = fb.run(l1, psm, sp)
                r["features_s"] = time.perf_counter() - t1
                r["predict_s_in_features"] = info["cache_stats"]["predict_s"]
                t1 = time.perf_counter()
                sc = model.apply(psm, feat)
                r["model_apply_s"] = time.perf_counter() - t1
                r["total_s"] = time.perf_counter() - t
                r["total_core_s"] = core_s() - c
                r["n_ms2"] = len(sp)
                r["run"] = sp.run
                runs.append(r)
                print(ds, json.dumps(r, default=float), flush=True)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()
        med = {k: float(np.median([r[k] for r in runs])) for k in runs[0] if k not in ("run",)}
        out["sets"][ds] = {"setup_s_unlocked": setup, "lock_wait_s": lock_wait, "runs": runs, "median": med,
                           "loadavg": os.getloadavg()}
        del gi
        torch.cuda.empty_cache()
    json.dump(out, open(os.path.join(C.RES, "bench.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
