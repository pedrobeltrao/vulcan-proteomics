#!/usr/bin/env python3
"""One-time predicted-library build for production scoring (per proteome x instrument class x NCE grid).

1. Copy the M3 prediction cache root (the cache the frozen reference features of D0/D1/D8 were computed from) to
   $VULCAN_SS_WORK/predlib.  The cache is append-only and first-copy-wins, so every value a reference feature
   used is preserved.
2. For every NCE of the frozen grid (vulcan.rescore.features.NCE_GRID) and each production index (human_contam_prod,
   pfu_contam_prod): PredCache.precompute (the frozen M3 code path: predicts MISSES only, batch 16384, fp32) for
   RT and MS2 at precursor charges 2, 3, 4 (z>=5, ~1% of PSMs, stays lazy).
3. Consolidate every key (vulcan.predict.library.consolidate) for memory-mapped loading.
Not a timed section (GPU shared with other agents; load recorded).
    /tmp/vulcan_env_*/pred/bin/python scripts/scoring_speed/build_library.py --out results/.../library_build.json
"""
import argparse
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_common as S  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--charges", default="2,3,4")
    ap.add_argument("--root", default=os.path.join(S.W, "predlib"))
    a = ap.parse_args()
    import torch
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.predict.library import consolidate, library_digest
    from vulcan.rescore.features import NCE_GRID
    from vulcan.hw import fingerprint
    out = {"root": a.root, "source": S.REF_CACHE, "charges": [int(z) for z in a.charges.split(",")], "steps": []}
    t0 = time.perf_counter()
    if not os.path.exists(os.path.join(a.root, ".copied")):
        os.makedirs(a.root, exist_ok=True)
        subprocess.run(["cp", "-a", S.REF_CACHE + "/.", a.root + "/"], check=True)
        open(os.path.join(a.root, ".copied"), "w").write(time.ctime())
    out["copy_s"] = time.perf_counter() - t0
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("VULCAN_SS_GPU_FRAC", "0.2")))  # shared GPU
    pred = APD()
    keys = []
    for fam in ("human_contam", "pfu_contam"):
        l1 = L1Index(S.C.index_dir(fam, "prod"))
        for nce in NCE_GRID:
            c = PredCache(a.root, pred, nce, "QE", S.C.STATIC, S.C.VAR)
            keys.append(c.key)
            t1 = time.perf_counter()
            c.precompute(l1, charges=tuple(out["charges"]))
            torch.cuda.synchronize()
            rec = {"family": fam, "nce": nce, "key": c.key, "n_rows": len(l1), "wall_s": time.perf_counter() - t1,
                   **{k: c.stats[k] for k in ("n_predicted_rt", "n_predicted_ms2", "predict_s", "library_load_s")},
                   "loadavg": os.getloadavg()}
            out["steps"].append(rec)
            print(json.dumps(rec), flush=True)
            del c
            torch.cuda.empty_cache()
            json.dump(out, open(a.out, "w"), indent=1, default=float)
    t2 = time.perf_counter()
    for k in sorted(set(keys)):
        consolidate(os.path.join(a.root, k))
    out["consolidate_s"] = time.perf_counter() - t2
    out["digest"] = library_digest(a.root, sorted(set(keys)))
    out["keys"] = sorted(set(keys))
    out["total_s"] = time.perf_counter() - t0
    out["bytes"] = sum(os.path.getsize(os.path.join(dp, f)) for k in set(keys)
                       for dp, _, fs in os.walk(os.path.join(a.root, k)) for f in fs)
    out["hw"] = fingerprint()
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
