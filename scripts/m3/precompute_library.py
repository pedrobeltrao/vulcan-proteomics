#!/usr/bin/env python3
"""Precompute the AlphaPeptDeep tryptic library (RT + MS2 at charges 2,3) for every peptidoform of an L1 index
(proteome x instrument class OT-HCD 'QE' x NCE), into the content-hashed cache; reports predictor throughput.

    flock /tmp/vulcan_bench.lock conda_envs/pred/bin/python scripts/m3/precompute_library.py --index DIR --nce 29 --out JSON
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import numpy as np  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--nce", type=float, default=29.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import torch
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.hw import fingerprint
    pred = APD()
    l1 = L1Index(a.index)
    # warm-up (per-length kernel setup) on 50k rows, timed separately
    rng = np.random.default_rng(11)
    w = rng.choice(len(l1), 50000, replace=False)
    t0 = time.perf_counter()
    pred.predict(l1.residues, l1.seq_off[w], l1.length[w], l1.mask[w].view(np.int64), np.full(len(w), 2), a.nce)
    torch.cuda.synchronize()
    warm = time.perf_counter() - t0
    # pure inference throughput on 1M rows (no cache I/O)
    s = rng.choice(len(l1), 1_000_000, replace=False)
    out = {}
    for kind, (ms2, rt) in {"rt": (False, True), "ms2": (True, False)}.items():
        torch.cuda.synchronize(); t0 = time.perf_counter()
        pred.predict(l1.residues, l1.seq_off[s], l1.length[s], l1.mask[s].view(np.int64), np.full(len(s), 2), a.nce,
                     want_ms2=ms2, want_rt=rt)
        torch.cuda.synchronize()
        out[f"{kind}_per_s"] = len(s) / (time.perf_counter() - t0)
    cache = PredCache(os.path.join(C.WORK, "predcache"), pred, a.nce, "QE", C.STATIC, C.VAR)
    t_lib = cache.precompute(l1, charges=(2, 3))
    out.update({"index": a.index, "index_content_hash": l1.meta["content_hash"], "n_peptidoforms": len(l1),
                "nce": a.nce, "instrument": "QE", "charges": [2, 3], "cache_key": cache.key,
                "library_wall_s": t_lib, "library_predict_s": cache.stats["predict_s"],
                "library_precursors_per_s_end_to_end": 2 * len(l1) / t_lib, "warmup_s": warm,
                "predictor": pred.version, "hw": fingerprint()})
    du = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(cache.dir) for f in fs)
    out["library_bytes"] = du
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "hw"}, indent=1, default=float))


if __name__ == "__main__":
    main()
