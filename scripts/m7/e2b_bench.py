#!/usr/bin/env python3
"""E2b: reranker throughput vs parameter count on the A100 (DESIGN §10 E2b).

For every size in vulcan.encoders.model.SIZES x dtype {fp16, bf16} x K {10, 20, 50}: synthetic but shape-realistic
inputs (63 peaks/spectrum, peptide lengths drawn from the D0 rank-1 length distribution 7..35, charge 2-3), full
split+joint forward in inference mode; spectra/s = S / wall (each spectrum's encoder runs once, its K candidates
through the peptide encoder and the joint stage).  Also the split stage alone (bi-encoder embedding rate:
spectra/s and peptides/s).  Each timed section holds /tmp/vulcan_bench.lock (fcntl) for a few seconds only;
loadavg + nvidia-smi utilisation of other processes recorded per section.  INT8: no GPU int8 path in this torch
build for these layers (torch.ao dynamic quant is CPU-only), so it is recorded as not measured.

    python scripts/m7/e2b_bench.py --out results/2026-10-28_M7_encoders/e2b.json
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import torch  # noqa: E402

from vulcan.encoders.model import SIZES, SplitJointEncoder, n_params  # noqa: E402

LOCK = "/tmp/vulcan_bench.lock"


def gpu_util():
    try:
        o = subprocess.run(["nvidia-smi", "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=10).stdout.strip()
        u, m = o.split(",")
        return {"util_pct": int(u), "mem_used_mib": int(m)}
    except Exception:  # noqa: BLE001
        return None


def make_batch(S, K, dev, rng, P=63):
    n = torch.full((S,), P, device=dev)
    mz = torch.sort(torch.rand(S, P, device=dev) * 1800 + 150, 1).values
    it = torch.rand(S, P, device=dev)
    ch = torch.randint(2, 4, (S,), device=dev)
    lens = rng.integers(7, 36, S * K)
    Lm = int(lens.max())
    tok = np.zeros((S * K, Lm), np.int64)
    for i, l in enumerate(lens):
        tok[i, :l] = rng.integers(1, 27, l)
    tok[tok == 2] = 1; tok[tok == 10] = 12; tok[tok == 15] = 16; tok[tok == 21] = 22; tok[tok == 24] = 25; tok[tok == 26] = 25
    return ({"mz": mz, "int": it, "n": n, "charge": ch},
            {"tokens": torch.from_numpy(tok).to(dev), "charge": ch.repeat_interleave(K)},
            torch.arange(S, device=dev).repeat_interleave(K))


def timed(fn, reps):
    torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    torch.cuda.synchronize()
    return (time.perf_counter() - t0) / reps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sizes", nargs="+", default=list(SIZES))
    ap.add_argument("--ks", nargs="+", type=int, default=[10, 20, 50])
    ap.add_argument("--dtypes", nargs="+", default=["fp16", "bf16"])
    ap.add_argument("--cand-budget", type=int, default=16384, help="candidates per forward batch")
    a = ap.parse_args()
    dev = torch.device("cuda")
    torch.backends.cuda.matmul.allow_tf32 = True
    rng = np.random.default_rng(11)
    out = {"rows": [], "int8": "not measured: no GPU int8 kernels for these layers in torch 2.x without TensorRT"}
    for size in a.sizes:
        cfg = SIZES[size]
        for dt in a.dtypes:
            dtype = torch.float16 if dt == "fp16" else torch.bfloat16
            torch.manual_seed(11)
            m = SplitJointEncoder(cfg).to(dev).to(dtype).eval()
            npar = n_params(m)
            for K in a.ks:
                S = max(64, a.cand_budget // K)
                spec, pep, sidx = make_batch(S, K, dev, rng)
                with torch.inference_mode():
                    def full():
                        m(spec, pep, sidx)

                    def split_spec():
                        m.encode_spectra(spec["mz"], spec["int"], spec["n"], spec["charge"])

                    def split_pep():
                        m.encode_peptides(pep["tokens"], pep["charge"])
                    for _ in range(3):
                        full()
                    torch.cuda.synchronize()
                    with open(LOCK, "w") as lk:
                        t_wait = time.perf_counter()
                        fcntl.flock(lk, fcntl.LOCK_EX)
                        waited = time.perf_counter() - t_wait
                        load = os.getloadavg()
                        gu = gpu_util()
                        tf = timed(full, 10)
                        ts = timed(split_spec, 10) if K == a.ks[0] else None
                        tp = timed(split_pep, 10)
                        fcntl.flock(lk, fcntl.LOCK_UN)
                mem = torch.cuda.max_memory_allocated() / 2 ** 30
                torch.cuda.reset_peak_memory_stats()
                r = {"size": size, "params": npar, "dtype": dt, "K": K, "batch_spectra": S, "batch_candidates": S * K,
                     "full_spectra_per_s": S / tf, "full_candidates_per_s": S * K / tf,
                     "split_spectra_per_s": (S / ts) if ts else None, "split_peptides_per_s": S * K / tp,
                     "peak_mem_gib": mem, "loadavg": load, "gpu_before": gu, "lock_wait_s": waited,
                     "config": cfg.to_dict()}
                out["rows"].append(r)
                print(json.dumps({k: (round(v, 1) if isinstance(v, float) else v) for k, v in r.items()
                                  if k not in ("config",)}), flush=True)
            del m
            torch.cuda.empty_cache()
    from vulcan import hw  # noqa: F401
    try:
        out["hwinfo"] = json.loads(subprocess.run([sys.executable, os.path.join(ROOT, "bench", "hwinfo.py")],
                                                  capture_output=True, text=True, timeout=120).stdout)
    except Exception as e:  # noqa: BLE001
        out["hwinfo"] = str(e)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
