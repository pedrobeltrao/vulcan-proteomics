#!/usr/bin/env python3
"""Predictor: reduced-precision (autocast fp16 / bf16) and batch-size check vs the frozen fp32 path.

Same 100k random peptidoforms of the human prod index (seed 1), z=2, NCE 26: predictions under fp32 (batch 16384,
the frozen setting), autocast fp16, autocast bf16, and fp32 with larger batches.  Reports exact-equality of the
stored values (MS2 is stored as fp16, RT as fp32), max |diff|, and throughput.  Run under the bench lock.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_common as S  # noqa: E402

import numpy as np  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import torch
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    l1 = L1Index(S.C.index_dir("human_contam", "prod"))
    pred = APD()
    rows = np.random.default_rng(1).choice(len(l1), 100_000, replace=False)

    def P(ms2, rt, batch=16384, dtype=None):
        torch.cuda.synchronize()
        t = time.perf_counter()
        if dtype is None:
            o = pred.predict(l1.residues, l1.seq_off[rows], l1.length[rows], l1.mask[rows].view(np.int64),
                             np.full(len(rows), 2), 26.0, batch=batch, want_ms2=ms2, want_rt=rt)
        else:
            with torch.autocast("cuda", dtype=dtype):
                o = pred.predict(l1.residues, l1.seq_off[rows], l1.length[rows], l1.mask[rows].view(np.int64),
                                 np.full(len(rows), 2), 26.0, batch=batch, want_ms2=ms2, want_rt=rt)
        torch.cuda.synchronize()
        return o, time.perf_counter() - t
    P(True, True)                       # warm-up
    out = {"n": len(rows), "charge": 2, "nce": 26.0}
    ref_m, tm = P(True, False)
    ref_r, tr = P(False, True)
    out["fp32_b16384"] = {"ms2_per_s": len(rows) / tm, "rt_per_s": len(rows) / tr}
    for name, kw in {"fp16_autocast": {"dtype": torch.float16}, "bf16_autocast": {"dtype": torch.bfloat16},
                     "fp32_b65536": {"batch": 65536}, "fp32_b262144": {"batch": 262144}}.items():
        try:
            m, tm = P(True, False, **kw)
            r, tr = P(False, True, **kw)
        except Exception as e:  # noqa: BLE001
            out[name] = {"error": repr(e)}
            continue
        dm = np.abs(m["ms2"].astype(np.float32) - ref_m["ms2"].astype(np.float32))
        dr = np.abs(r["rt"] - ref_r["rt"])
        out[name] = {"ms2_per_s": len(rows) / tm, "rt_per_s": len(rows) / tr,
                     "ms2_identical": bool(np.array_equal(m["ms2"], ref_m["ms2"])),
                     "ms2_frac_values_differ": float((dm > 0).mean()), "ms2_max_abs_diff": float(dm.max()),
                     "rt_identical": bool(np.array_equal(r["rt"], ref_r["rt"])),
                     "rt_frac_differ": float((dr > 0).mean()), "rt_max_abs_diff": float(dr.max())}
        print(name, out[name], flush=True)
    print("fp32", out["fp32_b16384"])
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
