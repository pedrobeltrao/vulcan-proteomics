#!/usr/bin/env python3
"""E2: K1 shortlist size K vs recall of the exhaustive best target / best decoy per spectrum, and K1/K2 GPU time.

    python scripts/m2/e2_topk_sweep.py [--datasets D0 D1] [--kinds prod E2001] [--ks 50 100 200 500 1000]

Recall is measured against K2 over every candidate in the precursor windows (exhaustive), for targets and for
decoys separately (labels used only for evaluation). Output: results/.../topk_sweep.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_common as C  # noqa: E402
from e2_correctness import _best_per_label, _rescored, labels_by_row  # noqa: E402

import torch  # noqa: E402

from vulcan.gpu import engine as E  # noqa: E402
from vulcan.index.build import L1Index  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1"])
    ap.add_argument("--kinds", nargs="+", default=["prod", f"E{C.ENTRAP_SEED}"])
    ap.add_argument("--ks", nargs="+", type=int, default=[50, 100, 200, 500, 1000])
    ap.add_argument("--out", default=os.path.join(C.RES, "topk_sweep.json"))
    a = ap.parse_args()
    res = json.load(open(a.out)) if os.path.exists(a.out) else {}
    for ds in a.datasets:
        raw = C.DATASETS[ds]["raw"][0]
        sp = E.load_store(C.store_path(raw))
        for kind in a.kinds:
            d = C.index_dir(ds, kind)
            l1 = L1Index(d)
            lab = labels_by_row(d)
            p0 = E.SearchParams()
            gi = E.GpuIndex(l1, p0, C.STATIC, C.VAR)
            n = len(sp)
            _, stx, dx = E.search_exhaustive(sp, gi, l1, p0)
            sx, hx, rx = _rescored(dx, l1, p0)
            btx, bdx = _best_per_label(sx, hx, lab["is_decoy"][rx], n)
            # spectra that matter for a 1% list: exhaustive best (any label) in the top 50% of scores
            best = np.maximum(btx, bdx)
            fin = np.isfinite(best)
            hi = fin & (best >= np.quantile(best[fin], 0.5))
            for K in a.ks:
                p = E.SearchParams(topk=K)
                E.search(sp, gi, l1, p)          # warm-up
                _, st, dk = E.search(sp, gi, l1, p)
                sk, hk, rk = _rescored(dk, l1, p)
                btk, bdk = _best_per_label(sk, hk, lab["is_decoy"][rk], n)
                ht, hd = np.isfinite(btx), np.isfinite(bdx)
                rt, rd = float(np.mean(btk[ht] == btx[ht])), float(np.mean(bdk[hd] == bdx[hd]))
                rec = {"K": K, "recall_target": rt, "recall_decoy": rd, "abs_delta": abs(rt - rd),
                       "recall_target_upper_half": float(np.mean((btk == btx)[ht & hi])),
                       "recall_decoy_upper_half": float(np.mean((bdk == bdx)[hd & hi])),
                       "n_rescored": int(st.counts["n_rescored"]), "n_window_candidates": int(st.counts["k1_window_candidates"]),
                       "gpu_ms": st.gpu_ms, "wall_s_total": st.wall_s["total"]}
                key = f"{ds}/{kind}/{C.run_name(raw)}/K{K}"
                res[key] = rec
                print(key, json.dumps(rec), flush=True)
            del gi
            torch.cuda.empty_cache()
            os.makedirs(os.path.dirname(a.out), exist_ok=True)
            with open(a.out, "w") as fh:
                json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
