#!/usr/bin/env python3
"""Validate family top-2 (vulcan/gpu/family.py): every (decoy seed j x entrapment seed s) replicate competition
reconstructed from ONE search over the multi-replicate index must equal a SEPARATE search over that replicate's
index (row subset T + E_s + dec_j(T) + dec_j(E_s), same peptidoforms), with the same per-family top-K rule.

Equality is exact per spectrum on the top-2 (pf_id, hyperscore bit pattern). Also reported, for information: agreement
of the reconstructed top-1 with a separate search in the default (global top-K) engine mode.

    python scripts/m2b/family_validate.py --ds D1 --out results/.../family_validation_D1.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m2"))

ENTRAP = [2001, 2002]
DECOYS = [1001, 1002]


def ensure_index(ds, work):
    import e2_common as C
    from vulcan.index.build import build_multi_decoy
    d = os.path.join(work, "index", f"{C.DATASETS[ds]['family']}_E{'-'.join(map(str, ENTRAP))}_D{'-'.join(map(str, DECOYS))}")
    if not os.path.exists(os.path.join(d, "meta.json")):
        t = time.perf_counter()
        build_multi_decoy([C.ensure_fasta(ds)], d, {"entrap_seeds": ENTRAP, "decoy_seeds": DECOYS})
        print(f"built {d} in {time.perf_counter() - t:.0f}s", flush=True)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ds", default="D1")
    ap.add_argument("--runs", type=int, default=1, help="number of the dataset's raw files to use")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=f"/tmp/vulcan_m2b_{os.getuid()}")
    a = ap.parse_args()
    import numpy as np
    import torch

    import e2_common as C
    from vulcan.gpu import engine as E
    from vulcan.gpu import family as F
    from vulcan.index.build import L1Index
    idx = ensure_index(a.ds, a.work)
    l1 = L1Index(idx)
    p = E.SearchParams()
    S, M = {"C": 57.021464}, {"M": 15.9949}
    t = time.perf_counter()
    gi = E.GpuIndex(l1, p, S, M)
    fi = F.FamilyIndex.from_l1(l1, gi)
    res = {"index": idx, "index_meta": {k: l1.meta[k] for k in ("n_peptidoforms", "n_peptidoforms_by_family",
                                                                   "decoy_sets", "entrapment", "timing_s")},
           "combined_gpu_index_build_s": time.perf_counter() - t, "combined_postings": gi.n_postings,
           "entrap_seeds": ENTRAP, "decoy_seeds": DECOYS, "runs": []}
    fwd_fam_np = fi.fwd_fam.cpu().numpy()
    side_np = fi.side.cpu().numpy()
    for raw in C.DATASETS[a.ds]["raw"][: a.runs]:
        sp = E.load_store(C.ensure_store(raw))
        torch.cuda.synchronize()
        t = time.perf_counter()
        tbl, st = F.search_family(sp, gi, l1, fi, p)
        torch.cuda.synchronize()
        rr = {"run": sp.run, "n_spectra": len(sp), "combined_search_s": time.perf_counter() - t,
              "combined_gpu_ms": st.gpu_ms, "combined_rows": tbl.num_rows, "combined_rescored": st.counts["n_rescored"],
              "replicates": []}
        t_sep = 0.0
        for j, dseed in enumerate(DECOYS):
            for si, eseed in enumerate(ENTRAP):
                keep = np.isin(fwd_fam_np, [0, 1 + si]) & np.isin(side_np, [-1, j])
                sub = F.SubIndex(l1, keep, j)
                t = time.perf_counter()
                gs = E.GpuIndex(sub, p, S, M)
                fs = F.FamilyIndex(fwd_fam_np[keep], np.where(side_np[keep] >= 0, 0, -1), [sub.partner], fi.F,
                                   gs.dev)
                ts, _ = F.search_family(sp, gs, sub, fs, p)
                torch.cuda.synchronize()
                t_sep += time.perf_counter() - t
                sep = F.reconstruct(ts, sorted(set(ts.column("fgroup").to_numpy().tolist())))
                rec = F.reconstruct(tbl, fi.replicate_groups(j, si))
                same_len = len(sep["spec_idx"]) == len(rec["spec_idx"])
                eq = {k: bool(same_len and np.array_equal(sep[k], rec[k])) for k in ("spec_idx", "rank", "pf_id")}
                eq["hyperscore_bits"] = bool(same_len and np.array_equal(sep["hyperscore"].view(np.uint64),
                                                                          rec["hyperscore"].view(np.uint64)))
                n_mis = None
                if not all(eq.values()):
                    ks = lambda d: set(zip(d["spec_idx"].tolist(), d["rank"].tolist(), d["pf_id"].tolist(),  # noqa
                                           d["hyperscore"].tolist()))
                    n_mis = len(ks(sep) ^ ks(rec))
                # informational: default engine mode (global top-K=500) separate search, top-1 agreement
                dt, _, _ = E.search(sp, gs, sub, p, report=1)
                d1 = dict(zip(dt.column("spec_idx").to_pylist(), dt.column("pf_id").to_pylist()))
                r1 = rec["rank"] == 1
                agree = np.mean([d1.get(s_) == pid for s_, pid in zip(rec["spec_idx"][r1].tolist(),
                                                                       rec["pf_id"][r1].tolist())])
                rr["replicates"].append({"decoy_seed": dseed, "entrap_seed": eseed, "n_rows_subset": int(keep.sum()),
                                         "n_top2_rows": int(len(rec["spec_idx"])), "exact": eq,
                                         "n_mismatched_rows": n_mis,
                                         "top1_agreement_vs_default_mode_separate": float(agree)})
                print(sp.run, dseed, eseed, eq, n_mis, f"default-mode top1 agreement {agree:.4f}", flush=True)
                del gs, fs
                torch.cuda.empty_cache()
        rr["separate_searches_total_s_incl_index_build"] = t_sep
        rr["all_exact"] = all(all(x["exact"].values()) for x in rr["replicates"])
        res["runs"].append(rr)
    res["all_exact"] = all(r["all_exact"] for r in res["runs"])
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print("ALL EXACT:", res["all_exact"])


if __name__ == "__main__":
    main()
