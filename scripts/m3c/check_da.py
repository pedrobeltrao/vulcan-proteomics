#!/usr/bin/env python3
"""Correctness of the Da fragment-tolerance mode (3c) on one D4 fraction (IT-CID, TMT):
  1. K1 index counts == brute-force counts (f32 Sage Tolerance::Da semantics) on every K1-appended candidate;
  2. top-K recall of the best target/decoy vs exhaustive K2 (K = 500);
  3. bin_ppm (K1 layout; results must be identical) 20 vs 200: identical PSMs, search time;
  4. Sage oracle: top-1 peptide agreement per scan with M1e's Sage D4 prod run (same tolerances; Sage searched TRFP
     mzML, we search our store, so a small disagreement from centroiding is expected).
    python scripts/m3c/check_da.py  -> results/2026-10-27_3c_itcid/check_da.json
"""
from __future__ import annotations

import json
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402


def main():
    import torch
    from vulcan.gpu import engine as E
    store = C.d4_stores()[0]
    fam = "human_tmt"
    l1 = C.load_l1(fam, "prod")
    st, var = C.mods(fam)
    sp = E.load_store(store)
    out = {"store": store}
    tabs = {}
    for bp in (20.0, 200.0):
        p = E.SearchParams(**{**C.SEARCH, "bin_ppm": bp, "report_psms": 1})
        gi = E.GpuIndex(l1, p, st, var)
        E.search(sp, gi, l1, p)                      # warm
        torch.cuda.synchronize()
        t0 = time.perf_counter()
        tbl, stt, aux = E.search(sp, gi, l1, p)
        out[f"search_s_bin{int(bp)}"] = time.perf_counter() - t0
        out[f"k1_gpu_ms_bin{int(bp)}"] = stt.gpu_ms.get("k1")
        tabs[bp] = tbl
        if bp == 200.0:
            # 1. K1 vs brute on appended candidates
            peaks, hyp = aux["peaks"], aux["hyp"]
            o_h, o_r, o_c, *_ = E.run_k1(gi, hyp, peaks, p, E.Stats())
            sel = torch.randperm(len(o_h), generator=torch.Generator().manual_seed(11))[:2_000_000].to(o_h.device)
            bru = E.k1_brute(gi, hyp, peaks, o_h[sel].long(), o_r[sel].long(), p)
            out["k1_vs_brute"] = {"n": int(len(sel)), "mismatch": int((bru != o_c[sel]).sum())}
            # 2. top-K recall vs exhaustive (first 3000 searched spectra)
            sub = np.unique(hyp["spec"])[:3000]
            ex, _, _ = E.search_exhaustive(sp, gi, l1, p, report=1, spec_subset=sub)
            from vulcan.index.build import read_labels
            lab = read_labels(C.index_dir(fam, "prod"))
            dec = lab.column("is_decoy").to_numpy()
            exs = {int(s): int(r) for s, r in zip(ex.column("spec_idx").to_numpy(), ex.column("pf_row").to_numpy())}
            p5 = E.SearchParams(**{**C.SEARCH, "report_psms": 1000})
            t5, _, _ = E.search(sp, gi, l1, p5)
            have = {}
            for s, r in zip(t5.column("spec_idx").to_numpy(), t5.column("pf_row").to_numpy()):
                have.setdefault(int(s), set()).add(int(r))
            rec = {"T": [0, 0], "D": [0, 0]}
            for s, r in exs.items():
                k = "D" if dec[r] else "T"
                rec[k][1] += 1
                rec[k][0] += r in have.get(s, ())
            out["topk_recall_best_exhaustive"] = {k: v[0] / max(v[1], 1) for k, v in rec.items()}
            out["topk_recall_n"] = {k: v[1] for k, v in rec.items()}
    out["bin20_vs_bin200_identical"] = bool(tabs[20.0].equals(tabs[200.0]))
    # 4. Sage oracle (M1e D4 prod, TRFP mzML, Sage 0.14.7, fragment +-0.5 Da)
    import csv
    run = os.path.basename(store)[:-len(".vspec.arrow")]
    sage = {}
    f = "/tmp/vulcan_m1e/cache/sage/D4/prod/results.sage.tsv"
    if os.path.exists(f):
        for row in csv.DictReader(open(f), delimiter="\t"):
            if row["filename"].startswith(run) and row["rank"] == "1":
                scan = int(row["scannr"].split("scan=")[-1])
                pep = re.sub(r"M\[\+15\.99\d*\]", "m", row["peptide"])
                pep = re.sub(r"\[[^\]]*\]-?", "", pep).replace("I", "L")
                sage[scan] = (pep, float(row["hyperscore"]), float(row["spectrum_q"]) <= 0.01 and row["label"] == "1")
        t = tabs[200.0]
        scans = t.column("scan").to_numpy()
        rows = t.column("pf_row").to_numpy()
        hs = t.column("hyperscore").to_numpy()
        n = agree = n_hi = agree_hi = 0
        dh = []
        for sc, r, h in zip(scans, rows, hs):
            if int(sc) not in sage:
                continue
            o = int(l1.seq_off[r])
            s = bytes(l1.residues[o:o + int(l1.length[r])]).decode()
            m = int(l1.mask[r]) & ((1 << 63) - 1)
            s = "".join(c.lower() if (m >> i) & 1 else c for i, c in enumerate(s)).replace("I", "L")
            n += 1
            agree += s == sage[int(sc)][0]
            if sage[int(sc)][2]:
                n_hi += 1
                agree_hi += s == sage[int(sc)][0]
            dh.append(h - sage[int(sc)][1])
        out["sage_oracle"] = {"n_common_scans": n, "same_top1_peptide": agree / max(n, 1),
                              "n_sage_target_psm_q_le_1pct": n_hi, "same_top1_peptide_sage_q_le_1pct": agree_hi / max(n_hi, 1),
                              "median_hyperscore_diff": float(np.median(dh)) if dh else None}
    os.makedirs(C.RES, exist_ok=True)
    json.dump(out, open(os.path.join(C.RES, "check_da.json"), "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
