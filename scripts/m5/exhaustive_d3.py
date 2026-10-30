#!/usr/bin/env python3
"""M5 non-cleavable search-space reduction on D3 (DSS): alpha-first (K1 top-K alpha, pair-closed, beta lookup) vs the
exhaustive pair space (every row pair whose masses + linker match the precursor at some isotope error).

    python scripts/m5/exhaustive_d3.py --outdir DIR [--topk 5 10 20 50] [--chunk 200]

Writes DIR/exhaustive.csm.parquet (exhaustive reference: every row taking part in a mass-matching pair scored by K2,
same pairing / ranking / linear-mono-loop competition / site localization) and DIR/alpha_k<K>.csm.parquet, plus
DIR/exhaustive_d3.json with search-space sizes (pairs and K2 candidates), wall times and label-free rank-1
agreement on the spectra where the exhaustive search has a CSM that beats its competitor. FDR / true FDP of each
table are computed afterwards by scripts/m5/run_m5.py (L4). No labels are read here.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--index", default="/tmp/vulcan_m5/index/D3")
    ap.add_argument("--store", default="/tmp/vulcan_m5/store/XLpeplib_Beveridge_QEx-HFX_DSS_R1.vspec.arrow")
    ap.add_argument("--topk", type=int, nargs="+", default=[5, 10, 20, 50])
    ap.add_argument("--chunk", type=int, default=200)
    a = ap.parse_args()
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch

    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ptm import search as PS
    from vulcan.xl import linkers as LK
    from vulcan.xl import search as XS
    from vulcan.xl import sites as S
    os.makedirs(a.outdir, exist_ok=True)
    STAT, VAR = {"C": 57.021464}, {"M": 15.9949}
    p = E.SearchParams(max_fragment_charge=3, max_peaks=200)
    l1 = L1Index(a.index)
    gi = E.GpuIndex(l1, p, STAT, VAR)
    sp = E.load_store(a.store)
    st = E.Stats()
    peaks = E.k0_preprocess(sp, p, gi.dev, st)
    xp0 = XS.XLParams()
    t0 = time.perf_counter()
    s_ok, n_pairs, n_rows = XS.exhaustive_pair_count(sp, gi, peaks, xp0, LK.DSS)
    t_count = time.perf_counter() - t0
    n_db = len(l1)
    out = {"n_index_rows": int(n_db), "n_xl_spectra": int(len(s_ok)),
           "exhaustive_pairs_total": int(n_pairs.sum()),
           "exhaustive_pairs_median_per_spectrum": float(np.median(n_pairs)),
           "exhaustive_rows_total": int(n_rows.sum()),
           "exhaustive_rows_median_per_spectrum": float(np.median(n_rows)),
           "all_unordered_pairs": float(n_db) * (n_db + 1) / 2, "count_time_s": t_count, "by_topk": {}}
    lin, _ = PS.search_modes(sp, gi, l1, p, PS.ModeParams(offsets=tuple(LK.DSS.mono_masses) + (LK.DSS.loop_mass,)),
                             report=1)

    def finish(tbl, name):
        m = (tbl.column("rank").to_numpy() == 1) & tbl.column("xl_wins").to_numpy(zero_copy_only=False)
        t1 = tbl.filter(pa.array(m))
        sa, sb, _, _ = S.localize_csms(t1, sp, l1, LK.DSS, STAT, VAR)
        t1 = t1.append_column("site_a", pa.array(sa)).append_column("site_b", pa.array(sb))
        pq.write_table(t1, os.path.join(a.outdir, f"{name}.csm.parquet"))
        return t1.to_pandas().set_index("spec_idx")

    # exhaustive reference in spectrum chunks (bounded GPU memory)
    t0 = time.perf_counter()
    parts, k2n = [], 0
    for i in range(0, len(s_ok), a.chunk):
        tb, est = XS.search_exhaustive_noncleavable(sp, gi, l1, p, xp0, LK.DSS, s_ok[i:i + a.chunk], lin=lin)
        parts.append(tb)
        k2n += est.counts["k2_scored"]
        torch.cuda.empty_cache()
    torch.cuda.synchronize()
    ex = pa.concat_tables(parts)
    out["exhaustive"] = {"k2_scored": int(k2n), "wall_s": time.perf_counter() - t0}
    exd = finish(ex, "exhaustive")
    print("exhaustive", out["exhaustive"], flush=True)
    for k in a.topk:
        xp = XS.XLParams(alpha_topk=k)
        t0 = time.perf_counter()
        tbl, _, st = XS.search(sp, gi, l1, p, xp, LK.DSS)
        torch.cuda.synchronize()
        w = time.perf_counter() - t0
        d1 = finish(tbl, f"alpha_k{k}")
        agree = sum(1 for s, e in exd.iterrows() if s in d1.index and
                    tuple(sorted((int(d1.loc[s].pf_id_a), int(d1.loc[s].pf_id_b)))) ==
                    tuple(sorted((int(e.pf_id_a), int(e.pf_id_b)))))
        # agreement among the top-1000 exhaustive CSMs by score (the region that can pass FDR)
        top = exd.sort_values("score", ascending=False).head(1000)
        agree_top = sum(1 for s, e in top.iterrows() if s in d1.index and
                        tuple(sorted((int(d1.loc[s].pf_id_a), int(d1.loc[s].pf_id_b)))) ==
                        tuple(sorted((int(e.pf_id_a), int(e.pf_id_b)))))
        out["by_topk"][str(k)] = {
            "alpha_shortlist": st.counts["alpha_shortlist"], "k2_candidates": st.counts["k2_window_candidates"],
            "pairs_considered": st.counts["n_pairs_considered"], "wall_s": w,
            "reduction_pairs_vs_exhaustive": float(n_pairs.sum()) / max(st.counts["n_pairs_considered"], 1),
            "reduction_k2_vs_exhaustive": float(k2n) / max(st.counts["k2_window_candidates"], 1),
            "rank1_agreement_all_exhaustive_csms": agree / max(len(exd), 1), "n_exhaustive_csms": int(len(exd)),
            "rank1_agreement_top1000_exhaustive": agree_top / max(len(top), 1),
            "wall_breakdown": st.wall_s}
        print(k, json.dumps({kk: v for kk, v in out["by_topk"][str(k)].items() if kk != "wall_breakdown"}), flush=True)
    from vulcan.hw import fingerprint
    out["hw"] = fingerprint()
    json.dump(out, open(os.path.join(a.outdir, "exhaustive_d3.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
