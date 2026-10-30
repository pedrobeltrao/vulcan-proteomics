#!/usr/bin/env python3
"""E7b cohort equivalence (dev only).

mode reference: a chimeric cohort run with --keep-stores (OUT, STAGE) vs the per-run REFERENCE pipeline on the same
  stores, built only from the pre-E7b components: vulcan.ms1 with the original kernels (exact_reference), the chimeric
  search with a fresh closed primary (engine.search, no worker-prepared hypotheses), FeatureBuilderV2.run_v2 +
  PredCache on the same library root, the frozen v2 model, chimeric_competition (q_sep) and picked protein groups.
  Compared per run: MS1 hypothesis table (cohort worker file vs reference), PSM table (all columns but `score`),
  `score` bitwise, the peptide table written by the post worker (pf_id, spec_idx, score, is_decoy, q), target IDs at
  q <= 0.1/1/5%, protein groups at 1%.
mode copies: every run of a large cohort OUT vs its base run in a reference cohort REF (name <base> or
  <base>__r<k>): PSM table incl. score, peptide table and protein table identical.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
ALPHAS = (0.001, 0.01, 0.05)


def pq_np(path, drop=()):
    import pyarrow.parquet as pq
    t = pq.read_table(path)
    return {c: t.column(c).to_numpy(zero_copy_only=False) for c in t.column_names if c not in drop}


def same(a: dict, b: dict) -> list:
    """Columns that differ (bitwise for numeric, exact for objects)."""
    bad = [c for c in set(a) ^ set(b)]
    for c in set(a) & set(b):
        x, y = a[c], b[c]
        if x.dtype != y.dtype or len(x) != len(y):
            bad.append(c)
        elif x.dtype.kind in "fiub":
            if x.tobytes() != y.tobytes():
                bad.append(c)
        elif not np.array_equal(x, y):
            bad.append(c)
    return sorted(bad)


def reference(a):
    from vulcan import hw
    hw.pin_thread_pools(8)
    import pyarrow.parquet as pq
    import torch

    from bench.t1.ours_chimeric import MODEL_V2
    from vulcan.chimeric import search as CS
    from vulcan.chimeric.compete import IonSets, chimeric_competition
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ms1.features import FeatureParams, find_features
    from vulcan.ms1.hypotheses import build_hypotheses
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore import post as P
    from vulcan.rescore.features_v2 import FeatureBuilderV2, model_psm
    from vulcan.rescore.static import StaticModel
    v2 = StaticModel.load(MODEL_V2)
    pred = APD()
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PredCache(a.lib, pred, nce, "QE", STATIC, VAR)
        return caches[nce]
    fb = FeatureBuilderV2(factory, STATIC, VAR)
    p = E.SearchParams(report_psms=5)
    l1 = L1Index(a.index)
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    ctx = P.IndexCtx(a.index)
    lab = ctx.lab
    out = {"mode": "reference", "cohort": a.out, "index": a.index, "lib": a.lib, "runs": {}}
    for f in sorted(os.listdir(a.stage)):
        if not f.endswith(".vspec.arrow"):
            continue
        run = f[:-len(".vspec.arrow")]
        t0 = time.time()
        store = os.path.join(a.stage, f)
        sp = E.load_store(store)
        fr = find_features(store, FeatureParams(), exact_reference=True)
        hyp, _ = build_hypotheses(store, fr)
        hyp_c = CS.read_hypotheses(os.path.join(a.out, run + ".hypotheses.arrow"))
        hyp_eq = same({c: hyp.column(c).to_numpy(zero_copy_only=False) for c in hyp.column_names},
                      {c: hyp_c.column(c).to_numpy(zero_copy_only=False) for c in hyp_c.column_names})
        tbl, _ = CS.search_chimeric(sp, gi, l1, p, CS.precursor_groups(hyp, sp))
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
        coh = pq_np(os.path.join(a.out, run + ".chim.psms.parquet"))
        psm_bad = same(psm, {k: v for k, v in coh.items() if k != "score"})
        ff = fr.features
        mf = {c: ff.column(c).to_numpy() for c in ("trace_cos", "n_isotopes", "rt_apex", "rt_start", "rt_end",
                                                  "intensity_apex")}
        feat, info = fb.run_v2(l1, psm, sp, mf)
        s = np.asarray(v2.apply(model_psm(psm), feat), np.float64)
        score_eq = s.tobytes() == coh["score"].tobytes()
        row = psm["pf_row"]
        cc = chimeric_competition(psm["spec_idx"].astype(np.int64), psm["group"], psm["chim_rank"], psm["pf_id"],
                                  lab["tdc_pair"][row], s, lab["is_decoy"][row], row, IonSets(l1), seed=11)

        class R:
            pass
        r = R()
        r.row, r.score, r.is_decoy, r.q = cc["win"], cc["res"].score, cc["res"].is_decoy, cc["q_sep"]
        res, prot = P.protein_groups(ctx, row, psm["pf_id"], r)
        pep_ref = {"pf_id": psm["pf_id"][r.row], "spec_idx": psm["spec_idx"][r.row], "score": r.score,
                   "is_decoy": r.is_decoy, "q": r.q}
        pep = pq_np(os.path.join(a.out, run + ".peptides.parquet"))
        pep_bad = same(pep_ref, {k: pep[k] for k in pep_ref})
        ids = {}
        for al in ALPHAS:
            x = set(pep_ref["pf_id"][(pep_ref["q"] <= al) & ~pep_ref["is_decoy"]].tolist())
            y = set(pep["pf_id"][(pep["q"] <= al) & ~pep["is_decoy"]].tolist())
            ids[str(al)] = {"ref": len(x), "cohort": len(y), "identical": x == y}
        pr = pq_np(os.path.join(a.out, run + ".proteins.parquet"))
        n_pc = int(((pr["q"] <= 0.01) & ~pr["is_decoy"] & ~pr["is_obsolete"]).sum())
        rec = {"hypotheses_differing_columns": hyp_eq, "psm_differing_columns": psm_bad, "score_bitwise": score_eq,
               "peptides_differing_columns": pep_bad, "ids": ids, "proteins_1pct": [prot[0.01]["n_groups"], n_pc],
               "nce": info["nce"], "wall_s": time.time() - t0}
        rec["ok"] = (not hyp_eq and not psm_bad and score_eq and not pep_bad and all(v["identical"] for v in ids.values())
                     and prot[0.01]["n_groups"] == n_pc)
        out["runs"][run] = rec
        print(run[:30], rec["ok"], hyp_eq, psm_bad, score_eq, pep_bad, ids["0.01"], rec["proteins_1pct"], flush=True)
        torch.cuda.empty_cache()
    out["all_ok"] = all(r["ok"] for r in out["runs"].values())
    return out


def copies(a):
    out = {"mode": "copies", "cohort": a.out, "reference_cohort": a.ref, "runs": {}}
    for f in sorted(os.listdir(a.out)):
        if not f.endswith(".chim.psms.parquet"):
            continue
        run = f[:-len(".chim.psms.parquet")]
        base = re.sub(r"_r\d\d$", "", run)
        rec = {"base": base}
        for suf in (".chim.psms.parquet", ".peptides.parquet", ".proteins.parquet"):
            x = pq_np(os.path.join(a.out, run + suf), drop=("run",))     # the run name differs by construction
            y = pq_np(os.path.join(a.ref, base + suf), drop=("run",))
            rec[suf] = same(x, y)
        rec["ok"] = not any(v for k, v in rec.items() if k.startswith("."))
        out["runs"][run] = rec
    out["n_runs"] = len(out["runs"])
    out["n_bases"] = len({r["base"] for r in out["runs"].values()})
    out["all_ok"] = all(r["ok"] for r in out["runs"].values())
    print("copies", out["n_runs"], "bases", out["n_bases"], "all_ok", out["all_ok"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["reference", "copies"])
    ap.add_argument("--out", required=True, help="cohort output dir")
    ap.add_argument("--stage", help="kept stores of the cohort (reference mode)")
    ap.add_argument("--index", default="/tmp/vulcan_e2/index/human_contam_prod")
    ap.add_argument("--lib", default="/tmp/vulcan_ss/predlib")
    ap.add_argument("--ref", help="reference cohort output dir (copies mode)")
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    out = reference(a) if a.mode == "reference" else copies(a)
    json.dump(out, open(a.json, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
