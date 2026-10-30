#!/usr/bin/env python3
"""BEFORE profile: the frozen production scoring path (M3e bench_prod configuration) with per-sub-stage timers.

Same path and cache state as scripts/m3e/bench_prod.py (persistent process; per-run label-blind NCE choice; fresh
"first-contact" cache root holding only the shipped NCE-29 human library; D0 run 1 warm-up excluded), with the
library code monkey-patched ONLY to add cuda-synchronized timers (no behaviour change).  Sub-stages:
  load_store, search, peaks (SpecPeaks), nce.{fragment_mz,observed,lookup,predict_ms2,fragment_features},
  main.{fragment_mz,observed,lookup_ms2,predict_ms2,fragment_features,lookup_rt,predict_rt,rt_calibrate},
  library_load (I/O), feature_matrix, gbm_decision, double_competition, protein_maps, picked_group_fdr, protein_glue.
    flock /tmp/vulcan_bench.lock /tmp/vulcan_env_*/pred/bin/python scripts/scoring_speed/profile_old.py --out JSON
"""
from __future__ import annotations

import time

T0 = time.perf_counter()
import argparse  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import sys  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_common as S  # noqa: E402

import numpy as np  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache-root", default=os.path.join(S.W, "predcache_firstcontact_old"))
    ap.add_argument("--only", default=None, help="comma list of tags to run (debug)")
    a = ap.parse_args()
    c0 = S.cpu()
    import torch
    import evaluate as EV
    import m3_data as MD
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.predict import apd as APDm
    from vulcan.predict import cache as PC
    from vulcan.rescore import features as FB
    from vulcan.rescore.static import StaticModel
    from vulcan.rescore import static as ST
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    one = {"python_imports": {"wall_s": time.perf_counter() - T0, "cpu_core_s": S.cpu() - c0}}
    acc = S.Acc(torch)
    phase = {"p": "main"}

    def wrap(mod, name, label):
        f = getattr(mod, name)

        def g(*x, **k):
            with acc(label() if callable(label) else label):
                return f(*x, **k)
        setattr(mod, name, g)
    # ---- timers (behaviour unchanged)
    wrap(FB, "fragment_mz", lambda: phase["p"] + ".fragment_mz")
    wrap(FB, "observed", lambda: phase["p"] + ".observed")
    wrap(FB, "fragment_features", lambda: phase["p"] + ".fragment_features")
    wrap(FB, "rt_calibrate", "main.rt_calibrate")
    wrap(FB.SpecPeaks, "__init__", "peaks")
    o_pred = APDm.APD.predict

    def pred(self, *x, **k):
        kind = "ms2" if k.get("want_ms2", True) else "rt"
        with acc(f"{phase['p']}.predict_{kind}"):
            return o_pred(self, *x, **k)
    APDm.APD.predict = pred
    o_lms2, o_lrt = PC.PredCache.lookup_ms2, PC.PredCache.lookup_rt

    def lms2(self, *x, **k):
        with acc(f"{phase['p']}.lookup_ms2_incl_predict"):
            return o_lms2(self, *x, **k)

    def lrt(self, *x, **k):
        with acc(f"{phase['p']}.lookup_rt_incl_predict"):
            return o_lrt(self, *x, **k)
    PC.PredCache.lookup_ms2, PC.PredCache.lookup_rt = lms2, lrt
    o_tab = PC.PredCache._tab

    def tab(self, kind):
        if kind in self.tables:
            return o_tab(self, kind)
        with acc("library_load"):
            return o_tab(self, kind)
    PC.PredCache._tab = tab
    o_choose = FB.FeatureBuilder.choose_nce

    def choose(self, *x, **k):
        phase["p"] = "nce"
        try:
            with acc("nce.total"):
                return o_choose(self, *x, **k)
        finally:
            phase["p"] = "main"
    FB.FeatureBuilder.choose_nce = choose
    wrap(ST, "feature_matrix", "feature_matrix")
    o_pm, o_pg = EV.protein_maps, EV.picked_group_fdr

    def pm(*x, **k):
        with acc("protein_maps"):
            return o_pm(*x, **k)

    def pg(*x, **k):
        with acc("picked_group_fdr"):
            return o_pg(*x, **k)
    EV.protein_maps, EV.picked_group_fdr = pm, pg

    root = a.cache_root
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root)
    shutil.copytree(os.path.join(S.REF_CACHE, S.LIB29_KEY), os.path.join(root, S.LIB29_KEY))
    with acc("model_load"):
        pred = APDm.APD()
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PC.PredCache(root, pred, nce, "QE", S.C.STATIC, S.C.VAR)
        return caches[nce]
    p = E.SearchParams(report_psms=S.C.REPORT)
    idx = {}
    with acc("index_load_upload"):
        for fam in ("human_contam",):
            d = S.C.index_dir(fam, "prod")
            l1 = L1Index(d)
            idx[fam] = (d, l1, E.GpuIndex(l1, p, S.C.STATIC, S.C.VAR), MD.labels(d))
    with acc("library_load_nce29_startup"):
        c29 = factory(29.0)
        for kind in ("rt", "ms2_z2", "ms2_z3"):
            c29._tab(kind)
    one.update(acc.rec)
    fb = FB.FeatureBuilder(factory, S.C.STATIC, S.C.VAR)
    model = StaticModel.load(S.MODEL)
    runs = {}
    only = set(a.only.split(",")) if a.only else None
    for tag, ds, path in S.plan():
        if only and tag not in only:
            continue
        fam = S.family(ds)
        if fam not in idx:
            with acc("index_load_upload_pfu"):
                d = S.C.index_dir(fam, "prod")
                l1 = L1Index(d)
                idx[fam] = (d, l1, E.GpuIndex(l1, p, S.C.STATIC, S.C.VAR), MD.labels(d))
            one["index_load_upload_pfu"] = acc.rec.pop("index_load_upload_pfu")
        d, l1, gi, lab = idx[fam]
        acc.reset()
        t0, cc0 = time.perf_counter(), S.cpu()
        with acc("load_store"):
            sp = E.load_store(path)
        with acc("search"):
            tbl, st, _ = E.search(sp, gi, l1, p)
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        ts, cs = time.perf_counter(), S.cpu()
        with acc("features_total"):
            feat, info = fb.run(l1, psm, sp)
        with acc("gbm_apply_total"):
            s = model.apply(psm, feat)
        with acc("double_competition"):
            spec = psm["spec_idx"].astype(np.int64)
            r = double_competition(spec, psm["pf_id"], lab["tdc_pair"][psm["pf_row"]], s, lab["is_decoy"][psm["pf_row"]],
                                   seed=11, tiebreak=hash_tiebreak(psm["pf_id"], 11))
        with acc("protein_total"):
            srch = type("S", (), {})()
            srch.psm, srch.lab, srch.index, srch.unit = psm, lab, d, psm["pf_id"]
            prot = EV.protein_level(srch, {"row": r.row, "score": r.score, "is_decoy": r.is_decoy, "q": r.q}, (0.01,))
        scoring = {"wall_s": time.perf_counter() - ts, "cpu_core_s": S.cpu() - cs}
        tot = {"wall_s": time.perf_counter() - t0, "cpu_core_s": S.cpu() - cc0}
        rec = dict(acc.rec)
        # exclusive sub-stages
        for ph in ("nce", "main"):
            for kind in ("ms2", "rt"):
                inc = rec.get(f"{ph}.lookup_{kind}_incl_predict")
                if inc:
                    pr = rec.get(f"{ph}.predict_{kind}", {"wall_s": 0.0, "cpu_core_s": 0.0})
                    rec[f"{ph}.lookup_{kind}_excl_predict"] = {m: inc[m] - pr[m] for m in ("wall_s", "cpu_core_s")}
        g = rec["gbm_apply_total"]
        rec["gbm_decision"] = {m: g[m] - rec["feature_matrix"][m] for m in ("wall_s", "cpu_core_s")}
        pt = rec["protein_total"]
        rec["protein_glue"] = {m: pt[m] - rec["protein_maps"][m] - rec["picked_group_fdr"][m] for m in ("wall_s", "cpu_core_s")}
        acc1 = (r.q <= 0.01) & ~r.is_decoy
        name = S.run_name(path)
        runs[name] = {"tag": tag, "ds": ds, "n_ms2": len(sp), "n_psms": tbl.num_rows, "nce": info["nce"],
                      "cache_stats": info["cache_stats"], "sub_stages": rec, "scoring_total": scoring,
                      "end_to_end_total": tot, "ids_1pct": int(acc1.sum()), "proteins_1pct": prot[0.01]["n_groups"]}
        print(tag, name[:30], len(sp), "NCE", info["nce"], "scoring", round(scoring["wall_s"], 2),
              round(scoring["cpu_core_s"], 2), flush=True)
    from vulcan.hw import fingerprint
    json.dump({"one_time": one, "runs": runs, "hw": fingerprint(), "loadavg_end": os.getloadavg(),
               "cache_root": root}, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
