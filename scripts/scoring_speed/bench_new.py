#!/usr/bin/env python3
"""AFTER: production scoring with vulcan.rescore.prod (fast path), timed per sub-stage, plus the equality record.

--mode bench     same loop as profile_old.py / scripts/m3e/bench_prod.py (persistent process; load store -> GPU
                 search -> features with the per-run label-blind NCE choice -> frozen GBM -> double competition ->
                 picked protein), on D0 run 1 (warm-up, excluded), 10 D8 runs, D0 run 2, D1.  --cache library: the
                 consolidated predicted library (build_library.py); --cache firstcontact: a fresh root holding only
                 the shipped NCE-29 library (the M3e state), other NCEs filled lazily.  Run under the bench lock.
--mode gpuprof   one pass with torch.profiler: CUDA kernel seconds per run (GPU-s), not wall-timed.
--mode equality  for every D0/D1/D8 run, on the REFERENCE PSM tables (M3 / M3e psms/*_prod): new features vs the
                 stored reference features (feat.parquet) and vs the old FeatureBuilder on the same library;
                 new GBM (GPU) vs sklearn; winners / q / q<=1% ID sets / protein groups, new vs reference.
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


def setup(a, acc):
    import torch
    import m3_data as MD
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.rescore.prod import ProdScorer
    root = a.lib_root
    if a.cache == "firstcontact":
        from vulcan.predict.library import consolidate
        root = os.path.join(S.W, "predcache_firstcontact_new")
        shutil.rmtree(root, ignore_errors=True)
        os.makedirs(root)
        shutil.copytree(os.path.join(S.REF_CACHE, S.LIB29_KEY), os.path.join(root, S.LIB29_KEY))
        consolidate(os.path.join(root, S.LIB29_KEY), log=lambda *x: None)   # shipped library = consolidated form
    with acc("predictor_and_model_load"):
        sc = ProdScorer(root, S.MODEL, "QE", S.C.STATIC, S.C.VAR, placement=a.placement)
    with acc("library_open_warm"):
        sc.warm() if a.cache == "library" else sc.warm(nces=[29.0], kinds=("rt", "ms2_z2", "ms2_z3"))
    p = E.SearchParams(report_psms=S.C.REPORT)
    idx = {}

    def index(fam):
        if fam not in idx:
            d = S.C.index_dir(fam, "prod")
            with acc(f"index_load_upload_{fam}"):
                l1 = L1Index(d)
                idx[fam] = (d, l1, E.GpuIndex(l1, p, S.C.STATIC, S.C.VAR), sc.index_ctx(d))
        return idx[fam]
    index("human_contam")
    return torch, E, MD, sc, p, index, root


def bench(a):
    acc = S.Acc()
    c0 = S.cpu()
    import torch
    acc.torch = torch
    one = {"python_imports": {"wall_s": time.perf_counter() - T0, "cpu_core_s": S.cpu() - c0}}
    torch, E, MD, sc, p, index, root = setup(a, acc)
    from vulcan.rescore import post as P
    one.update(acc.rec)
    runs = {}
    for tag, ds, path in S.plan():
        d, l1, gi, ctx = index(S.family(ds))
        if f"index_load_upload_{S.family(ds)}" in acc.rec and S.family(ds) != "human_contam":
            one[f"index_load_upload_{S.family(ds)}"] = acc.rec.pop(f"index_load_upload_{S.family(ds)}")
        acc.reset()
        t0, cc0 = time.perf_counter(), S.cpu()
        with acc("load_store"):
            sp = E.load_store(path)
        with acc("search"):
            tbl, st, _ = E.search(sp, gi, l1, p)
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        ts, cs = time.perf_counter(), S.cpu()
        with acc("features_total"):
            feat, info = sc.fb.run(l1, psm, sp)
        with acc("gbm_apply_total"):
            s = sc.model.apply(psm, feat)
        with acc("double_competition"):
            r = P.compete(ctx, psm["pf_row"], psm["pf_id"], psm["spec_idx"], s)
        with acc("protein_total"):
            _, prot = P.protein_groups(ctx, psm["pf_row"], psm["pf_id"], r)
        scoring = {"wall_s": time.perf_counter() - ts, "cpu_core_s": S.cpu() - cs}
        tot = {"wall_s": time.perf_counter() - t0, "cpu_core_s": S.cpu() - cc0}
        rec = dict(acc.rec)
        for k, v in info["timing_s"].items():
            rec[f"features.{k}"] = {"wall_s": v}
        name = S.run_name(path)
        runs[name] = {"tag": tag, "ds": ds, "n_ms2": len(sp), "n_psms": tbl.num_rows, "nce": info["nce"],
                      "cache_stats": info["cache_stats"], "sub_stages": rec, "scoring_total": scoring,
                      "end_to_end_total": tot, "ids_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum()),
                      "proteins_1pct": prot[0.01]["n_groups"]}
        print(tag, name[:30], len(sp), "NCE", info["nce"], "scoring", round(scoring["wall_s"], 3),
              round(scoring["cpu_core_s"], 3), "ids", runs[name]["ids_1pct"], flush=True)
    from vulcan.hw import fingerprint
    from vulcan.predict.library import library_digest
    import resource
    json.dump({"one_time": one, "runs": runs, "hw": fingerprint(), "loadavg_end": os.getloadavg(), "cache": a.cache,
               "library_placement": a.placement, "library_resident_bytes": getattr(sc, "library_bytes", None),
               "max_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
               "gpu_max_mem_bytes": torch.cuda.max_memory_allocated(),
               "cache_root": root, "library_digest": library_digest(root, os.listdir(root))},
              open(a.out, "w"), indent=1, default=float)


def gpuprof(a):
    acc = S.Acc()
    torch, E, MD, sc, p, index, root = setup(a, acc)
    from torch.profiler import ProfilerActivity, profile
    fb, model = sc.fb, sc.model
    if a.impl == "frozen":                     # BEFORE: frozen components on a fresh first-contact cache
        from vulcan.predict.cache import PredCache
        from vulcan.rescore.features import FeatureBuilder
        froot = os.path.join(S.W, "predcache_firstcontact_gpuprof")
        shutil.rmtree(froot, ignore_errors=True)
        os.makedirs(froot)
        shutil.copytree(os.path.join(S.REF_CACHE, S.LIB29_KEY), os.path.join(froot, S.LIB29_KEY))
        oc = {}

        def of(nce):
            if nce not in oc:
                oc[nce] = PredCache(froot, sc.pred, nce, "QE", S.C.STATIC, S.C.VAR)
            return oc[nce]
        fb, model = FeatureBuilder(of, S.C.STATIC, S.C.VAR), sc.model_static
    out = {}
    for tag, ds, path in S.plan():
        d, l1, gi, ctx = index(S.family(ds))
        sp = E.load_store(path)
        res = {}
        with profile(activities=[ProfilerActivity.CUDA]) as pr:
            tbl, st, _ = E.search(sp, gi, l1, p)
            torch.cuda.synchronize()
        res["search_gpu_s"] = _kernel_s(pr)
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        with profile(activities=[ProfilerActivity.CUDA]) as pr:
            feat, info = fb.run(l1, psm, sp)
            torch.cuda.synchronize()
        res["features_gpu_s"] = _kernel_s(pr)
        with profile(activities=[ProfilerActivity.CUDA]) as pr:
            model.apply(psm, feat)
            torch.cuda.synchronize()
        res["gbm_gpu_s"] = _kernel_s(pr)
        res["scoring_gpu_s"] = res["features_gpu_s"] + res["gbm_gpu_s"]
        res["tag"] = tag
        res["predicted_ms2"] = info["cache_stats"]["n_predicted_ms2"]
        out[S.run_name(path)] = res
        print(tag, S.run_name(path)[:30], {k: round(v, 4) for k, v in res.items() if isinstance(v, float)}, flush=True)
    json.dump({"runs": out, "cache": a.cache, "impl": a.impl}, open(a.out, "w"), indent=1)


def _kernel_s(pr):
    tot = 0.0
    for e in pr.events():
        if e.device_type.name == "CUDA":
            tot += e.device_time_total if hasattr(e, "device_time_total") else e.cuda_time_total
    return tot / 1e6


def equality(a):
    import pyarrow.parquet as pq
    acc = S.Acc()
    torch, E, MD, sc, p, index, root = setup(a, acc)
    import evaluate as EV
    from vulcan.predict.cache import PredCache
    from vulcan.rescore import post as P
    from vulcan.rescore.features import FeatureBuilder
    old_caches = {}

    def old_factory(nce):
        if nce not in old_caches:
            old_caches[nce] = PredCache(root, sc.pred, nce, "QE", S.C.STATIC, S.C.VAR)
        return old_caches[nce]
    old_fb = FeatureBuilder(old_factory, S.C.STATIC, S.C.VAR)
    seen = set()
    out = {}
    for tag, ds, path in S.plan():
        name = S.run_name(path)
        if name in seen:
            continue
        seen.add(name)
        d, l1, gi, ctx = index(S.family(ds))
        refd = S.ref_psm_dir(ds)
        t = pq.read_table(os.path.join(refd, name + ".psms.parquet"))
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        ft = pq.read_table(os.path.join(refd, name + ".feat.parquet"))
        ref = {c: ft.column(c).to_numpy() for c in ft.column_names}
        sp = E.load_store(path)
        # search reproducibility (same store, index, params) vs the reference PSM table
        tbl, _, _ = E.search(sp, gi, l1, p)
        srch = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS if c != "run"}
        search_eq = all(np.array_equal(srch[c], psm[c]) for c in srch)
        new, info = sc.fb.run(l1, psm, sp)
        old, info_old = old_fb.run(l1, psm, sp)
        rec = {"ds": ds, "n_psms": len(psm["pf_row"]), "search_reproduces_reference_psms": bool(search_eq),
               "nce_new": info["nce"], "nce_old_code": info_old["nce"],
               "nce_reference": json.load(open(os.path.join(refd, "feat_info.json")))["runs"][name]["nce"],
               "new_cache_stats": info["cache_stats"], "features": {}}
        for k in ref:
            x, y, z = new[k], ref[k], old[k]
            dif = np.abs(x.astype(np.float64) - y.astype(np.float64))
            rec["features"][k] = {"bitwise_eq_reference": bool(np.array_equal(x.view(np.uint32), y.view(np.uint32))),
                                  "n_differ_reference": int((x.view(np.uint32) != y.view(np.uint32)).sum()),
                                  "max_abs_diff_reference": float(np.nanmax(dif)) if len(dif) else 0.0,
                                  "bitwise_eq_old_code_same_library": bool(np.array_equal(x.view(np.uint32),
                                                                                         z.view(np.uint32)))}
        s_ref = sc.model_static.apply(psm, ref)                # frozen path: sklearn on the reference features
        s_new = sc.model.apply(psm, new)                       # fast path: GPU GBM on the new features
        s_gpu_ref = sc.model.apply(psm, ref)                   # GPU GBM on the reference features
        rec["gbm_gpu_vs_sklearn_same_features_bitwise"] = bool(np.array_equal(s_gpu_ref, s_ref))
        rec["score_bitwise_eq_reference"] = bool(np.array_equal(s_new, s_ref))
        rec["score_max_abs_diff_reference"] = float(np.abs(s_new - s_ref).max())
        r_ref = P.compete(ctx, psm["pf_row"], psm["pf_id"], psm["spec_idx"], s_ref)
        r_new = P.compete(ctx, psm["pf_row"], psm["pf_id"], psm["spec_idx"], s_new)
        for al in (0.001, 0.01, 0.05):
            A = set(psm["pf_id"][r_ref.row[(r_ref.q <= al) & ~r_ref.is_decoy]].tolist())
            B = set(psm["pf_id"][r_new.row[(r_new.q <= al) & ~r_new.is_decoy]].tolist())
            rec[f"ids_q{al}"] = {"reference": len(A), "new": len(B), "identical": A == B,
                                 "only_ref": len(A - B), "only_new": len(B - A)}
        rec["winners_identical"] = bool(np.array_equal(r_ref.row, r_new.row) and np.array_equal(r_ref.q, r_new.q))
        # protein glue: new post.protein_groups vs the frozen evaluate.protein_level on the same winners
        sr = type("S", (), {})()
        sr.psm, sr.lab, sr.index, sr.unit = psm, ctx.lab, d, psm["pf_id"]
        w = {"row": r_ref.row, "score": r_ref.score, "is_decoy": r_ref.is_decoy, "q": r_ref.q}
        old_prot = EV.protein_level(sr, w, (0.01,))
        res_new, new_prot = P.protein_groups(ctx, psm["pf_row"], psm["pf_id"], r_ref)
        rec["protein_glue_identical_n_groups_1pct"] = old_prot[0.01] == new_prot[0.01]
        res_new2, new_prot2 = P.protein_groups(ctx, psm["pf_row"], psm["pf_id"], r_new)
        rec["proteins_1pct"] = {"reference": old_prot[0.01]["n_groups"], "new": new_prot2[0.01]["n_groups"]}
        out[name] = rec
        nfeat_eq = sum(v["bitwise_eq_reference"] for v in rec["features"].values())
        print(ds, name[:28], "search_eq", search_eq, "feat_eq_ref", f"{nfeat_eq}/{len(rec['features'])}",
              "old_code_eq", all(v["bitwise_eq_old_code_same_library"] for v in rec["features"].values()),
              "score_eq", rec["score_bitwise_eq_reference"], "gbm_eq", rec["gbm_gpu_vs_sklearn_same_features_bitwise"],
              "ids1%", rec["ids_q0.01"], "prot", rec["proteins_1pct"], rec["protein_glue_identical_n_groups_1pct"],
              flush=True)
    json.dump({"runs": out, "library_root": root}, open(a.out, "w"), indent=1, default=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["bench", "gpuprof", "equality"], required=True)
    ap.add_argument("--cache", choices=["library", "firstcontact"], default="library")
    ap.add_argument("--lib-root", default=os.path.join(S.W, "predlib"))
    ap.add_argument("--impl", choices=["fast", "frozen"], default="fast", help="gpuprof only")
    ap.add_argument("--placement", choices=["gpu", "host", "mmap"], default="gpu", help="predicted-library placement")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    {"bench": bench, "gpuprof": gpuprof, "equality": equality}[a.mode](a)


if __name__ == "__main__":
    main()
