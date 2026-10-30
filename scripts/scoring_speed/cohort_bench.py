#!/usr/bin/env python3
"""End-to-end cohort throughput RAW -> PSMs -> scored peptides -> picked protein groups (vulcan.cohort).

Configs (one `python -m vulcan.cohort` process each, every one under `flock /tmp/vulcan_bench.lock`):
  noscore   search only (the M2b driver, for reference)
  frozen    + scoring with the frozen pre-optimization components (FeatureBuilder + PredCache in the M3e
            "first-contact" state, sklearn GBM, competition/proteins inline in the GPU thread)
  fast      + production scoring (vulcan.rescore.cohort_scoring: resident consolidated library, GPU features + GBM,
            competition/proteins in post worker processes)
RAWs: the 10 D8 holdout runs + the 2 D0 runs (human_contam_prod index), cycled to --repeat-to runs; D1 (pfu index) as a
separate small cohort (fast only).  Writes cohort_<config>.json (the cohort stats) and a summary JSON.
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_common as S  # noqa: E402

PY = sys.executable


def raws_human():
    # D0 RAWs copied once to local NVMe (as the D8 ones are): the Lustre path made 3.4 GB decodes time out under load
    d0 = []
    for r in S.C.DEV["D0"]["raw"]:
        src = os.path.join(S.ROOT, "data", "raw", r)
        dst = os.path.join(S.W, "raw_d0", os.path.basename(r))
        if not os.path.exists(dst) or os.path.getsize(dst) != os.path.getsize(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(src, dst)
        d0.append(dst)
    return sorted(glob.glob(os.path.join(S.M3E_W, "raw_d8", "*.raw"))) + d0


def run(cfg, raws, index, repeat_to, workers, lib, post_workers, outdir, log, search_threads=1):
    out = os.path.join(S.W, "cohort_out", cfg)
    stage = os.path.join(S.W, "cohort_stage", cfg)
    shutil.rmtree(out, ignore_errors=True)
    shutil.rmtree(stage, ignore_errors=True)
    lst = os.path.join(S.W, f"rawlist_{cfg}.txt")
    open(lst, "w").write("\n".join(raws) + "\n")
    cmd = ["flock", "/tmp/vulcan_bench.lock", PY, "-m", "vulcan.cohort", "--index", index, "--raw-list", lst,
           "--out", out, "--stage", stage, "--workers", str(workers), "--stats", os.path.join(outdir, f"cohort_{cfg}.json"),
           "--params", json.dumps({"report_psms": S.C.REPORT}),         # the M3 search configuration
           "--search-threads", str(search_threads)]
    if repeat_to:
        cmd += ["--repeat-to", str(repeat_to)]
    if cfg.startswith("frozen"):
        root = os.path.join(S.W, "predcache_firstcontact_cohort")
        shutil.rmtree(root, ignore_errors=True)
        os.makedirs(root)
        shutil.copytree(os.path.join(S.REF_CACHE, S.LIB29_KEY), os.path.join(root, S.LIB29_KEY))
        cmd += ["--score-lib", root, "--score-impl", "frozen"]
    elif cfg.startswith("fast"):
        cmd += ["--score-lib", lib, "--score-post-workers", str(post_workers),
                "--score-placement", "host" if "_host" in cfg else ("mmap" if "_mmap" in cfg else "gpu")]
    env = dict(os.environ, PYTHONPATH=S.ROOT, PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True",
               VULCAN_COHORT_TIMEOUT_SCALE=os.environ.get("VULCAN_COHORT_TIMEOUT_SCALE", "3"))   # shared node
    t0 = time.time()
    with open(os.path.join(S.W, f"cohort_{cfg}.log"), "w") as fh:
        subprocess.run(cmd, check=True, cwd=S.ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT)
    st = json.load(open(os.path.join(outdir, f"cohort_{cfg}.json")))
    rec = {k: st[k] for k in ("n_runs", "workers", "makespan_s", "runs_per_hour", "runs_per_hour_steady", "bottleneck",
                              "stage_busy_fraction", "core_s", "gpu_s", "setup")}
    rec["wall_incl_process_start_s"] = time.time() - t0
    rec["process_age_at_end_s"] = st["process"]["age_at_end_s"]
    rec["runs_per_hour_incl_startup"] = st["n_runs"] / st["process"]["age_at_end_s"] * 3600
    rec["loadavg_end"] = os.getloadavg()
    if "scoring" in st:
        rec["scoring"] = st["scoring"]
        per = {}
        for r in st["runs"]:
            base = r["run"].rsplit("_r", 1)[0] if r["run"][-4:-2] == "_r" else r["run"]
            per.setdefault(base, set()).add((r["post"]["ids_1pct"], r["post"]["proteins_1pct"]))
        rec["ids_proteins_by_base_run"] = {k: sorted(v) for k, v in per.items()}
        g = [r["score_gpu_stage_wall_s"] for r in st["runs"]]
        rec["score_gpu_stage_wall_s_median"] = sorted(g)[len(g) // 2]
    rec["vs_reference"] = check_outputs(out, cfg)
    log(cfg, json.dumps({k: rec[k] for k in ("n_runs", "makespan_s", "runs_per_hour", "runs_per_hour_steady",
                                             "bottleneck")}, default=float))
    return rec


def check_outputs(out, cfg):
    """First copy of every base run: PSM table vs the reference PSM table (M3/M3e, same index + params), and
    q<=1% peptide set / protein count vs the equality record."""
    import numpy as np
    import pyarrow.parquet as pq
    import m3_data as MD
    from vulcan.rescore import post as P
    from vulcan.rescore.static import StaticModel
    eq = json.load(open(os.path.join(S.RES, "equality.json")))["runs"]
    model = StaticModel.load(S.MODEL)
    ctxs = {}
    res = {}
    for name, r in eq.items():
        f = os.path.join(out, name + ".psms.parquet")
        if not os.path.exists(f):
            continue
        a = pq.read_table(f)
        b = pq.read_table(os.path.join(S.ref_psm_dir(r["ds"]), name + ".psms.parquet"))
        cols = [c for c in b.column_names if c != "run"]
        same = a.num_rows == b.num_rows and all(np.array_equal(a.column(c).to_numpy(zero_copy_only=False),
                                                               b.column(c).to_numpy(zero_copy_only=False)) for c in cols)
        rec = {"psm_table_identical_to_reference": bool(same)}
        pf = os.path.join(out, name + ".peptides.parquet")
        if os.path.exists(pf):
            w = pq.read_table(pf)
            acc = (w.column("q").to_numpy() <= 0.01) & ~w.column("is_decoy").to_numpy()
            rec["ids_1pct"] = int(acc.sum())
            rec["ids_1pct_equal_reference"] = rec["ids_1pct"] == r["ids_q0.01"]["new"]
            pr = pq.read_table(os.path.join(out, name + ".proteins.parquet"))
            pa_ = (pr.column("q").to_numpy() <= 0.01) & ~pr.column("is_decoy").to_numpy() & ~pr.column("is_obsolete").to_numpy()
            rec["proteins_1pct"] = int(pa_.sum())
            rec["proteins_1pct_equal_reference"] = rec["proteins_1pct"] == r["proteins_1pct"]["new"]
            if "score" in a.column_names and same:
                d = S.ref_psm_dir(r["ds"])
                psm = {c: b.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
                ft = pq.read_table(os.path.join(d, name + ".feat.parquet"))
                s_ref = model.apply(psm, {c: ft.column(c).to_numpy() for c in ft.column_names})
                rec["score_bitwise_equal_frozen_reference"] = bool(np.array_equal(a.column("score").to_numpy(), s_ref))
                fam = S.family(r["ds"])
                if fam not in ctxs:
                    ctxs[fam] = P.IndexCtx(S.C.index_dir(fam, "prod"))
                rr = P.compete(ctxs[fam], psm["pf_row"], psm["pf_id"], psm["spec_idx"], s_ref)
                A = set(psm["pf_id"][rr.row[(rr.q <= 0.01) & ~rr.is_decoy]].tolist())
                B = set(w.column("pf_id").to_numpy()[acc].tolist())
                rec["ids_1pct_set_identical_to_reference"] = A == B
        res[name] = rec
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="noscore,frozen,fast")
    ap.add_argument("--repeat-to", type=int, default=48)
    ap.add_argument("--repeat-to-frozen", type=int, default=24, help="shorter: the frozen config is GPU-thread bound")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--post-workers", type=int, default=3)
    ap.add_argument("--lib", default=os.path.join(S.W, "predlib"))
    ap.add_argument("--outdir", default=S.RES)
    ap.add_argument("--d1", action="store_true")
    a = ap.parse_args()
    res = {}
    summ = os.path.join(a.outdir, "cohort_summary.json")
    if os.path.exists(summ):
        res = json.load(open(summ))
    human = S.C.index_dir("human_contam", "prod")
    for cfg in a.configs.split(","):
        rt = a.repeat_to_frozen if cfg.startswith("frozen") else a.repeat_to
        st = 2 if cfg.endswith("_st2") else 1
        res[cfg] = run(cfg, raws_human(), human, rt, a.workers, a.lib, a.post_workers, a.outdir, print, st)
        json.dump(res, open(summ, "w"), indent=1, default=float)
    if a.d1:
        raws = [os.path.join(S.ROOT, "data", "raw", r) for r in S.C.DEV["D1"]["raw"]]
        res["fast_D1"] = run("fast_D1", raws, S.C.index_dir("pfu_contam", "prod"), 4, 4, a.lib, 1, a.outdir, print)
        json.dump(res, open(summ, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
