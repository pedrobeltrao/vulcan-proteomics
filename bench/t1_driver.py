#!/usr/bin/env python3
"""Pre-registered T1 evaluation driver (DESIGN §1 "no loss", §6.1-6.3, §6.0 Amendment A1; review B1/M4/M7).

One command per dataset:
    python bench/t1_driver.py --dataset D1                 # dev dry run (entrapment seeds 2001-2003)
    python bench/t1_driver.py --dataset D0
    python3 bench/sealed_run.py --dataset T1 --purpose "..." -- python bench/t1_driver.py --dataset T1   # sealed

Pipeline (each stage is idempotent: finished outputs under --work are reused, so a re-run resumes):
  db           one multi-replicate index (T, arm-C foreign F, E_s for 3 entrapment seeds, shuffled decoys of the 5
               DECOY_SEEDS) -> 25 replicate DBs (5 production, 15 arm P, 5 arm C) as index dirs + peptide FASTAs,
               §6.1 uniqueness verified (bench/t1/dbs.py)
  ours         vulcan.cohort ingest + production search, engine search of every replicate, L3 features + the frozen
               static_gbm_pred (A1.1) (bench/t1/ours.py)
  comparators  Sage (LDA, predict_rt on), Sage+mokapot, Sage+AlphaPeptDeep features+mokapot, Comet+mokapot on the
               identical replicate FASTAs (bench/t1/comparators.py)
  evaluate     IDs (peptidoform / sequence / picked protein), §6.3 verdicts under both A1.6 CI conventions, §1
               non-inferiority ratios with run x decoy-shuffle bootstrap CIs at own q and at matched FDP
               (bench/t1/evaluate.py)
Outputs: <out>/{metrics.json, tables.md, config.yaml, timing.json, git_sha}.

Sealed mode (--dataset T1): raw files come ONLY from $SEALED_FILES (set by bench/sealed_run.py), entrapment seeds are
ENTRAP_SEEDS_SEALED (2004-2006), and staged copies / stores / mzML (and, with A3, MS1 tables) of sealed runs are deleted
at the end unless --keep-sealed-derivatives. Dev datasets refuse sealed seeds and sealed paths.

Amendment A3 (DESIGN §6.0c; bench/t1/a3.py): `--a3` adds the pre-registered SECONDARY arm vulcan_v2_sep (stage `chim`)
and every A3.3.3 comparator not in the A1 set (sage_chim2, sage_chim2+pred+mokapot), computes the ratios for both
numerators (v1 primary first, unchanged) and writes `a3_secondary` (LOO backstop, per-merged-stratum verdicts, void
flags, Sage-rank split).  For the sealed T1 look A3 is forced (A3.3.8) and the driver refuses to run unless the tree is
clean and its frozen code paths equal config/a3.yaml's freeze (A3.1.3; vulcan.gates.a3_config.require_freeze_commit).
Without --a3 the driver does exactly what it did before A3 (dev).
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import resource
import shutil
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

DATASETS = {
    "D0": {"fasta": ["data/fasta/human_sp_canonical.fasta", "data/fasta/contaminants_frankenfield2022.fasta"],
           "foreign": ["data/fasta/arath_UP000006548.fasta"],
           "raws": ["data/raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw",
                    "data/raw/PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw"],
           "instrument": "QE", "prec_ppm": 10.0, "frag_ppm": 20.0, "role": "dev"},
    "D1": {"fasta": ["data/fasta/pfu_UP000001013.fasta", "data/fasta/contaminants_frankenfield2022.fasta"],
           "foreign": ["data/fasta/human_sp_canonical.fasta"],
           "raws": ["data/raw/PXD001077/Velos005137.raw"],
           "instrument": "QE", "prec_ppm": 10.0, "frag_ppm": 20.0, "role": "dev"},
    # D8h (E7/v2, 2026-10-28): the A1.1 D8 holdout = first 10 human runs of data/manifests/e9_d8.tsv (manifest order),
    # fetched + sha-verified by scripts/m3e/fetch_d8.py to node-local disk ($VULCAN_D8_RAW, default
    # /tmp/vulcan_m3e/raw_d8).  Dev role; never used for any selection.
    "D8h": {"fasta": ["data/fasta/human_sp_canonical.fasta", "data/fasta/contaminants_frankenfield2022.fasta"],
            "foreign": ["data/fasta/arath_UP000006548.fasta"],
            "raws": sorted(os.path.join(os.environ.get("VULCAN_D8_RAW", "/tmp/vulcan_m3e/raw_d8"), f)
                           for f in (os.listdir(os.environ.get("VULCAN_D8_RAW", "/tmp/vulcan_m3e/raw_d8"))
                                     if os.path.isdir(os.environ.get("VULCAN_D8_RAW", "/tmp/vulcan_m3e/raw_d8")) else [])
                           if f.endswith(".raw")),
            "instrument": "QE", "prec_ppm": 10.0, "frag_ppm": 20.0, "role": "dev"},
    # T1 = ProteoBench DDA-QExactive (Q Exactive HF-X, HYE).  FASTA: the dev copy of the ProteoBench HYE zip, which is
    # byte-identical to T1's (DATA.md D9, sha256 cae44ed2e3cb9ce5...), includes Cont_ contaminants.  Arm C:
    # Arabidopsis (absent from the sample).  Raws: $SEALED_FILES only.
    "T1": {"fasta": ["data/truth/ProteoBench-DDA-Astral/ProteoBenchFASTA_MixedSpecies_HYE.zip"],
           "foreign": ["data/fasta/arath_UP000006548.fasta"],
           "fasta_sha256": {"data/truth/ProteoBench-DDA-Astral/ProteoBenchFASTA_MixedSpecies_HYE.zip":
                            "cae44ed2e3cb9ce5972c28bd7303c8c9e6c1182d7fa30af06f722f5ac0613363"},   # == T1's (dev_sealed.sha256)
           "raws": None, "instrument": "QE", "prec_ppm": 10.0, "frag_ppm": 20.0, "role": "sealed"},
}
ALL_METHODS = ["vulcan", "sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot"]
CHIM_METHODS = {"vulcan_chim_v1": "score_v1", "vulcan_v2": "score_v2",       # E7/v2 (bench/t1/ours_chimeric.py), dev only
                "vulcan_chim_v1_sep": "score_v1", "vulcan_v2_sep": "score_v2"}   # amendment E7-A1: per-stratum TDC


def sage_chim_r(m):
    """E7b fair comparator (results/2026-10-28_E7b_fair_speed/config.yaml): 'sage_chim<R>' (Sage's own LDA, all ranks)
    and 'sage_chim<R>+pred+mokapot' -> R (Sage chimera: true, report_psms R); None for every other method."""
    import re
    mm = re.fullmatch(r"sage_chim(\d+)(\+pred\+mokapot)?", m)
    return int(mm.group(1)) if mm else None


def staged_python():
    """Re-exec under the node-local staged pred env (scripts/e9/stage_pred_env.sh; Lustre imports take minutes)."""
    if os.environ.get("VULCAN_T1_STAGED") == "1":
        return
    stager = os.path.join(ROOT, "scripts", "e9", "stage_pred_env.sh")
    if os.path.exists(stager):
        py = subprocess.run(["bash", stager], check=True, capture_output=True, text=True).stdout.strip().splitlines()[-1]
    else:                       # Lustre env (slow imports, same packages)
        py = os.path.join(ROOT, "conda_envs", "pred", "bin", "python")
    dest = os.path.dirname(os.path.dirname(os.path.dirname(py)))
    env = {**os.environ, "VULCAN_T1_STAGED": "1", "PYTHONPATH": ROOT}
    if os.path.exists(os.path.join(dest, "vulcanraw", "VulcanRaw.dll")):
        env.update(VULCAN_DOTNET_ROOT=os.path.join(dest, "dotnet"), VULCAN_TRFP_BIN=os.path.join(dest, "trfp"),
                   VULCAN_RAW_DLL=os.path.join(dest, "vulcanraw", "VulcanRaw.dll"))
    os.execve(py, [py, os.path.abspath(__file__)] + sys.argv[1:], env)


LOOK_TAGS = ("prereg-v0.4.3", "prereg-v0.4.4")
LOOK_STAGED_VARS = ("VULCAN_DOTNET_ROOT", "VULCAN_TRFP_BIN", "VULCAN_RAW_DLL")


def look_mode() -> bool:
    """A5 look (DESIGN §6.0d A5.11.2): the driver runs under `bench/sealed_run.py --look a5` (SEALED_LOOK_ID set)."""
    return bool(os.environ.get("SEALED_LOOK_ID"))


def look_guard(a, ap) -> None:
    """A5.11.2 items 1-3 (defence in depth; sealed_run.py checks the same before unlocking).  Refuses non-default
    analysis arguments, --out/--work outside $SEALED_LOOK_ROOT or non-empty (no resume), a look id that is not the open
    ledger line of this dataset, missing prereg tags, and any VULCAN_* variable that is not config/a5.yaml look.env's
    value (except VULCAN_T1_STAGED and the staged-tool paths the stager sets, which must point into the staged env)."""
    import yaml
    sys.path.insert(0, os.path.join(ROOT, "bench"))
    import sealed_run as SR
    lid, lroot, ds = os.environ["SEALED_LOOK_ID"], os.environ.get("SEALED_LOOK_ROOT", ""), os.environ.get("SEALED_DATASET")
    bad = []
    for k in ("stages", "methods", "ours", "mzml_cache", "keep_sealed_derivatives"):
        if getattr(a, k) != ap.get_default(k):
            bad.append(f"--{k.replace('_', '-')} is not allowed in the look")
    if not lroot or not os.path.isabs(lroot):
        bad.append("SEALED_LOOK_ROOT unset")
    for k in ("out", "work"):
        v = getattr(a, k)
        if not v:
            bad.append(f"--{k} must be given in the look")
            continue
        rv, rr = os.path.realpath(v), os.path.realpath(lroot or "/nonexistent")
        scr = os.environ.get("SEALED_LOOK_SCRATCH")
        rs = os.path.realpath(os.path.join(scr, ds or "")) if (scr and k == "work") else None   # node-local scratch
        if not rv.startswith(rr + os.sep) and not (rs and rv.startswith(rs + os.sep)):
            bad.append(f"--{k} is not under $SEALED_LOOK_ROOT")
        if os.path.exists(v) and (not os.path.isdir(v) or os.listdir(v)):
            bad.append(f"--{k} exists and is not empty (no resume in the look)")
    _, looks = SR.look_state(ROOT)
    L = looks.get(lid)
    if not L or ds not in L["datasets"] or ds in L["results"]:
        bad.append("SEALED_LOOK_ID is not the open ledger line of this dataset")
    head = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    for t in LOOK_TAGS:
        if not SR.is_ancestor(ROOT, t, head):
            bad.append(f"tag {t} is not an ancestor of HEAD")
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config", "a5.yaml")))
    lk = cfg.get("look") if isinstance(cfg, dict) else None
    if not isinstance(lk, dict) or not isinstance(lk.get("env"), dict):
        bad.append("config/a5.yaml look.env missing")
        lk = {"env": {}}
    want = {str(k): str(v) for k, v in lk["env"].items()}
    staged_dest = os.path.dirname(os.path.dirname(os.path.dirname(sys.executable)))
    for k, v in os.environ.items():
        if not k.startswith("VULCAN_") or k == "VULCAN_T1_STAGED":
            continue
        if k in LOOK_STAGED_VARS and os.path.realpath(v).startswith(os.path.realpath(staged_dest) + os.sep):
            continue
        if want.get(k) != v:
            bad.append(f"environment variable {k} is not in config/a5.yaml look.env with this value")
    for k, v in want.items():
        if k.startswith("VULCAN_") and k not in LOOK_STAGED_VARS and os.environ.get(k) != v:
            bad.append(f"environment variable {k} differs from config/a5.yaml look.env")
    if bad:
        sys.exit("t1_driver: REFUSED in the look: " + "; ".join(bad))


def write_manifest(out: str) -> str:
    """<out>/manifest.sha256: sha256 of every file under out (sorted), A5.5.3."""
    import hashlib
    lines = []
    for dp, dns, fns in os.walk(out):
        dns.sort()
        for f in sorted(fns):
            p = os.path.join(dp, f)
            if os.path.relpath(p, out) == "manifest.sha256":
                continue
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 22), b""):
                    h.update(b)
            lines.append(f"{h.hexdigest()}  {os.path.relpath(p, out)}")
    path = os.path.join(out, "manifest.sha256")
    open(path, "w").write("\n".join(lines) + "\n")
    return path


def core_s():
    a = resource.getrusage(resource.RUSAGE_SELF)
    b = resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


class Clock:
    def __init__(self):
        self.t = {}

    def __call__(self, name):
        clock = self

        class _C:
            def __enter__(self):
                self.t0, self.c0, self.l0 = time.time(), core_s(), os.getloadavg()[0]

            def __exit__(self, *a):
                clock.t[name] = {"wall_s": time.time() - self.t0, "core_s_self_and_children": core_s() - self.c0,
                                 "loadavg_start": self.l0, "loadavg_end": os.getloadavg()[0]}
        return _C()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    ap.add_argument("--out", default=None)
    ap.add_argument("--work", default=None, help="local-disk work dir (default /tmp/vulcan_t1_<uid>/<dataset>)")
    ap.add_argument("--stages", nargs="+", default=["db", "ours", "comparators", "evaluate"])
    ap.add_argument("--methods", nargs="+", default=ALL_METHODS)
    ap.add_argument("--threads", type=int, default=None, help="Sage threads (default vulcan.hw.n_threads())")
    ap.add_argument("--comet-jobs", type=int, default=6, help="concurrent Comet searches (each 4 threads)")
    ap.add_argument("--mokapot-jobs", type=int, default=6)
    ap.add_argument("--cohort-workers", type=int, default=None)
    ap.add_argument("--mzml-cache", default=None, help="reuse TRFP mzML (same flags) from this dir (dev only)")
    ap.add_argument("--rt-off-control", action="store_true", help="extra Sage run with predict_rt=false on prod_<j1>")
    ap.add_argument("--keep-sealed-derivatives", action="store_true")
    ap.add_argument("--ours", default="vulcan", help="method in the numerator of the §1 ratios (default: the frozen T1 "
                    "scorer 'vulcan'; E7/v2 dev runs use vulcan_v2)")
    ap.add_argument("--a3", action="store_true", help="Amendment A3 secondary arm (bench/t1/a3.py); forced for T1")
    a = ap.parse_args()
    staged_python()
    if look_mode():                         # A5.11.2: sealed hardening (asserted after the re-exec)
        look_guard(a, ap)
        a.rt_off_control = True             # A5.11.2.3: in the pre-registered 3f command
    from bench.t1 import a3 as A3
    from bench.t1 import a5 as A5
    if DATASETS[a.dataset]["role"] == "sealed":
        a.a3 = True                         # A3.3.8: v2_sep and every A3.3.3 comparator run in the look, unconditionally
    if a.a3:
        a.methods = list(dict.fromkeys(list(a.methods) + [m for m in A3.A3_METHODS if m not in a.methods]))
        if any(s in a.stages for s in ("ours",)) and "chim" not in a.stages:
            a.stages = [s for s in a.stages] + ["chim"]
        a.ours = A3.A3_NUMERATORS

    from vulcan import hw
    from vulcan.stats._util import DECOY_SEEDS, ENTRAP_SEEDS_DEV, ENTRAP_SEEDS_SEALED

    from bench.t1 import comparators as CMP
    from bench.t1 import dbs as DB
    from bench.t1 import evaluate as EV
    from bench.t1 import ours as OURS

    cfg = DATASETS[a.dataset]
    sealed = cfg["role"] == "sealed"
    a3_freeze = None
    if sealed:
        from vulcan.gates import a3_config
        a3_freeze = a3_config.require_freeze_commit()          # A3.1.3: freeze commit, clean tree, tree hashes
        a3_freeze["models"] = a3_config.check_models()          # v1, v2_sep, ce_s sha256 as recorded
        assert set(A3.A3_METHODS) <= set(a.methods) and a.ours == A3.A3_NUMERATORS and "chim" in a.stages
    if sealed:
        if os.environ.get("SEALED_DATASET") != a.dataset or not os.environ.get("SEALED_FILES"):
            sys.exit("t1_driver: T1 runs only under bench/sealed_run.py --dataset T1 (SEALED_FILES unset)")
        raws = sorted(os.environ["SEALED_FILES"].split())
        entrap = list(ENTRAP_SEEDS_SEALED)
        assert not a.mzml_cache, "no mzML cache for sealed data"
    else:
        raws = [os.path.join(ROOT, r) for r in cfg["raws"]]
        assert not any("data/sealed" in r for r in raws)
        entrap = list(ENTRAP_SEEDS_DEV)
    for f, h in cfg.get("fasta_sha256", {}).items():
        got = DB.sha256(os.path.join(ROOT, f))
        if got != h:
            sys.exit(f"t1_driver: FASTA {f} sha256 {got} != {h}")
    decoys = list(DECOY_SEEDS)
    assert not (set(entrap) & set(ENTRAP_SEEDS_SEALED)) or sealed, "sealed entrapment seeds on dev data"
    today = datetime.date.today().isoformat()
    out = a.out or os.path.join(ROOT, "results", f"{today}_T1driver_dryrun" if not sealed else f"{today}_T1_sealed",
                                a.dataset)
    work = a.work or os.path.join(f"/tmp/vulcan_t1_{os.getuid()}", a.dataset)
    os.makedirs(out, exist_ok=True)
    os.makedirs(work, exist_ok=True)
    threads = a.threads or hw.n_threads()
    # A3: the secondary arm and its comparators predict through their OWN cache (work/predcache_a3), so the A1 path
    # (v1 and Sage+pred+mokapot) sees exactly the prediction-cache state of a run without A3.  AlphaPeptDeep is not
    # batch-composition invariant in the last bits (3i, E7b), so a shared cache filled first by the chimeric stage
    # changed a few hundred predicted RTs of Sage+pred peptides (A3 dry run, results/2026-10-28_A3_freeze/NOTES.md).
    pc_a3 = os.path.join(work, "predcache_a3" if a.a3 else "predcache")
    logf = open(os.path.join(out, "driver.log"), "a")

    def log(msg):
        line = f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    clock = Clock()
    t_start, c_start = time.time(), core_s()
    sha = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True).stdout.strip())
    log(f"t1_driver {a.dataset} sha={sha[:12]} dirty_tracked={dirty} work={work} out={out}")
    reps = DB.replicate_names(decoys, entrap)
    first_prod = f"prod_{decoys[0]}"
    runs = [os.path.basename(r).rsplit(".", 1)[0] for r in raws]
    meta = {"dataset": a.dataset, "role": cfg["role"], "git_sha": sha, "git_dirty_tracked": dirty, "raws": raws,
            "runs": runs, "decoy_seeds": decoys, "entrap_seeds": entrap, "replicates": reps, "work": work,
            "hw": hw.fingerprint() if hasattr(hw, "fingerprint") else None, "threads": threads}

    # ---------------------------------------------------------------- db
    rep_info = None
    if "db" in a.stages or any(s in a.stages for s in ("ours", "chim", "comparators", "evaluate")):
        with clock("db"):
            full = DB.build_full(cfg, work, decoys, entrap, log)
            rep_info = DB.materialize(full, work, decoys, entrap, log)
            vf = os.path.join(work, "index", "verify.json")
            if not os.path.exists(vf):
                json.dump(DB.verify(full), open(vf, "w"), indent=1)
            meta["db_verify"] = json.load(open(vf))
            meta["index_meta"] = {k: v for k, v in json.load(open(os.path.join(full, "meta.json"))).items()
                                  if k not in ("fasta",)}
            assert meta["db_verify"]["ok"], meta["db_verify"]
            log(f"db verify: {meta['db_verify']['violations']}")

    # ---------------------------------------------------------------- ours
    scorer = None
    stage = os.path.join(work, "stage")
    cohort_dir = os.path.join(work, "ours", "cohort")
    if "ours" in a.stages and "vulcan" in a.methods:
        with clock("ours_cohort_ingest_and_prod_search"):
            st = OURS.run_cohort(raws, rep_info[first_prod]["index"], cohort_dir, stage,
                                 a.cohort_workers or max(1, min(len(raws), threads - 4)), log)
            meta["cohort"] = {k: st.get(k) for k in ("makespan_s", "runs_per_hour", "n_runs", "setup") if k in st}
        stores = OURS.stores_of(stage, runs)
        with clock("ours_replicate_searches"):
            meta["ours_search"] = OURS.search_replicates(stores, rep_info, os.path.join(work, "ours"), cohort_dir,
                                                         first_prod, log)
        with clock("ours_features_and_scoring"):
            scorer = OURS.Scorer(os.path.join(work, "predcache"), cfg["instrument"], a.dataset)
            meta["ours_score"] = OURS.score_replicates(scorer, stores, rep_info, os.path.join(work, "ours"), log)

    # ---------------------------------------------------------------- chimeric (E7/v2, dev only)
    if "chim" in a.stages and any(m in CHIM_METHODS for m in a.methods):
        assert not sealed or [m for m in a.methods if m in CHIM_METHODS] == [A3.SECONDARY], \
            "E7/v2 chimeric methods are dev-only (A3: only vulcan_v2_sep runs in the T1 look)"
        from bench.t1 import ours_chimeric as OC
        stores = OURS.stores_of(stage, runs)
        ms1_dir = os.path.join(work, "ms1")
        with clock("chim_ms1"):
            meta["chim_ms1"] = OC.run_ms1(stores, ms1_dir, log)
        with clock("chim_replicate_searches"):
            meta["chim_search"] = OC.search_replicates(stores, rep_info, os.path.join(work, "ours"),
                                                       os.path.join(work, "ours_chim"), ms1_dir, log)
        with clock("chim_features_and_scoring"):
            sc2 = OC.ScorerV2(pc_a3, cfg["instrument"], a.dataset)
            meta["chim_score"] = OC.score_replicates(sc2, stores, rep_info, os.path.join(work, "ours_chim"), ms1_dir,
                                                     log, force_v2=True)
            del sc2

    # ---------------------------------------------------------------- comparators
    cmp_root = os.path.join(work, "cmp")
    cm = [m for m in a.methods if m != "vulcan" and m not in CHIM_METHODS]
    if "comparators" in a.stages and cm:
        from concurrent.futures import ThreadPoolExecutor
        CMP.ensure_mokapot(log)
        with clock("cmp_mzml"):
            mz = CMP.ensure_mzml(raws, work, a.mzml_cache, log=log)
            meta["mzml"] = mz["info"]
        mzmls = [mz["files"][r] for r in runs]
        mk_pool = ThreadPoolExecutor(a.mokapot_jobs)
        futs = []
        if any(m.startswith("comet") for m in cm):
            comet_pool = ThreadPoolExecutor(a.comet_jobs)
            cfuts = {rep: [comet_pool.submit(CMP.run_comet, rep, rep_info[rep]["fasta"], f, cmp_root, cfg, 4)
                           for f in mzmls] for rep in reps}
        with clock("cmp_sage_searches"):
            if any(m.startswith("sage") and not sage_chim_r(m) for m in cm):
                for rep in reps:
                    CMP.run_sage(rep, rep_info[rep]["fasta"], mzmls, cmp_root, cfg, threads)
                    if "sage+mokapot" in cm:
                        futs.append(mk_pool.submit(CMP.run_mokapot, [os.path.join(cmp_root, "sage", rep, "results.sage.pin")],
                                                   os.path.join(cmp_root, "sage_mokapot", rep)))
                    log(f"sage {rep} done")
                if a.rt_off_control:
                    CMP.run_sage(first_prod, rep_info[first_prod]["fasta"], mzmls, cmp_root, cfg, threads,
                                 predict_rt=False, tag="sage_rt_off")
        chim_rs = sorted({sage_chim_r(m) for m in cm if sage_chim_r(m)})
        if chim_rs:
            with clock("cmp_sage_chim_searches"):
                for R in chim_rs:
                    for rep in reps:
                        CMP.run_sage(rep, rep_info[rep]["fasta"], mzmls, cmp_root, cfg, threads, tag=f"sage_chim{R}",
                                     chimera=True, report_psms=R)
                    log(f"sage chimera R={R} done")
        chim_pred = sorted({sage_chim_r(m) for m in cm if sage_chim_r(m) and m.endswith("+pred+mokapot")})
        if chim_pred:
            with clock("cmp_sage_chim_pred_features"):
                if a.a3:
                    scorer_chim = OURS.Scorer(pc_a3, cfg["instrument"], a.dataset)
                else:
                    if scorer is None:
                        scorer = OURS.Scorer(os.path.join(work, "predcache"), cfg["instrument"], a.dataset)
                    scorer_chim = scorer
                from vulcan.gpu.engine import load_store
                sps = {}
                for s in OURS.stores_of(stage, runs):
                    sp = load_store(s)
                    sps[sp.run] = sp
                for R in chim_pred:
                    meta[f"sage_chim{R}_pred_pin"] = {}
                    for rep in reps:
                        d = os.path.join(cmp_root, f"sage_chim{R}_pred_mokapot", rep)
                        meta[f"sage_chim{R}_pred_pin"][rep] = CMP.sage_pred_pin(
                            scorer_chim, os.path.join(cmp_root, f"sage_chim{R}", rep), rep_info[rep]["index"], sps, d, log,
                            chimeric=True)
                        futs.append(mk_pool.submit(CMP.run_mokapot, [os.path.join(d, "sage_pred.pin")], d))
                del sps, scorer_chim
        if "sage+pred+mokapot" in cm:
            with clock("cmp_sage_pred_features"):
                if scorer is None:
                    scorer = OURS.Scorer(os.path.join(work, "predcache"), cfg["instrument"], a.dataset)
                from vulcan.gpu.engine import load_store
                stores = OURS.stores_of(stage, runs)
                sps = {}
                for s in stores:
                    sp = load_store(s)
                    sps[sp.run] = sp
                meta["sage_pred_pin"] = {}
                for rep in reps:
                    d = os.path.join(cmp_root, "sage_pred_mokapot", rep)
                    meta["sage_pred_pin"][rep] = CMP.sage_pred_pin(scorer, os.path.join(cmp_root, "sage", rep),
                                                                   rep_info[rep]["index"], sps, d, log)
                    futs.append(mk_pool.submit(CMP.run_mokapot, [os.path.join(d, "sage_pred.pin")], d))
                del sps
        with clock("cmp_comet_searches_wait"):
            if any(m.startswith("comet") for m in cm):
                for rep in reps:
                    [f.result() for f in cfuts[rep]]
                    if "comet+mokapot" in cm:
                        pins = [os.path.join(cmp_root, "comet", rep, r + ".pin") for r in runs]
                        futs.append(mk_pool.submit(CMP.run_mokapot, pins, os.path.join(cmp_root, "comet_mokapot", rep)))
                comet_pool.shutdown()
        with clock("cmp_mokapot_wait"):
            [f.result() for f in futs]
            mk_pool.shutdown()
        meta["sage_rt_check"] = CMP.rt_check(os.path.join(cmp_root, "sage", first_prod))
        if a.rt_off_control:
            meta["sage_rt_off_check"] = CMP.rt_check(os.path.join(cmp_root, "sage_rt_off", first_prod))
        log(f"sage RT check: {meta['sage_rt_check']}")

    # ---------------------------------------------------------------- evaluate
    if "evaluate" in a.stages:
        with clock("evaluate"):
            def loader(m, rep, rl):
                try:
                    if m == "vulcan":
                        return EV.ours_units(os.path.join(work, "ours", rep), rl, runs)
                    if m in CHIM_METHODS:
                        from bench.t1 import ours_chimeric as OC
                        return OC.chim_units(os.path.join(work, "ours_chim", rep), rl, runs, CHIM_METHODS[m],
                                            qcol="q_sep" if m.endswith("_sep") else "q_gw")
                    if m == "sage":
                        return EV.sage_units(os.path.join(cmp_root, "sage", rep), rl, runs)
                    if m == "sage+mokapot":
                        return EV.mokapot_units(os.path.join(cmp_root, "sage_mokapot", rep), rl, runs, "sage")
                    if m == "sage+pred+mokapot":
                        return EV.mokapot_units(os.path.join(cmp_root, "sage_pred_mokapot", rep), rl, runs, "sage")
                    if m == "comet+mokapot":
                        return EV.mokapot_units(os.path.join(cmp_root, "comet_mokapot", rep), rl, runs, "comet")
                    if sage_chim_r(m):
                        R = sage_chim_r(m)
                        if m.endswith("+pred+mokapot"):
                            return EV.mokapot_units(os.path.join(cmp_root, f"sage_chim{R}_pred_mokapot", rep), rl,
                                                    runs, "sage")
                        return EV.sage_units(os.path.join(cmp_root, f"sage_chim{R}", rep), rl, runs, all_ranks=True)
                except (FileNotFoundError, KeyError) as e:
                    log(f"missing output {m} {rep}: {e}")
                    return None
            n_missing = []
            _loader = loader

            def loader(m, rep, rl):             # A5.5.4 check 1: count every "missing output" event
                r = _loader(m, rep, rl)
                if r is None:
                    n_missing.append((m, rep))
                return r
            metrics = EV.evaluate(work, runs, decoys, entrap, a.methods, loader, log, ours=a.ours)
            if a.a3:
                with clock("evaluate_a3_secondary"):
                    if "chim_search" not in meta:
                        p = os.path.join(work, "ours_chim")
                        meta["chim_search"] = {r: json.load(open(os.path.join(p, r, "search_done.json"))) for r in reps
                                               if os.path.exists(os.path.join(p, r, "search_done.json"))}
                    metrics["a3_secondary"] = A3.secondary(work, runs, decoys, entrap, metrics, EV.LAST_CURVES,
                                                           a.methods, meta, log)
                with clock("evaluate_a5_additions"):         # A5 (DESIGN §6.0d) A5.11.2.5; additive, nothing re-searched
                    metrics["a5"] = A5.additions(work, runs, decoys, entrap, metrics, EV.LAST_CURVES, a.methods, meta,
                                                 log, loader=loader)
            if os.path.exists(os.path.join(cmp_root, "sage", first_prod, "results.sage.tsv")):
                meta["sage_rt_check"] = CMP.rt_check(os.path.join(cmp_root, "sage", first_prod))
            if os.path.exists(os.path.join(cmp_root, "sage_rt_off", first_prod, "results.sage.tsv")):
                meta["sage_rt_off_check"] = CMP.rt_check(os.path.join(cmp_root, "sage_rt_off", first_prod))
            if os.path.exists(os.path.join(cmp_root, "sage_rt_off", first_prod, "results.sage.tsv")):
                rl = EV.RepLabels(os.path.join(work, "index", first_prod))
                on = EV.sage_units(os.path.join(cmp_root, "sage", first_prod), rl, runs)
                off = EV.sage_units(os.path.join(cmp_root, "sage_rt_off", first_prod), rl, runs)
                metrics["sage_predict_rt_control"] = {
                    lv: {str(al): {"rt_on": EV.ids_dataset(f(on), lv, al), "rt_off": EV.ids_dataset(f(off), lv, al)}
                         for al in (0.001, 0.01)}
                    for lv, f in (("peptidoform", lambda x: x), ("sequence", EV.seq_units))}
    t_total = {"wall_s": time.time() - t_start, "core_s_self_and_children": core_s() - c_start,
               "stages": clock.t}
    meta["driver_timing"] = t_total
    with open(os.path.join(out, "timing.jsonl"), "a") as fh:   # one line per invocation (stages resume from cache)
        fh.write(json.dumps({"start": datetime.datetime.fromtimestamp(t_start).isoformat(timespec="seconds"),
                             "argv": sys.argv[1:], "git_sha": sha, **t_total}) + "\n")
    log(f"driver total {t_total['wall_s']:.0f}s")
    json.dump(meta, open(os.path.join(out, "run_meta.json"), "w"), indent=1, default=str)
    if "evaluate" in a.stages:
        metrics["driver_timing"] = t_total
        metrics["meta"] = {k: meta[k] for k in ("dataset", "role", "git_sha", "git_dirty_tracked", "runs",
                                                "decoy_seeds", "entrap_seeds") if k in meta}
        if a.a3:
            metrics["meta"]["a3"] = {"methods": a.methods, "numerators": a.ours, "freeze_status": a3_freeze}
        metrics["sage_rt_check"] = {"predict_rt_true": meta.get("sage_rt_check"),
                                    "predict_rt_false": meta.get("sage_rt_off_check")}
        json.dump(metrics, open(os.path.join(out, "metrics.json"), "w"), indent=1, default=float)
        from bench.t1.report import tables
        open(os.path.join(out, "tables.md"), "w").write(tables(metrics))
        log("wrote metrics.json and tables.md")
    if sealed and not a.keep_sealed_derivatives and not look_mode():
        for d in (stage, os.path.join(work, "mzml"), os.path.join(work, "ms1")):
            shutil.rmtree(d, ignore_errors=True)
        log("deleted staged stores and mzML of sealed runs")
    if look_mode():
        # A5.5.3: delete the whole work dir except the replicate indexes (built from FASTAs, not sealed data); then the
        # output manifest; then the completeness rule (A5.5.4 check 1: non-zero exit on any missing output)
        for e in os.listdir(work):
            if e != "index":
                p = os.path.join(work, e)
                shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) and not os.path.islink(p) else os.remove(p)
        log("deleted work except replicate indexes")
        from bench.t1 import completeness as COMPL
        problems = []
        if "evaluate" not in a.stages:
            problems.append("evaluate stage not run")
        else:
            problems = [f"missing output {m} {r}" for m, r in n_missing] + COMPL.check(metrics, runs=runs)
        if problems:
            open(os.path.join(out, "INCOMPLETE"), "w").write(f"{len(problems)} completeness violations\n")
        logf.close()
        write_manifest(out)
        print("manifest written", flush=True)
        if problems:
            print(f"look failed: incomplete ({len(problems)} violations)", flush=True)
            sys.exit(3)
        print("completeness rule: pass", flush=True)


if __name__ == "__main__":
    main()
