#!/usr/bin/env python3
"""Arm M (DESIGN §6.3 (M), A1.5) for v1 and v2_sep: the parametrised module of Amendment A5 (DESIGN §6.0d, A5.6 and
A5.11.3), shared by dev mode (D7, dev seeds) and sealed mode (T5 = PXD052597, both runs, in the A5 look).

    dev:     python bench/t5_armM.py --mode dev --runs data/raw/PXD031616/LM2486_DDA.raw ... \
                 --entrap-seeds 2001 2002 2003 --decoy-seeds 1001 --out <dir> --work <dir> [--stores DIR]
                 [--targets FASTA ...] [--db-recipe dev_foreign]
    sealed:  python3 bench/sealed_run.py --look a5   (config/a5.yaml look.commands.T5 runs)
             python bench/t5_armM.py --mode sealed --out $SEALED_LOOK_ROOT/out --work $SEALED_LOOK_ROOT/work

Design (A5.6; carried from scripts/m3e/armM_prepare.py + scripts/m3e/armM.py (v1) and scripts/a3/dev_armM.py (v2_sep),
which stay unchanged with their committed outputs):
  databases   targets = yeast_UP000002311 + contaminants_frankenfield2022 (sealed: fixed; no spiked protein, no sample
              sheet read).  Entrapment = human Swiss-Prot peptides neighbour-filtered against the targets exactly as arm C
              (scripts/entrapment/build_dbs.py --arms C: I/L-identical, isobaric N<->GG / Q<->AG, single substitution,
              decoy collisions for DECOY_SEEDS 1001-1005), restricted to the effective search space (in_space) and
              subsampled to r = 1 (|E_eff| = |T_eff|, numpy default_rng(entrapment seed)), one draw per entrapment seed.
              Decoys:
                recipe "shuffle" (sealed; forced): one vulcan.index.build.build_multi_decoy index per entrapment seed
                  (T + E_s + C-term-fixed hash shuffles of every DECOY_SEED), materialised by bench/t1/dbs.materialize as
                  replicate M_<j>_<s> = T + E_s + dec_j(T) + dec_j(E_s)  -> 5 x 3 = 15 replicates (§6.3 milestone);
                recipe "dev_foreign" (dev only; the dev design): one vulcan.index.foreign.build_foreign index per seed
                  (pseudo-reverse decoys, collisions shuffled with decoy seed 1001), as armM_prepare.py.
  v1          closed engine search (report_psms 5) of every replicate (bench/t1/ours.search_replicates), L3 features +
              frozen static_gbm_pred through bench/t1/ours.py `Scorer` (the T1 path), paired double competition
              (vulcan.stats.tdc.double_competition, seed 11, hash tiebreak) over all runs as one experiment.
  v2_sep      vulcan.ms1 hypotheses, chimeric search (bench/t1/ours_chimeric.search_replicates), v2 features + frozen v2
              (ours_chimeric.ScorerV2), vulcan.chimeric.compete.chimeric_competition (seed 11), q = q_sep.
              Lists: overall; each MERGED stratum (frozen label-blind >= 200-pair merge; decides voids); each raw
              chimeric stratum {1, 2, >=3} (reported).
  reference   hyperscore through each pipeline (hyperscore_v1: closed; hyperscore_v2: chimeric) -- reported only.
  estimator   combined FDP_M = N_E (1 + 1/r) / (N_T + N_E), r = peptide-level effective r (1.0 by construction);
              units peptidoform (decides) and sequence (reported); alpha 0.1% (reported), 1%, 5%;
              vulcan.stats.entrapment.evaluate_arm, CI convention (a) decides, (b) reported.
  caches      fresh prediction caches under --work (predcache_t5_v1, predcache_t5_v2), never shared with T1.
Outputs: <out>/armM.json (schema below), <out>/manifest.sha256.  Nothing is written outside --out and --work.
Sealed mode: raws only from $SEALED_FILES (SEALED_DATASET == T5); seeds fixed (ENTRAP_SEEDS_SEALED x DECOY_SEEDS);
--out/--work must lie under $SEALED_LOOK_ROOT and be absent or empty (no resume); dev-only options refused; at the end
work/ is deleted except the replicate indexes (db_*/), which are built from FASTAs only.
Completeness: every replicate x run present in every list, else armM.json has complete=false and the exit code is 3.

armM.json = {"mode", "runs", "entrap_seeds", "decoy_seeds", "n_replicates", "complete", "missing",
             "v1": {"overall": L}, "v2_sep": {"overall": L, "merged": {block: L}, "raw": {"1"|"2"|">=3": L}},
             "hyperscore_v1": {"overall": L}, "hyperscore_v2": {"overall": L}, + provenance keys}
L = {unit: {alpha: {"a": {verdict, mean_fdp, ci95, mean_lb, lb_ci95, mean_n_t}, "b": {...}, "counts": {rep: {N_T, N_E}}}}}
Merged block keys: the block label ("1", "2_3", ...) when every replicate has the same merge map (else
"contains_<k>"); "merge_maps" records each replicate's map.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

TARGETS = ["data/fasta/yeast_UP000002311.fasta", "data/fasta/contaminants_frankenfield2022.fasta"]
FOREIGN = ["data/fasta/human_sp_canonical.fasta"]
INSTRUMENT = "QE"                 # Q Exactive family (T5: Q Exactive HF-X; D7 dev used "QE")
UNITS = ("peptidoform", "sequence")
ALPHAS = (0.001, 0.01, 0.05)
GATED = (0.01, 0.05)
SEED = 11
STRATA = (1, 2, 3)
ARM_KEYS = ("verdict", "mean_fdp", "ci95", "mean_lb", "lb_ci95", "mean_n_t")
SEALED_DATASET = "T5"
# sealed console whitelist (config/a5.yaml look.console; sealed_run filters again): data-free lines only
CONSOLE_RE = __import__("re").compile(r"^t5_armM: (stage [a-z0-9_]+( done)?|complete|wrote armM\.json and manifest\.sha256)$")


# ------------------------------------------------------------------------------------------------ arguments

def _dev_look_standins(env) -> bool:
    """A5.11.8 dev look dry run: SEALED_FILES (set by sealed_run --look a5 --dev) may hold stand-in files only, never a
    file listed as sealed in data/manifests/dev_sealed.tsv; the dev raws themselves come from --runs."""
    import csv
    if not env.get("SEALED_LOOK_ID"):
        return False
    real = set()
    p = os.path.join(ROOT, "data", "manifests", "dev_sealed.tsv")
    if os.path.exists(p):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("role") == "sealed":
                    real.add(r["filename"])
    return not any(os.path.basename(f) in real for f in env["SEALED_FILES"].split())

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["dev", "sealed"])
    ap.add_argument("--runs", nargs="*", default=None, help="RAW files (dev); sealed: omitted or '{files}' "
                    "expanded by sealed_run (must equal $SEALED_FILES)")
    ap.add_argument("--entrap-seeds", nargs="+", type=int, default=None)
    ap.add_argument("--decoy-seeds", nargs="+", type=int, default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--stores", default=None, help="dev only: existing <run>.vspec.arrow stores (skip ingest)")
    ap.add_argument("--targets", nargs="+", default=None, help="dev only: target FASTAs (default yeast + contaminants)")
    ap.add_argument("--db-recipe", choices=["shuffle", "dev_foreign"], default=None,
                    help="decoy recipe (sealed: shuffle, forced; dev default: dev_foreign if one decoy seed)")
    ap.add_argument("--workers", type=int, default=None, help="cohort decode workers")
    return ap.parse_args(argv)


def _under(path, root):
    p, r = os.path.realpath(path), os.path.realpath(root)
    return p == r or p.startswith(r.rstrip("/") + "/")


def _fresh(path):
    return not os.path.exists(path) or (os.path.isdir(path) and not os.listdir(path))


def check_args(a, env=None) -> dict:
    """Mode rules (A5.11.3.1, A5.11.2.1).  Returns the resolved settings; SystemExit on any refusal."""
    from vulcan.stats._util import DECOY_SEEDS, ENTRAP_SEEDS_DEV, ENTRAP_SEEDS_SEALED
    env = os.environ if env is None else env
    if a.mode == "sealed":
        if env.get("SEALED_DATASET") != SEALED_DATASET or not env.get("SEALED_FILES"):
            sys.exit("t5_armM: REFUSED, sealed mode runs only under bench/sealed_run.py --look a5 (SEALED_DATASET=T5)")
        files = sorted(env["SEALED_FILES"].split())
        if a.runs and sorted(a.runs) != files:
            sys.exit("t5_armM: REFUSED, --runs differs from $SEALED_FILES")
        if a.entrap_seeds is not None and list(a.entrap_seeds) != list(ENTRAP_SEEDS_SEALED):
            sys.exit("t5_armM: REFUSED, sealed entrapment seeds are ENTRAP_SEEDS_SEALED")
        if a.decoy_seeds is not None and list(a.decoy_seeds) != list(DECOY_SEEDS):
            sys.exit("t5_armM: REFUSED, sealed decoy seeds are DECOY_SEEDS (all 5)")
        for opt in ("stores", "targets"):
            if getattr(a, opt) is not None:
                sys.exit(f"t5_armM: REFUSED, --{opt.replace('_', '-')} is dev only")
        if a.db_recipe not in (None, "shuffle"):
            sys.exit("t5_armM: REFUSED, sealed mode uses the shuffle decoy recipe")
        lr = env.get("SEALED_LOOK_ROOT")
        if not lr:
            sys.exit("t5_armM: REFUSED, $SEALED_LOOK_ROOT unset")
        sc = env.get("SEALED_LOOK_SCRATCH")
        if not _under(a.out, lr):
            sys.exit(f"t5_armM: REFUSED, --out {a.out} is not under $SEALED_LOOK_ROOT")
        if not (_under(a.work, lr) or (sc and _under(a.work, os.path.join(sc, SEALED_DATASET)))):
            sys.exit(f"t5_armM: REFUSED, --work {a.work} is not under $SEALED_LOOK_ROOT or $SEALED_LOOK_SCRATCH/T5")
        for nm, p in (("--out", a.out), ("--work", a.work)):
            if not _fresh(p):
                sys.exit(f"t5_armM: REFUSED, {nm} {p} exists and is not empty (no resume in sealed mode)")
        return {"raws": files, "entrap": list(ENTRAP_SEEDS_SEALED), "decoys": list(DECOY_SEEDS),
                "targets": [os.path.join(ROOT, t) for t in TARGETS], "recipe": "shuffle", "stores": None,
                "dataset": SEALED_DATASET}
    raws = [os.path.abspath(r) for r in (a.runs or [])]
    if not raws:
        sys.exit("t5_armM: give --runs in dev mode")
    entrap = list(a.entrap_seeds or ENTRAP_SEEDS_DEV)
    decoys = list(a.decoy_seeds or DECOY_SEEDS)
    if set(entrap) & set(ENTRAP_SEEDS_SEALED):
        sys.exit("t5_armM: REFUSED, sealed entrapment seeds on dev data")
    bad = [p for p in [*raws, a.out, a.work, a.stores or "", *(a.targets or [])] if "data/sealed" in os.path.realpath(p)
           or "data/sealed" in p]
    if bad or (env.get("SEALED_FILES") and not _dev_look_standins(env)):
        sys.exit("t5_armM: REFUSED, sealed paths in dev mode")
    recipe = a.db_recipe or ("dev_foreign" if len(decoys) == 1 else "shuffle")
    return {"raws": raws, "entrap": entrap, "decoys": decoys,
            "targets": [os.path.abspath(t) for t in (a.targets or [os.path.join(ROOT, t) for t in TARGETS])],
            "recipe": recipe, "stores": a.stores, "dataset": "D7"}


def staged_python():
    """Re-exec under the node-local staged pred env (as bench/t1_driver.staged_python; scripts/e9/stage_pred_env.sh)."""
    if os.environ.get("VULCAN_T1_STAGED") == "1":
        return
    stager = os.path.join(ROOT, "scripts", "e9", "stage_pred_env.sh")
    if os.path.exists(stager):
        py = subprocess.run(["bash", stager], check=True, capture_output=True, text=True).stdout.strip().splitlines()[-1]
    else:
        py = os.path.join(ROOT, "conda_envs", "pred", "bin", "python")
    dest = os.path.dirname(os.path.dirname(os.path.dirname(py)))
    env = {**os.environ, "VULCAN_T1_STAGED": "1", "PYTHONPATH": ROOT}
    if os.path.exists(os.path.join(dest, "vulcanraw", "VulcanRaw.dll")):
        env.update(VULCAN_DOTNET_ROOT=os.path.join(dest, "dotnet"), VULCAN_TRFP_BIN=os.path.join(dest, "trfp"),
                   VULCAN_RAW_DLL=os.path.join(dest, "vulcanraw", "VulcanRaw.dll"))
    os.execve(py, [py, os.path.abspath(__file__)] + sys.argv[1:], env)


# ------------------------------------------------------------------------------------------------ databases
def rep_name(j, s):
    return f"M_{j}_{s}"


def target_fasta(paths, work):
    out = os.path.join(work, "db", "targets.fasta")
    if os.path.exists(out):
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out + ".partial", "w") as oh:
        for p in paths:
            txt = open(p).read()
            oh.write(txt if txt.endswith("\n") else txt + "\n")
    os.replace(out + ".partial", out)
    return out


def armc_peptides(tfa, work, log):
    """Neighbour-filtered human peptides (scripts/entrapment/build_dbs.py --arms C, unchanged) -> (T, E) lists."""
    d = os.path.join(work, "db", "armC")
    tsv = os.path.join(d, "t5.C.peptides.tsv")
    if not os.path.exists(tsv):
        os.makedirs(d, exist_ok=True)
        cmd = [sys.executable, os.path.join(ROOT, "scripts", "entrapment", "build_dbs.py"), "--name", "t5",
               "--target", tfa, "--foreign", *[os.path.join(ROOT, f) for f in FOREIGN], "--arms", "C", "--out", d]
        log("stage db_armc")
        with open(os.path.join(d, "build.log"), "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
    T, E = [], []
    with open(tsv) as fh:
        next(fh)
        for line in fh:
            seq, kind = line.split("\t", 2)[:2]
            (T if kind == "T" else E if kind == "E" else []).append(seq)
    return T, E


def subsample(T, E, seed):
    """armM_prepare.py: effective search space, r = 1 draw with numpy default_rng(seed)."""
    p = os.path.join(ROOT, "scripts", "entrapment")
    if p not in sys.path:
        sys.path.insert(0, p)
    from peptides import in_space  # noqa: E402
    Teff = [t for t in T if in_space(t)]
    Eeff = sorted(e for e in E if in_space(e))
    pick = np.sort(np.random.default_rng(seed).choice(len(Eeff), size=len(Teff), replace=False))
    sub = [Eeff[i] for i in pick]
    return sub, {"n_target": len(T), "n_target_eff": len(Teff), "n_foreign_filtered": len(E),
                 "n_foreign_eff": len(Eeff), "r_peptide_eff": len(sub) / len(Teff)}


def prepare_db(st, work, log) -> tuple[dict, dict]:
    """-> ({rep: {"index", "decoy_seed", "entrap_seed"}}, db info)."""
    tfa = target_fasta(st["targets"], work)
    T, E = armc_peptides(tfa, work, log)
    reps, info = {}, {"recipe": st["recipe"], "target_fasta_sha256": _sha(tfa), "seeds": {}}
    for s in st["entrap"]:
        sub, si = subsample(T, E, s)
        d = os.path.join(work, f"db_{s}")
        os.makedirs(d, exist_ok=True)
        fl = os.path.join(d, "foreign.txt")
        if not os.path.exists(fl):
            with open(fl + ".partial", "w") as fh:
                fh.write("\n".join(sub) + "\n")
            os.replace(fl + ".partial", fl)
        si["foreign_sha256"] = _sha(fl)
        if st["recipe"] == "dev_foreign":
            assert len(st["decoys"]) == 1, "dev_foreign recipe: one decoy seed"
            j = st["decoys"][0]
            idx = os.path.join(d, "index", rep_name(j, s))
            if not os.path.exists(os.path.join(idx, "meta.json")):
                from vulcan.index.foreign import build_foreign
                build_foreign([tfa], idx, fl, params={"decoy_seed": j})
            reps[rep_name(j, s)] = {"index": idx, "decoy_seed": j, "entrap_seed": s}
        else:
            from bench.t1 import dbs as DB
            from vulcan.index.build import build_multi_decoy
            full = os.path.join(d, "index", "full")
            if not os.path.exists(os.path.join(full, "meta.json")):
                meta = build_multi_decoy([tfa], full, {"entrap_seeds": [], "decoy_seeds": list(st["decoys"]),
                                                       "foreign_sha256": si["foreign_sha256"]}, foreign=sub)
                json.dump(meta, open(os.path.join(full, "meta.json"), "w"), indent=1, default=str)
            mat = DB.materialize(full, d, st["decoys"], [], lambda m: None)
            for j in st["decoys"]:
                reps[rep_name(j, s)] = {"index": mat[f"C_{j}"]["index"], "decoy_seed": j, "entrap_seed": s}
        nb = json.load(open(os.path.join(reps[rep_name(st["decoys"][0], s)]["index"], "meta.json")))
        si["n_peptidoforms_by_family"] = nb.get("n_peptidoforms_by_family")
        info["seeds"][str(s)] = si
        log(f"stage db_seed_{s}")
    return reps, info


# ------------------------------------------------------------------------------------------------ pipelines
def stores_for(st, work, reps, log):
    """Kept stores: dev may reuse --stores; otherwise vulcan.cohort ingest (the T1 path, bench/t1/ours.run_cohort)
    with the production search on the first replicate (checked bit-for-bit against the engine search)."""
    from bench.t1 import ours as OURS
    runs = sorted(os.path.basename(r).rsplit(".", 1)[0] for r in st["raws"])
    if st["stores"]:
        return [os.path.join(st["stores"], r + ".vspec.arrow") for r in runs], runs, None
    first = next(iter(reps))
    cohort = os.path.join(work, "cohort")
    stage = os.path.join(work, "stage")
    OURS.run_cohort(st["raws"], reps[first]["index"], cohort, stage, st.get("workers") or max(1, min(len(runs), 8)),
                    log)
    log("stage ingest")
    return OURS.stores_of(stage, runs), runs, first


def run_pipelines(st, work, reps, stores, first, log) -> dict:
    from bench.t1 import ours as OURS
    from bench.t1 import ours_chimeric as OC
    info = {}
    ours = os.path.join(work, "ours")
    info["v1_search"] = OURS.search_replicates(stores, reps, ours, os.path.join(work, "cohort"), first, log)
    log("stage v1_search")
    sc = OURS.Scorer(os.path.join(work, "predcache_t5_v1"), INSTRUMENT, st["dataset"])
    info["v1_score"] = OURS.score_replicates(sc, stores, reps, ours, log)
    del sc
    log("stage v1_score")
    ms1 = os.path.join(work, "ms1")
    info["ms1"] = OC.run_ms1(stores, ms1, log)
    log("stage ms1")
    chim = os.path.join(work, "ours_chim")
    info["v2_search"] = OC.search_replicates(stores, reps, ours, chim, ms1, log)
    log("stage v2_search")
    sc2 = OC.ScorerV2(os.path.join(work, "predcache_t5_v2"), INSTRUMENT, st["dataset"])
    info["v2_sha256"] = sc2.v2_sha
    info["v2_score"] = OC.score_replicates(sc2, stores, reps, chim, ms1, log, force_v2=True)
    del sc2
    log("stage v2_score")
    return info


# ------------------------------------------------------------------------------------------------ evaluation
def _labels(index_dir):
    import pyarrow as pa
    with pa.OSFile(os.path.join(index_dir, "labels.arrow"), "rb") as src:
        t = pa.ipc.open_file(src).read_all()
    return {k: t.column(k).to_numpy(zero_copy_only=False) for k in ("is_decoy", "is_entrapment", "tdc_pair")}


def _seqs(index_dir, rows):
    import pyarrow as pa

    from bench.t1.dbs import seq_bytes
    with pa.OSFile(os.path.join(index_dir, "peptidoforms.arrow"), "rb") as src:
        pf = pa.ipc.open_file(src).read_all()
    res = np.load(os.path.join(index_dir, "residues.npy"), mmap_mode="r")
    rows = np.asarray(rows, np.int64)
    off = pf.column("seq_off").to_numpy()[rows]
    ln = pf.column("length").to_numpy()[rows]
    return seq_bytes(res, off, ln)


def _read_runs(d, runs, suffix, missing, tag):
    import pandas as pd
    import pyarrow.parquet as pq
    parts = []
    for i, r in enumerate(runs):
        p = os.path.join(d, f"{r}{suffix}")
        if not os.path.exists(p):
            missing.append(f"{tag}:{r}")
            continue
        t = pq.read_table(p).to_pandas()
        if len(t) == 0:
            missing.append(f"{tag}:{r}:empty")
        t["ri"] = i
        parts.append(t)
    return pd.concat(parts, ignore_index=True) if parts else None


def winners_closed(rep, info, runs, work, missing) -> dict | None:
    """v1 / hyperscore through the closed pipeline: paired double competition over all runs (as scripts/m3e/armM.py)."""
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    psm = _read_runs(os.path.join(work, "ours", rep), runs, ".scored.parquet", missing, f"v1/{rep}")
    if psm is None:
        return None
    lab = _labels(info["index"])
    rows = psm.pf_row.values
    spec = np.array([f"{runs[a]}:{b}" for a, b in zip(psm.ri.values, psm.spec_idx.values)])
    unit = psm.pf_id.values
    tb = hash_tiebreak(unit, SEED)
    out = {}
    for name, col in (("v1", "score"), ("hyperscore_v1", "hyperscore")):
        r = double_competition(spec, unit, lab["tdc_pair"][rows], psm[col].values.astype(float),
                               lab["is_decoy"][rows], seed=SEED, tiebreak=tb)
        wr = rows[r.row]
        out[name] = {"is_decoy": np.asarray(r.is_decoy, bool), "q": np.asarray(r.q), "ent": lab["is_entrapment"][wr],
                     "unit": unit[r.row], "seq": _seqs(info["index"], wr), "lists": {"overall": None}}
    return out


def winners_chim(rep, info, runs, work, missing, merge_maps) -> dict | None:
    """v2_sep / hyperscore through the chimeric competition (as scripts/a3/dev_armM.py), q = q_sep."""
    from vulcan.chimeric.compete import IonSets, chimeric_competition
    from vulcan.index.build import L1Index
    psm = _read_runs(os.path.join(work, "ours_chim", rep), runs, ".chim.scored.parquet", missing, f"v2/{rep}")
    if psm is None:
        return None
    lab = _labels(info["index"])
    ions = IonSets(L1Index(info["index"]))
    row = psm.pf_row.values
    spec = psm.ri.values.astype(np.int64) * (1 << 26) + psm.spec_idx.values.astype(np.int64)
    out = {}
    for name, col in (("v2_sep", "score_v2"), ("hyperscore_v2", "hyperscore")):
        cc = chimeric_competition(spec, psm.group.values, psm.chim_rank.values, psm.pf_id.values, lab["tdc_pair"][row],
                                  psm[col].values.astype(np.float64), lab["is_decoy"][row], row, ions, seed=SEED)
        res = cc["res"]
        wr = row[cc["win"]]
        merged = np.asarray(cc["merged"])
        stratum = np.asarray(cc["stratum"])
        lists = {"overall": None}
        if name == "v2_sep":
            merge_maps[rep] = {str(k): list(v) for k, v in cc["merge"].items()}
            for k in STRATA:
                lists[("raw", k)] = stratum == k
                lists[("merged", k)] = merged == "_".join(map(str, cc["merge"][k]))
        out[name] = {"is_decoy": np.asarray(res.is_decoy, bool), "q": np.asarray(cc["q_sep"]), "ent": lab["is_entrapment"][wr],
                     "unit": psm.pf_id.values[cc["win"]], "seq": _seqs(info["index"], wr), "lists": lists}
    return out


def _counts(x, list_key, unit, a):
    from vulcan.stats import entrapment as ENT
    ts = ~x["is_decoy"]
    if list_key is not None:
        ts = ts & x["lists"][list_key]
    acc = x["q"] <= a
    u = x["unit"] if unit == "peptidoform" else x["seq"]
    return (ENT.combined_counts(u[ts], x["ent"][ts], acc[ts]),
            {"N_T": int((acc & ts & ~x["ent"]).sum()), "N_E": int((acc & ts & x["ent"]).sum())})


def _arm_eval(job):
    """One evaluate_arm call (convention (a)); convention (b) is read off the SAME bootstrap: evaluate_arm returns
    ci95_b / lb_ci95_b, and its verdict for (b) is ENT.verdict on those bounds -- identical to a second call with
    ci_convention="b" (same seed, same draws; checked bit-for-bit against scripts/a3/dev_armM.py's two calls on D7)."""
    from vulcan.stats import entrapment as ENT
    rs, a, r_eff = job
    e = ENT.evaluate_arm(rs, a, "combined", r=r_eff, ci_convention="a")
    out = {"a": {k: e[k] for k in ARM_KEYS}}
    b = {"mean_fdp": e["mean_fdp"], "ci95": e["ci95_b"], "mean_lb": e["mean_lb"], "lb_ci95": e["lb_ci95_b"],
         "mean_n_t": e["mean_n_t"]}
    b["verdict"] = ENT.verdict(a, e["mean_fdp"], e["ci95_b"][1], e["lb_ci95_b"][0], e["mean_n_t"])
    out["b"] = {k: b[k] for k in ARM_KEYS}
    return out


def arm_tables(specs: dict, reps: dict, r_eff: float, workers: int = 1) -> dict:
    """specs = {name: (x_by_rep, list_key)} -> {name: L}.  L = {unit: {alpha: {"a", "b", "counts"}}} (combined
    estimator, one Replicate per replicate DB).  The bootstrap calls run in a fork pool (deterministic: each call
    seeds its own generator, so the result does not depend on the pool)."""
    from vulcan.stats import entrapment as ENT
    jobs, keys, out = [], [], {}
    for name, (x_by_rep, list_key) in specs.items():
        out[name] = {}
        for unit in UNITS:
            out[name][unit] = {}
            for a in ALPHAS:
                rs, counts = [], {}
                for rep, x in x_by_rep.items():
                    c, n = _counts(x, list_key, unit, a)
                    rs.append(ENT.Replicate(c, reps[rep]["entrap_seed"], reps[rep]["decoy_seed"]))
                    counts[rep] = n
                out[name][unit][str(a)] = {"counts": counts}
                jobs.append((rs, a, r_eff))
                keys.append((name, unit, str(a)))
    if workers > 1 and len(jobs) > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(workers, len(jobs)), mp_context=mp.get_context("fork")) as ex:
            res = list(ex.map(_arm_eval, jobs))
    else:
        res = [_arm_eval(j) for j in jobs]
    for (name, unit, a), r in zip(keys, res):
        out[name][unit][a].update(r)
    return out


def arm_m(x_by_rep: dict, reps: dict, list_key, r_eff: float, workers: int = 1) -> dict:
    """L for one list (see arm_tables)."""
    return arm_tables({"_": (x_by_rep, list_key)}, reps, r_eff, workers)["_"]


def eval_workers() -> int:
    try:
        n = int(os.environ.get("VULCAN_THREADS", "0"))
    except ValueError:
        n = 0
    return max(1, min(16, n or (os.cpu_count() or 2) // 2))


def evaluate(st, reps, runs, work, dbinfo, log) -> dict:
    missing = []
    r_effs = {v["r_peptide_eff"] for v in dbinfo["seeds"].values()}
    assert len(r_effs) == 1, r_effs
    r_eff = float(r_effs.pop())
    per = {k: {} for k in ("v1", "hyperscore_v1", "v2_sep", "hyperscore_v2")}
    merge_maps = {}
    for rep, info in reps.items():
        c = winners_closed(rep, info, runs, work, missing)
        h = winners_chim(rep, info, runs, work, missing, merge_maps)
        for d in (c, h):
            for k, v in (d or {}).items():
                per[k][rep] = v
    for k, d in per.items():
        for rep in reps:
            if rep not in d:
                missing.append(f"{k}/{rep}:replicate")
    res = {"mode": st["mode"], "runs": runs, "entrap_seeds": st["entrap"], "decoy_seeds": st["decoys"],
           "n_replicates": len(reps), "complete": False, "missing": sorted(set(missing)), "r_eff": r_eff}
    if res["missing"]:
        return res
    specs = {(k, "overall"): (per[k], None) for k in ("v1", "hyperscore_v1", "hyperscore_v2", "v2_sep")}
    maps = {json.dumps(m, sort_keys=True) for m in merge_maps.values()}
    same = len(maps) == 1
    seen = set()
    for k in STRATA:
        kk = str(k) if k < 3 else ">=3"
        specs[("raw", kk)] = (per["v2_sep"], ("raw", k))
        if same:
            m = next(iter(merge_maps.values()))
            lab = "_".join(map(str, m[str(k)]))
            if lab in seen:
                continue
            seen.add(lab)
        else:
            lab = f"contains_{kk}"
        specs[("merged", lab)] = (per["v2_sep"], ("merged", k))
    tabs = arm_tables(specs, reps, r_eff, eval_workers())
    for k in ("v1", "hyperscore_v1", "hyperscore_v2"):
        res[k] = {"overall": tabs[(k, "overall")]}
    v2 = {"overall": tabs[("v2_sep", "overall")],
          "merged": {lab: t for (kind, lab), t in tabs.items() if kind == "merged"},
          "raw": {kk: t for (kind, kk), t in tabs.items() if kind == "raw"}}
    res["v2_sep"] = v2
    res["merge_maps"] = merge_maps
    res["merge_map_identical_across_replicates"] = same
    res["complete"] = True
    return res


def violated(res: dict) -> dict:
    """Void inputs (A5.7.1.2 / A5.7.2): gated 'violated' under convention (a); peptidoform; v1 overall; v2_sep overall +
    merged blocks.  Reported for the claim script; it decides nothing here."""
    hits = {"v1": [], "v2_sep": []}
    if not res.get("complete"):
        return hits
    lists = {"v1": {"overall": res["v1"]["overall"]},
             "v2_sep": {"overall": res["v2_sep"]["overall"], **{f"merged:{k}": v for k, v in res["v2_sep"]["merged"].items()}}}
    for sc, ls in lists.items():
        for ln, L in ls.items():
            for a in GATED:
                if L["peptidoform"][str(a)]["a"]["verdict"] == "violated":
                    hits[sc].append(f"{ln}@{a}")
    return hits


# ------------------------------------------------------------------------------------------------ outputs
def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def write_manifest(out):
    lines = []
    for dp, _, fs in os.walk(out):
        for f in fs:
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, out)
            if rel == "manifest.sha256":
                continue
            lines.append(f"{_sha(p)}  {rel}")
    lines.sort(key=lambda l: l.split("  ", 1)[1])
    with open(os.path.join(out, "manifest.sha256"), "w") as fh:
        fh.write("\n".join(lines) + ("\n" if lines else ""))
    return os.path.join(out, "manifest.sha256")


def clean_work(work):
    """Sealed: delete every sealed-derived file under work/ (stores, mzML, MS1, PSM tables, caches); keep db_*/ (built
    from FASTAs only) and db/ (target FASTA + arm-C peptide lists)."""
    for n in os.listdir(work):
        if n == "db" or n.startswith("db_"):
            continue
        p = os.path.join(work, n)
        shutil.rmtree(p) if os.path.isdir(p) and not os.path.islink(p) else os.remove(p)


def main(argv=None):
    a = parse_args(argv)
    st = check_args(a)
    st["mode"] = a.mode
    st["workers"] = a.workers
    staged_python()
    sealed = a.mode == "sealed"
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    logf = open(os.path.join(a.work, "t5_armM.log"), "a")

    def log(msg):
        line = f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}"
        logf.write(line + "\n")
        logf.flush()
        if not sealed:
            print(line, flush=True)
        elif CONSOLE_RE.match(f"t5_armM: {msg}"):        # sealed console: data-free whitelisted lines only
            print(f"t5_armM: {msg}", flush=True)

    t0 = time.time()
    sha = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    reps, dbinfo = prepare_db(st, a.work, log)
    stores, runs, first = stores_for(st, a.work, reps, log)
    pinfo = run_pipelines(st, a.work, reps, stores, first, log)
    res = evaluate(st, reps, runs, a.work, dbinfo, log)
    log("stage evaluate")
    res.update({"git_sha": sha, "db": dbinfo, "recipe": st["recipe"], "instrument": INSTRUMENT,
                "violated_gated_convention_a": violated(res), "v2_sha256": pinfo.get("v2_sha256"),
                "v1_cohort_identity": {r: v.get("identical_to_cohort_output") for r, v in pinfo["v1_search"].items()
                                       if "identical_to_cohort_output" in v},
                "wall_s": time.time() - t0, "estimator": "combined FDP_M = N_E (1 + 1/r) / (N_T + N_E)",
                "deciding_convention": "a"})
    json.dump(res, open(os.path.join(a.out, "armM.json"), "w"), indent=1, default=float)
    if sealed:
        clean_work(a.work)
        log("stage cleanup")
    write_manifest(a.out)
    if not res["complete"]:
        log("stage incomplete")
        sys.exit(3)
    log("wrote armM.json and manifest.sha256")
    log("complete")


if __name__ == "__main__":
    main()
