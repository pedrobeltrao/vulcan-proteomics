#!/usr/bin/env python3
"""E7/v2 training data (pre-registered rule: results/2026-10-28_E7_chimeric_v2/config.yaml `training`).

  --freeze   write results/2026-10-28_E7_chimeric_v2/train_v2_folds.tsv: data/manifests/train.tsv (M3, 13 runs) + every
             data/manifests/train_m7.tsv row with source == m7_new whose RAW is complete (size == expected_size) NOW;
             folds: M3 projects keep their fold, new projects round-robin over 0..2 in manifest order; raw sha256.
  --prep     per frozen run: RAW copied to NVMe -> store (vulcan.io.convert) -> OT-HCD check (>= 0.95, as M3) ->
             vulcan.ms1 features + hypotheses -> chimeric search (species prod index, report 5) -> v2 features
             (incl. v1 features) -> $W/train/<run>.chim.{psms,feat}.parquet.  Labels are not read here.

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/e7/prep_train.py --freeze | --prep
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_E7_WORK", "/tmp/vulcan_e7")
RES = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2")
FOLDS = os.path.join(RES, "train_v2_folds.tsv")
MAN = os.path.join(ROOT, "data", "manifests")
INDEX = {"human": "/tmp/vulcan_e2/index/human_contam_prod", "mouse": "/tmp/vulcan_m3/index/mouse_contam_prod",
         "ecoli": "/tmp/vulcan_m3/index/ecoli_contam_prod"}
M3_STORE = "/tmp/vulcan_m3/store_train"
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}


def raw_path(r):
    return os.path.join(ROOT, "data", "raw", "train", r["accession"], r["filename"])


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def freeze():
    assert not os.path.exists(FOLDS), "training set already frozen"
    m3 = list(csv.DictReader(open(os.path.join(MAN, "train.tsv")), delimiter="\t"))
    m7 = [r for r in csv.DictReader(open(os.path.join(MAN, "train_m7.tsv")), delimiter="\t") if r["source"] == "m7_new"]
    rows = []
    for r in m3:
        rows.append({**{k: r[k] for k in ("accession", "filename", "species", "instrument")}, "fold": int(r["fold"]),
                     "source": "m3_train"})
    k = 0
    skipped = []
    for r in m7:
        p = raw_path(r)
        if os.path.exists(p) and os.path.getsize(p) == int(r["expected_size"]):
            rows.append({**{x: r[x] for x in ("accession", "filename", "species", "instrument")}, "fold": k % 3,
                         "source": "m7_new"})
            k += 1
        else:
            skipped.append(r["accession"])
    for r in rows:
        r["raw_sha256"] = sha256(raw_path(r))
    os.makedirs(RES, exist_ok=True)
    with open(FOLDS, "w") as fh:
        w = csv.DictWriter(fh, ["accession", "fold", "species", "instrument", "source", "filename", "raw_sha256"],
                           delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({x: r[x] for x in w.fieldnames})
    json.dump({"frozen_at": time.strftime("%Y-%m-%d %H:%M:%S"), "n_runs": len(rows),
               "m7_not_complete_at_freeze": skipped}, open(FOLDS + ".log.json", "w"), indent=1)
    print(f"frozen {len(rows)} runs ({k} new); not complete: {len(skipped)}")


def frozen():
    return list(csv.DictReader(open(FOLDS), delimiter="\t"))


def store_of(r):
    name = r["filename"][:-4]
    p = os.path.join(M3_STORE, name + ".vspec.arrow")
    if r["source"] == "m3_train" and os.path.exists(p):
        return p
    out = os.path.join(W, "store_train", name + ".vspec.arrow")
    if not os.path.exists(out):
        os.makedirs(os.path.join(W, "raw_train"), exist_ok=True)
        loc = os.path.join(W, "raw_train", r["filename"])
        if not os.path.exists(loc):
            shutil.copyfile(raw_path(r), loc + ".partial")
            os.replace(loc + ".partial", loc)
        subprocess.run([sys.executable, "-m", "vulcan.io.convert", loc, os.path.dirname(out), "--threads", "4"],
                       check=True, cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT}, stdout=subprocess.DEVNULL)
        os.remove(loc)
    return out


def is_ot_hcd(store):
    import pyarrow.ipc as ipc
    t = ipc.open_file(store).read_all().select(["ms_level", "activation", "analyzer"])
    m = t.column("ms_level").to_numpy() == 2
    act = t.column("activation").to_numpy(zero_copy_only=False)[m]
    ana = t.column("analyzer").to_numpy(zero_copy_only=False)[m]
    return float(np.mean((act == 2) & (ana == 1)))


def prep():
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch

    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features_v2 import FeatureBuilderV2, read_ms1_features
    out = os.path.join(W, "train")
    ms1 = os.path.join(W, "ms1", "train")
    os.makedirs(out, exist_ok=True)
    os.makedirs(ms1, exist_ok=True)
    cache_root = os.path.join(W, "predcache_train")
    if not os.path.exists(cache_root):          # M3's cache (values of every M3 training peptide preserved)
        shutil.copytree("/tmp/vulcan_m3/predcache", cache_root + ".partial")
        os.replace(cache_root + ".partial", cache_root)
    pred = APD()
    caches = {}

    def fac(nce):
        if nce not in caches:
            caches[nce] = PredCache(cache_root, pred, nce, "QE", STATIC, VAR)
        return caches[nce]
    fb = FeatureBuilderV2(fac, STATIC, VAR)
    p = E.SearchParams(report_psms=5)
    log = {}
    lp = os.path.join(out, "prep_log.json")
    if os.path.exists(lp):
        log = json.load(open(lp))
    rows = frozen()
    for sp_name in ("human", "mouse", "ecoli"):
        todo = [r for r in rows if r["species"] == sp_name and not os.path.exists(
            os.path.join(out, r["filename"][:-4] + ".chim.feat.parquet")) and r["filename"][:-4] not in log]
        if not todo:
            continue
        l1 = L1Index(INDEX[sp_name])
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        for r in todo:
            run = r["filename"][:-4]
            t0 = time.perf_counter()
            try:
                st = store_of(r)
            except subprocess.CalledProcessError as e:     # technical, label-independent exclusion (logged)
                log[run] = {"excluded": f"RAW not decodable by RawFileReader (vulcan.io.convert exit {e.returncode})"}
                json.dump(log, open(lp, "w"), indent=1)
                print(run, log[run], flush=True)
                continue
            frac = is_ot_hcd(st)
            if frac < 0.95:
                log[run] = {"excluded": f"not OT-HCD (FTMS+HCD fraction of MS2 = {frac:.2f})"}
                json.dump(log, open(lp, "w"), indent=1)
                print(run, log[run], flush=True)
                continue
            if not os.path.exists(os.path.join(ms1, run + ".hypotheses.arrow")):
                subprocess.run([sys.executable, "-m", "vulcan.ms1", st, ms1, "--threads", "8"], check=True, cwd=ROOT,
                               env={**os.environ, "PYTHONPATH": ROOT}, stdout=subprocess.DEVNULL)
            spx = E.load_store(st)
            g = CS.precursor_groups(CS.read_hypotheses(os.path.join(ms1, run + ".hypotheses.arrow")), spx)
            tbl, sst = CS.search_chimeric(spx, gi, l1, p, g)
            pq.write_table(tbl, os.path.join(out, run + ".chim.psms.parquet"))
            psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
            feats, finfo = fb.run_v2(l1, psm, spx, read_ms1_features(os.path.join(ms1, run + ".features.arrow")))
            pq.write_table(pa.table({k: np.asarray(v, np.float32) for k, v in feats.items()}),
                           os.path.join(out, run + ".chim.feat.parquet"))
            log[run] = {"species": sp_name, "accession": r["accession"], "fold": int(r["fold"]), "ot_hcd": frac,
                        "n_spectra": len(spx), "n_psms": tbl.num_rows, "n_secondary_groups": sst["n_secondary_groups"],
                        "nce": finfo["nce"], "wall_s": time.perf_counter() - t0}
            json.dump(log, open(lp, "w"), indent=1, default=float)
            print(run, log[run], flush=True)
        del gi, l1
        torch.cuda.empty_cache()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--freeze", action="store_true")
    ap.add_argument("--prep", action="store_true")
    a = ap.parse_args()
    if a.freeze:
        freeze()
    if a.prep:
        prep()
