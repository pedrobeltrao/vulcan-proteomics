#!/usr/bin/env python3
"""M7 training-data preparation (label side is L4-only; the encoder never sees labels.arrow).

Stages (idempotent, resumable; work dir $VULCAN_M7_WORK, default /tmp/vulcan_m7):
  --ingest    new M7 runs: RAW copied to local NVMe, vulcan.io.convert -> store (LESSONS 4); M3 runs reuse
              /tmp/vulcan_m3/store_train.  Runs whose MS2 are < 95% FTMS+HCD are excluded (OT-HCD class only).
  --search    per species prod index (T + pseudo-reverse decoys), report_psms = 20 (hard negatives come from ranks 2..20)
  --features  L3 predicted features (vulcan.rescore.features, AlphaPeptDeep, per-run label-blind NCE choice)
  --label     frozen static_gbm_pred score -> paired double competition per run; positives = target winners at
              q <= 0.1% (spec: "high-confidence target PSMs (q <= 0.1%) from our pipeline").  Writes
              <run>.label.parquet (spec_idx, pf_row of the positive, q).

    VULCAN_M7_WORK=/tmp/vulcan_m7 python scripts/m7/prep_train.py --ingest --search --features --label
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M7_WORK", "/tmp/vulcan_m7")
REPORT = 20
Q_POS = 0.001
FAM = {"human": "human_contam", "mouse": "mouse_contam", "ecoli": "ecoli_contam"}


def rows():
    return list(csv.DictReader(open(os.path.join(ROOT, "data", "manifests", "train_m7.tsv")), delimiter="\t"))


def run_name(r):
    return r["filename"][:-4]


def store_of(r):
    if r["source"] == "m3_train":
        return C.train_store(r)
    return os.path.join(W, "store", run_name(r) + ".vspec.arrow")


def psm_dir(sp):
    return os.path.join(W, "psms", f"train_{sp}")


def ot_hcd_fraction(store):
    import pyarrow.ipc as ipc
    t = ipc.open_file(store).read_all().select(["ms_level", "activation", "analyzer"])
    m = t.column("ms_level").to_numpy() == 2
    act = t.column("activation").to_numpy(zero_copy_only=False)[m]
    ana = t.column("analyzer").to_numpy(zero_copy_only=False)[m]
    return float(np.mean((act == 2) & (ana == 1)))


def ingest_one(r):
    out = store_of(r)
    if os.path.exists(out):
        return f"[skip] {r['id']}"
    raw = os.path.join(ROOT, "data", "raw", "train", r["accession"], r["filename"])
    if not os.path.exists(raw):
        return f"[missing raw] {r['id']}"
    bad = os.path.join(W, "store", run_name(r) + ".FAILED")
    if os.path.exists(bad):
        return f"[failed earlier] {r['id']}"
    tmp = os.path.join(W, "rawcopy", r["filename"])
    os.makedirs(os.path.dirname(tmp), exist_ok=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    shutil.copyfile(raw, tmp)
    t0 = time.time()
    try:
        subprocess.run([sys.executable, "-m", "vulcan.io.convert", tmp, os.path.dirname(out), "--threads", "4"],
                       check=True, cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT}, timeout=1800,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:  # noqa: BLE001  (unreadable RAW: recorded, run excluded)
        open(bad, "w").write(repr(e))
        return f"[FAILED] {r['id']} {e!r}"[:300]
    finally:
        os.remove(tmp)
    return f"[ok] {r['id']} {time.time() - t0:.0f}s"


def ingest(jobs):
    from concurrent.futures import ThreadPoolExecutor
    todo = [r for r in rows() if r["source"] == "m7_new"]
    with ThreadPoolExecutor(jobs) as ex:
        for m in ex.map(ingest_one, todo):
            print(m, flush=True)


def class_ok():
    path = os.path.join(W, "class_check.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for r in rows():
        s = store_of(r)
        if run_name(r) not in res and os.path.exists(s):
            res[run_name(r)] = ot_hcd_fraction(s)
    json.dump(res, open(path, "w"), indent=1)
    return {k for k, v in res.items() if v >= 0.95}, res


def search():
    import pyarrow.parquet as pq
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    ok, _ = class_ok()
    p = E.SearchParams(report_psms=REPORT)
    for sp, fam in FAM.items():
        todo = [r for r in rows() if r["species"] == sp and run_name(r) in ok and os.path.exists(store_of(r))
                and not os.path.exists(os.path.join(psm_dir(sp), run_name(r) + ".psms.parquet"))]
        if not todo:
            continue
        os.makedirs(psm_dir(sp), exist_ok=True)
        l1 = L1Index(C.index_dir(fam, "prod"))
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        for r in todo:
            spc = E.load_store(store_of(r))
            tbl, st, _ = E.search(spc, gi, l1, p)
            pq.write_table(tbl, os.path.join(psm_dir(sp), run_name(r) + ".psms.parquet"))
            print(sp, run_name(r), len(spc), tbl.num_rows, round(st.wall_s["total"], 2), flush=True)
        del gi
        torch.cuda.empty_cache()


def features():
    import pyarrow as pa
    import pyarrow.parquet as pq
    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    pred = APD()
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PredCache(os.path.join(W, "predcache"), pred, nce, "QE", C.STATIC, C.VAR)
        return caches[nce]
    fb = FeatureBuilder(factory, C.STATIC, C.VAR)
    by_run = {run_name(r): r for r in rows()}
    for sp, fam in FAM.items():
        l1 = None
        for f in sorted(glob.glob(os.path.join(psm_dir(sp), "*.psms.parquet"))):
            run = os.path.basename(f)[:-len(".psms.parquet")]
            out = os.path.join(psm_dir(sp), run + ".feat.parquet")
            if os.path.exists(out):
                continue
            l1 = l1 or L1Index(C.index_dir(fam, "prod"))
            t = pq.read_table(f)
            psm = {k: t.column(k).to_numpy() for k in ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt")}
            feats, info = fb.run(l1, psm, load_store(store_of(by_run[run])))
            pq.write_table(pa.table(feats), out + ".partial")
            os.replace(out + ".partial", out)
            print(sp, run, info["nce"], json.dumps(info["cache_stats"], default=float)[:200], flush=True)


def label():
    import m3_data as MD
    import pyarrow as pa
    import pyarrow.parquet as pq
    from vulcan.rescore.static import StaticModel
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
    summ = {}
    for sp, fam in FAM.items():
        if not os.path.isdir(psm_dir(sp)):
            continue
        lab = MD.labels(C.index_dir(fam, "prod"))
        for run, psm, feat in MD.load_runs(psm_dir(sp)):
            out = os.path.join(psm_dir(sp), run + ".label.parquet")
            if os.path.exists(out):
                continue
            rows_ = psm["pf_row"]
            unit = psm["pf_id"]
            s = model.apply(psm, feat)
            r = double_competition(psm["spec_idx"].astype(np.int64), unit, lab["tdc_pair"][rows_], s,
                                   lab["is_decoy"][rows_], seed=11, tiebreak=hash_tiebreak(unit, 11))
            # positives: PSM-level winners (best candidate of the spectrum, any label) that are targets, whose
            # peptidoform is accepted at peptide-level q <= 0.1% and whose PSM score >= that list's score threshold
            pos = (~r.is_decoy) & (r.q <= Q_POS)
            if pos.any():
                thr = float(r.score[pos].min())
                acc_units = set(unit[r.row[pos]].tolist())
                w1 = r.psm_winners
                keep = (~lab["is_decoy"][rows_[w1]]) & (s[w1] >= thr) & np.array([u in acc_units for u in unit[w1]], bool)
                ri = w1[keep]
            else:
                ri = np.zeros(0, np.int64)
            pq.write_table(pa.table({"spec_idx": psm["spec_idx"][ri].astype(np.int64), "row": ri.astype(np.int64),
                                     "pf_row": rows_[ri].astype(np.int64), "score": s[ri].astype(np.float64)}), out)
            summ[run] = {"species": sp, "n_spectra_with_psm": int(np.unique(psm["spec_idx"]).size),
                         "n_peptidoforms_q0.001": int(pos.sum()), "n_pos_psms": int(len(ri)),
                         "n_target_peptidoforms_q0.01": int(((~r.is_decoy) & (r.q <= 0.01)).sum())}
            print(run, summ[run], flush=True)
    p = os.path.join(W, "label_summary.json")
    old = json.load(open(p)) if os.path.exists(p) else {}
    old.update(summ)
    json.dump(old, open(p, "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", action="store_true")
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--features", action="store_true")
    ap.add_argument("--label", action="store_true")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    if a.ingest:
        ingest(a.jobs)
    if a.search:
        search()
    if a.features:
        features()
    if a.label:
        label()


if __name__ == "__main__":
    main()
