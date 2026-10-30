#!/usr/bin/env python3
"""E9 stream-and-delete driver (DESIGN §7 E9, tree 3d).  Stdlib only (this process never imports torch).

    python3 scripts/e9/stream.py --start 0 --end 100 [--streams 32] [--batch 50] [--workers 24]

Pipeline (all on local NVMe under --work, default /tmp/vulcan_e9; keep < 400 GB):
  1. downloader: --streams parallel curl transfers of E9/D8 manifest rows (data/manifests/e9_d8.tsv) to work/raw/,
     bounded by --max-pending-gb of downloaded-but-unprocessed RAW; size + PRIDE SHA-1 verified; sha256 appended to
     data/manifests/e9_d8.sha256 (one line per run, `sha256  e9/<id>__<file>`).
  2. batcher: every --batch ready runs -> one `python -m vulcan.cohort` process (staged env, --keep-stores,
     --skip-failed, report_psms=5) against the arm-P augmented index (T + E_2001 + decoys), then the RAWs are deleted,
     then `python -m vulcan.gpu.run` on the kept stores against the arm-C index (T + Arabidopsis foreign + decoys).
  3. scorer: a persistent `scripts/e9/score_worker.py` (pred env) takes each finished batch from work/queue/, builds
     L3 features, applies the frozen default scorer (static_gbm_pred, D7), writes compact per-run tables to
     work/out/{P,C}/ and deletes the stores.
Every stage appends timing records to work/timeline.jsonl.  Re-running skips runs already in work/out/P (resume).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests", "e9_d8.tsv")
SHA = os.path.join(ROOT, "data", "manifests", "e9_d8.sha256")
IDX_P = "/tmp/vulcan_e2/index/human_contam_E2001"
IDX_C = "/tmp/vulcan_m3/index/human_contam_C"
PARAMS = {"report_psms": 5}

_lock = threading.Lock()


def tlog(work, rec):
    rec["t"] = time.time()
    with _lock, open(os.path.join(work, "timeline.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def run_name(row):
    return f"{row['id']}__{safe(row['filename'])[:-4]}"


class Budget:
    def __init__(self, cap):
        self.cap, self.used, self.cv = cap, 0, threading.Condition()

    def take(self, n):
        with self.cv:
            while self.used > 0 and self.used + n > self.cap:
                self.cv.wait()
            self.used += n

    def give(self, n):
        with self.cv:
            self.used -= n
            self.cv.notify_all()


def _fetch(row, part, size):
    """curl with resume; returns (ok, err)."""
    err = ""
    for attempt in range(5):
        r = subprocess.run(["curl", "-s", "-S", "--retry", "3", "--connect-timeout", "60", "--speed-time", "120",
                            "--speed-limit", "20000", "-C", "-", "-o", part, row["url"]], capture_output=True, text=True)
        if os.path.exists(part) and os.path.getsize(part) == size:
            return True, ""
        err = (r.stderr or "")[-300:]
        if os.path.exists(part) and os.path.getsize(part) > size:
            os.remove(part)
        time.sleep(5 * (attempt + 1))
    return False, err


def _hash(path):
    h256, h1 = hashlib.sha256(), hashlib.sha1()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h256.update(b)
            h1.update(b)
    return h256.hexdigest(), h1.hexdigest()


def download(row, work, budget, ready, failed):
    """Size AND PRIDE SHA-1 verified download.  A SHA-1 mismatch (seen on ~5% of resumed transfers through the cluster
    proxy, and the cause of most 'ingest failed' runs) deletes the file and re-downloads from scratch (<= 3 times)."""
    size = int(row["expected_size"])
    budget.take(size)
    path = os.path.join(work, "raw", run_name(row) + ".raw")
    part = path + ".part"
    if os.path.exists(path) and os.path.getsize(path) == size:     # already downloaded + verified (resume)
        tlog(work, {"ev": "download_reuse", "id": row["id"], "bytes": size})
        ready.put((row, path, size))
        return
    t0 = time.time()
    for full in range(3):
        ok, err = _fetch(row, part, size)
        if not ok:
            break
        t1 = time.time()
        s256, s1 = _hash(part)
        sha1_ok = (not row["pride_sha1"]) or s1 == row["pride_sha1"]
        if sha1_ok:
            break
        tlog(work, {"ev": "sha1_mismatch", "id": row["id"], "attempt": full})
        os.remove(part)
        ok, err = False, "PRIDE SHA-1 mismatch"
    t_dl = time.time() - t0
    if not ok:
        failed.append({"stage": "download", "id": row["id"], "error": err})
        if os.path.exists(part):
            os.remove(part)
        budget.give(size)
        tlog(work, {"ev": "download_failed", "id": row["id"], "wall_s": t_dl, "error": err})
        return
    os.replace(part, path)
    with _lock, open(SHA, "a") as fh:
        fh.write(f"{s256}  e9/{run_name(row)}.raw\n")
    tlog(work, {"ev": "download", "id": row["id"], "bytes": size, "wall_s": t_dl, "hash_s": time.time() - t1,
                "pride_sha1_match": True, "pride_sha1_present": bool(row["pride_sha1"])})
    ready.put((row, path, size))


def sh(cmd, log):
    with open(log, "a") as fh:
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode


def process_batch(k, items, work, a, budget):
    bdir = os.path.join(work, "batches", f"b{k:03d}")
    os.makedirs(bdir, exist_ok=True)
    stage = os.path.join(bdir, "stage")
    lst = os.path.join(bdir, "raws.txt")
    with open(lst, "w") as fh:
        fh.write("\n".join(p for _, p, _ in items) + "\n")
    vpy = os.path.join(ROOT, "scripts", "m2b", "vpy.sh")
    t0 = time.time()
    for attempt in range(4):   # a whole-cohort failure (e.g. GPU OOM at index build on a shared GPU) is retried
        rc = sh([vpy, "-m", "vulcan.cohort", "--index", a.index_p, "--raw-list", lst, "--out", os.path.join(bdir, "P"),
                 "--stage", stage, "--workers", str(a.workers), "--keep-stores", "--skip-failed",
                 "--params", json.dumps(PARAMS), "--stats", os.path.join(bdir, "cohort_stats.json")],
                os.path.join(bdir, "cohort.log"))
        if rc == 0:
            break
        tlog(work, {"ev": "cohort_retry", "batch": k, "attempt": attempt, "rc": rc})
        time.sleep(60 * (attempt + 1))
    t_cohort = time.time() - t0
    for _, p, size in items:
        if os.path.exists(p):
            os.remove(p)
        budget.give(size)
    stores = sorted(os.path.join(stage, f) for f in os.listdir(stage) if f.endswith(".vspec.arrow")) \
        if os.path.isdir(stage) else []
    # stage dir also holds staged RAW copies only transiently; drop anything that is not a store
    t1 = time.time()
    rc_c = None
    if stores and not a.no_c:
        rc_c = sh([vpy, "-m", "vulcan.gpu.run", "--index", a.index_c, "--stores", *stores, "--out",
                   os.path.join(bdir, "C"), "--params", json.dumps(PARAMS)], os.path.join(bdir, "csearch.log"))
    t_c = time.time() - t1
    tlog(work, {"ev": "batch_search", "batch": k, "n": len(items), "n_stores": len(stores), "cohort_rc": rc,
                "cohort_wall_s": t_cohort, "c_search_rc": rc_c, "c_search_wall_s": t_c})
    q = os.path.join(work, "queue")
    os.makedirs(q, exist_ok=True)
    with open(os.path.join(q, f"b{k:03d}.json.tmp"), "w") as fh:
        json.dump({"batch": k, "dir": bdir, "stores": stores, "ids": [r["id"] for r, _, _ in items]}, fh)
    os.replace(os.path.join(q, f"b{k:03d}.json.tmp"), os.path.join(q, f"b{k:03d}.json"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=100)
    ap.add_argument("--work", default="/tmp/vulcan_e9")
    ap.add_argument("--streams", type=int, default=32)
    ap.add_argument("--batch", type=int, default=50)
    ap.add_argument("--first-batch", type=int, default=None, help="size of the first batch (fast first result)")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--max-pending-gb", type=float, default=150)
    ap.add_argument("--index-p", default=IDX_P)
    ap.add_argument("--index-c", default=IDX_C)
    ap.add_argument("--no-c", action="store_true")
    ap.add_argument("--max-unscored", type=int, default=3, help="max searched batches waiting for the scorer")
    a = ap.parse_args()
    for d in ("raw", "batches", "queue", "out"):
        os.makedirs(os.path.join(a.work, d), exist_ok=True)
    rows = list(csv.DictReader(open(MAN), delimiter="\t"))[a.start:a.end]
    def scored(arm):
        od = os.path.join(a.work, "out", arm)
        return {f.split("__")[0] for f in os.listdir(od) if f.endswith(".scored.parquet")} if os.path.isdir(od) else set()
    done = scored("P") & (scored("C") if not a.no_c else scored("P"))     # resume: a run is done when every arm is
    queued = set()
    todo = [r for r in rows if r["id"] not in done]
    kb = len([d for d in os.listdir(os.path.join(a.work, "batches")) if d.startswith("b")])
    tlog(a.work, {"ev": "stream_start", "start": a.start, "end": a.end, "todo": len(todo), "streams": a.streams,
                  "batch": a.batch, "workers": a.workers})
    budget = Budget(a.max_pending_gb * 1e9)
    ready: queue.Queue = queue.Queue()
    failed: list = []
    t0 = time.time()
    ex = ThreadPoolExecutor(a.streams)
    futs = [ex.submit(download, r, a.work, budget, ready, failed) for r in todo]
    pending = []
    first = a.first_batch or a.batch
    while True:
        try:
            pending.append(ready.get(timeout=5))
        except queue.Empty:
            pass
        all_dl = all(f.done() for f in futs) and ready.empty()
        if all_dl and not pending:
            break
        want = first if not queued else a.batch
        if len(pending) >= want or (all_dl and pending):
            # backpressure: stores of searched-but-unscored batches stay on NVMe; cap the backlog
            qd = os.path.join(a.work, "queue")
            while sum(not os.path.exists(os.path.join(qd, f[:-5] + ".done"))
                      for f in os.listdir(qd) if f.endswith(".json")) >= a.max_unscored:
                time.sleep(5)
            batch, pending = pending[:want], pending[want:]
            process_batch(kb, batch, a.work, a, budget)
            queued.add(kb)
            kb += 1
    ex.shutdown()
    tlog(a.work, {"ev": "stream_end", "wall_s": time.time() - t0, "failed": failed})
    print(json.dumps({"wall_s": time.time() - t0, "failed": failed}))


if __name__ == "__main__":
    sys.exit(main())
