#!/usr/bin/env python3
"""M6d: run the unchanged M6 front end (ingest -> MS1 features -> GPU search -> frozen static_gbm_pred -> pooled paired
TDC) on an extra dev set (D11 / D12) and register it for scripts/m6d (dataset.json).  Work dir /tmp/vulcan_m6d/<name>.
    python3 scripts/m6d/run_set.py --name d11 --raws R1.raw ... --A run ... --B run ... --expected '{"HUMAN":0,...}'
           [--C run ...] [--fasta F] [--index DIR] [--sel run run run] [--jobs 6]
Runs listed under --C are extra runs of the same experiment (e.g. single-proteome MBR acceptors in D12)."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
VPY = os.path.join(ROOT, "scripts", "m2b", "vpy.sh")
PE = f"/tmp/vulcan_env_{os.getuid()}/pe/bin/python"
PRED_SITE = "/tmp/vulcan_m6/pred_site"


def sh(cmd, env=None, log=None):
    t0 = time.perf_counter()
    with open(log, "w") if log else open(os.devnull, "w") as fh:
        rc = subprocess.call(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if rc:
        raise RuntimeError(f"failed ({rc}): {' '.join(cmd)}  log {log}")
    return time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--raws", nargs="+", required=True)
    ap.add_argument("--A", nargs="+", required=True)
    ap.add_argument("--B", nargs="+", default=[])
    ap.add_argument("--C", nargs="*", default=[])
    ap.add_argument("--expected", required=True)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--index", default="/tmp/vulcan_m6/index/hye_prod")
    ap.add_argument("--sel", nargs="*", default=[])
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    W = f"/tmp/vulcan_m6d/{a.name}"
    for d in ("store", "ms1", "psms", "logs", "rawstage", "ids"):
        os.makedirs(f"{W}/{d}", exist_ok=True)
    T = {}
    env = dict(os.environ, VULCAN_M6_WORK=W)

    def conv(raw):
        b = os.path.basename(raw)[:-4]
        if os.path.exists(f"{W}/store/{b}.vspec.arrow"):
            return 0.0
        st = f"{W}/rawstage/{b}.raw"
        shutil.copy(raw, st)            # never decode a RAW in place on Lustre (LESSONS)
        t = sh([VPY, "-m", "vulcan.io.convert", st, f"{W}/store/", "--threads", "2"], log=f"{W}/logs/convert_{b}.log")
        os.remove(st)
        return t

    def ms1(raw):
        b = os.path.basename(raw)[:-4]
        if os.path.exists(f"{W}/ms1/{b}.features.arrow"):
            return 0.0
        return sh([VPY, "-m", "vulcan.ms1", f"{W}/store/{b}.vspec.arrow", f"{W}/ms1", "--threads", "4"],
                  log=f"{W}/logs/ms1_{b}.log")
    with ThreadPoolExecutor(a.jobs) as ex:
        T["convert_s"] = list(ex.map(conv, a.raws))
    with ThreadPoolExecutor(a.jobs) as ex:
        T["ms1_s"] = list(ex.map(ms1, a.raws))
    index = a.index
    if a.fasta:
        index = f"{W}/index"
        if not os.path.exists(f"{index}/meta.json"):
            code = f"import sys; sys.path.insert(0,'{ROOT}'); from vulcan.index.build import build; print(build(['{a.fasta}'], '{index}', {{}}))"
            T["index_s"] = sh([VPY, "-c", code], log=f"{W}/logs/index.log")
    stores = [f"{W}/store/{os.path.basename(r)[:-4]}.vspec.arrow" for r in a.raws]
    if not os.path.exists(f"{W}/psms/done"):
        T["search_s"] = sh([VPY, "-m", "vulcan.gpu.run", "--index", index, "--stores", *stores, "--out", f"{W}/psms",
                            "--params", json.dumps({"report_psms": 5})], env=env, log=f"{W}/logs/search.log")
        open(f"{W}/psms/done", "w").close()
    if not os.path.exists(f"{W}/psms/score_info.json"):
        envp = dict(env, PYTHONPATH=f"{ROOT}:{PRED_SITE}")
        T["score_s"] = sh([PE, f"{ROOT}/scripts/m6/score.py", "--psms", f"{W}/psms", "--index", index,
                           "--dataset", a.name.upper()], env=envp, log=f"{W}/logs/score.log")
    runs = [os.path.basename(r)[:-4] for r in a.raws]
    if not os.path.exists(f"{W}/ids/all/ids.parquet"):
        T["ids_all_s"] = sh([VPY, f"{ROOT}/scripts/m6/ids.py", "--psms", f"{W}/psms", "--index", index, "--out",
                             f"{W}/ids/all", "--runs", *runs], env=env, log=f"{W}/logs/ids_all.log")
    ab = a.A + a.B
    if a.C and not os.path.exists(f"{W}/ids/ab/ids.parquet"):
        T["ids_ab_s"] = sh([VPY, f"{ROOT}/scripts/m6/ids.py", "--psms", f"{W}/psms", "--index", index, "--out",
                            f"{W}/ids/ab", "--runs", *ab], env=env, log=f"{W}/logs/ids_ab.log")
    ds = {"name": a.name, "A": a.A, "B": a.B, "expected": json.loads(a.expected),
          "ids": f"{W}/ids/ab" if a.C else f"{W}/ids/all", "ms1": f"{W}/ms1", "store_dir": f"{W}/store",
          "sel": a.sel, "extra": {"C": a.C, "ids_all": f"{W}/ids/all", "index": index}}
    json.dump(ds, open(f"{W}/dataset.json", "w"), indent=1)
    json.dump(T, open(f"{W}/frontend_timing.json", "w"), indent=1)
    print(json.dumps(T))


if __name__ == "__main__":
    main()
