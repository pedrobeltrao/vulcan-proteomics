#!/usr/bin/env python3
"""Ingest RAW -> spectral store for M3e (D7 remaining runs, DN null runs, D8 confirmation runs).  RAWs on Lustre are
first copied to local NVMe (LESSONS: never let RawFileReader mmap Lustre), converted with vulcan.io.convert, copy removed.
    python scripts/m3e/ingest.py --set d7|dn|d8 [--jobs 4]"""
import argparse, csv, glob, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
M3NULL = "/tmp/vulcan_m3/store_null"


def todo(which):
    if which == "d7":
        return [(os.path.join(ROOT, "data", "raw", "PXD031616", f"LM{n}_DDA.raw"), M3NULL) for n in (2487, 2490, 2492)]
    if which == "dn":
        rows = csv.DictReader(open(os.path.join(ROOT, "data", "manifests", "g2null.tsv")), delimiter="\t")
        return [(os.path.join(ROOT, "data", "raw", r["accession"], r["filename"]), os.path.join(W, "store_dn")) for r in rows]
    if which == "d8":
        return [(p, os.path.join(W, "store_d8")) for p in sorted(glob.glob(os.path.join(W, "raw_d8", "*.raw")))]
    raise SystemExit(which)


def conv(job):
    raw, outdir = job
    os.makedirs(outdir, exist_ok=True)
    base = os.path.basename(raw)[:-4]
    out = os.path.join(outdir, base + ".vspec.arrow")
    if os.path.exists(out) and os.path.exists(out + ".json"):
        return out, "exists"
    src = raw
    if not raw.startswith("/tmp/"):
        os.makedirs(os.path.join(W, "rawcopy"), exist_ok=True)
        src = os.path.join(W, "rawcopy", os.path.basename(raw))
        shutil.copyfile(raw, src)
    r = subprocess.run([sys.executable, "-m", "vulcan.io.convert", src, outdir, "--threads", "4"], cwd=ROOT,
                       env={**os.environ, "PYTHONPATH": ROOT}, capture_output=True, text=True, timeout=3600)
    if src != raw:
        os.remove(src)
    if r.returncode:
        return out, "FAILED " + r.stderr[-500:]
    return out, "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    jobs = [j for j in todo(a.set) if os.path.exists(j[0])]
    with ThreadPoolExecutor(a.jobs) as ex:
        for out, st in ex.map(conv, jobs):
            print(out, st, flush=True)


if __name__ == "__main__":
    main()
