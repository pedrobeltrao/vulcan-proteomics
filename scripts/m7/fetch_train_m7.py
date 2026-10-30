#!/usr/bin/env python3
"""Fetch the NEW M7 training runs (data/manifests/train_m7.tsv, source=m7_new) to data/raw/train/<acc>/.
Resumable (HTTP Range), size AND PRIDE SHA-1 verified (LESSONS 6; a mismatch deletes and re-fetches once),
idempotent.  Writes data/manifests/train_m7.sha256 (sha256 of every M7 training raw, M3 rows included).

    python3 scripts/m7/fetch_train_m7.py [--jobs 6]
"""
import argparse, csv, hashlib, os, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
import fetch_train as FT  # noqa: E402

DATA = os.path.join(ROOT, "data")


def hfile(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def fetch(r):
    path = FT.dest(r)
    for attempt in range(2):
        msg = FT.download(r)
        if not r.get("pride_sha1"):
            return msg + " (no PRIDE sha1)"
        if hfile(path, "sha1") == r["pride_sha1"]:
            return msg + " sha1 OK"
        print(f"  sha1 mismatch {r['id']} attempt {attempt}; deleting", file=sys.stderr)
        os.remove(path)
    raise RuntimeError(f"{r['id']} sha1 mismatch twice")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(DATA, "manifests", "train_m7.tsv")), delimiter="\t"))
    new = [r for r in rows if r["source"] == "m7_new"]
    failed = []
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(fetch, r): r for r in new}
        for f, r in futs.items():
            try:
                print(f.result(), flush=True)
            except Exception as e:  # noqa: BLE001
                failed.append(r["id"]); print("FAILED", r["id"], e, flush=True)
    lines = [f"{hfile(FT.dest(r), 'sha256')}  {os.path.relpath(FT.dest(r), DATA)}" for r in rows
             if os.path.exists(FT.dest(r))]
    open(os.path.join(DATA, "manifests", "train_m7.sha256"), "w").write("\n".join(lines) + "\n")
    print("wrote train_m7.sha256;", len(lines), "files; failed:", failed)


if __name__ == "__main__":
    main()
