#!/usr/bin/env python3
"""A1.1 confirmation set: the first N (=10) runs of the E9/D8 manifest (data/manifests/e9_d8.tsv, manifest order) are
downloaded independently of the E9 streamer (its files are never touched) to $VULCAN_M3E_WORK/raw_d8/ (node-local,
re-fetchable), size-checked, and their sha256 compared with the E9 record (data/manifests/e9_d8.sha256).
    python3 scripts/m3e/fetch_d8.py [--n 10] [--jobs 5]"""
import argparse, csv, hashlib, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WORK = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
SAFE = __import__("re").compile(r"[^A-Za-z0-9._-]")


def run_name(r):
    return f"{r['id']}__{SAFE.sub('_', r['filename'])[:-4]}"


def get(r):
    out = os.path.join(WORK, "raw_d8", run_name(r) + ".raw")
    if not (os.path.exists(out) and os.path.getsize(out) == int(r["expected_size"])):
        subprocess.run(["curl", "-sS", "-L", "--retry", "5", "-C", "-", "-o", out, r["url"]], check=True)
    assert os.path.getsize(out) == int(r["expected_size"]), out
    h = hashlib.sha256()
    with open(out, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return run_name(r), h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=5)
    a = ap.parse_args()
    os.makedirs(os.path.join(WORK, "raw_d8"), exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "manifests", "e9_d8.tsv")), delimiter="\t"))[:a.n]
    ref = {}
    sp = os.path.join(ROOT, "data", "manifests", "e9_d8.sha256")
    for line in open(sp):
        h, p = line.split()
        ref[os.path.basename(p)[:-4]] = h
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(get, rows))
    rep = {n: {"sha256": h, "e9_sha256": ref.get(n), "match": ref.get(n) == h} for n, h in res}
    json.dump(rep, open(os.path.join(WORK, "raw_d8", "verify.json"), "w"), indent=1)
    print(json.dumps(rep, indent=1))
    assert all(v["match"] for v in rep.values()), "sha256 mismatch vs E9 record"


if __name__ == "__main__":
    main()
