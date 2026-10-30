#!/usr/bin/env python3
"""Fetch the static-scorer TRAINING runs (data/manifests/train.tsv, role=train) to data/raw/train/<acc>/ and the
mouse Swiss-Prot FASTA to data/fasta/.  Stdlib only; resumable (HTTP Range on .part), size-verified, idempotent.
Writes data/manifests/train.sha256 (sha256sum -c compatible from within data/).

    python3 scripts/m3/fetch_train.py [--jobs 6] [--check]
"""
import argparse, csv, hashlib, os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "data")
UA = {"User-Agent": "vulcan-m3-train-fetch (research)"}
FASTAS = {"mouse_sp_canonical.fasta":
          "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:10090+AND+reviewed:true&format=fasta"}


def dest(r):
    return os.path.join(DATA, "raw", "train", r["accession"], r["filename"])


def download(r):
    path, size = dest(r), int(r["expected_size"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and os.path.getsize(path) == size:
        return f"[skip] {r['id']}"
    part = path + ".part"
    for attempt in range(8):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        try:
            req = urllib.request.Request(r["url"], headers={**UA, **({"Range": f"bytes={have}-"} if have else {})})
            with urllib.request.urlopen(req, timeout=120) as resp, open(part, "ab" if have else "wb") as fh:
                if have and resp.status != 206:
                    fh.truncate(0)
                while True:
                    b = resp.read(1 << 22)
                    if not b:
                        break
                    fh.write(b)
            if os.path.getsize(part) == size:
                os.replace(part, path)
                return f"[ok] {r['id']} {size/1e9:.2f} GB"
        except Exception as e:  # noqa: BLE001
            print(f"  retry {r['id']} {attempt}: {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{r['id']} failed")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(DATA, "manifests", "train.tsv")), delimiter="\t"))
    if a.check:
        for r in rows:
            p = dest(r)
            print("OK  " if os.path.exists(p) and os.path.getsize(p) == int(r["expected_size"]) else "MISS", r["id"])
        return
    for name, url in FASTAS.items():
        p = os.path.join(DATA, "fasta", name)
        if not os.path.exists(p):
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=600) as resp:
                open(p + ".part", "wb").write(resp.read())
            os.replace(p + ".part", p)
    with ThreadPoolExecutor(a.jobs) as ex:
        for m in ex.map(download, rows):
            print(m, flush=True)
    lines = [f"{sha(dest(r))}  {os.path.relpath(dest(r), DATA)}" for r in rows]
    lines += [f"{sha(os.path.join(DATA, 'fasta', n))}  fasta/{n}" for n in FASTAS]
    open(os.path.join(DATA, "manifests", "train.sha256"), "w").write("\n".join(lines) + "\n")
    print("wrote train.sha256")


if __name__ == "__main__":
    main()
