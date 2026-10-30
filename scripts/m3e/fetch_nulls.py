#!/usr/bin/env python3
"""Fetch the A1.3 null-source data (dev set DN): RAW runs of data/manifests/g2null.tsv via scripts/data/fetch_tier_s.py
(-> data/raw/<acc>/, sha256 -> data/manifests/g2null.sha256) and the homolog-filter proteomes (UniProt REST, all entries
of each proteome) -> data/fasta/dn_<upid>.fasta, sha256 -> data/manifests/g2null_fasta.sha256.
    python3 scripts/m3e/fetch_nulls.py [--fasta-only] [--jobs 4]"""
import argparse, hashlib, os, runpy, subprocess, sys, urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
UPIDS = ["UP000008524", "UP000001939", "UP000008816", "UP000001688", "UP000001417"]


def fasta(upid):
    out = os.path.join(ROOT, "data", "fasta", f"dn_{upid}.fasta")
    if not os.path.exists(out):
        url = (f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=false&query=proteome%3A{upid}" if upid != "UP000001688" else "https://rest.uniprot.org/uniparc/stream?format=fasta&query=upid%3AUP000001688")
        with urllib.request.urlopen(url, timeout=300) as r, open(out + ".partial", "wb") as fh:
            fh.write(r.read())
        os.replace(out + ".partial", out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-only", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    with open(os.path.join(ROOT, "data", "manifests", "g2null_fasta.sha256"), "w") as fh:
        for u in UPIDS:
            p = fasta(u)
            n = sum(1 for line in open(p) if line.startswith(">"))
            fh.write(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  fasta/{os.path.basename(p)}  # {n} entries\n")
            print(u, n, flush=True)
    if not a.fasta_only:
        sys.argv = [sys.argv[0], "--manifest", "g2null", "--roles", "dev", "--jobs", str(a.jobs)]
        runpy.run_path(os.path.join(ROOT, "scripts", "data", "fetch_tier_s.py"), run_name="__main__")


if __name__ == "__main__":
    main()
