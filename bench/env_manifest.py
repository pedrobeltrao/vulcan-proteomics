#!/usr/bin/env python3
"""Environment manifests for the A5 look (DESIGN §6.0d; reviews/amendment_A5_draft.md A5.11.5 `hashes`).

A manifest is a sorted list of "<sha256>  <relative path>" lines for one tool or environment; its digest is the sha256
of the manifest text.  `config/a5.yaml env_manifests` records each digest at the freeze, and `bench/sealed_run.py
--look a5` recomputes every one before writing a ledger line and refuses on any mismatch.

What is hashed (implementation choice, recorded in scripts/a3/IMPLEMENTATION_NOTES.md):
  - a file (the Sage binary): its bytes;
  - a plain directory (Comet, TRFP): every regular file under it, by bytes;
  - a Python environment (a directory with conda-meta/ or a lib/python*/site-packages/): the files that pin every
    installed file -- conda-meta/*.json (conda's per-package file lists with sha256), every
    site-packages/*.dist-info/RECORD (pip's per-file sha256 lists), pyvenv.cfg, and bin/python*'s link targets -- plus
    the bytes of every top-level package's __init__.py.  Hashing all bytes of torch/CUDA libraries (several GB on
    Lustre) at look time would take tens of minutes; a package change rewrites its RECORD / conda-meta entry.
Symlinks are recorded by their target string, not followed.

    python bench/env_manifest.py <name>=<path> ...          # prints "<name> <digest>" per entry
"""
from __future__ import annotations

import glob
import hashlib
import os
import sys


def _sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _is_pyenv(d: str) -> bool:
    return os.path.isdir(os.path.join(d, "conda-meta")) or bool(glob.glob(os.path.join(d, "lib", "python*", "site-packages")))


def _walk(d: str):
    for dp, dns, fns in os.walk(d, followlinks=False):
        dns.sort()
        for f in sorted(fns):
            yield os.path.join(dp, f)


def manifest_lines(path: str) -> list:
    path = os.path.abspath(path)
    if not os.path.exists(path):
        return [f"MISSING  {path}"]
    if os.path.isfile(path):
        return [f"{_sha(path)}  {os.path.basename(path)}"]
    files = []
    if _is_pyenv(path):
        files += glob.glob(os.path.join(path, "conda-meta", "*.json"))
        files += glob.glob(os.path.join(path, "pyvenv.cfg"))
        files += glob.glob(os.path.join(path, "bin", "python*"))
        for sp in glob.glob(os.path.join(path, "lib", "python*", "site-packages")):
            files += glob.glob(os.path.join(sp, "*.dist-info", "RECORD"))
            files += glob.glob(os.path.join(sp, "*", "__init__.py"))
    else:
        files = list(_walk(path))
    out = []
    for f in sorted(set(files)):
        rel = os.path.relpath(f, path)
        if os.path.islink(f):
            out.append(f"link:{os.readlink(f)}  {rel}")
        elif os.path.isfile(f):
            out.append(f"{_sha(f)}  {rel}")
    return out


def digest(path: str) -> str:
    return hashlib.sha256(("\n".join(manifest_lines(path)) + "\n").encode()).hexdigest()


def main(argv):
    for a in argv:
        name, p = a.split("=", 1)
        print(name, digest(p), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
