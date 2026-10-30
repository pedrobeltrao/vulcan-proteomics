#!/usr/bin/env python3
"""Fill config/a5.yaml `hashes` and `env_manifests` (A5.11.5) from the committed tree at the freeze.

    python scripts/a5/fill_hashes.py [--check]

hashes: sha256 of bench/sealed_run.py, bench/t1_driver.py, every bench/t1/*.py, bench/t5_armM.py, bench/env_manifest.py
and the claim-script review (reviews/claims_a5_review.md), read with `git show HEAD:<path>` (so the values are the
committed bytes).  env_manifests: bench/env_manifest.py digests of the Sage binary, the Comet env + binary, TRFP, the
mokapot env, and the pe / pred envs that scripts/e9/stage_pred_env.sh stages for the look.  The two top-level blocks are
rewritten in place (text replacement; the rest of the file, with its comments, is untouched).  --check only compares.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import env_manifest as EM  # noqa: E402

FILES = ["bench/sealed_run.py", "bench/t1_driver.py", "bench/t5_armM.py", "bench/env_manifest.py",
         "bench/t1/claims_a5.py", "reviews/claims_a5_review.md"]
ENVS = {"sage": "tools/bin/sage", "comet_env": "conda_envs/comet", "trfp": "conda_envs/trfp",
        "mokapot_env": "/tmp/vulcan_m1e/envs/mokapot", "pe_env": "conda_envs/pe", "pred_env": "conda_envs/pred"}


def git_bytes(path):
    r = subprocess.run(["git", "-C", ROOT, "show", f"HEAD:{path}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def compute():
    files = list(FILES) + sorted(p for p in subprocess.run(["git", "-C", ROOT, "ls-files", "bench/t1/*.py"],
                                                           capture_output=True, text=True).stdout.split())
    hashes = {}
    for f in dict.fromkeys(files):
        b = git_bytes(f)
        if b is None:
            raise SystemExit(f"fill_hashes: {f} is not committed")
        hashes[f] = hashlib.sha256(b).hexdigest()
    envs = {}
    for k, p in ENVS.items():
        ap = p if os.path.isabs(p) else os.path.join(ROOT, p)
        envs[k] = {"path": ap, "digest": EM.digest(ap)}
    return hashes, envs


def render(hashes, envs):
    h = "hashes:                     # A5.11.5: sha256 of the committed bytes (scripts/a5/fill_hashes.py); checked by sealed_run --look a5\n"
    h += "".join(f"  {k}: {v}\n" for k, v in hashes.items())
    e = "env_manifests:              # bench/env_manifest.py digests; recomputed and compared by sealed_run --look a5\n"
    e += "".join(f"  {k}: {{path: {v['path']}, digest: {v['digest']}}}\n" for k, v in envs.items())
    return h, e


def replace_block(text, key, new):
    m = re.search(rf"^{key}:.*\n((?:[ #].*\n|\n)*)", text, flags=re.M)
    if not m:
        raise SystemExit(f"fill_hashes: no top-level `{key}:` block")
    return text[:m.start()] + new + text[m.end():]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    hashes, envs = compute()
    p = os.path.join(ROOT, "config", "a5.yaml")
    import yaml
    cur = yaml.safe_load(open(p))
    if a.check:
        ok = cur.get("hashes") == hashes and cur.get("env_manifests") == envs
        print("hashes match" if ok else "MISMATCH")
        sys.exit(0 if ok else 1)
    h, e = render(hashes, envs)
    t = open(p).read()
    t = replace_block(t, "hashes", h)
    t = replace_block(t, "env_manifests", e)
    open(p, "w").write(t)
    new = yaml.safe_load(open(p))
    assert new["hashes"] == hashes and new["env_manifests"] == envs
    print(f"wrote {len(hashes)} file hashes and {len(envs)} env manifests")


if __name__ == "__main__":
    main()
