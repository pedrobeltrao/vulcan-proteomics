#!/usr/bin/env python3
"""Mechanical re-freeze when ONLY the T2 driver files changed inside a frozen tree after the A3/A5 freeze code commit.

The A3.1.3 freeze hashes the whole `bench/` tree, and the T2 driver (bench/t2_driver.py, bench/t2/) lives there and is
still being committed by its own agent.  Those files are imported by nothing on the T1, T5, DN2 or v2_sep paths (the
v1 closure guard and the A5 hash list do not contain them).  This script:
  1. refuses unless HEAD is clean and a descendant of the recorded freeze commit, and every difference between the
     recorded freeze commit and HEAD inside the frozen trees (vulcan.gates.a3_config.FROZEN_TREES) is under
     bench/t2/ or is bench/t2_driver.py;
  2. refuses unless the v1 closure guard passes at HEAD and config/a5.yaml `hashes` still match the committed bytes;
  3. rewrites config/a3.yaml freeze.commit to HEAD and freeze.tree_hashes.bench to HEAD's bench tree, appending a
     dated comment that names the old commit and the T2-only diff.
The new config is committed by the caller (`git commit -- config/a3.yaml`); every other freeze record stays valid
because no frozen file outside bench/t2* changed.  Any other difference -> a full re-freeze (G0, D1 dry run, guard).

    python scripts/a5/refreeze_t2_only.py [--dry-run]
"""
import argparse
import datetime
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def git(*a):
    return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True, check=True).stdout.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    from vulcan.gates import a3_config as A
    cfg = A.load()
    old = cfg["freeze"]["commit"]
    head = git("rev-parse", "HEAD")
    if git("status", "--porcelain", "--untracked-files=no"):
        sys.exit("refreeze: tracked files modified; commit first")
    if subprocess.run(["git", "-C", ROOT, "merge-base", "--is-ancestor", old, head]).returncode:
        sys.exit("refreeze: the recorded freeze commit is not an ancestor of HEAD")
    changed = git("diff", "--name-only", old, head, "--", *A.FROZEN_TREES).split()
    bad = [f for f in changed if not (f.startswith("bench/t2/") or f == "bench/t2_driver.py")]
    if bad:
        sys.exit("refreeze: frozen files other than the T2 driver changed -> full re-freeze needed:\n  " + "\n  ".join(bad))
    g = subprocess.run([sys.executable, os.path.join(ROOT, "scripts/a3/v1_closure_guard.py"), "--commit", head],
                       capture_output=True, text=True)
    if g.returncode:
        sys.exit("refreeze: v1 closure guard fails at HEAD:\n" + g.stdout[-2000:])
    c = subprocess.run([sys.executable, os.path.join(ROOT, "scripts/a5/fill_hashes.py"), "--check"], capture_output=True,
                       text=True)
    if c.returncode:
        sys.exit("refreeze: config/a5.yaml hashes do not match the committed bytes")
    new_bench = git("rev-parse", f"{head}:bench")
    print(f"T2-only changes since {old[:8]}: {len(changed)} files; new bench tree {new_bench}")
    if a.dry_run:
        return
    p = os.path.join(ROOT, "config", "a3.yaml")
    s = open(p).read()
    s, n1 = re.subn(rf"^  commit: {old}\b.*$", f"  commit: {head}   # re-frozen {datetime.date.today()} "
                    f"(scripts/a5/refreeze_t2_only.py): only bench/t2* changed since {old[:8]}", s, count=1, flags=re.M)
    s, n2 = re.subn(r"^    bench: [0-9a-f]{40}$", f"    bench: {new_bench}", s, count=1, flags=re.M)
    if n1 != 1 or n2 != 1:
        sys.exit("refreeze: could not locate freeze.commit / tree_hashes.bench in config/a3.yaml")
    open(p, "w").write(s)
    st = A.freeze_status(A.load())
    print("freeze_status after rewrite (before commit):", {k: st.get(k) for k in ("tree_hashes_match", "head_tree_hashes_match")})


if __name__ == "__main__":
    main()
