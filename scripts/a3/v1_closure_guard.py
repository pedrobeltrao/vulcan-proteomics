#!/usr/bin/env python3
"""A3.0 Q1 / A3.7: v1 isolation guard.  Computes the transitive import closure of the v1 step-5 path
(scripts/a2/gates_a2.py) and of bench/t1_driver.py (the whole driver: stricter than "its v1 path", since any module
it imports can change v1's T1 numbers), at the tag prereg-v0.4.2 AND at the target commit (union), and diffs every
file in it against prereg-v0.4.2.

Rule (A3.0 Q1): every changed file must be listed in scripts/a3/v1_closure_allow.yaml with the sha256 of its EXACT diff
(`git diff prereg-v0.4.2 <commit> -- <path>`, so any further edit invalidates the entry) and one of
    kind: deviation        a walk-agent dated deviation note (field `note`: where it is recorded)
    kind: bit_identity     committed bit-identity evidence on the v1 path (field `evidence`: committed path(s))
Anything else -> FAIL (exit 1).  An entry whose evidence path is not committed -> FAIL.  With --require-clean, every
closure file must also be unmodified in the working tree (step 5 runs on a clean tree).

Closure: AST of every reachable file (module-level AND function-level `import` / `from ... import`, relative imports
resolved inside packages), resolved against the repository roots the entry points put on sys.path (ROOT, scripts/a2,
scripts/m3, scripts/m3e, scripts/m2, the importing file's directory); third-party modules are outside the repository
and ignored.  Package __init__.py files on the way are included.  Non-Python inputs of the v1 path are added
explicitly (EXTRA): the CUDA kernel source and its build script (vulcan/gpu/_lib.py rebuilds the .so from it), the A2
config and tables.  Files are read with `git show <commit>:<path>`, so the closure is a property of the commit.

    python scripts/a3/v1_closure_guard.py [--commit HEAD] [--require-clean] [--out report.json]
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TAG = "prereg-v0.4.2"
ENTRY = ["scripts/a2/gates_a2.py", "bench/t1_driver.py"]
ROOTS = ["", "scripts/a2", "scripts/m3", "scripts/m3e", "scripts/m2",
         "scripts/e7", "scripts/e7b"]          # A3 freeze: bench/t1/a3.py imports scripts/e7/strata, e7b/sage_rank_strata
EXTRA = ["vulcan/gpu/kernels.cu", "scripts/m2/build_kernels.sh", "config/a2.yaml", "config/a2_panels.tsv",
         "config/a2_excluded_accessions.tsv"]
ALLOW = os.path.join(ROOT, "scripts", "a3", "v1_closure_allow.yaml")


def git(*args, check=True) -> str:
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=check).stdout


def tree_files(commit: str) -> set:
    return set(git("ls-tree", "-r", "--name-only", commit).split())


def show(commit: str, path: str) -> str | None:
    r = subprocess.run(["git", "-C", ROOT, "show", f"{commit}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def _candidates(mod: str, base_dirs):
    parts = mod.split(".")
    for b in base_dirs:
        stem = os.path.join(b, *parts) if b else os.path.join(*parts)
        yield stem + ".py"
        yield os.path.join(stem, "__init__.py")


def _pkg_inits(path: str, files: set):
    d = os.path.dirname(path)
    out = []
    while d:
        f = os.path.join(d, "__init__.py")
        if f in files:
            out.append(f)
        d = os.path.dirname(d)
    return out


def imports_of(src: str, path: str):
    """Module names imported anywhere in the file; relative imports returned as absolute dotted names.  Also modules run
    as a subprocess with `python -m <module>` (a literal "-m" followed by a literal module name in a list/tuple, e.g.
    bench/t1/ours.py -> vulcan.cohort, bench/t1/ours_chimeric.py -> vulcan.ms1), plus their `<module>.__main__`: they
    run on the same path although no import statement names them (A3 freeze: stricter closure)."""
    tree = ast.parse(src)
    pkg = os.path.dirname(path).replace(os.sep, ".")
    names = []
    for n in ast.walk(tree):
        if isinstance(n, (ast.List, ast.Tuple)):
            el = n.elts
            for x, y in zip(el, el[1:]):
                if isinstance(x, ast.Constant) and x.value == "-m" and isinstance(y, ast.Constant) and \
                        isinstance(y.value, str) and all(p.isidentifier() for p in y.value.split(".")):
                    names.extend([y.value, f"{y.value}.__main__"])
        if isinstance(n, ast.Import):
            names.extend(a.name for a in n.names)
        elif isinstance(n, ast.ImportFrom):
            if n.level:
                base = pkg.split(".")
                base = base[:len(base) - (n.level - 1)] if n.level > 1 else base
                mod = ".".join([x for x in base if x] + ([n.module] if n.module else []))
            else:
                mod = n.module or ""
            names.append(mod)
            names.extend(f"{mod}.{a.name}" for a in n.names)     # `from pkg import module`
    return names


def closure(commit: str, entries=ENTRY) -> set:
    files = tree_files(commit)
    todo = [e for e in entries if e in files]
    seen = set()
    while todo:
        f = todo.pop()
        if f in seen:
            continue
        seen.add(f)
        for init in _pkg_inits(f, files):
            if init not in seen:
                todo.append(init)
        if not f.endswith(".py"):
            continue
        src = show(commit, f)
        if src is None:
            continue
        bases = ROOTS + [os.path.dirname(f)]
        for mod in imports_of(src, f):
            if not mod:
                continue
            for c in _candidates(mod, bases):
                c = os.path.normpath(c)
                if c in files:
                    if c not in seen:
                        todo.append(c)
                    break
    return seen | {e for e in EXTRA if e in files}


def diff_sha(commit: str, path: str) -> tuple[str, int]:
    d = git("diff", TAG, commit, "--", path)
    return (hashlib.sha256(d.encode()).hexdigest() if d else ""), d.count("\n@@")


def load_allow():
    import yaml
    if not os.path.exists(ALLOW):
        return []
    return yaml.safe_load(open(ALLOW)) or []


def committed(path: str, commit: str) -> bool:
    return show(commit, path) is not None


def check(commit: str = "HEAD", require_clean: bool = False) -> dict:
    full = git("rev-parse", commit).strip()
    cl = sorted(closure(TAG) | closure(full))
    allow = {(a["path"], a["diff_sha256"]): a for a in load_allow()}
    per_entry = {e: closure(TAG, [e]) | closure(full, [e]) for e in ENTRY}
    rep = {"tag": TAG, "commit": full, "n_files": len(cl), "files": cl, "changed": [], "unallowed": [],
           "dirty": [], "by_entry": {e: {"n_files": len(v)} for e, v in per_entry.items()}}
    for f in cl:
        sha, nh = diff_sha(full, f)
        if not sha:
            continue
        a = allow.get((f, sha))
        ok = a is not None and a.get("kind") in ("deviation", "bit_identity")
        if ok and a["kind"] == "bit_identity":
            ev = a.get("evidence") or []
            ok = bool(ev) and all(committed(e, full) for e in ([ev] if isinstance(ev, str) else ev))
        if ok and a["kind"] == "deviation":
            ok = bool(a.get("note"))
        rec = {"path": f, "diff_sha256": sha, "hunks": nh + (1 if sha else 0), "allowed": ok,
               "in_closure_of": [e for e, v in per_entry.items() if f in v],
               "entry": a, "commits": git("log", "--format=%h %s", f"{TAG}..{full}", "--", f).strip().splitlines()}
        rep["changed"].append(rec)
        if not ok:
            rep["unallowed"].append(f)
    if require_clean:
        st = git("status", "--porcelain", "--", *cl)
        rep["dirty"] = [x[3:] for x in st.splitlines()]
    rep["pass"] = not rep["unallowed"] and not rep["dirty"]
    for e in ENTRY:
        rep["by_entry"][e]["unallowed"] = [c["path"] for c in rep["changed"] if not c["allowed"] and e in c["in_closure_of"]]
        rep["by_entry"][e]["pass"] = not rep["by_entry"][e]["unallowed"] and not [d for d in rep["dirty"]
                                                                                   if d in per_entry[e]]
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", default="HEAD")
    ap.add_argument("--require-clean", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    rep = check(a.commit, a.require_clean)
    if a.out:
        json.dump(rep, open(a.out, "w"), indent=1)
    for c in rep["changed"]:
        print(("OK   " if c["allowed"] else "FAIL ") + c["path"], c["diff_sha256"][:12], "|", "; ".join(c["commits"])[:150])
    for d in rep["dirty"]:
        print("DIRTY", d)
    for e, v in rep["by_entry"].items():
        print(f"  {e}: {v['n_files']} files, unallowed {len(v['unallowed'])} -> {'PASS' if v['pass'] else 'FAIL'}")
    print(f"v1 closure: {rep['n_files']} files, {len(rep['changed'])} changed since {TAG}, "
          f"{len(rep['unallowed'])} without allowed evidence -> {'PASS' if rep['pass'] else 'FAIL'}")
    sys.exit(0 if rep["pass"] else 1)


if __name__ == "__main__":
    main()
