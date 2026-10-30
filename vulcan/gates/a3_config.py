"""Loader and guards for config/a3.yaml (Amendment A3: freeze, model hashes, A-syn, V rows, tiers)."""
from __future__ import annotations

import os
import subprocess

from vulcan.gates.a2_config import sha256_file

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PATH = os.path.join(ROOT, "config", "a3.yaml")
TODO = "TODO"


def load(path: str = PATH) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


def check_models(cfg: dict | None = None) -> dict:
    """Verify every scorer file hash recorded in config/a3.yaml; raises on mismatch."""
    cfg = cfg or load()
    out = {}
    for name, s in cfg["scorers"].items():
        for fk, hk in (("model", "sha256"), ("preregistration", "preregistration_sha256")):
            if fk in s:
                got = sha256_file(os.path.join(ROOT, s[fk]))
                if got != s[hk]:
                    raise RuntimeError(f"A3 config: {name} {s[fk]} sha256 {got} != recorded {s[hk]}")
                out[f"{name}.{fk}"] = got
    return out


def check_offsets(cfg: dict | None = None) -> str:
    cfg = cfg or load()
    got = sha256_file(os.path.join(ROOT, cfg["asyn"]["offsets_table"]))
    if got != cfg["asyn"]["offsets_sha256"]:
        raise RuntimeError(f"A3 config: offsets table sha256 {got} != recorded {cfg['asyn']['offsets_sha256']}")
    return got


def _git(*args) -> str:
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=True).stdout.strip()


def tree_hashes(commit: str, paths) -> dict:
    """`git ls-tree` object hash of each path at commit."""
    out = {}
    for p in paths:
        line = _git("ls-tree", commit, p.rstrip("/"))
        out[p] = line.split()[2] if line else None
    return out


# A3.1.3 frozen trees, extended by A5.11.6 (DESIGN §6.0d; stricter): T5 and the comparators import from scripts/m3,
# m3e, entrapment, e7b; scripts/e9 stages the env; scripts/baselines builds the mokapot env.
FROZEN_TREES = ("vulcan", "bench", "scripts/a3", "scripts/e7", "scripts/e7b", "scripts/m3", "scripts/m3e",
                "scripts/entrapment", "scripts/e9", "scripts/baselines")


def freeze_status(cfg: dict | None = None) -> dict:
    """Is the working tree at the recorded freeze commit, clean, with the recorded tree hashes?"""
    cfg = cfg or load()
    fz = cfg["freeze"]
    st = {"freeze_commit": fz["commit"], "recorded": fz["commit"] != TODO}
    if not st["recorded"]:
        st["ok"] = False
        st["reason"] = "freeze commit not recorded yet (config/a3.yaml freeze.commit is TODO)"
        return st
    missing = [t for t in FROZEN_TREES if t not in (fz.get("tree_hashes") or {})]
    if missing:
        st["ok"] = False
        st["reason"] = f"config/a3.yaml freeze.tree_hashes lacks frozen trees {missing} (A3.1.3 + A5.11.6)"
        return st
    head = _git("rev-parse", "HEAD")
    st["head"] = head
    # clean = no modified tracked file anywhere AND no untracked file under the code / config paths (an untracked
    # module could be imported); untracked result files elsewhere are allowed
    st["clean"] = (_git("status", "--porcelain", "--untracked-files=no") == "" and
                   _git("status", "--porcelain", "--untracked-files=all", "--", "vulcan", "bench", "scripts",
                        "config", "tests") == "")
    full = _git("rev-parse", fz["commit"])
    st["at_freeze"] = head == full
    # config/a3.yaml cannot name the commit that contains it, so it is committed AFTER the freeze commit.  "Running
    # from the freeze commit" therefore means: HEAD is the freeze commit or a descendant of it, AND every frozen code
    # path (tree_hashes) is byte-identical at HEAD to the freeze commit (git tree object hashes), AND the recorded
    # hashes are the freeze commit's.  Any change to a frozen path after the freeze makes the HEAD hashes differ.
    st["freeze_is_ancestor"] = subprocess.run(["git", "-C", ROOT, "merge-base", "--is-ancestor", full, head],
                                              capture_output=True).returncode == 0
    at_fz = tree_hashes(full, fz["tree_hashes"].keys())
    at_head = tree_hashes(head, fz["tree_hashes"].keys())
    st["tree_hashes_match"] = all(at_fz[k] == v for k, v in fz["tree_hashes"].items())
    st["head_tree_hashes_match"] = all(at_head[k] == v for k, v in fz["tree_hashes"].items())
    st["code_at_freeze"] = bool(st["freeze_is_ancestor"] and st["tree_hashes_match"] and st["head_tree_hashes_match"])
    st["ok"] = bool(st["clean"] and st["code_at_freeze"])
    if not st["ok"]:
        st["reason"] = ("not at (a descendant of) the freeze commit with identical frozen code paths, tree not clean, "
                        "or tree hashes differ")
    return st


def require_freeze_commit(cfg: dict | None = None) -> dict:
    """A3.1.3: step 4b, the A3 part of step 5 and the T1 v2_sep numerator run from the freeze commit on a clean tree."""
    st = freeze_status(cfg)
    if not st["ok"]:
        raise SystemExit(f"A3: refusing to run: {st.get('reason')} ({st})")
    return st
