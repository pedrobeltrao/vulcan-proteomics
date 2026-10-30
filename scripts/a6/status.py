#!/usr/bin/env python3
"""A6.2 status file: the ONLY channel from A6 to the A5 look.

`results/2026-10-28_A6_hn1/a6_status.json` holds a closed-enum `status`, the step-5 commit sha (or null) and the
sha256 of the gate output (or null).  `read_status` is what the A5 claim script calls: if the file is absent,
unparseable, carries a value outside the enum, is not committed, or its gate output does not hash to the recorded
sha256, the state is "stopped: node time limit".  No one sets the state by hand: `write_status` is called only by
scripts/a6/gates_hn1.py (step 5) or by the label-free stop paths (not tagged / label-free stop).
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
STATUS_PATH = os.path.join(ROOT, "results", "2026-10-28_A6_hn1", "a6_status.json")
ENUM = ("PASS", "FAILED", "not evaluable", "invalid pool", "not completed", "label-free stop", "not tagged",
        "stopped: node time limit")
DEFAULT = "stopped: node time limit"
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def write_status(status: str, step5_commit: str | None = None, gate_output: str | None = None, note: str = "",
                 path: str = STATUS_PATH) -> dict:
    """Write the status file.  Refuses any status outside the closed enum and malformed sha values.
    gate_output: path of the step-5 gate JSON (its sha256 is recorded) or None."""
    if status not in ENUM:
        raise ValueError(f"a6 status {status!r} is not in the closed enum {ENUM}")
    if step5_commit is not None and not SHA_RE.match(step5_commit):
        raise ValueError(f"step5_commit {step5_commit!r} is not a 40-hex git sha")
    if status in ("PASS", "FAILED") and (step5_commit is None or gate_output is None):
        raise ValueError("PASS/FAILED need the step-5 commit sha and the gate output")
    rec = {"amendment": "A6", "status": status, "step5_commit": step5_commit,
           "gate_output": os.path.relpath(gate_output, ROOT) if gate_output else None,
           "gate_output_sha256": sha256_file(gate_output) if gate_output else None, "note": note}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(rec, fh, indent=1)
        fh.write("\n")
    return rec


def _committed_unmodified(path: str) -> bool:
    rel = os.path.relpath(path, ROOT)
    ls = subprocess.run(["git", "-C", ROOT, "ls-files", "--error-unmatch", rel], capture_output=True)
    st = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", rel], capture_output=True, text=True)
    return ls.returncode == 0 and not st.stdout.strip()


def read_status(path: str = STATUS_PATH, require_committed: bool = True) -> dict:
    """-> {"status": <enum value>, "verified": bool, "reason": str, "record": dict|None}.
    Anything absent or unverifiable -> "stopped: node time limit" (A6.2)."""
    def stop(reason, rec=None):
        return {"status": DEFAULT, "verified": False, "reason": reason, "record": rec}
    if not os.path.exists(path):
        return stop("absent")
    try:
        rec = json.load(open(path))
    except (OSError, ValueError):
        return stop("unparseable")
    if not isinstance(rec, dict) or rec.get("status") not in ENUM:
        return stop("status outside the closed enum", rec)
    c = rec.get("step5_commit")
    if c is not None and not (isinstance(c, str) and SHA_RE.match(c)):
        return stop("malformed step5_commit", rec)
    if require_committed and not _committed_unmodified(path):
        return stop("status file not committed or modified", rec)
    if rec["status"] in ("PASS", "FAILED") and (c is None or not rec.get("gate_output")):
        return stop("PASS/FAILED without step-5 commit or gate output", rec)
    go = rec.get("gate_output")
    if go is not None:
        gp = os.path.join(ROOT, go)
        if not os.path.exists(gp) or not HEX64.match(str(rec.get("gate_output_sha256") or "")) \
                or sha256_file(gp) != rec["gate_output_sha256"]:
            return stop("gate output missing or sha256 mismatch", rec)
    return {"status": rec["status"], "verified": True, "reason": "ok", "record": rec}


if __name__ == "__main__":  # pragma: no cover
    print(json.dumps(read_status(), indent=1))
