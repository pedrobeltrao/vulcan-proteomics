"""Loader for config/a2.yaml (Amendment A2 thresholds, seeds and fixed implementation choices)."""
from __future__ import annotations

import hashlib
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PATH = os.path.join(ROOT, "config", "a2.yaml")


def load(path: str = PATH) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def check_panel_files(cfg: dict | None = None) -> dict:
    """Verify panel / cRAP / Frankenfield FASTA hashes against the config; raises on any mismatch."""
    cfg = cfg or load()
    p = cfg["panels"]
    out = {}
    for k in ("fasta", "crap_fasta", "frankenfield_fasta"):
        key = {"fasta": "fasta_sha256", "crap_fasta": "crap_sha256", "frankenfield_fasta": "frankenfield_sha256"}[k]
        got = sha256_file(os.path.join(ROOT, p[k]))
        if got != p[key]:
            raise RuntimeError(f"A2 config: {p[k]} sha256 {got} != recorded {p[key]}")
        out[k] = got
    return out


def check_seqclf(cfg: dict | None = None) -> str:
    cfg = cfg or load()
    s = cfg["gates"]["seq_classifier"]
    got = sha256_file(s["path"])
    if got != s["sha256"]:
        raise RuntimeError(f"A2 config: seq classifier sha256 {got} != recorded {s['sha256']}")
    return got
