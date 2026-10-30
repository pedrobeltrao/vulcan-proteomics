"""Shared E2 settings and idempotent preparation (stores, L1 indexes, Sage inputs) on local disk.

Work dir: $VULCAN_E2_WORK (default /tmp/vulcan_e2). Everything there is rebuilt by `python scripts/m2/e2_prepare.py`.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
WORK = os.environ.get("VULCAN_E2_WORK", "/tmp/vulcan_e2")
RES = os.path.join(ROOT, "results", "2026-10-27_E2_gpu_search")
STATIC = {"C": 57.021464}
VAR = {"M": 15.9949}
MAX_VAR = 2
ENTRAP_SEED = 2001          # ENTRAP_SEEDS_DEV[0] (vulcan.stats._util)

DATASETS = {
    "D0": {"accession": "PXD001468", "fasta": ["human_sp_canonical.fasta", "contaminants_frankenfield2022.fasta"],
           "family": "human_contam",
           "raw": ["PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw",
                   "PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw"]},
    "D1": {"accession": "PXD001077", "fasta": ["pfu_UP000001013.fasta", "contaminants_frankenfield2022.fasta"],
           "family": "pfu_contam", "raw": ["PXD001077/Velos005137.raw"]},
}


def run_name(raw: str) -> str:
    return os.path.basename(raw)[:-4]


def store_path(raw: str) -> str:
    return os.path.join(WORK, "store", run_name(raw) + ".vspec.arrow")


def mgf_path(raw: str) -> str:
    return os.path.join(WORK, "mgf", run_name(raw) + ".mgf")


def fasta_path(ds: str) -> str:
    return os.path.join(WORK, DATASETS[ds]["family"] + ".fasta")


def index_dir(ds: str, kind: str) -> str:
    """kind: 'prod' (T + dec(T)) or 'E<seed>' (T + E + dec(T) + dec(E))."""
    return os.path.join(WORK, "index", f"{DATASETS[ds]['family']}_{kind}")


def sage_fasta(ds: str, kind: str) -> str:
    return os.path.join(WORK, "sage", f"{DATASETS[ds]['family']}_{kind}.peptides.fasta")


def ensure_store(raw: str):
    out = store_path(raw)
    if os.path.exists(out):
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    subprocess.run([sys.executable, "-m", "vulcan.io.convert", os.path.join(ROOT, "data", "raw", raw),
                    os.path.dirname(out), "--threads", "4"], check=True, cwd=ROOT,
                   env={**os.environ, "PYTHONPATH": ROOT})
    return out


def ensure_fasta(ds: str):
    out = fasta_path(ds)
    if not os.path.exists(out):
        os.makedirs(WORK, exist_ok=True)
        with open(out + ".partial", "w") as fh:
            for f in DATASETS[ds]["fasta"]:
                fh.write(open(os.path.join(ROOT, "data", "fasta", f)).read())
        os.replace(out + ".partial", out)
    return out


def ensure_index(ds: str, kind: str):
    from vulcan.index.build import build
    d = index_dir(ds, kind)
    if os.path.exists(os.path.join(d, "meta.json")):
        return d
    seeds = [] if kind == "prod" else [int(kind[1:])]
    build([ensure_fasta(ds)], d, {"entrap_seeds": seeds})
    return d


def ensure_mgf(raw: str):
    from vulcan.gpu.engine import load_store
    from vulcan.gpu.sage_bridge import write_mgf
    out = mgf_path(raw)
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        write_mgf(load_store(ensure_store(raw)), out)
    return out


def ensure_sage_fasta(ds: str, kind: str):
    from vulcan.gpu.sage_bridge import write_peptide_fasta
    out = sage_fasta(ds, kind)
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        info = write_peptide_fasta(ensure_index(ds, kind), out)
        with open(out + ".json", "w") as fh:
            json.dump(info, fh)
    return out


def hwinfo() -> dict:
    try:
        from vulcan.hw import fingerprint
        return fingerprint()
    except Exception as e:  # pragma: no cover
        out = {"error": repr(e)}
        for k, cmd in {"nvidia_smi": ["nvidia-smi", "-q"], "lscpu": ["lscpu"]}.items():
            try:
                out[k] = subprocess.run(cmd, capture_output=True, text=True).stdout
            except Exception:
                pass
        return out
