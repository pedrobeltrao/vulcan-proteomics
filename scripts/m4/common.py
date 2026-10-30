"""Shared M4 settings and idempotent preparation (FASTAs, stores, L1 indexes) on node-local disk.

Work dir: $VULCAN_M4_WORK (default /tmp/vulcan_m4). Rebuilt by `python scripts/m4/prepare.py`.
Dev data only: nothing under data/sealed/ or data/truth/sealed/ is read here (DESIGN §6.6, M4 rule).
"""
from __future__ import annotations

import csv
import glob
import json
import os
import re
import subprocess
import sys
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
WORK = os.environ.get("VULCAN_M4_WORK", "/tmp/vulcan_m4")
E2WORK = os.environ.get("VULCAN_E2_WORK", "/tmp/vulcan_e2")
RES = os.path.join(ROOT, "results", "2026-10-27_M4_ptm")

STATIC = {"C": 57.021464}
PHOSPHO = 79.966331
VAR_CLOSED = {"M": 15.9949}
VAR_PHOSPHO = {"S": PHOSPHO, "T": PHOSPHO, "Y": PHOSPHO, "M": 15.994915, "^": 42.010565}
MAX_VAR_PHOSPHO = 3

for p in ("/tmp/vulcan_m4/store",):
    pass

D0_STORES = [os.path.join(E2WORK, "store", f"{r}.vspec.arrow") for r in
             ("b1906_293T_proteinID_01A_QE3_122212", "b1937_293T_proteinID_01B_QE3_122212")]
D0_INDEX = os.path.join(E2WORK, "index", "human_contam_prod")       # E2 closed index (T + pseudo-reverse)
D0_MGF = [os.path.join(E2WORK, "mgf", os.path.basename(s).replace(".vspec.arrow", ".mgf")) for s in D0_STORES]

D2_RAW = {
    "PXD000138": sorted(glob.glob(os.path.join(ROOT, "data/raw/PXD000138/*.raw"))),
    "PXD009449": sorted(glob.glob(os.path.join(ROOT, "data/raw/PXD009449/*Ymod_Phospho*3xHCD*.raw"))),
    "PXD007058": sorted(glob.glob(os.path.join(ROOT, "data/raw/PXD007058/*.raw"))),
}
POOL_WELL = {"58": "E10", "60": "E12", "66": "F6", "72": "F12"}     # DATA.md D2 dev pools
SEALED_T2_WELLS = {"C1", "G4", "B12", "F1"}                       # PXD000138 sealed pools 25, 76, 24, 61


def run_name(raw: str) -> str:
    return os.path.basename(raw)[:-4]


def store_path(raw: str) -> str:
    return os.path.join(WORK, "store", run_name(raw) + ".vspec.arrow")


def d2_stores() -> dict:
    return {acc: [store_path(r) for r in raws] for acc, raws in D2_RAW.items()}


def ensure_store(raw: str) -> str:
    out = store_path(raw)
    if os.path.exists(out):
        return out
    os.makedirs(os.path.join(WORK, "rawstage"), exist_ok=True)
    loc = os.path.join(WORK, "rawstage", os.path.basename(raw))
    subprocess.run(["cp", raw, loc], check=True)                  # never decode on Lustre (LESSONS)
    subprocess.run([os.path.join(ROOT, "scripts/m2b/vpy.sh"), "-m", "vulcan.io.convert", loc,
                    os.path.dirname(out), "--threads", "4"], check=True, cwd=ROOT)
    os.remove(loc)
    return out


# ------------------------------------------------------------------------------------------------ library FASTAs
def ymod_sequences() -> list[str]:
    r = csv.DictReader(open(os.path.join(ROOT, "data/truth/PXD009449/Ymod_Phospho_evidence.txt")), delimiter="\t")
    return sorted({x["Sequence"] for x in r})


def ferries_tables() -> dict:
    """pool -> {'peptides': rows of txt/peptides.txt, 'msms': rows of txt/msms.txt} (submitter MaxQuant)."""
    out = {}
    for z in sorted(glob.glob(os.path.join(ROOT, "data/truth/PXD007058/*.zip"))):
        pool = re.search(r"pool(\d)", os.path.basename(z)).group(1)
        zf = zipfile.ZipFile(z)
        t = {}
        for k in ("peptides", "msms"):
            lines = zf.read(f"txt/{k}.txt").decode().splitlines()
            t[k] = list(csv.DictReader(lines, delimiter="\t"))
        out[pool] = t
    return out


def ferries_sequences() -> list[str]:
    return sorted({r["Sequence"] for t in ferries_tables().values() for r in t["peptides"]})


def d2_fasta() -> str:
    """human SwissProt canonical + contaminants + D2 dev libraries (PXD000138 dev pool entries, PXD009449 Ymod
    peptides, PXD007058 submitter peptides; one entry per library peptide for the latter two)."""
    out = os.path.join(WORK, "d2_phospho.fasta")
    if os.path.exists(out):
        return out
    os.makedirs(WORK, exist_ok=True)
    with open(out + ".partial", "w") as fh:
        for f in ("human_sp_canonical.fasta", "contaminants_frankenfield2022.fasta"):
            fh.write(open(os.path.join(ROOT, "data/fasta", f)).read().rstrip("\n") + "\n")
        for f in sorted(glob.glob(os.path.join(ROOT, "data/truth/PXD000138/pool_*.fasta"))):
            fh.write(open(f).read().rstrip("\n") + "\n")
        for i, s in enumerate(ymod_sequences()):
            fh.write(f">lib_PXD009449_Ymod_{i}\n{s}\n")
        for i, s in enumerate(ferries_sequences()):
            fh.write(f">lib_PXD007058_{i}\n{s}\n")
    os.replace(out + ".partial", out)
    return out


def index_dir(name: str) -> str:
    return os.path.join(WORK, "index", name)


def ensure_index(name: str, fasta: str, var: dict, max_var: int, entrap_seeds=()) -> str:
    from vulcan.index.build import build
    d = index_dir(name)
    if os.path.exists(os.path.join(d, "meta.json")):
        return d
    build([fasta], d, {"variable_mods": var, "max_variable_mods": max_var, "static_mods": STATIC,
                       "entrap_seeds": list(entrap_seeds)})
    return d


def hwinfo() -> dict:
    from vulcan.hw import fingerprint
    return fingerprint()


def dump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
