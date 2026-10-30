"""M3 (scoring / E4) shared settings and idempotent preparation on local disk.

Work dir: $VULCAN_M3_WORK (default /tmp/vulcan_m3); dev stores/indexes of E2 are reused from $VULCAN_E2_WORK
(/tmp/vulcan_e2, rebuilt by scripts/m2/e2_prepare.py).  Everything here is rebuilt by scripts/m3/prepare.py.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m2"))
import e2_common as E2  # noqa: E402

WORK = os.environ.get("VULCAN_M3_WORK", "/tmp/vulcan_m3")
RES = os.path.join(ROOT, "results", "2026-10-27_M3_scoring")
M1E_DB = os.environ.get("VULCAN_M1E_DB", "/tmp/vulcan_m1e/entrapment/db")
STATIC, VAR = E2.STATIC, E2.VAR
REPORT = 5                      # candidates per spectrum kept for rescoring (rank by hyperscore)
ENTRAP_SEEDS = (2001, 2002, 2003)   # ENTRAP_SEEDS_DEV (§6.1)
SEARCH = {"report_psms": REPORT}

DEV = {  # dataset -> (family, fasta files, raw files, M1e C-arm name)
    "D0": {"family": "human_contam", "raw": E2.DATASETS["D0"]["raw"], "c_db": "human_contam"},
    "D1": {"family": "pfu_contam", "raw": E2.DATASETS["D1"]["raw"], "c_db": "pfu_contam"},
}
TRAIN_FASTA = {"human": ["human_sp_canonical.fasta", "contaminants_frankenfield2022.fasta"],
               "mouse": ["mouse_sp_canonical.fasta", "contaminants_frankenfield2022.fasta"],
               "ecoli": ["ecoli_K12_UP000000625.fasta", "contaminants_frankenfield2022.fasta"]}


def train_rows():
    with open(os.path.join(ROOT, "data", "manifests", "train.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def train_raw(r):
    return os.path.join(ROOT, "data", "raw", "train", r["accession"], r["filename"])


def train_store(r):
    return os.path.join(WORK, "store_train", r["filename"][:-4] + ".vspec.arrow")


def index_dir(family: str, kind: str) -> str:
    """kind: prod | E<seed> | C.  D0/D1 prod and E2001 come from E2's work dir."""
    if family in ("human_contam", "pfu_contam") and kind in ("prod", "E2001"):
        return os.path.join(E2.WORK, "index", f"{family}_{kind}")
    return os.path.join(WORK, "index", f"{family}_{kind}")


def fasta_path(family: str) -> str:
    if family in ("human_contam", "pfu_contam"):
        return E2.fasta_path({"human_contam": "D0", "pfu_contam": "D1"}[family])
    sp = family.split("_")[0]
    out = os.path.join(WORK, family + ".fasta")
    if not os.path.exists(out):
        os.makedirs(WORK, exist_ok=True)
        with open(out + ".partial", "w") as fh:
            for f in TRAIN_FASTA[sp]:
                fh.write(open(os.path.join(ROOT, "data", "fasta", f)).read())
        os.replace(out + ".partial", out)
    return out


def foreign_list(family: str) -> str:
    """Arm-C foreign peptides (kind E) from the M1e C database (neighbour-filtered, DESIGN §6.3)."""
    out = os.path.join(WORK, "foreign", f"{family}.C.foreign.txt")
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        src = os.path.join(M1E_DB, f"{family}.C.peptides.tsv")
        with open(src) as fh, open(out + ".partial", "w") as oh:
            next(fh)
            for line in fh:
                seq, kind = line.split("\t", 2)[:2]
                if kind == "E":
                    oh.write(seq + "\n")
        os.replace(out + ".partial", out)
    return out


def c_arm_r(family: str) -> float:
    with open(os.path.join(M1E_DB, f"{family}.C.json")) as fh:
        return float(json.load(fh)["r"])


def ensure_index(family: str, kind: str):
    d = index_dir(family, kind)
    if os.path.exists(os.path.join(d, "meta.json")):
        return d
    if kind == "C":
        from vulcan.index.foreign import build_foreign
        build_foreign([fasta_path(family)], d, foreign_list(family))
    else:
        from vulcan.index.build import build
        build([fasta_path(family)], d, {"entrap_seeds": [] if kind == "prod" else [int(kind[1:])]})
    return d


def ensure_train_store(r):
    out = train_store(r)
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        subprocess.run([sys.executable, "-m", "vulcan.io.convert", train_raw(r), os.path.dirname(out),
                        "--threads", "4"], check=True, cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT})
    return out


def psm_dir(tag: str) -> str:
    return os.path.join(WORK, "psms", tag)
