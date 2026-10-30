"""M6 (quant) shared settings.  Work dir $VULCAN_M6_WORK (default /tmp/vulcan_m6); results in
results/2026-10-27_M6_quant.  Dev data only: D9 (ProteoBench DDA-Astral HYE + human-only DDA-Astral runs) and D4."""
from __future__ import annotations

import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
WORK = os.environ.get("VULCAN_M6_WORK", "/tmp/vulcan_m6")
RES = os.path.join(ROOT, "results", "2026-10-27_M6_quant")
STATIC = {"C": 57.021464}
VAR = {"M": 15.9949}
TMT = 229.162932
SEED = 11

D9_FASTA = os.path.join(WORK, "fasta", "ProteoBenchFASTA_MixedSpecies_HYE.fasta")
D9_INDEX = os.path.join(WORK, "index", "hye_prod")
D4_FASTA = "/tmp/vulcan_e2/human_contam.fasta"
D4_INDEX = os.path.join(WORK, "index", "human_contam_tmtK_prod")
# ProteoBench HYE convention: A/B log2 fold change per species (A = 65% H / 30% Y / 5% E, B = 65 / 15 / 20)
EXPECTED_LOG2FC = {"HUMAN": 0.0, "YEAST": 1.0, "ECOLI": -2.0}


def d9_raws():
    a = sorted(glob.glob(os.path.join(ROOT, "data", "raw", "ProteoBench-DDA-Astral", "*.raw")))
    h = sorted(glob.glob(os.path.join(ROOT, "data", "raw", "ProteoBench-DDA-Astral-entrapment", "*.raw")))
    return a, h


def run_name(raw):
    return os.path.basename(raw).rsplit(".", 1)[0]


def condition(run):
    """'A' / 'B' for the HYE runs, 'H' for the human-only MBR control runs."""
    if "_Condition_A_" in run:
        return "A"
    if "_Condition_B_" in run:
        return "B"
    if "_Human_" in run:
        return "H"
    raise ValueError(run)


def store(run):
    return os.path.join(WORK, "store", run + ".vspec.arrow")


def d4_raws():
    return sorted(glob.glob(os.path.join(ROOT, "data", "raw", "PXD016766", "*.raw")))
