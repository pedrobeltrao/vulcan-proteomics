"""M5 LiP (D7) + CF-MS (D10cf, E10) shared settings.  Work dir $VULCAN_M5LC_WORK (default /tmp/vulcan_m5lc);
results in results/2026-10-27_M5_lip_cf.  Dev data only."""
from __future__ import annotations

import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "baselines"))
from timed import parse_time_v  # noqa: E402,F401

WORK = os.environ.get("VULCAN_M5LC_WORK", "/tmp/vulcan_m5lc")
RES = os.path.join(ROOT, "results", "2026-10-27_M5_lip_cf")
STATIC = {"C": 57.021464}
VAR = {"M": 15.9949}
SEED = 11

# ---------------------------------------------------------------- D7 (PXD031616), LiP/TrP design (DATA.md D7)
LIP_FASTA = os.path.join(WORK, "fasta", "yeast_contam_snca.fasta")
LIP_INDEX = os.path.join(WORK, "index", "yeast_semi")
LIP_INDEX_SPECIFIC = os.path.join(WORK, "index", "yeast_specific")
# run -> (digest, condition, replicate); mapping from the submitter MSstats reports (DATA.md D7)
D7_DESIGN = {
    "LM2481_DDA": ("LiP", "M1", 1), "LM2492_DDA": ("LiP", "M1", 2), "LM2498_DDA": ("LiP", "M1", 3),
    "LM2480_DDA": ("LiP", "F1", 1), "LM2494_DDA": ("LiP", "F1", 2), "LM2497_DDA": ("LiP", "F1", 3),
    "LM2486_DDA": ("TrP", "M1", 1), "LM2490_DDA": ("TrP", "M1", 2), "LM2500_DDA": ("TrP", "M1", 3),
    "LM2487_DDA": ("TrP", "F1", 1), "LM2489_DDA": ("TrP", "F1", 2), "LM2502_DDA": ("TrP", "F1", 3),
}
SNCA = "P37840"


def d7_raws():
    return [os.path.join(ROOT, "data", "raw", "PXD031616", r + ".raw") for r in sorted(D7_DESIGN)]


# ---------------------------------------------------------------- D10cf (PXD015406 control SEC series, fractions 10-38)
CF_FASTA = os.path.join(WORK, "fasta", "human_contam.fasta")
CF_INDEX = os.path.join(WORK, "index", "human_entrap2001")
CF_ENTRAP_SEED = 2001


def cf_raws():
    return sorted(glob.glob(os.path.join(ROOT, "data", "raw", "PXD015406", "Control_HEK_SEC_020118_*.raw")),
                  key=cf_fraction)


def cf_fraction(path_or_run) -> int:
    b = os.path.basename(path_or_run)
    return int(b.split("_")[4].rstrip("bc.raw").rstrip("bc"))


def store(run):
    return os.path.join(WORK, "store", run + ".vspec.arrow")
