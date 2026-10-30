"""3c (IT-CID instrument-class static-first scorer) shared settings.  Work dir $VULCAN_3C_WORK (default /tmp/vulcan_3c);
results in results/2026-10-27_3c_itcid.  Dev data only (D4 = PXD016766, D10cf = PXD015406); never data/sealed.

Indexes (L1, vulcan.index.build; targets + paired pseudo-reverse decoys):
  family human     human Swiss-Prot + contaminants (E2/M3 FASTA), CAM static, Met-ox variable  (reuses E2/M3 indexes)
  family mouse     mouse + contaminants (M3 FASTA)                                              (reuses M3 prod index)
  family human_tmt / mouse_tmt   same FASTAs, CAM + TMT-K static; TMT peptide-N-term static added at load time
                   (mass32 + TMT, mask bit 63 + kernel N-term delta = Sage's static `^`, as M6 D4)
  kinds: prod | E2001..E2003 (arm P, ENTRAP_SEEDS_DEV) | C (arm C, M1e neighbour-filtered Arabidopsis peptides)
Search (IT-CID): fragment tolerance +-0.5 Da (fragment_unit "da"; Sage's IT setting, as M1e's Sage D4 config),
precursor +-10 ppm, isotope errors 0..2, top-5 candidates per spectrum (report_psms 5), K1 bins 200 ppm.
"""
from __future__ import annotations

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

WORK = os.environ.get("VULCAN_3C_WORK", "/tmp/vulcan_3c")
RES = os.path.join(ROOT, "results", "2026-10-27_3c_itcid")
TMT = 229.162932
CAM = {"C": 57.021464}
OX = {"M": 15.9949}
ENTRAP_SEEDS = (2001, 2002, 2003)
KINDS = ("prod", "E2001", "E2002", "E2003", "C")
SEARCH = {"report_psms": 5, "fragment_unit": "da", "fragment_ppm": (-0.5, 0.5), "precursor_ppm": (-10.0, 10.0),
          "bin_ppm": 200.0}
FASTA = {"human": "/tmp/vulcan_e2/human_contam.fasta", "mouse": "/tmp/vulcan_m3/mouse_contam.fasta"}
FOREIGN_C = "/tmp/vulcan_m3/foreign/human_contam.C.foreign.txt"   # M3 arm-C list (Arabidopsis, neighbour-filtered)
M1E_DB = "/tmp/vulcan_m1e/entrapment/db"


def is_tmt(family: str) -> bool:
    return family.endswith("_tmt")


def mods(family: str):
    """(static, var) as the GPU engine takes them (var '^' = peptide-N-term delta on mask bit 63)."""
    if is_tmt(family):
        return dict(CAM, K=TMT), dict(OX, **{"^": TMT})
    return dict(CAM), dict(OX)


def index_dir(family: str, kind: str) -> str:
    if family == "human":
        if kind in ("prod", "E2001"):
            return f"/tmp/vulcan_e2/index/human_contam_{kind}"
        return f"/tmp/vulcan_m3/index/human_contam_{kind}"
    if family == "mouse" and kind == "prod":
        return "/tmp/vulcan_m3/index/mouse_contam_prod"
    return os.path.join(WORK, "index", f"{family}_{kind}")


def ensure_index(family: str, kind: str) -> str:
    d = index_dir(family, kind)
    if os.path.exists(os.path.join(d, "meta.json")):
        return d
    assert is_tmt(family), d
    sp = family.split("_")[0]
    par = {"static_mods": dict(CAM, K=TMT)}
    if kind == "C":
        from vulcan.index.foreign import build_foreign
        build_foreign([FASTA[sp]], d, FOREIGN_C, par)
    else:
        from vulcan.index.build import build
        build([FASTA[sp]], d, dict(par, entrap_seeds=[] if kind == "prod" else [int(kind[1:])]))
    return d


def load_l1(family: str, kind: str):
    import numpy as np
    from vulcan.index.build import L1Index
    l1 = L1Index(index_dir(family, kind))
    if is_tmt(family):
        l1.mass32 = (l1.mass32.astype(np.float32) + np.float32(TMT)).astype(np.float32)
        assert (np.diff(l1.mass32) >= 0).all()
        l1.mask = l1.mask | np.uint64(1 << 63)
    return l1


def c_arm_r() -> float:
    import json
    return float(json.load(open(os.path.join(M1E_DB, "human_contam.C.json")))["r"])


def train_rows():
    with open(os.path.join(ROOT, "data", "manifests", "train_itcid.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def train_store(r):
    return os.path.join(WORK, "store_train", os.path.splitext(r["filename"])[0] + ".vspec.arrow")


def train_family(r):
    return r["species"] + ("_tmt" if r["stratum"] == "TMT" else "")


def d4_stores():
    import glob
    return sorted(glob.glob("/tmp/vulcan_m6/d4/store/*.vspec.arrow"))


def d10_stores():
    import glob
    return sorted(glob.glob("/tmp/vulcan_m5lc/store/Control_HEK_SEC_020118_*.vspec.arrow"))


DEV = {"D4": {"family": "human_tmt", "stores": d4_stores},
       "D10cf": {"family": "human", "stores": d10_stores}}


def psm_dir(tag: str) -> str:
    return os.path.join(WORK, "psms", tag)
