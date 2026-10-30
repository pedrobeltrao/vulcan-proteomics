"""M5 LiP: semi-specific digest classes, semi index G0 (label swap / row permutation), half-tryptic calls."""
import os

import numpy as np
import pytest

from vulcan.lip import halftryptic as H
from vulcan.lip.digest import digest, sites

SNCA = ("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAP"
        "QEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA")
P2 = "MSTRPEPTIDEKAAGGLLWPKRNNNSTVQEDKLLMMSGHAAGGFWER"


@pytest.fixture(scope="module")
def fasta(tmp_path_factory):
    d = tmp_path_factory.mktemp("lip")
    p = d / "t.fasta"
    p.write_text(f">sp|P37840|SYUA_HUMAN\n{SNCA}\n>sp|X00001|TEST\n{P2}\n")
    return str(p)


def test_sites_trypsin_p():
    s = sites("AAKPAARAAK")
    assert s == [0, 7, 10]          # K-P not cleaved


def test_semi_classes(fasta):
    d = digest([fasta], semi=True, workers=1)
    cls = dict(zip(d["keys"], d["cls"].tolist()))
    assert cls["EGVVAAAEK"] == 0                    # specific
    assert cls["GVVAAAEK"] == 1                     # N-semi (N-end not at a site)
    assert cls["EGVVAAAE"] == 2                     # C-semi
    spec = digest([fasta], semi=False, workers=1)
    assert set(spec["keys"]) <= set(d["keys"])
    assert all(c == 0 for c in spec["cls"])


def test_semi_index_g0(fasta, tmp_path):
    pytest.importorskip("pyarrow")
    from vulcan.lip.index import build_semi
    a = build_semi([fasta], str(tmp_path / "a"), workers=1)
    build_semi([fasta], str(tmp_path / "b"), workers=1, shuffle_rows_seed=7, swap_labels=True)
    for f in ("peptidoforms.arrow", "residues.npy", "lip_class.parquet"):
        assert open(tmp_path / "a" / f, "rb").read() == open(tmp_path / "b" / f, "rb").read(), f
    n = a["n_peptidoforms_by_family"]
    assert n["T"] == n["DT"]


def test_halftryptic_calls():
    db = H.ProteinDB(["sp|P37840|SYUA_HUMAN", "sp|X00001|TEST"], [SNCA, P2],
                     {"X00001": {"n": {P2.index("NNSTVQ")}, "c": set()}})
    il = lambda s: s.replace("I", "L")
    c = H.call("EGVVAAAEK", [0], db)
    assert c["call"] == "tryptic" and c["start"] == 13
    c = H.call("GVVAAAEK", [0], db)                 # N-side PK cut after E13
    assert c["call"] == "half_tryptic" and c["cut_side"] == "N" and c["cut_site"] == 13
    c = H.call("EGVVAAAE", [0], db)                 # C-side cut after E20
    assert c["call"] == "half_tryptic" and c["cut_side"] == "C" and c["cut_site"] == 20
    c = H.call("DVFMK", [0], db)                     # Met removal: excluded
    assert c["call"] == "bio_terminal" and "MetRemoval" in c["bio_reason"]
    c = H.call(il("NNSTVQEDK"), [1], db)             # annotated mature-chain start
    assert c["call"] == "bio_terminal"
    c = H.call(il("NSTVQEDK"), [1], db)
    assert c["call"] == "half_tryptic" and c["cut_site"] == P2.index("NSTVQ")


def test_insource_flag():
    import pandas as pd
    prec = pd.DataFrame({"run": ["r"] * 3, "target_seq": ["EGVVAAAEK", "GVVAAAEK", "VVAAAEK"],
                         "call": ["tryptic", "half_tryptic", "half_tryptic"], "rt_apex": [10.0, 10.05, 12.0]})
    f = H.insource_flags(prec)
    assert f.tolist() == [False, True, False]
