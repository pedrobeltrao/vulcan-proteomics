"""Unit tests for the T1 driver helpers (bench/t1): peptide parsing, replicate naming, matched-FDP interpolation,
vectorized paired counts == vulcan.stats.entrapment.paired_counts."""
import numpy as np
import pandas as pd

from bench.t1 import comparators as CMP
from bench.t1 import dbs as DB
from bench.t1 import evaluate as EV


def test_parse_peptide():
    assert CMP.parse_peptide("PEPM[+15.9949]IK") == ("PEPMLK", 1 << 3)
    assert CMP.parse_peptide("-.AC[57.0215]M[15.9949]R.-") == ("ACMR", 1 << 2)
    assert CMP.parse_peptide("K.LMMK.R") == ("LMMK", 0)


def test_replicates():
    r = DB.replicate_names([1001, 1002], [2001, 2002, 2003])
    assert len(r) == 10 and r[0] == "prod_1001" and r[4] == "C_1001"
    assert DB.parse_rep("P_1002_2003") == ("P", 1002, 2003)


def test_matched_ids():
    fdp = np.interp(EV.GRID, [0, 0.25], [0, 0.25])          # FDP == q
    ids = EV.GRID * 1000
    v, q, ok = EV.matched_ids(fdp, ids, 0.01)
    assert ok and abs(q - 0.01) < 1e-9 and abs(v - 10) < 1e-6


def test_paired_grid_matches_official():
    rng = np.random.default_rng(1)
    n = 3000
    tseq = [f"S{i}" for i in rng.integers(0, 1500, n)]
    fam = rng.choice([0, 1], n, p=[0.8, 0.2])
    score = rng.normal(size=n)
    order = np.argsort(-score)
    q = np.empty(n)
    q[order] = np.linspace(0.0, 0.2, n)               # q non-increasing in score
    df = pd.DataFrame({"seq": [f"x{i}" for i in range(n)], "mask": np.zeros(n, np.uint64), "is_decoy": False,
                       "q": q, "score": score, "runbits": rng.integers(1, 4, n), "fam": fam, "target_seq": tseq,
                       "nvar": 0})
    g = EV.paired_grid(df, "sequence", 2, check=True)   # raises if != official at 1% and 5%
    assert g.shape == (2, len(EV.GRID), 4)
