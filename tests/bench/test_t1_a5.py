"""bench/t1/a5.py (Amendment A5, DESIGN §6.0d A5.11.2.5): synthetic unit tests for the claim-level arm-C quantity and
its first-crossing threshold (A5.7.1 a-f), the stratum-1 restriction (A5.5.5), the A1.8 rank-change fractions and
the D18 ladder / taxon screen."""
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from bench.t1 import a5 as A5
from bench.t1 import evaluate as EV


def _curve(fdp_by_grid, ids_by_grid, R=1, J=5, S=3):
    """cnt such that fdp_paired = fdp exactly: n_t = 1000, n_e = fdp*1000, n_eo = n_egt = 0 -> fdp = n_e/(n_t+n_e)."""
    G = len(EV.GRID)
    cnt = np.zeros((R, J, S, G, 4))
    # fdp = n_e/(n_t+n_e)  ->  n_e = fdp*n_t/(1-fdp)
    nt = 1000.0
    cnt[..., 0] = nt / R
    cnt[..., 1] = (fdp_by_grid * nt / (1 - fdp_by_grid))[None, None, None, :] / R
    ids = np.broadcast_to(ids_by_grid / R, (R, J, G)).copy()
    return {"ids": {lv: ids for lv in EV.LEVELS}, "cnt": {lv: cnt for lv in EV.LEVELS}}


def test_q_at_match_first_crossing_and_equals_evaluate():
    fdp = EV.GRID.copy()                        # FDP == q
    c = _curve(fdp, EV.GRID * 1000)
    q, reached = A5.q_at_match(c, "peptidoform", 0.01)
    assert reached and abs(q - 0.01) < 1e-9
    # non-monotone curve: the FIRST crossing counts
    f2 = np.where(EV.GRID < 0.02, EV.GRID * 0.5, 0.0)
    f2[np.argmin(np.abs(EV.GRID - 0.003))] = 0.02
    q2, _ = A5.q_at_match(_curve(f2, EV.GRID * 1000), "peptidoform", 0.01)
    assert q2 < 0.003 + 1e-12
    # never reached -> grid end, reached False
    q3, r3 = A5.q_at_match(_curve(np.full(len(EV.GRID), 0.001), EV.GRID), "peptidoform", 0.05)
    assert not r3 and q3 == EV.GRID[-1]
    # identical to evaluate.noninferiority's ours_q_at_match
    curves = {"a": _curve(fdp, EV.GRID * 1000), "b": _curve(fdp * 0.5, EV.GRID * 900)}
    out = {"methods": {"b": {"verdicts": {lv: {} for lv in EV.LEVELS}}}}
    ni = EV.noninferiority(out, curves, ["a", "b"], "a")
    for f in (0.01, 0.05):
        assert np.isclose(ni["b"]["peptidoform"]["matched_fdp"][str(f)]["ours_q_at_match"],
                          A5.q_at_match(curves["a"], "peptidoform", f)[0], rtol=1e-12)


def _combined_tab(n_t, n_e, q_t, q_e):
    rows = [(f"T{i}", 0, False, q_t) for i in range(n_t)] + [(f"F{i}", 0, False, q_e) for i in range(n_e)]
    df = pd.DataFrame(rows, columns=["seq", "mask", "is_decoy", "q"])
    df["fam"] = [0] * n_t + [-1] * n_e
    df["mask"] = df["mask"].astype(np.uint64)
    return df


def test_armc_lb_threshold_and_quantity():
    # 5 replicates: 990 T at q 0.001, 10 F at q 0.001, 100 F at q 0.04 -> LB at q*=0.01 is 10/1000
    tabs = [(j, _combined_tab(990, 10, 0.001, 0.001)) for j in range(1001, 1006)]
    extra = [(j, pd.concat([t, _combined_tab(0, 100, 0.0, 0.04).assign(seq=[f"G{i}" for i in range(100)])],
                           ignore_index=True)) for j, t in tabs]
    d = A5.armc_lb(extra, "peptidoform", 0.01, r=1.0, f=0.01)
    assert np.allclose(d["lb_reps"], 0.01) and np.isclose(d["mean_lb"], 0.01)
    assert d["r"] == 1.0 and np.isclose(d["mean_fdp_c"], 10 * 2 / 1000)
    # identical replicates: no between-seed variance; CI (a) contains the mean
    assert d["lb_ci95"][0] <= 0.01 <= d["lb_ci95"][1]
    # at q*=0.05 the 100 extra F count: LB = 110/1100 = 0.1, lower CI far above f=0.05 -> void boundary fires
    d5 = A5.armc_lb(extra, "peptidoform", 0.05, r=1.0, f=0.05)
    assert np.isclose(d5["mean_lb"], 0.1) and d5["lb_ci95"][0] > 0.05


def test_stratum1_view_restricts_rows():
    tab = pd.DataFrame({"seq": list("abcdef"), "q": np.linspace(0, 0.05, 6)})
    cc = {"merge": {1: [1], 2: [2, 3], 3: [2, 3]}, "merged": ["1", "2_3", "1", "2_3", "1", "2_3"]}
    s1, lab = A5._s1_view(tab, cc)
    assert lab == "1" and list(s1.seq) == ["a", "c", "e"]
    cc2 = {"merge": {1: [1, 2], 2: [1, 2], 3: [3]}, "merged": ["1_2", "1_2", "3", "1_2", "3", "3"]}
    s1b, lab2 = A5._s1_view(tab, cc2)
    assert lab2 == "1_2" and list(s1b.seq) == ["a", "b", "d"]


def test_stratum1_paired_counts_only_stratum1_pairs():
    """An entrapment in stratum 1 whose target partner is in stratum 2 is an E-only entry of the stratum-1 list."""
    tab = pd.DataFrame({"seq": ["T1", "E1", "T2", "E2"], "mask": np.zeros(4, np.uint64), "is_decoy": False,
                        "q": [0.001] * 4, "score": [5.0, 6.0, 5.0, 4.0], "runbits": 1, "fam": [0, 1, 0, 1],
                        "target_seq": ["T1", "T1", "T2", "T2"], "nvar": 0})
    cc = {"merge": {1: [1], 2: [2], 3: [3]}, "merged": ["2", "1", "1", "1"]}
    s1, _ = A5._s1_view(tab, cc)
    g = EV.paired_grid(s1, "peptidoform", 1, check=True)
    k = int(np.flatnonzero(np.isclose(EV.GRID, 0.01))[0])
    n_t, n_e, n_eo, n_egt = g[0, k]
    assert (n_t, n_e, n_eo, n_egt) == (1, 2, 1, 0)        # T2 vs E2 (E2 lower), E1 alone (T1 not in stratum 1)


def test_rank_v1_fractions(tmp_path):
    class RL:                                   # minimal RepLabels: 4 pf rows, rows 2/3 are decoys of 0/1
        is_decoy_row = np.array([False, False, True, True])
        tdc_pair = np.array([2, 3, 0, 1])
    # 3 spectra; spectrum 0: rank-2 PSM scores best; spectrum 1: rank 1 best; spectrum 2: decoy rank 1 best
    psm = pd.DataFrame({"run": ["r"] * 6, "spec_idx": [0, 0, 1, 1, 2, 2], "rank": [1, 2, 1, 2, 1, 2],
                        "pf_row": [0, 1, 1, 0, 2, 3], "pf_id": [0, 1, 1, 0, 2, 3],
                        "score": [1.0, 9.0, 8.0, 0.5, 7.0, 0.1]})
    d = tmp_path / "prod"
    os.makedirs(d)
    pq.write_table(pa.Table.from_pandas(psm), d / "r.scored.parquet")
    out = A5.rank_v1(str(d), RL(), ["r"])
    assert out["n_spectra"] == 3 and np.isclose(out["frac_spectra_rank_gt1"], 1 / 3)
    assert out["frac_accepted_units_rank_gt1"] is None or 0.0 <= out["frac_accepted_units_rank_gt1"] <= 1.0


def test_ladder_screen():
    t = ["PEPTIDEKAAR", "SAMPLEMASSK"]
    f = ["PEPTWWWWWK",          # shares prefix PEPT but mass far from any target -> kept
         "AAAAGGGWDKAAR",       # suffix KAAR with... mass far -> kept
         "PEPTLDEKAAR",         # I/L-collapsed identical to the target: same mass, shares prefix -> excluded
         "WWWWSAMPLE"]
    hit = A5.ladder_screen(t, f)
    assert list(hit) == [False, False, True, False]
    # same mass and same suffix (permuted middle) -> excluded
    hit2 = A5.ladder_screen(["PEPTIDEKAAR"], ["PETPIDEKAAR"])
    assert hit2[0]
    # same mass, neither terminal 4-mer shared -> kept
    hit3 = A5.ladder_screen(["PEPTIDEKAAR"], ["EPPTIDEKARA"])
    assert not hit3[0]


def test_taxon_pattern():
    assert A5._TAXON.search("sp|P10795|RBS1A_ARATH Ribulose bisphosphate carboxylase small subunit 1A, chloroplastic")
    assert A5._TAXON.search("sp|X|Y_ARATH Something OS=Arabidopsis GN=rbcL PE=1")
    assert A5._TAXON.search("sp|X|PSBA Photosystem II protein D1")
    assert not A5._TAXON.search("sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens GN=HBB")
