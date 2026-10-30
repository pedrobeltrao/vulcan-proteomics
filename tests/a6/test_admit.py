"""A6.4.6 unit tests (synthetic labels): symmetric I, the |I| count assertion, the label guard, and the reduced index
row subset.  No real data."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "a6"))
sys.path.insert(0, ROOT)

import admit_hn1 as AH  # noqa: E402
from reduced_index import ReducedL1, pair_close  # noqa: E402
from vulcan.stats.tdc import double_competition  # noqa: E402


def _toy(thr_gap=0.0):
    """Units 0..2*n-1; pair k = (2k, 2k+1); even = target, odd = decoy.  One spectrum per unit (rank-1 = itself).
    Pair 0: decoy wins with 30 (above thr 20).  Pair 1: target wins with 19.999 (just below).  Plus filler targets
    so that thr is defined: many targets at 50 (accepted)."""
    units, scores = [], []
    units += [0, 1]
    scores += [10.0, 30.0]            # pair 0: decoy (unit 1) wins, 30 >= thr
    units += [2, 3]
    scores += [20.0 - 1e-3, 5.0]      # pair 1: target wins just below thr
    for k in range(2, 202):           # filler pairs, target 50 >> decoy 1
        units += [2 * k, 2 * k + 1]
        scores += [50.0, 1.0]
    units = np.array(units)
    spec = np.arange(units.size)
    hs = np.array(scores)
    pair = units // 2
    dec = (units % 2).astype(bool)
    return spec, units, pair, hs, dec


def test_symmetric_I_decoy_won_in_target_below_out_and_swap_invariant():
    spec, units, pair, hs, dec = _toy()
    tb = units.astype(np.uint64) * np.uint64(2654435761)
    thr = 20.0
    dc = double_competition(spec, units, pair, hs, dec, tiebreak=tb)
    pk = pair                         # label-free pair key
    I = AH.identified_pairs(dc, pk, thr)
    assert 0 in I                     # decoy-won pair above threshold is in I
    assert 1 not in I                 # target-won pair just below threshold is not
    AH.check_I_count(I, spec, hs, pk, tb, thr)
    dc_sw = double_competition(spec, units, pair, hs, ~dec, tiebreak=tb)
    I_sw = AH.identified_pairs(dc_sw, pk, thr)
    assert np.array_equal(I, I_sw)    # swapping labels leaves I unchanged


def test_count_assertion_catches_target_only_I():
    """A pilot-style I built from accepted TARGETS only (drops the decoy-won pair) must fail the count assertion."""
    spec, units, pair, hs, dec = _toy()
    tb = units.astype(np.uint64) * np.uint64(2654435761)
    dc = double_competition(spec, units, pair, hs, dec, tiebreak=tb)
    I_target_only = np.unique(pair[dc.row[(dc.score >= 20.0) & ~dc.is_decoy]])
    with pytest.raises(AssertionError):
        AH.check_I_count(I_target_only, spec, hs, pair, tb, 20.0)


def test_label_view_refuses_non_host_rows():
    lv = AH.LabelView(np.array([False, True, False]), np.array([0, 0, 1]))
    d, p = lv.at([0, 1], "host_full_search")
    assert d.tolist() == [False, True]
    for origin in ("reduced_index_search", "null_copy", "mixture"):
        with pytest.raises(PermissionError):
            lv.at([0], origin)
    assert sum(x.size for x in lv.rows_read) == 2


class _L1:
    def __init__(self, n=10):
        self.pf_id = np.arange(n, dtype=np.uint64)
        self.mass32 = np.linspace(500, 600, n).astype(np.float32)
        self.seq_off = np.arange(n, dtype=np.int64) * 7
        self.length = np.full(n, 7, np.uint8)
        self.mask = np.zeros(n, np.uint64)
        self.n_var = np.zeros(n, np.uint8)
        self.missed = np.zeros(n, np.uint8)
        self.partner = np.array([(i + 1 if i % 2 == 0 else i - 1) for i in range(n)], np.int64)
        self.residues = np.frombuffer(b"PEPTIDE" * n, np.uint8)
        self.meta = {}

    def __len__(self):
        return len(self.pf_id)


def test_reduced_index_subset_and_pair_closure():
    l1 = _L1()
    flag = np.zeros(10, bool)
    flag[2] = True
    with pytest.raises(ValueError):
        ReducedL1(l1, flag)           # not pair-closed
    R = pair_close(l1, flag)
    assert R[2] and R[3] and R.sum() == 2
    red = ReducedL1(l1, R)
    assert len(red) == 8
    assert np.array_equal(red.row_map, np.array([0, 1, 4, 5, 6, 7, 8, 9]))
    # partner remapped inside the reduced index and maps back to the full partner
    assert np.array_equal(red.row_map[red.partner.astype(np.int64)], l1.partner[red.row_map])


def test_no_label_argument_after_certification():
    """project_removal and null_stage take no label array (signature check)."""
    import inspect
    for f in (AH.project_removal, AH.null_stage, AH.freeze):
        names = set(inspect.signature(f).parameters)
        assert not names & {"labels", "lab", "is_decoy", "tdc_pair", "dec"}, f.__name__
