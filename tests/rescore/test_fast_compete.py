"""E7b: fast chimeric competition (numba ion bins + neighbour filter) equals the reference exactly."""
import numpy as np

from vulcan.chimeric.compete import IonSets, chimeric_competition
from vulcan.chimeric.fast_compete import chimeric_competition_fast, ion_bins_of, neighbour_filter_fast
from vulcan.ptm.conga import neighbour_filter

AA = np.frombuffer(b"ACDEFGHIKLMNPQRSTVWY", np.uint8)


class FakeL1:
    def __init__(self, n, seed):
        rng = np.random.default_rng(seed)
        self.length = rng.integers(7, 51, n).astype(np.uint8)
        self.seq_off = np.concatenate([[0], np.cumsum(self.length.astype(np.int64))[:-1]])
        self.residues = AA[rng.integers(0, len(AA), int(self.length.astype(np.int64).sum()))]
        m = np.zeros(n, np.uint64)
        for i in range(n):        # Met-ox bits on a few M positions
            o, L = self.seq_off[i], int(self.length[i])
            for k in np.flatnonzero(self.residues[o:o + L] == ord("M"))[:2]:
                if rng.random() < 0.5:
                    m[i] |= np.uint64(1) << np.uint64(k)
        self.mask = m


def test_ion_bins_match_reference():
    l1 = FakeL1(3000, 1)
    ions = IonSets(l1)
    rows = np.arange(3000)
    off, b = ion_bins_of(rows, ions)
    for i in rows:
        assert np.array_equal(np.array(sorted(ions.get(int(i))), np.int64), b[off[i]:off[i + 1]])


def test_competition_matches_reference():
    rng = np.random.default_rng(3)
    n_rows, n = 2000, 20000
    l1 = FakeL1(n_rows, 2)
    # shared ion sets between some rows: copy a few peptides so the neighbour filter has work to do
    spec = rng.integers(0, 3000, n).astype(np.int64)
    group = rng.integers(0, 4, n)
    rows = rng.integers(0, n_rows, n)
    pair = rows // 2
    dec = (rows % 2).astype(bool)
    score = np.round(rng.normal(size=n), 2)          # ties on purpose
    keep_ref = neighbour_filter(spec, score, IonSets(l1).sets(rows), 0.05)
    keep_new = neighbour_filter_fast(spec, score, rows, IonSets(l1), 0.05)
    assert np.array_equal(keep_ref, keep_new)
    a = chimeric_competition(spec, group, group + 1, rows, pair, score, dec, rows, IonSets(l1), seed=11, min_size=50)
    b = chimeric_competition_fast(spec, group, group + 1, rows, pair, score, dec, rows, IonSets(l1), seed=11,
                                  min_size=50, group_walk=True)
    for k in ("win", "survivors", "stratum", "merged", "q_sep", "q_pooled", "q_gw"):
        assert np.array_equal(a[k], b[k]), k
