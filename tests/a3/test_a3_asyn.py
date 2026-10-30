"""A3.2.2 A-syn construction: deterministic draw, envelope exclusion, screens, and the P1 swap test (the construction
depends only on the unordered T∪D set)."""
import os
import zlib

import numpy as np

from vulcan.gates import a3_asyn as AS
from vulcan.gates.a2_identity import GateDB, peptide_mass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def test_offsets_table_loads_with_recorded_hash():
    from vulcan.gates import a3_config
    cfg = a3_config.load()
    o = AS.load_offsets(os.path.join(ROOT, cfg["asyn"]["offsets_table"]), cfg["asyn"]["offsets_sha256"])
    assert len(o) == cfg["asyn"]["offsets_rows"]
    assert np.all(np.diff(o.delta_mz) >= 0)                       # committed order = empirical CDF order
    assert set(np.unique(o.charge)) <= {2, 3, 4, 5, 6}
    try:
        AS.load_offsets(os.path.join(ROOT, cfg["asyn"]["offsets_table"]), "0" * 64)
        raise AssertionError("hash mismatch not detected")
    except RuntimeError:
        pass


def test_draw_is_deterministic_and_uniform():
    rows = [AS.draw_row(f"run|{i}", 2, 0, 1000) for i in range(20000)]
    assert rows == [AS.draw_row(f"run|{i}", 2, 0, 1000) for i in range(20000)]
    h = np.bincount(rows, minlength=1000)
    assert h.min() > 0 and abs(h.mean() - 20) < 1e-9
    assert AS.draw_row("k", 2, 0, 10**6) != AS.draw_row("k", 3, 0, 10**6)          # ranks are independent draws
    assert AS.draw_row("k", 2, 0, 10**6) != AS.draw_row("k", 2, 1, 10**6)          # the hash chain advances
    assert AS.draw_row("k", 2, 0, 10**6, seed=47) != AS.draw_row("k", 2, 0, 10**6, seed=48)


def test_envelope_exclusion():
    mz, z = 650.33, 2
    M = (mz - AS.PROTON) * z
    for k in (-1, 0, 1, 2, 3):
        for zz in (1, 2, 3, 6):
            assert AS.envelope_excluded(mz, z, (M + k * AS.NEUTRON) / zz + AS.PROTON, zz)
    # searched isotope 1 of a heavier draw lands on the envelope
    assert AS.envelope_excluded(mz, z, (M + 2 * AS.NEUTRON) / 3 + AS.PROTON, 3)
    assert AS.envelope_excluded(mz, z, mz + 5e-6 * mz, 2)              # within 20 ppm
    assert not AS.envelope_excluded(mz, z, mz + 0.37, 2)
    assert not AS.envelope_excluded(mz, z, mz - 1.21, 3)


def _db(T, D):
    return GateDB.from_lists(T, D)


def _family_screen(db, peaks_seed, mz_s, z_s, hs_p):
    """Label-free stand-in for the exhaustive family screen: a deterministic pseudo-hyperscore of every candidate in
    W' (function of the spectrum and the candidate's sequence only), max >= hyperscore(p) -> excluded."""
    best = -np.inf
    for r in db.window_rows(mz_s, z_s, 10.0, (0, 1, 2)):
        h = zlib.crc32(f"{peaks_seed}|{db.seq(r)}|{int(db.mask[r])}".encode()) % 1000 / 100.0
        best = max(best, h)
    return best >= hs_p


def _peptides(n, seed):
    g = np.random.default_rng(seed)
    aa = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    return ["".join(g.choice(aa, size=int(g.integers(7, 16)))) + "K" for _ in range(n)]


def test_construct_and_p1_swap():
    T = [(p, 0) for p in _peptides(3000, 1)]
    D = [(p[::-1][1:] + "K", 0) for p, _ in T]
    peps = _peptides(40, 2)
    keys = [f"runX|{i}" for i in range(40)]
    z = np.array([2] * 20 + [3] * 20)
    mz = np.array([peptide_mass(p) / zz + AS.PROTON for p, zz in zip(peps, z)])
    o = AS.Offsets(np.linspace(-2.5, 2.5, 2001), np.array([2, 3, 4] * 667), "test")
    out = []
    for db in (_db(T, D), _db(D, T)):                      # P1: the same construction whichever list is "T"
        idf = AS.identity_screen_fn(db, np.array(peps, dtype=object), np.zeros(40, np.int64))

        def screen(i, mz_s, z_s, db=db, idf=idf):
            return idf(i, mz_s, z_s), _family_screen(db, i, mz_s, z_s, 9.5)
        out.append(AS.construct(keys, mz, z, o, screen))
    assert AS.table_hash(out[0]) == AS.table_hash(out[1])
    t = out[0]
    assert set(np.unique(t["rank"])) <= {2, 3}
    assert t["counts"]["accepted"] == len(t["a_idx"])
    assert t["counts"]["accepted"] + t["counts"]["none_after_max_tries"] == 80
    for i, r, m, zz in zip(t["a_idx"], t["rank"], t["mz"], t["charge"]):
        assert not AS.envelope_excluded(mz[i], z[i], m, zz)


def test_construct_exclusion_order_and_exhaustion():
    o = AS.Offsets(np.array([0.0, 0.37]), np.array([2, 2]), "t")
    keys = ["a|1"]
    mz, z = np.array([500.0]), np.array([2])
    # row 0 (delta 0) is always envelope-excluded; row 1 passes env; screens exclude everything -> no hypothesis
    t = AS.construct(keys, mz, z, o, lambda i, m, zz: (True, False))
    assert t["counts"]["accepted"] == 0 and t["counts"]["none_after_max_tries"] == 2
    assert t["counts"]["env"] + t["counts"]["b"] == 2 * AS.MAX_TRIES
    t = AS.construct(keys, mz, z, o, lambda i, m, zz: (False, False))
    assert t["counts"]["accepted"] == 2 and np.allclose(t["mz"], 500.37)
