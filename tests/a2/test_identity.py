"""A2.2(b) identity screen pieces: fragment ladders, f(c), windows, isobaric Hamming (vs brute force)."""
import random

import numpy as np

from vulcan.gates import a2_identity as I


def test_sage_masses_match_engine_table():
    import re
    src = open(I.__file__.replace("gates/a2_identity.py", "gpu/engine.py")).read()
    block = src[src.index("SAGE_MONO = np.array(["):]
    block = block[:block.index("np.float32)")]
    vals = [float(x) for x in re.findall(r"[0-9]+\.[0-9]+", block)]
    assert np.allclose(np.array(vals, np.float32).astype(float), I.SAGE_MONO)


def test_ion_ladder_and_mass():
    seq = "PEPTIDEK"
    r = I.residue_masses(seq)
    ions = I.ion_ladder(seq, first=1)
    L = len(seq)
    assert ions.size == 2 * (L - 1)
    assert np.isclose(ions[0], r[0] + I.PROTON)
    assert np.isclose(ions[L - 1], r[-1] + I.H2O + I.PROTON)          # y1
    assert np.isclose(I.peptide_mass(seq), r.sum() + I.H2O)
    # variable Met-ox only where the mask bit is set; static CAM on C
    assert np.isclose(I.residue_masses("MCM", 0b100)[2] - I.residue_masses("MCM", 0)[2], 15.9949)
    assert np.isclose(I.residue_masses("C")[0], I.SAGE_MONO[2] + 57.021464)


def test_f_overlap_self_and_reverse():
    p = "ELVISLIVESK"
    pi = I.ion_ladder(p, first=2)
    assert I.f_overlap(pi, I.ion_ladder(p, first=1)) == 1.0
    rev = p[-2::-1] + p[-1]
    f = I.f_overlap(pi, I.ion_ladder(rev, first=1))
    assert 0.0 <= f < 0.5


def test_window_and_max_f():
    T = [("ELVISLIVESK", 0), ("AAAAGGGGK", 0), ("PEPTIDEKR", 0)]
    D = [(s[-2::-1] + s[-1], m) for s, m in T]
    db = I.GateDB.from_lists(T, D)
    p = "ELVISLIVESK"
    mz = (I.peptide_mass(p) + 2 * I.PROTON) / 2
    rows = db.window_rows(mz, 2)
    assert {db.seq(r) for r in rows} == {"ELVISLIVESK", "SEVILSIVLEK"}
    assert I.max_f_in_window(db, p, 0, mz, 2) == 1.0
    # isotope error 1: a peptide 1 neutron lighter than the observed mono is in W
    mz_iso = mz + I.NEUTRON / 2
    assert I.max_f_in_window(db, p, 0, mz_iso, 2) == 1.0
    assert I.max_f_in_window(db, p, 0, mz + 5.0, 2) == 0.0


def test_gate_db_label_swap_invariant():
    T = [("ELVISLIVESK", 0), ("MAMMALK", 0b1), ("PEPTIDEKR", 0)]
    D = [(s[-2::-1] + s[-1], m) for s, m in T]
    a, b = I.GateDB.from_lists(T, D), I.GateDB.from_lists(D, T)
    assert a.fingerprint() == b.fingerprint()
    assert np.array_equal(a.mass, b.mass) and np.array_equal(a.residues, b.residues)


def test_isobaric_variants():
    v = I.isobaric_variants("ANK")
    assert set(v) == {"ANK", "AGGK"}
    v = I.isobaric_variants("AGGK")
    assert "ANK" in v and "QGK" in v and "AGGK" in v           # AG->Q and GG->N
    assert I.isobaric_variants("NNNNNNNNNNNNN", cap=100) is None


def test_hamming_isobaric_matches_bruteforce():
    rng = random.Random(5)
    AA = "ACDEFGHIKLMNPQRSTVWY"
    peps = ["".join(rng.choice(AA) for _ in range(rng.randint(7, 14))) for _ in range(400)]
    ph = I.PeptideHamming(peps)
    queries = []
    for p in peps[:120]:
        q = list(p)
        for _ in range(rng.randint(0, 3)):
            q[rng.randrange(len(q))] = rng.choice(AA)
        queries.append("".join(q).replace("GG", "N", 1))
    queries += ["".join(rng.choice(AA) for _ in range(rng.randint(7, 12))) for _ in range(200)]
    for q in queries:
        assert I.hamming_isobaric_near(ph, q) == I.hamming_isobaric_bruteforce(q, peps), q


def test_hamming_il_collapsed():
    ph = I.PeptideHamming(["LLLLLLLK"])
    assert ph.near("IIIIIIIK")
    assert I.hamming_isobaric_near(ph, "IIIAAIIK")          # 2 substitutions
    assert not I.hamming_isobaric_near(ph, "IIAAAIIK")      # 3
