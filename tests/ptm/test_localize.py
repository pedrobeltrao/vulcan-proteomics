"""L3 localization operator: exhaustive configurations, true-site recovery on synthetic spectra, and the decoy-AA
G0 cases (scoring blind to residue identity given masses; DESIGN §3 label blindness, §4 L3, §6.5)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from vulcan.ptm import localize as LZ  # noqa: E402

DEV = "cpu"


def _spectrum(seq, sites, static, table=None, nl=True, z=2, seed=0, noise=40):
    rng = np.random.default_rng(seed)
    m = LZ.residue_masses(seq, static, table=table)
    m[list(sites)] += LZ.PHOSPHO
    pre = np.cumsum(m)[:-1]
    tot = m.sum() + LZ.H2O
    ions = []
    for j, b in enumerate(pre):
        y = tot - b
        ions += [b + LZ.PROTON, y + LZ.PROTON]
        nb = sum(1 for s in sites if s <= j)
        if nl and nb:
            ions.append(b - LZ.H3PO4 + LZ.PROTON)
        if nl and len(sites) - nb:
            ions.append(y - LZ.H3PO4 + LZ.PROTON)
    mz = np.concatenate([np.array(ions), rng.uniform(150, 1500, noise)])
    it = np.concatenate([rng.uniform(50, 100, len(ions)), rng.uniform(1, 60, noise)])
    o = np.argsort(mz)
    return mz[o], it[o]


def _run(seqs, sites_true, acceptors, static=None, table=None, k=1):
    static = static or {"C": 57.021464}
    p = LZ.LocParams()
    mzs, its, rm, cp = [], [], [], []
    for i, (s, st) in enumerate(zip(seqs, sites_true)):
        mz, it = _spectrum(s, st, static, table=table, seed=i)
        mzs.append(mz); its.append(it)
        rm.append(LZ.residue_masses(s, static, table=table))
        cp.append(LZ.candidate_positions(s, acceptors))
    peaks = LZ.prepare_peaks(mzs, its, p)
    n = len(seqs)
    sc = LZ.score_configs(rm, np.zeros(n), cp, np.full(n, k), np.full(n, 2), np.arange(n), peaks, p, device=DEV)
    return sc, LZ.summarize(sc, cp)


def test_enumeration_is_exhaustive():
    assert LZ.enumerate_configs(5, 2).shape == (10, 2)
    assert LZ.enumerate_configs(7, 3).shape == (35, 3)


def test_recovers_true_site():
    seqs = ["AGSPTLYKR", "VSSAPTEGR", "LLQSAYTPK", "ASDTSPEAYR"]
    sites = [(2,), (1,), (5,), (4,)]
    sc, sm = _run(seqs, sites, "STYA")
    for s, st in zip(sm, sites):
        assert tuple(s["best_cfg"]) == st
        assert s["site_delta"][0] > 0


def test_decoy_aa_residue_blind():
    """G0 decoy-AA case: with a mass table where A weighs exactly what S weighs, a peptide with A at position i and
    the same peptide with S at position i give bit-identical configuration scores (no letter-specific ion, loss or
    rule), whichever residue is called the decoy."""
    tab = dict(LZ.RESIDUE)
    tab["A"] = tab["S"]
    a_seq, s_seq = "VSAPTEGLR", "VSSPTEGLR"          # position 2: A vs S, same mass under `tab`
    sc_a, _ = _run([a_seq], [(2,)], "STYA", table=tab)
    sc_s, _ = _run([s_seq], [(2,)], "STYA", table=tab)
    assert np.array_equal(sc_a[0]["cfg"], sc_s[0]["cfg"])
    assert np.array_equal(sc_a[0]["score"], sc_s[0]["score"])
    # and with the true (different) masses the decoy-A configuration still gets the same ion types: the neutral loss
    # applies to the A-placed mod exactly as to S/T/Y (scores differ only through masses)
    tab2 = dict(LZ.RESIDUE)
    tab2["Y"] = tab2["A"]
    sc_y, _ = _run(["VSYPTEGLR"], [(2,)], "STYA", table=tab2)
    assert np.array_equal(sc_y[0]["score"], _run(["VSAPTEGLR"], [(2,)], "STYA", table=tab2)[0][0]["score"])


def test_decoy_aa_permutation_of_letters():
    """Letters of candidate positions never reach the scorer: permuting S/T/Y/A among candidate positions with
    their masses carried along (here: all acceptor masses equal) leaves scores identical."""
    tab = dict(LZ.RESIDUE)
    for c in "TYA":
        tab[c] = tab["S"]
    s1, s2 = "ASTPYEGK", "TYAPSEGK"
    sc1, _ = _run([s1], [(1,)], "STYA", table=tab)
    sc2, _ = _run([s2], [(1,)], "STYA", table=tab)
    assert np.array_equal(sc1[0]["score"], sc2[0]["score"])


def test_multi_site():
    sc, sm = _run(["ASPTSEGSYKR"], [(1, 8)], "STY", k=2)
    assert tuple(sm[0]["best_cfg"]) == (1, 8)
