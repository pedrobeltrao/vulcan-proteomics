"""A2.4 fragment-shifted nulls and A2.2(e) single-precursor hypotheses (synthetic Spectra; torch-free host code)."""
import numpy as np

from vulcan.gates.a2_fs import FS_DELTA_TH, a2_hypotheses, certified_precursor, fragment_shift
from vulcan.gpu.host import SearchParams, Spectra, make_hypotheses


def spectra(n=6, peaks=20, seed=0):
    rng = np.random.default_rng(seed)
    off = np.arange(n + 1, dtype=np.int64) * peaks
    return Spectra(run="syn", scan=np.arange(n, dtype=np.int64) + 1, index=np.arange(n, dtype=np.int64),
                   rt=np.linspace(1, 60, n), precursor_mz=rng.uniform(400, 900, n), charge=np.array([2, 0, 3, 2, 2, 0], np.int8)[:n],
                   mz=rng.uniform(150, 1500, n * peaks).astype(np.float32),
                   intensity=rng.uniform(1, 100, n * peaks).astype(np.float32), offsets=off)


def test_fragment_shift():
    sp = spectra()
    idx = np.array([0, 3])
    mono = np.array([500.25, 610.5])
    z = np.array([2, 3])
    fs = fragment_shift(sp, idx, mono, z)
    assert abs(FS_DELTA_TH - 2 / 7) < 1e-15
    for j, i in enumerate(idx):
        a, b = sp.offsets[i], sp.offsets[i + 1]
        c, d = fs.offsets[j], fs.offsets[j + 1]
        assert np.array_equal(fs.mz[c:d], (sp.mz[a:b] + np.float32(2 / 7)).astype(np.float32))
        assert np.array_equal(fs.intensity[c:d], sp.intensity[a:b])
    assert np.array_equal(fs.precursor_mz, mono) and fs.charge.tolist() == [2, 3]
    assert np.array_equal(fs.rt, sp.rt[idx]) and not set(fs.scan) & set(sp.scan)


def test_certified_precursor_copy():
    sp = spectra()
    c = certified_precursor(sp, np.array([1, 4]), np.array([450.0, 700.0]), np.array([2, 2]))
    assert c.charge.tolist() == [2, 2] and c.precursor_mz.tolist() == [450.0, 700.0]
    assert np.array_equal(c.mz[:20], sp.mz[20:40])


def test_single_hypothesis_for_A_and_FS():
    sp = spectra()
    mass32 = np.sort(np.random.default_rng(1).uniform(500, 3000, 5000)).astype(np.float32)
    p = SearchParams()
    valid = np.ones(len(sp), bool)
    single = np.array([True, False, False, True, False, False])
    h = a2_hypotheses(sp, valid, mass32, p, single)
    ref = make_hypotheses(sp, valid, mass32, p)
    for s in range(len(sp)):
        n = int((h["spec"] == s).sum())
        if single[s]:
            assert n == 1 and h["iso"][h["spec"] == s][0] == 0 and h["z"][h["spec"] == s][0] == sp.charge[s]
        else:
            assert n == int((ref["spec"] == s).sum())
    o = np.lexsort((h["iso"], h["z"], h["spec"]))
    assert np.array_equal(o, np.arange(o.size))


def test_single_needs_known_charge():
    import pytest
    sp = spectra()
    with pytest.raises(ValueError):
        a2_hypotheses(sp, np.ones(len(sp), bool), np.linspace(500, 3000, 100).astype(np.float32), SearchParams(),
                      np.array([False, True, False, False, False, False]))
