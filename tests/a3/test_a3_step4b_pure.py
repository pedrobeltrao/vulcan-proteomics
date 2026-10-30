"""Pure helpers of scripts/a3/step4b_chimeric (label-free)."""
import numpy as np

import step4b_chimeric as S4
from vulcan.gpu.host import Spectra


def _sp(n=5):
    lens = np.array([3, 1, 4, 2, 5])[:n]
    off = np.r_[0, np.cumsum(lens)].astype(np.int64)
    return Spectra("r", np.arange(n) + 100, np.arange(n) + 7, np.linspace(1, 2, n), np.linspace(400, 800, n),
                   np.array([2, 3, 2, 0, 4], np.int8)[:n], np.arange(off[-1], dtype=np.float32) + 100,
                   np.ones(off[-1], np.float32), off)


def test_subset_spectra_peaks_and_precursor():
    sp = _sp()
    v = S4.subset_spectra(sp, np.array([2, 0]), pmz=np.array([1.0, 2.0]), charge=np.array([5, 6]))
    assert v.offsets.tolist() == [0, 4, 7]
    assert v.mz.tolist() == sp.mz[4:8].tolist() + sp.mz[0:3].tolist()
    assert v.precursor_mz.tolist() == [1.0, 2.0] and v.charge.tolist() == [5, 6] and v.index.tolist() == [9, 7]
    w = S4.subset_spectra(sp, np.array([1]))
    assert w.precursor_mz[0] == sp.precursor_mz[1] and w.charge[0] == sp.charge[1]


def test_offset_delta_mass_recovers_recipe():
    host_mz, z = np.array([500.0, 700.0]), np.array([2, 3])
    delta = np.array([1.53, -0.47])
    mix_mz = host_mz + delta / z
    assert np.allclose(S4.offset_delta_mass(mix_mz, host_mz, z), delta)


def test_keys_and_unit_counts_are_label_free():
    keys = S4.a_spec_keys({0: "runA", 1: "runB"}, np.array([1, 0]), np.array([5, 9]))
    assert keys.tolist() == ["runB|5", "runA|9"]
    prov = np.array([0, 20, 21, 21])
    spec = np.array([0, 0, 1, 2, 2, 3, 3, 3])
    grp = np.array([0, 0, 0, 0, 1, 0, 1, 2])
    rank = grp + 1
    c = S4.unit_counts(prov, rank, spec, grp)
    assert c["host"] == {"1": 1, "2": 0, "3": 0} and c["A"] == {"1": 1, "2": 0, "3": 0}
    assert c["B"] == {"1": 2, "2": 2, "3": 1}


def test_neutral_rt_equals_gates_a1():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "m3e"))
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "m3"))
    from gates_a1 import RT_FEATS, neutral_rt
    g = np.random.default_rng(0)
    f = {k: g.normal(size=50).astype(np.float32) for k in ("rt_abs_resid", "rt_resid_z", "pred_sa")}
    a, b = S4.neutral_rt(f), neutral_rt(f)
    assert tuple(RT_FEATS) == S4.RT_FEATS
    assert all(np.array_equal(a[k], b[k]) for k in f)
