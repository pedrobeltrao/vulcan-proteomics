"""E7 chimeric search (vulcan.chimeric): closed search unchanged, every co-isolated precursor gets its own candidate
list, G0 label-swap invariance, and the neighbour-filter ion convention.  Skipped without CUDA."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
torch = pytest.importorskip("torch")
pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")

FASTA = os.path.join(ROOT, "data", "fasta", "contaminants_frankenfield2022.fasta")
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}


def _frags(l1, r, rng, keep=0.7):
    from vulcan.index._core import core
    seq, mask = l1.sequence(r), int(l1.mask[r])
    ref = [x for x in core.sage_reference([seq], [("C", 57.021464)], [("M", 15.9949)], 2) if x[2] == mask][0]
    f = np.array(ref[4] + ref[5], np.float32) + np.float32(1.0072764)
    return f[rng.random(len(f)) < keep]


def _chimeric_spectra(l1, n=200, seed=7):
    """Each spectrum = fragments of peptide A (primary precursor) + peptide B (co-isolated) + noise."""
    from vulcan.gpu.engine import Spectra
    rng = np.random.default_rng(seed)
    rows = rng.choice(len(l1), 2 * n, replace=False)
    A, B = rows[:n], rows[n:]
    mz, it, off, pmz = [], [], [0], []
    for a, b in zip(A, B):
        m = np.sort(np.concatenate([_frags(l1, a, rng), _frags(l1, b, rng), rng.uniform(150, 2000, 60).astype(np.float32)]))
        mz.append(m)
        it.append(rng.uniform(1, 100, len(m)).astype(np.float32))
        off.append(off[-1] + len(m))
        pmz.append((float(l1.mass32[a]) + 2 * 1.0072764) / 2)
    sp = Spectra("synthetic", np.arange(n, dtype=np.uint32), np.arange(n, dtype=np.uint32), np.zeros(n),
                 np.array(pmz), np.full(n, 2, np.int8), np.concatenate(mz), np.concatenate(it), np.array(off, np.int64))
    # one secondary group per spectrum: peptide B at z = 3 (feature mono m/z)
    sec = {"spec": np.arange(n, dtype=np.int64), "group": np.ones(n, np.uint16),
           "mz": (l1.mass32[B].astype(np.float64) + 3 * 1.0072764) / 3, "charge": np.full(n, 3, np.int8),
           "feature_id": np.arange(n, dtype=np.int64), "isolation_fraction": np.full(n, 0.3, np.float32),
           "feature_intensity": np.ones(n, np.float32), "iso_cos": np.ones(n, np.float32),
           "rt_apex": np.zeros(n, np.float32), "nfw": np.full(n, 2, np.uint16)}
    prim = {"feature_id": np.full(n, -1, np.int64), "isolation_fraction": np.zeros(n, np.float32),
            "feature_intensity": np.zeros(n, np.float32), "iso_cos": np.zeros(n, np.float32),
            "rt_apex": np.full(n, np.nan, np.float32), "nfw": np.full(n, 2, np.uint16)}
    return sp, {"prim": prim, "sec": sec, "n_groups": np.full(n, 2, np.uint16)}, A, B


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    from vulcan.index.build import build
    d = tmp_path_factory.mktemp("idx")
    build([FASTA], str(d / "a"), {"entrap_seeds": [2001]})
    build([FASTA], str(d / "b"), {"entrap_seeds": [2001]}, shuffle_rows_seed=3, swap_labels=True)
    return str(d / "a"), str(d / "b")


def test_closed_unchanged_and_secondary_found(built):
    import pyarrow as pa
    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    l1 = L1Index(built[0])
    p = E.SearchParams(min_peaks=5, report_psms=5)
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    sp, groups, A, B = _chimeric_spectra(l1)
    closed, _, _ = E.search(sp, gi, l1, p)
    tbl, st = CS.search_chimeric(sp, gi, l1, p, groups)
    g = tbl.column("group").to_numpy()
    prim = tbl.filter(pa.array(g == 0)).select(closed.column_names)
    assert prim.equals(closed)                                # chimeric mode off == group 0 == closed, bit-identical
    t1 = tbl.filter(pa.array((g == 1) & (tbl.column("rank").to_numpy() == 1)))
    top_b = dict(zip(t1.column("spec_idx").to_pylist(), t1.column("pf_row").to_pylist()))
    assert np.mean([top_b.get(i) == b for i, b in enumerate(B)]) > 0.85      # each precursor: its own candidate list
    c1 = closed.filter(pa.array(closed.column("rank").to_numpy() == 1))
    top_a = dict(zip(c1.column("spec_idx").to_pylist(), c1.column("pf_row").to_pylist()))
    assert np.mean([top_a.get(i) == a for i, a in enumerate(A)]) > 0.85
    assert set(np.unique(tbl.column("chim_rank").to_numpy())) == {1, 2}
    # no secondary groups -> exactly the closed table (+ columns)
    empty = {k: v[:0] for k, v in groups["sec"].items()}
    t0, _ = CS.search_chimeric(sp, gi, l1, p, {**groups, "sec": empty, "n_groups": np.ones(len(sp), np.uint16)})
    assert t0.select(closed.column_names).equals(closed)


def test_label_swap_invariance(built):
    """G0 for the chimeric search: role-swapped, row-shuffled index -> same peptidoforms and scores."""
    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    outs = []
    for d in built:
        l1 = L1Index(d)
        p = E.SearchParams(min_peaks=5, report_psms=5)
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        if not outs:
            sp, groups, _, _ = _chimeric_spectra(l1)
        tbl, _ = CS.search_chimeric(sp, gi, l1, p, groups)
        r1 = tbl.column("rank").to_numpy() == 1
        key = sorted(zip(tbl.column("spec_idx").to_numpy()[r1].tolist(), tbl.column("group").to_numpy()[r1].tolist(),
                         tbl.column("pf_id").to_numpy()[r1].tolist(), tbl.column("hyperscore").to_numpy()[r1].tolist()))
        outs.append(key)
    assert outs[0] == outs[1]


def test_neighbour_ion_convention():
    from vulcan.chimeric.compete import IonSets
    from vulcan.ptm.conga import neighbour_filter

    class V:   # two tryptic peptides ending in K with no shared ion index >= 3
        residues = np.frombuffer(b"PEPTIDEKAGSWYNMK", np.uint8)
        seq_off = np.array([0, 8])
        length = np.array([8, 8])
        mask = np.array([0, 0], np.uint64)
    ions = IonSets(V())
    a, b = ions.get(0), ions.get(1)
    assert len(a) == 2 * (8 - 3)                              # b3..b7, y3..y7
    keep = neighbour_filter(np.array([0, 0]), np.array([2.0, 1.0]), [a, b], 0.05)
    assert keep.tolist() == [True, True]                      # shared y1/y2-type ions do not make neighbours
    keep = neighbour_filter(np.array([0, 0]), np.array([2.0, 1.0]), [a, a], 0.05)
    assert keep.tolist() == [True, False]
