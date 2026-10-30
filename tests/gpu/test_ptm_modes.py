"""M4: variable-mod (STY phospho, M-ox, peptide-N-term acetyl) index and the mass-offset / open modes.
* masses, masks and b/y of the N-term-acetyl / phospho peptidoforms equal Sage's (sage_reference);
* K1 == f32 brute force on the var-mod index (fragments start at the N-term delta);
* G0: role swap + pair-list permutation leave the index side and the closed, mass-offset and open PSM tables
  bit-identical (label blindness of the new modes)."""
import os
import sys

import numpy as np
import pyarrow as pa
import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
torch = pytest.importorskip("torch")
pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")

FASTA = os.path.join(ROOT, "data", "fasta", "contaminants_frankenfield2022.fasta")
STAT = {"C": 57.021464}
VAR = {"S": 79.966331, "T": 79.966331, "Y": 79.966331, "M": 15.994915, "^": 42.010565}
PARAMS = {"variable_mods": VAR, "max_variable_mods": 2, "entrap_seeds": [2001]}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    from vulcan.index.build import build
    d = tmp_path_factory.mktemp("ptmidx")
    build([FASTA], str(d / "a"), PARAMS)
    build([FASTA], str(d / "b"), PARAMS, shuffle_rows_seed=5, swap_labels=True)
    return str(d / "a"), str(d / "b")


def _ref(l1, r):
    from vulcan.index._core import core
    seq, mask = l1.sequence(r), int(l1.mask[r])
    refs = core.sage_reference([seq], sorted(STAT.items()), sorted(VAR.items()), 2)
    return [x for x in refs if x[2] == mask]


def _spectra(l1, n=300, seed=7, shift=0.0):
    from vulcan.gpu.engine import Spectra
    rng = np.random.default_rng(seed)
    rows = rng.choice(len(l1), n, replace=False)
    mz, it, off, pmz = [], [], [0], []
    for r in rows:
        ref = _ref(l1, r)[0]
        f = np.array(ref[4] + ref[5], np.float32) + np.float32(1.0072764)
        f = f[rng.random(len(f)) < 0.7]
        m = np.sort(np.concatenate([f, rng.uniform(150, 2000, 60).astype(np.float32)]))
        mz.append(m)
        it.append(rng.uniform(1, 100, len(m)).astype(np.float32))
        off.append(off[-1] + len(m))
        pmz.append((float(l1.mass32[r]) + shift + 2 * 1.0072764) / 2)
    n = len(rows)
    return Spectra("synthetic", np.arange(n, dtype=np.uint32), np.arange(n, dtype=np.uint32), np.zeros(n),
                   np.array(pmz), np.full(n, 2, np.int8), np.concatenate(mz), np.concatenate(it),
                   np.array(off, np.int64)), rows


def test_masses_and_masks_match_sage(built):
    from vulcan.index.build import L1Index
    l1 = L1Index(built[0])
    nt = (l1.mask.astype(np.uint64) >> np.uint64(63)) & np.uint64(1)
    assert nt.any() and (l1.n_var > 0).any()
    rng = np.random.default_rng(1)
    cand = np.concatenate([np.flatnonzero(nt == 1)[:200], rng.choice(len(l1), 200, replace=False)])
    for r in cand:
        ref = _ref(l1, int(r))
        assert len(ref) == 1, (l1.sequence(r), hex(int(l1.mask[r])))
        assert np.float32(ref[0][3]) == l1.mass32[r]


def test_k1_bruteforce_varmod(built):
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    l1 = L1Index(built[0])
    p = E.SearchParams(min_peaks=5)
    gi = E.GpuIndex(l1, p, STAT, VAR)
    sp, rows = _spectra(l1)
    st = E.Stats()
    pk = E.k0_preprocess(sp, p, gi.dev, st)
    hyp = E.make_hypotheses(sp, pk["valid"], gi.mass32_np, p)
    o_h, o_r, o_c, *_ = E.run_k1(gi, hyp, pk, p, st)
    h, r = E.window_candidates(hyp, gi.dev)
    bc = E.k1_brute(gi, hyp, pk, h, r, p).long()
    a = sorted(zip(((o_h.long() << 32) | o_r.long()).tolist(), o_c.tolist()))
    nz = bc > 0
    assert a == sorted(zip(((h[nz] << 32) | r[nz]).tolist(), bc[nz].tolist()))
    tbl, _, _ = E.search(sp, gi, l1, p)
    top = dict(zip(tbl.column("spec_idx").to_pylist(), tbl.column("pf_row").to_pylist()))
    assert np.mean([top.get(i) == r0 for i, r0 in enumerate(rows)]) > 0.85


@pytest.mark.parametrize("mode", ["closed", "offset", "open"])
def test_g0_modes(built, mode):
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.ptm import search as PS
    a, b = L1Index(built[0]), L1Index(built[1])
    for k in ("pf_id", "mass32", "length", "mask", "partner", "residues"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k
    p = E.SearchParams(min_peaks=5)
    mp = PS.ModeParams(offsets=(15.994915, 42.010565) if mode == "offset" else (),
                       open_da=(-150.0, 500.0) if mode == "open" else None, spectra_per_batch=64)
    sp, _ = _spectra(a, shift=15.994915 if mode != "closed" else 0.0)
    outs = []
    for l1 in (a, b):
        gi = E.GpuIndex(l1, p, STAT, VAR)
        t = PS.search_modes(sp, gi, l1, p, mp, report=5)[0]
        # hyp_offset is NaN for open hypotheses (NaN != NaN in Table.equals): compare its bit pattern instead
        ho = t.column("hyp_offset").to_numpy().view(np.uint32)
        outs.append(t.set_column(t.column_names.index("hyp_offset"), "hyp_offset", pa.array(ho)))
    assert outs[0].equals(outs[1])
    if mode != "closed":   # the shifted spectra are found through the offset / open hypotheses
        m = np.asarray(outs[0].column("mode").to_numpy())
        r1 = outs[0].column("rank").to_numpy() == 1
        assert (m[r1] > 0).mean() > 0.5
