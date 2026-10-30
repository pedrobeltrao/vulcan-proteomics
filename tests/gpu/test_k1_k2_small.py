"""Small end-to-end GPU checks (skipped without CUDA): K1 counts == f32 brute force; G0 role swap invariance."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
torch = pytest.importorskip("torch")
pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")

FASTA = os.path.join(ROOT, "data", "fasta", "contaminants_frankenfield2022.fasta")


def _spectra(l1, n=300, seed=5):
    """Synthetic spectra from random index peptidoforms (b/y + noise), so every code path gets real matches."""
    from vulcan.gpu.engine import Spectra
    from vulcan.index._core import core
    rng = np.random.default_rng(seed)
    rows = rng.choice(len(l1), n, replace=False)
    mz, it, off, pmz = [], [], [0], []
    for r in rows:
        seq, mask = l1.sequence(r), int(l1.mask[r])
        ref = [x for x in core.sage_reference([seq], [("C", 57.021464)], [("M", 15.9949)], 2) if x[2] == mask][0]
        f = np.array(ref[4] + ref[5], np.float32) + np.float32(1.0072764)
        f = f[rng.random(len(f)) < 0.7]
        noise = rng.uniform(150, 2000, 60).astype(np.float32)
        m = np.sort(np.concatenate([f, noise]))
        mz.append(m)
        it.append(rng.uniform(1, 100, len(m)).astype(np.float32))
        off.append(off[-1] + len(m))
        pmz.append((float(l1.mass32[r]) + 2 * 1.0072764) / 2)
    n = len(rows)
    return Spectra("synthetic", np.arange(n, dtype=np.uint32), np.arange(n, dtype=np.uint32), np.zeros(n),
                   np.array(pmz), np.full(n, 2, np.int8), np.concatenate(mz), np.concatenate(it),
                   np.array(off, np.int64)), rows


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    from vulcan.index.build import build
    d = tmp_path_factory.mktemp("idx")
    build([FASTA], str(d / "a"), {"entrap_seeds": [2001]})
    build([FASTA], str(d / "b"), {"entrap_seeds": [2001]}, shuffle_rows_seed=3, swap_labels=True)
    return str(d / "a"), str(d / "b")


def test_k1_equals_bruteforce_and_finds_truth(built):
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    l1 = L1Index(built[0])
    p = E.SearchParams(min_peaks=5)
    gi = E.GpuIndex(l1, p, {"C": 57.021464}, {"M": 15.9949})
    sp, rows = _spectra(l1)
    st = E.Stats()
    pk = E.k0_preprocess(sp, p, gi.dev, st)
    hyp = E.make_hypotheses(sp, pk["valid"], gi.mass32_np, p)
    o_h, o_r, o_c, *_ = E.run_k1(gi, hyp, pk, p, st)
    h, r = E.window_candidates(hyp, gi.dev)
    bc = E.k1_brute(gi, hyp, pk, h, r, p).long()
    a = sorted(zip(((o_h.long() << 32) | o_r.long()).tolist(), o_c.tolist()))
    nz = bc > 0
    b = sorted(zip(((h[nz] << 32) | r[nz]).tolist(), bc[nz].tolist()))
    assert a == b
    tbl, _, _ = E.search(sp, gi, l1, p)
    top = dict(zip(tbl.column("spec_idx").to_pylist(), tbl.column("pf_row").to_pylist()))
    hit = np.mean([top.get(i) == r0 for i, r0 in enumerate(rows)])
    assert hit > 0.9


def test_g0_role_swap_and_permutation(built):
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    a, b = L1Index(built[0]), L1Index(built[1])
    for k in ("pf_id", "mass32", "length", "mask", "partner", "residues"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k
    p = E.SearchParams(min_peaks=5)
    sp, _ = _spectra(a)
    outs = []
    for l1 in (a, b):
        gi = E.GpuIndex(l1, p, {"C": 57.021464}, {"M": 15.9949})
        outs.append(E.search(sp, gi, l1, p, report=5)[0])
    assert outs[0].equals(outs[1])
