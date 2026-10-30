"""Family top-2 (vulcan/gpu/family.py) on a small multi-replicate index (skipped without CUDA):
every (decoy seed x entrapment seed) replicate reconstructed from one search == a separate search of the replicate."""
import importlib.util
import os
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
torch = pytest.importorskip("torch")
pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")

FASTA = os.path.join(ROOT, "data", "fasta", "contaminants_frankenfield2022.fasta")
_spec = importlib.util.spec_from_file_location("k12", os.path.join(os.path.dirname(__file__), "test_k1_k2_small.py"))
_k12 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_k12)


@pytest.fixture(scope="module")
def multi(tmp_path_factory):
    from vulcan.index.build import L1Index, build_multi_decoy, read_labels
    d = str(tmp_path_factory.mktemp("idx") / "m")
    build_multi_decoy([FASTA], d, {"entrap_seeds": [2001, 2002], "decoy_seeds": [1001, 1002]})
    return L1Index(d), read_labels(d)


def test_multi_decoy_index_structure(multi):
    l1, lab = multi
    assert len(l1.partners) == 2
    side = lab.column("decoy_set").to_numpy()
    n = np.arange(len(l1))
    for j, pj in enumerate(l1.partners):
        fwd = side < 0
        assert (side[pj[fwd]] == j).all()                       # forward -> its set-j decoy
        assert (pj[pj[fwd]] == n[fwd]).all()                    # symmetric
        other = (side >= 0) & (side != j)
        assert (pj[other] == n[other]).all()                    # other sets: self
    assert len(np.unique(l1.pf_id)) == len(l1)
    fam = lab.column("family").to_numpy()
    assert np.bincount(fam, minlength=4)[2] == 2 * np.bincount(fam, minlength=4)[0]


def test_family_top2_equals_separate_searches(multi):
    from vulcan.gpu import engine as E
    from vulcan.gpu import family as F
    l1, _ = multi
    p = E.SearchParams(min_peaks=5)
    gi = E.GpuIndex(l1, p, {"C": 57.021464}, {"M": 15.9949})
    fi = F.FamilyIndex.from_l1(l1, gi)
    sp, _ = _k12._spectra(l1, n=400, seed=7)
    tbl, _ = F.search_family(sp, gi, l1, fi, p)
    ff, side = fi.fwd_fam.cpu().numpy(), fi.side.cpu().numpy()
    for j in range(2):
        for si in range(2):
            keep = np.isin(ff, [0, 1 + si]) & np.isin(side, [-1, j])
            sub = F.SubIndex(l1, keep, j)
            gs = E.GpuIndex(sub, p, {"C": 57.021464}, {"M": 15.9949})
            fs = F.FamilyIndex(ff[keep], np.where(side[keep] >= 0, 0, -1), [sub.partner], fi.F, gs.dev)
            ts, _ = F.search_family(sp, gs, sub, fs, p)
            a = F.reconstruct(ts, sorted(set(ts.column("fgroup").to_numpy().tolist())))
            b = F.reconstruct(tbl, fi.replicate_groups(j, si))
            assert len(a["spec_idx"]) > 200
            for k in ("spec_idx", "rank", "pf_id"):
                assert np.array_equal(a[k], b[k]), (j, si, k)
            assert np.array_equal(a["hyperscore"].view(np.uint64), b["hyperscore"].view(np.uint64))
