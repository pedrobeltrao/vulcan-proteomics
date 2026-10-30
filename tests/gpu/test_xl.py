"""M5 XL search (vulcan.xl): G0 label blindness + recovery of synthetic cross-links.

* G0 (DESIGN §3): role swap of every target/decoy pair + permutation of the pair list leave the index side and the
  XL CSM tables (cleavable DSSO path, non-cleavable DSS alpha-first path, incl. the linear/mono/loop competitor
  columns, homo-family flag and link sites) bit-identical.
* Synthetic spectra built from two peptides (b/y with the partner + linker shift on the ions that carry the linked
  residue; DSSO adds stub ions and the alkene/thiol doublets) are recovered as the rank-1 pair.
* Doublet detection finds the planted doublet masses.
"""
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
VAR = {"M": 15.9949}
PARAMS = {"enzyme": {"cleave_at": "KR", "restrict": "P", "c_terminal": True, "missed_cleavages": 3,
                     "min_len": 5, "max_len": 40}, "peptide_min_mass": 400.0, "peptide_max_mass": 6000.0,
          "max_variable_mods": 1}
PROTON = 1.0072764


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    from vulcan.index.build import build
    d = tmp_path_factory.mktemp("xlidx")
    build([FASTA], str(d / "a"), PARAMS)
    build([FASTA], str(d / "b"), PARAMS, shuffle_rows_seed=11, swap_labels=True)
    return str(d / "a"), str(d / "b")


def _ions(l1, r, site, shift, max_c=2):
    from vulcan.xl.sites import residue_masses
    seq = l1.sequence(int(r))
    rm = residue_masses(seq, int(l1.mask[r]), STAT, VAR)
    cb = np.cumsum(rm)[:-1]
    cy = rm.sum() + 18.010565 - cb
    k = np.arange(len(seq) - 1)
    b = np.where(k >= site, cb + shift, cb)
    y = np.where(site > k, cy + shift, cy)
    out = []
    for c in range(1, max_c + 1):
        out.append((np.concatenate([b, y]) + c * PROTON) / c)
    return np.concatenate(out)


def _spectra(l1, linker, n=150, seed=3, keep_doublets=False):
    """Synthetic z=4 cross-link spectra of random K-containing peptide pairs (forward or decoy rows alike). 80% of
    fragment peaks are kept at random; keep_doublets=True keeps every precursor-level doublet peak (else 80% too)."""
    from vulcan.gpu.engine import Spectra
    rng = np.random.default_rng(seed)
    cand = np.flatnonzero([("K" in l1.sequence(i)[:-1]) and 7 <= l1.length[i] <= 20 and l1.n_var[i] == 0
                           for i in range(len(l1))])
    mz, it, off, pmz, truth = [], [], [0], [], []
    for _ in range(n):
        a, b = rng.choice(cand, 2, replace=False)
        sa = l1.sequence(int(a))[:-1].index("K")
        sb = l1.sequence(int(b))[:-1].index("K")
        ma, mb = float(l1.mass32[a]), float(l1.mass32[b])
        M = ma + mb + linker.mass
        peaks = [_ions(l1, a, sa, mb + linker.mass), _ions(l1, b, sb, ma + linker.mass)]
        dbl = []
        if linker.cleavable:
            for r, s_, m in ((a, sa, ma), (b, sb, mb)):
                for stub in linker.doublet:
                    peaks.append(_ions(l1, r, s_, stub, max_c=1))
                    for c in (1, 2):
                        (dbl if keep_doublets else peaks).append(np.array([(m + stub + c * PROTON) / c,
                                                                           (m + stub + 1.00335 + c * PROTON) / c]))
        f = np.concatenate(peaks)
        f = f[(rng.random(len(f)) < 0.8) & (f > 150) & (f < 2000)]
        if dbl:
            g = np.concatenate(dbl)
            f = np.concatenate([f, g[(g > 150) & (g < 2000)]])
        m = np.sort(np.concatenate([f, rng.uniform(150, 2000, 40)])).astype(np.float32)
        mz.append(m)
        it.append(rng.uniform(10, 100, len(m)).astype(np.float32))
        off.append(off[-1] + len(m))
        pmz.append((M + 4 * PROTON) / 4)
        truth.append((min(l1.pf_id[a], l1.pf_id[b]), max(l1.pf_id[a], l1.pf_id[b])))
    sp = Spectra("synthetic", np.arange(n, dtype=np.uint32), np.arange(n, dtype=np.uint32), np.zeros(n),
                 np.array(pmz), np.full(n, 4, np.int8), np.concatenate(mz), np.concatenate(it),
                 np.array(off, np.int64))
    return sp, truth


def _run(l1, sp, linker):
    from vulcan.gpu import engine as E
    from vulcan.xl import search as XS
    from vulcan.xl import sites as S
    p = E.SearchParams(min_peaks=5, max_fragment_charge=3, max_peaks=200)
    gi = E.GpuIndex(l1, p, STAT, VAR)
    xp = XS.XLParams(isotope_errors=(0, 1))
    tbl, lin, st = XS.search(sp, gi, l1, p, xp, linker)
    sa, sb, na, nb = S.localize_csms(tbl, sp, l1, linker, STAT, VAR)
    tbl = tbl.append_column("site_a", pa.array(sa)).append_column("site_b", pa.array(sb))
    return tbl.drop_columns(["pf_row_a", "pf_row_b", "competitor_pf_row"]), tbl


@pytest.mark.parametrize("lk", ["DSSO", "DSS"])
def test_g0_and_recovery(built, lk):
    from vulcan.index.build import L1Index
    from vulcan.xl.linkers import LINKERS
    linker = LINKERS[lk]
    a, b = L1Index(built[0]), L1Index(built[1])
    for k in ("pf_id", "mass32", "length", "mask", "partner", "residues"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k
    sp, truth = _spectra(a, linker)
    ta, full = _run(a, sp, linker)
    tb, _ = _run(b, sp, linker)
    assert ta.equals(tb)
    # pf_row is index-specific but the index side is identical here, so rows agree too
    r1 = full.filter(pa.array(full.column("rank").to_numpy() == 1))
    got = {s: (min(x, y), max(x, y)) for s, x, y in zip(r1.column("spec_idx").to_pylist(),
                                                         r1.column("pf_id_a").to_pylist(),
                                                         r1.column("pf_id_b").to_pylist())}
    rec = np.mean([got.get(i) == t for i, t in enumerate(truth)])
    assert rec > 0.8, rec


def test_doublets_found(built):
    from vulcan.index.build import L1Index
    from vulcan.xl import doublets as DB
    from vulcan.xl.linkers import DSSO
    l1 = L1Index(built[0])
    sp, truth = _spectra(l1, DSSO, n=60, seed=5, keep_doublets=True)
    d = DB.detect(sp, DSSO)
    idx = {int(pf): i for i, pf in enumerate(l1.pf_id)}
    hit = 0
    for s, (x, y) in enumerate(truth):
        ms = d["mass"][d["spec"] == s]
        for pf in (x, y):
            m = float(l1.mass32[idx[int(pf)]])
            hit += bool(np.any(np.abs(ms - m) < 0.02))
    assert hit / (2 * len(truth)) > 0.9
