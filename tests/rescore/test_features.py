"""L3 predicted-feature assembly: fragment masses, peak matching, similarity features, RT calibration (no peptdeep)."""
import numpy as np
import pytest

torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("needs CUDA", allow_module_level=True)

from vulcan.rescore import features as F  # noqa: E402


class _L1:
    def __init__(self, seqs, masks):
        self.residues = np.frombuffer("".join(seqs).encode(), np.uint8).copy()
        self.length = np.array([len(s) for s in seqs], np.uint8)
        self.seq_off = np.r_[0, np.cumsum(self.length)[:-1]].astype(np.int64)
        self.mask = np.array(masks, np.uint64)


def _ref_frag(seq, mask, static, var):
    m = [F.MONO[a] + static.get(a, 0) + (var.get(a, 0) if (mask >> i) & 1 else 0) for i, a in enumerate(seq)]
    out = []
    for f in range(len(seq) - 1):
        b = sum(m[:f + 1])
        y = sum(m[f + 1:]) + F.H2O
        out.append([b + F.PROTON, (b + 2 * F.PROTON) / 2, y + F.PROTON, (y + 2 * F.PROTON) / 2])
    return np.array(out)


def test_fragment_mz_matches_reference():
    seqs, masks = ["PEPTCIDEK", "MLMAAR"], [0, 0b100]
    l1 = _L1(seqs, masks)
    st, vr = {"C": 57.021464}, {"M": 15.9949}
    mz, valid = F.fragment_mz(l1, np.arange(2), torch.device("cuda"), st, vr, 9)
    for i, (s, m) in enumerate(zip(seqs, masks)):
        ref = _ref_frag(s, m, st, vr)
        got = mz[i, :len(s) - 1].cpu().numpy()
        assert np.allclose(got, ref, atol=1e-9)
        assert valid[i].sum().item() == len(s) - 1


def test_observed_and_features():
    from vulcan.gpu.engine import Spectra
    # one spectrum with peaks exactly at two of the fragments (+5 ppm) and one noise peak
    l1 = _L1(["PEPTIDEK"], [0])
    ref = _ref_frag("PEPTIDEK", 0, {}, {})
    peaks = np.sort(np.array([ref[2, 0] * (1 + 5e-6), ref[4, 2], 333.3], np.float32))
    sp = Spectra("r", np.array([1]), np.array([0]), np.array([1.0]), np.array([500.0]), np.array([2], np.int8),
                 peaks, np.array([10.0, 20.0, 5.0], np.float32)[np.argsort([ref[2, 0], ref[4, 2], 333.3])],
                 np.array([0, 3]))
    pk = F.SpecPeaks(sp, torch.device("cuda"))
    mz, valid = F.fragment_mz(l1, np.arange(1), torch.device("cuda"), {}, {}, 8)
    obs = F.observed(pk, torch.zeros(1, dtype=torch.long, device="cuda"), mz, 20.0)
    o = obs[0].cpu().numpy()
    assert (o > 0).sum() == 2 and o[2, 0] > 0 and o[4, 2] > 0
    # perfect prediction of the observed vector -> SA = 1, entropy similarity = 1
    p = obs.reshape(1, -1).clone()
    f = F.fragment_features(p, obs.reshape(1, -1), torch.ones_like(p), pk.tic[:1])
    assert abs(f["pred_sa"].item() - 1) < 1e-6 and abs(f["pred_entropy"].item() - 1) < 1e-5
    assert abs(f["pred_obs_tic"].item() - 30 / 35) < 1e-6


def test_rt_calibration_recovers_monotone_map():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1, 5000)
    y = 20 + 100 * x ** 1.3 + rng.normal(0, 1, x.size)
    y[:200] = rng.uniform(0, 150, 200)                  # outliers
    fit, mad = F.rt_calibrate(y, x, np.ones(x.size, bool))
    assert np.median(np.abs(fit[200:] - (20 + 100 * x[200:] ** 1.3))) < 0.5
    assert 0.7 < mad < 1.5
