"""Fast production scoring components reproduce the frozen path bit for bit (results/2026-10-27_scoring_speed)."""
import os

import numpy as np
import pytest

torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("CUDA required", allow_module_level=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl")


def test_gpu_hgb_bit_identical_random_model():
    from sklearn.ensemble import HistGradientBoostingClassifier
    from vulcan.rescore.fast import GpuHGB
    rng = np.random.default_rng(0)
    X = rng.normal(size=(20000, 7))
    X[rng.random(X.shape) < 0.05] = np.nan
    y = (np.nan_to_num(X[:, 0]) + 0.5 * np.nan_to_num(X[:, 1]) ** 2 + rng.normal(size=len(X)) > 0.3).astype(int)
    m = HistGradientBoostingClassifier(max_iter=60, max_depth=5, random_state=0).fit(X, y)
    Xt = rng.normal(size=(30001, 7))
    Xt[rng.random(Xt.shape) < 0.05] = np.nan
    Xt[:50, 2] = m._predictors[0][0].nodes["num_threshold"][0]       # values exactly on a threshold
    assert np.array_equal(GpuHGB(m).decision_function(Xt, chunk=7000), m.decision_function(Xt))


@pytest.mark.skipif(not os.path.exists(MODEL), reason="frozen model not present")
def test_gpu_hgb_bit_identical_frozen_model():
    from vulcan.rescore.fast import GpuHGB
    from vulcan.rescore.static import StaticModel
    sm = StaticModel.load(MODEL)
    rng = np.random.default_rng(1)
    X = rng.normal(size=(50000, len(sm.names))) * 5
    thr = np.concatenate([t[0].nodes["num_threshold"] for t in sm.model._predictors])
    X[:2000] = rng.choice(thr, size=(2000, X.shape[1]))
    assert np.array_equal(GpuHGB(sm.model).decision_function(X), sm.model.decision_function(X))


def _spectra(rng, n=300):
    from vulcan.gpu.host import Spectra
    ln = rng.integers(0, 80, n)
    off = np.zeros(n + 1, np.int64)
    np.cumsum(ln, out=off[1:])
    mz = rng.uniform(100, 2000, off[-1]).astype(np.float32)
    mz[rng.random(len(mz)) < 0.01] = np.nan
    mz[5:9] = mz[4]                                     # duplicates within a spectrum
    it = rng.exponential(1e5, off[-1]).astype(np.float32)
    it[rng.random(len(it)) < 0.05] = 0
    z = np.zeros(n)
    return Spectra("x", np.arange(n), np.arange(n), z, z, z.astype(np.int8), mz, it, off)


def test_gpu_peaks_identical():
    from vulcan.rescore.fast import GpuPeaks
    from vulcan.rescore.features import SpecPeaks
    sp = _spectra(np.random.default_rng(2))
    a, b = SpecPeaks(sp, torch.device("cuda")), GpuPeaks(sp, torch.device("cuda"))
    for k in ("key", "it", "tic", "mzmin", "mzmax"):
        x, y = getattr(a, k).cpu().numpy(), getattr(b, k).cpu().numpy()
        assert x.dtype == y.dtype and np.array_equal(x, y), k


class _FakePred:
    """Deterministic stand-in predictor; `salt` changes values (simulates batch-composition differences)."""
    version = {"model_zip_sha256": "fake"}

    def __init__(self, salt=0):
        self.salt = salt

    def predict(self, residues, seq_off, length, mask, charge, nce, instrument="QE", want_ms2=True, want_rt=True,
                batch=16384):
        n = len(length)
        fo = np.zeros(n + 1, np.int64)
        np.cumsum((np.asarray(length) - 1) * 4, out=fo[1:])
        rng = np.random.default_rng(int(seq_off.sum()) + self.salt)
        return {"rt": rng.random(n).astype(np.float32), "ms2": rng.random(fo[-1]).astype(np.float16), "ms2_off": fo}


class _L1:
    def __init__(self, n, rng):
        self.length = rng.integers(7, 30, n).astype(np.int32)
        self.seq_off = np.r_[0, np.cumsum(self.length)[:-1]].astype(np.int64)
        self.residues = rng.integers(65, 90, int(self.length.sum())).astype(np.uint8)
        self.mask = np.zeros(n, np.int64).view(np.uint64)
        self.pf_id = rng.permutation(np.arange(10 * n, dtype=np.uint64))[:n]

    def __len__(self):
        return len(self.length)


def test_libcache_matches_predcache(tmp_path):
    from vulcan.predict.cache import PredCache
    from vulcan.predict.library import LibCache, consolidate
    rng = np.random.default_rng(3)
    l1 = _L1(3000, rng)
    root = str(tmp_path)
    # three writers with different values for overlapping rows (first copy must win everywhere)
    for salt, rows in ((1, np.arange(0, 1500)), (2, np.arange(1000, 2500)), (3, np.arange(200, 1200))):
        c = PredCache(root, _FakePred(salt), 26.0, "QE", {"C": 57.021464}, {"M": 15.9949})
        c.lookup_rt(l1, rows)
        for z in (2, 3):
            c.lookup_ms2(l1, rows, np.full(len(rows), z))
    key = c.key
    consolidate(os.path.join(root, key))
    # a shard written after consolidation
    c = PredCache(root, _FakePred(4), 26.0, "QE", {"C": 57.021464}, {"M": 15.9949})
    c.lookup_ms2(l1, np.arange(2400, 2700), np.full(300, 2))
    q = rng.choice(2500, 4000)
    zq = rng.choice([2, 3], len(q))
    old = PredCache(root, _FakePred(9), 26.0, "QE", {"C": 57.021464}, {"M": 15.9949})
    new = LibCache(root, _FakePred(9), 26.0, "QE", {"C": 57.021464}, {"M": 15.9949})
    assert np.array_equal(old.lookup_rt(l1, q), new.lookup_rt(l1, q))
    o1, f1 = old.lookup_ms2(l1, q, zq)
    o2, f2 = new.lookup_ms2(l1, q, zq)
    assert np.array_equal(o1, o2) and np.array_equal(f1.view(np.uint16), f2.view(np.uint16))
    # misses (incl. a charge never seen) are predicted, appended and then served as hits
    q2 = np.arange(2600, 2900)
    o3, f3 = new.lookup_ms2(l1, q2, np.r_[np.full(150, 2), np.full(150, 4)])
    o4, f4 = new.lookup_ms2(l1, q2, np.r_[np.full(150, 2), np.full(150, 4)])
    assert np.array_equal(f3.view(np.uint16), f4.view(np.uint16)) and new.stats["ms2_miss"] == 200   # z2: rows 2700-2749; z4: 150


def test_libcache_device_mirror_identical(tmp_path):
    from vulcan.predict.cache import PredCache
    from vulcan.predict.library import LibCache, consolidate
    rng = np.random.default_rng(5)
    l1 = _L1(4000, rng)
    l1.pf_id[:5] = np.array([2 ** 64 - 1, 2 ** 63, 2 ** 63 - 1, 0, 1], np.uint64)   # sign-bit edge cases
    root = str(tmp_path)
    for salt, rows in ((1, np.arange(0, 2500)), (2, np.arange(2000, 3500))):
        c = PredCache(root, _FakePred(salt), 29.0, "QE", {"C": 57.021464}, {"M": 15.9949})
        c.lookup_rt(l1, rows)
        for z in (2, 3):
            c.lookup_ms2(l1, rows, np.full(len(rows), z))
    consolidate(os.path.join(root, c.key))
    host = LibCache(root, _FakePred(9), 29.0, "QE", {"C": 57.021464}, {"M": 15.9949}, resident=True)
    gpu = LibCache(root, _FakePred(9), 29.0, "QE", {"C": 57.021464}, {"M": 15.9949}, device="cuda")
    q = np.r_[np.arange(5), rng.choice(3500, 5000)]
    zq = rng.choice([2, 3], len(q))
    assert np.array_equal(host.lookup_rt(l1, q), gpu.lookup_rt_dev(l1, q).cpu().numpy())
    o1, f1 = host.lookup_ms2(l1, q, zq)
    o2, f2 = gpu.lookup_ms2_dev(l1, q, zq)
    assert np.array_equal(o1, o2) and np.array_equal(f1.view(np.uint16), f2.cpu().numpy().view(np.uint16))
    # misses (rows never predicted, a new charge) go through the host path and are then served from the device
    q2 = np.r_[np.arange(3400, 3700), np.arange(10, 20)]
    z2 = np.r_[np.full(300, 2), np.full(10, 5)]
    o3, f3 = gpu.lookup_ms2_dev(l1, q2, z2)
    o4, f4 = gpu.lookup_ms2_dev(l1, q2, z2)
    assert np.array_equal(f3.cpu().numpy().view(np.uint16), f4.cpu().numpy().view(np.uint16))
    assert gpu.stats["ms2_miss"] == 210
