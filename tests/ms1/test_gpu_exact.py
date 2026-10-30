"""E11: the opt-in GPU MS1 backend (vulcan.ms1.gpu) reproduces the CPU feature finder bitwise.

Checks every Hills array, every group_isotopes array, the feature table (+ schema metadata), peak_feature/peak_iso and
the hypothesis table.  Skipped without a usable CUDA device."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from synth import make_run  # noqa: E402

from vulcan.ms1 import isotopes as IS  # noqa: E402
from vulcan.ms1.features import FeatureParams, find_features  # noqa: E402


def _gpu_ok():
    try:
        from vulcan.ms1 import gpu
        return gpu.available()
    except Exception:  # noqa: BLE001
        return False


pytestmark = pytest.mark.skipif(not _gpu_ok(), reason="needs CUDA + nvcc (vulcan.ms1.gpu)")

HILL_FIELDS = ("off", "p_scan", "p_mz", "p_int", "p_src", "h_mz", "h_start", "h_end", "h_apex_scan", "h_apex_int",
               "h_area", "h_len")


def _same(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def assert_equal_results(fc, fg):
    for k in HILL_FIELDS:
        assert _same(getattr(fc.hills, k), getattr(fg.hills, k)), k
    assert fc.features.equals(fg.features)
    assert fc.features.schema.metadata == fg.features.schema.metadata
    assert _same(fc.peak_feature, fg.peak_feature)
    assert _same(fc.peak_iso, fg.peak_iso)
    assert fc.params == fg.params
    assert fg.timings["backend"] == "gpu"


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_gpu_backend_bitwise(seed):
    ms1, _ = make_run(n_feat=400, n_scans=250, seed=seed, noise_per_scan=600)
    fc = find_features(ms1, FeatureParams())
    fg = find_features(ms1, FeatureParams(), backend="gpu")
    assert fc.features.num_rows > 100
    assert_equal_results(fc, fg)


def test_gpu_hill_peaks_off_same_tables():
    """Production mode (hill peak arrays left on the device): identical tables and peak maps."""
    ms1, _ = make_run(n_feat=300, n_scans=200, seed=8, noise_per_scan=500)
    fc = find_features(ms1, FeatureParams())
    fg = find_features(ms1, FeatureParams(), backend="gpu", hill_peaks=False)
    assert fg.hills.p_src is None and fg.hills.p_int is None
    for k in ("off", "h_mz", "h_start", "h_end", "h_apex_scan", "h_apex_int", "h_area", "h_len"):
        assert _same(getattr(fc.hills, k), getattr(fg.hills, k)), k
    assert fc.features.equals(fg.features) and fc.features.schema.metadata == fg.features.schema.metadata
    assert _same(fc.peak_feature, fg.peak_feature) and _same(fc.peak_iso, fg.peak_iso)


def test_cpu_default_untouched():
    """backend defaults to cpu; an unknown backend is an error (the frozen path is never switched implicitly)."""
    import inspect
    assert inspect.signature(find_features).parameters["backend"].default == "cpu"
    ms1, _ = make_run(n_feat=50, n_scans=80, seed=1, noise_per_scan=100)
    assert "backend" not in find_features(ms1, FeatureParams()).timings
    with pytest.raises(ValueError):
        find_features(ms1, FeatureParams(), backend="cuda")


def test_gpu_gated_early_release_same_tables(monkeypatch):
    """Cohort mode (cross-process gate set): the device is released after the isotope arrays and the greedy pass +
    feature table (CPU peak map) run on the host -- identical output, hill_peaks on or off."""
    import threading

    from vulcan.ms1 import gpu
    monkeypatch.setattr(gpu, "_gate", threading.Semaphore(1))
    ms1, _ = make_run(n_feat=300, n_scans=200, seed=12, noise_per_scan=500)
    fc = find_features(ms1, FeatureParams())
    for hp in (True, False):
        fg = find_features(ms1, FeatureParams(), backend="gpu", hill_peaks=hp)
        assert fg.timings["backend"] == "gpu" and "gpu_held" in fg.timings
        assert fc.features.equals(fg.features) and fc.features.schema.metadata == fg.features.schema.metadata
        assert _same(fc.peak_feature, fg.peak_feature) and _same(fc.peak_iso, fg.peak_iso)
        if hp:
            for k in HILL_FIELDS:
                assert _same(getattr(fc.hills, k), getattr(fg.hills, k)), k
        else:
            assert fg.hills.p_src is None


def test_gpu_memory_fallback_same_tables(monkeypatch):
    """Too little free device memory -> the run takes the CPU path (recorded), with identical output."""
    from vulcan.ms1 import gpu
    monkeypatch.setattr(gpu.Session, "bytes_needed", staticmethod(lambda n: 1 << 62))
    ms1, _ = make_run(n_feat=100, n_scans=120, seed=9, noise_per_scan=200)
    fc = find_features(ms1, FeatureParams())
    fg = find_features(ms1, FeatureParams(), backend="gpu")
    assert fg.timings["backend"] == "cpu_fallback" and "gpu_mem" in fg.timings["fallback_reason"]
    assert fc.features.equals(fg.features) and fc.features.schema.metadata == fg.features.schema.metadata


def test_gpu_isotope_arrays_bitwise():
    """The group_isotopes dict itself (incl. n_candidates and the -1/0 padding)."""
    from vulcan.ms1 import gpu
    from vulcan.ms1.hills import detect_hills
    ms1, _ = make_run(n_feat=300, n_scans=200, seed=7, noise_per_scan=500)
    p = FeatureParams()
    lens = np.diff(ms1.off)
    scan = np.repeat(np.arange(len(lens), dtype=np.int32), lens)
    hc = detect_hills(scan, ms1.mz, ms1.intensity, ms1.ms1_rt, htol_ppm=p.htol_ppm, minlh=p.minlh, hvf=p.hvf)
    gc = IS.group_isotopes(hc)
    with gpu.session_lock:
        S = gpu.session()
        hg = S.hills(scan, ms1.mz, ms1.intensity, ms1.ms1_rt, None, htol_ppm=p.htol_ppm, minlh=p.minlh, hvf=p.hvf)
        gg = S.isotopes(hg, charges=p.charges, itol_ppm=p.itol_ppm, max_iso=p.max_iso, min_trace_cos=p.min_trace_cos,
                        min_overlap_for_cos=p.min_overlap_for_cos, min_avg_cos=p.min_avg_cos,
                        min_explained=p.min_explained, use_charge_ban=True, scan_bin=IS.SCAN_BIN)
        S.release()
    for k in HILL_FIELDS:
        assert _same(getattr(hc, k), getattr(hg, k)), k
    assert gc["n_candidates"] == gg["n_candidates"] and gg["fast_path"]
    for k, v in gc.items():
        if isinstance(v, np.ndarray):
            assert _same(v, gg[k]), k


def test_gpu_small_scan_bin(monkeypatch):
    """Hills spanning many scan blocks exercise the first-common-block de-duplication on the GPU too."""
    monkeypatch.setattr(IS, "SCAN_BIN", 2)
    ms1, _ = make_run(n_feat=300, n_scans=200, seed=5, noise_per_scan=400)
    assert_equal_results(find_features(ms1, FeatureParams()), find_features(ms1, FeatureParams(), backend="gpu"))


def test_gpu_param_variants():
    """mz window (peak subset -> p_src remap), min_intensity, no charge ban, shorter chains, other tolerances."""
    ms1, _ = make_run(n_feat=300, n_scans=200, seed=11, noise_per_scan=500)
    for kw in (dict(mz_min=500.0, mz_max=1200.0), dict(min_intensity=2e4), dict(charge_ban=False),
               dict(max_iso=4, min_avg_cos=0.7), dict(htol_ppm=5.0, itol_ppm=12.0, minlh=3, hvf=1.5),
               dict(min_overlap_for_cos=3), dict(charges=(2, 3, 4))):
        p = FeatureParams(**kw)
        assert_equal_results(find_features(ms1, p), find_features(ms1, p, backend="gpu"))


def test_gpu_unsupported_falls_back_to_cpu():
    ms1, _ = make_run(n_feat=100, n_scans=120, seed=3, noise_per_scan=200)
    p = FeatureParams(gap=1)
    fc = find_features(ms1, p)
    fg = find_features(ms1, p, backend="gpu")
    assert "backend" not in fg.timings and fc.features.equals(fg.features)


def test_gpu_empty_run():
    from vulcan.ms1.features import MS1Data
    ms1 = MS1Data(1, np.arange(3), np.arange(3), np.arange(3) / 60.0, np.empty(0, np.float32),
                  np.empty(0, np.float32), np.zeros(4, np.int64))
    fc = find_features(ms1, FeatureParams())
    fg = find_features(ms1, FeatureParams(), backend="gpu")
    assert fg.features.num_rows == 0 and fc.features.equals(fg.features)


def test_gpu_hypotheses_equal(tmp_path):
    """Hypothesis table built from the GPU feature result equals the CPU one (needs a store; synthetic MS1 only)."""
    ms1, _ = make_run(n_feat=200, n_scans=150, seed=4, noise_per_scan=300)
    fc = find_features(ms1, FeatureParams())
    fg = find_features(ms1, FeatureParams(), backend="gpu")
    from vulcan.ms1.hypotheses import window_features
    s_ord = np.arange(0, len(ms1.off) - 1, 3, dtype=np.int64)
    lo = np.full(len(s_ord), 600.0)
    hi = np.full(len(s_ord), 602.0)
    a = window_features(ms1, fc.peak_feature, s_ord, lo, hi)
    b = window_features(ms1, fg.peak_feature, s_ord, lo, hi)
    for x, y in zip(a, b):
        assert _same(np.asarray(x), np.asarray(y))
