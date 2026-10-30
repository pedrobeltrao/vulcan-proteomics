"""E7b: the fast MS1 paths (fused hill passes, RT-binned isotope lookups) are bitwise equal to the reference kernels."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from synth import make_run  # noqa: E402

from vulcan.ms1 import hills as H  # noqa: E402
from vulcan.ms1 import isotopes as IS  # noqa: E402
from vulcan.ms1.features import FeatureParams, find_features  # noqa: E402

HILL_FIELDS = ("off", "p_scan", "p_mz", "p_int", "p_src", "h_mz", "h_start", "h_end", "h_apex_scan", "h_apex_int",
               "h_area", "h_len")


def _hills(ms1, reference):
    p = FeatureParams()
    lens = np.diff(ms1.off)
    scan = np.repeat(np.arange(len(lens), dtype=np.int32), lens)
    return H.detect_hills(scan, ms1.mz, ms1.intensity, ms1.ms1_rt, htol_ppm=p.htol_ppm, minlh=p.minlh, hvf=p.hvf,
                          reference=reference)


def test_fast_paths_bitwise_equal():
    for seed in (0, 1, 2):
        ms1, _ = make_run(n_feat=400, n_scans=250, seed=seed, noise_per_scan=600)
        hf, hr = _hills(ms1, False), _hills(ms1, True)
        for k in HILL_FIELDS:
            a, b = getattr(hf, k), getattr(hr, k)
            assert a.dtype == b.dtype and a.tobytes() == b.tobytes(), (seed, k)
        gf = IS.group_isotopes(hf)
        gr = IS.group_isotopes(hr, exact_reference=True)
        assert gf["fast_path"] and not gr["fast_path"]
        assert gf["n_candidates"] == gr["n_candidates"]
        for k, v in gf.items():
            if isinstance(v, np.ndarray):
                assert v.dtype == gr[k].dtype and v.tobytes() == gr[k].tobytes(), (seed, k)
        tf = find_features(ms1, FeatureParams()).features
        tr = find_features(ms1, FeatureParams(), exact_reference=True).features
        assert tf.equals(tr) and tf.schema.metadata == tr.schema.metadata


def test_small_scan_bin_still_equal(monkeypatch):
    """Many hills spanning several scan blocks exercise the first-common-block de-duplication."""
    monkeypatch.setattr(IS, "SCAN_BIN", 2)
    ms1, _ = make_run(n_feat=300, n_scans=200, seed=5, noise_per_scan=400)
    h = _hills(ms1, False)
    gf = IS.group_isotopes(h)
    gr = IS.group_isotopes(h, exact_reference=True)
    for k, v in gf.items():
        if isinstance(v, np.ndarray):
            assert v.tobytes() == gr[k].tobytes(), k
