"""Integration tests on real vendor files (skipped when data/tools are absent).

Thermo: direct (RawFileReader in-process) on a scan range vs TRFP-mzML + pyteomics on the same scans.
Bruker: vulcan_tims native frames == timsrust_pyo3 frames on the first frames; invariants.
Set VULCAN_SLOW=1 to also run the whole-file direct vs TRFP-mzML store comparison.
"""
import glob
import os
import subprocess

import numpy as np
import pytest

from vulcan.io import check
from vulcan.io import store as St

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RAW = os.path.join(ROOT, "data/raw/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw")
TRFP = os.path.join(ROOT, "conda_envs/trfp/bin/ThermoRawFileParser")
DLL = os.path.join(ROOT, "tools/vulcanraw/VulcanRaw.dll")
BRUKER_D = sorted(glob.glob("/tmp/vulcan_stage/raw/*.d")) + sorted(glob.glob(os.path.join(ROOT, "data/cache/*.d")))

need_thermo = pytest.mark.skipif(not (os.path.exists(RAW) and os.path.exists(DLL) and os.path.exists(TRFP)),
                                 reason="Thermo RAW / VulcanRaw.dll / TRFP not available")
need_bruker = pytest.mark.skipif(not BRUKER_D, reason="no unzipped Bruker .d in /tmp/vulcan_stage/raw or data/cache")


@pytest.fixture(scope="module")
def trfp_mzml(tmp_path_factory):
    d = tmp_path_factory.mktemp("trfp")
    subprocess.run([TRFP, "-i", RAW, "-o", str(d), "-f", "2", "-l", "4"], check=True)
    return str(d / (os.path.basename(RAW)[:-4] + ".mzML"))


@need_thermo
def test_thermo_direct_scan_range_vs_pyteomics(tmp_path, trfp_mzml):
    from vulcan.io import thermo_direct
    out = str(tmp_path / "d.vspec.arrow")
    meta = thermo_direct.convert(RAW, out, first_scan=2200, last_scan=2700, batch_scans=128)
    assert meta["n_spectra"] == 501 and meta["n_ms2"] > 0
    t = St.read_table(out)
    assert t.column("scan").to_numpy()[0] == 2200
    res = check.check_against_pyteomics(out, trfp_mzml, n=60)
    assert res["n_bad"] == 0, res
    lvl = t.column("ms_level").to_numpy()
    assert t.column("isolation_target_mz").null_count == int((lvl == 1).sum())
    assert t.column("filter_string").null_count == 0


@need_thermo
@pytest.mark.skipif(os.environ.get("VULCAN_SLOW") != "1", reason="set VULCAN_SLOW=1")
def test_thermo_direct_equals_mzml_path(tmp_path, trfp_mzml):
    from vulcan.io import mzml, thermo_direct
    a, b = str(tmp_path / "a.vspec.arrow"), str(tmp_path / "b.vspec.arrow")
    thermo_direct.convert(RAW, a)
    mzml.convert(trfp_mzml, b)
    res = check.compare_stores(a, b)
    assert res["ok"], res


@need_bruker
def test_bruker_native_equals_pyo3(tmp_path):
    from vulcan.io import bruker
    if bruker._native() is None:
        pytest.skip("vulcan_tims not built (scripts/ingest/build_vulcan_tims.sh)")
    d = BRUKER_D[0]
    a, b = str(tmp_path / "n.vspec.arrow"), str(tmp_path / "p.vspec.arrow")
    ma = bruker.convert(d, a, mode="frames", max_frames=300)
    mb = bruker.convert(d, b, mode="frames_pyo3", max_frames=300)
    assert ma["n_spectra"] == mb["n_spectra"] and ma["n_peaks"] == mb["n_peaks"]
    assert ma["events_outside_windows"] == mb["events_outside_windows"]
    ta, tb = St.read_table(a), St.read_table(b)
    pa_, pb = St.peaks_numpy(ta, True), St.peaks_numpy(tb, True)
    assert np.array_equal(pa_.offsets, pb.offsets) and np.array_equal(pa_.mz, pb.mz)
    assert np.array_equal(pa_.intensity, pb.intensity) and np.array_equal(pa_.im_index, pb.im_index)
    for c in ("scan", "ms_level", "isolation_target_mz", "window_index", "im_scan_begin", "im_scan_end", "nce"):
        assert ta.column(c).equals(tb.column(c)), c
    # invariants: im_index within the row's scan range for window rows
    lvl = ta.column("ms_level").to_numpy()
    sb = St.column_numpy(ta, "im_scan_begin")
    se = St.column_numpy(ta, "im_scan_end")
    for i in np.nonzero(lvl == 2)[0][:50]:
        seg = pa_.im_index[pa_.offsets[i]:pa_.offsets[i + 1]]
        assert len(seg) == 0 or (seg.min() >= sb[i] and seg.max() < se[i])
    meta = St.run_meta(a)
    assert meta["acquisition"] in ("DIA-PASEF", "DDA-PASEF", "MS1") and len(meta["im_lut"]) > 0
