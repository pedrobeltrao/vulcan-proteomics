"""vulcan.cohort.load_ms2 == vulcan.gpu.engine.load_store (same MS2 arrays), on a synthetic multi-batch store."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
pytest.importorskip("torch")   # engine imports torch


def test_load_ms2_matches_load_store(tmp_path):
    from vulcan.cohort import load_ms2
    from vulcan.gpu.engine import load_store
    from vulcan.io.store import StoreWriter
    rng = np.random.default_rng(3)
    w = StoreWriter(str(tmp_path / "x.vspec.arrow"),
                    dict(run_name="x", vendor="test", instrument_model="m", source_path="s", source_size=0,
                         converter="test", converter_version="0", acquisition="DDA", peak_order="mz_ascending"))
    scan0 = 1
    for nb in (37, 4096, 5):                      # several record batches, incl. a tiny one
        lvl = rng.choice([1, 2, 2, 2, 3], nb).astype(np.uint8)
        lens = rng.integers(0, 40, nb)
        off = np.zeros(nb + 1, np.int64)
        np.cumsum(lens, out=off[1:])
        pmz = np.where(lvl >= 2, rng.uniform(300, 1500, nb), np.nan)
        z = np.where(lvl >= 2, rng.integers(0, 5, nb), 0).astype(np.int16)
        w.write({"scan": np.arange(scan0, scan0 + nb), "ms_level": lvl, "rt_min": rng.uniform(0, 90, nb),
                 "precursor_mz": pmz, "precursor_charge": z, "precursor_charge_mask": z <= 0,
                 "peak_offsets": off, "mz": np.sort(rng.uniform(100, 2000, off[-1])).astype(np.float32),
                 "intensity": rng.uniform(1, 1e6, off[-1]).astype(np.float32)})
        scan0 += nb
    w.close()
    a, b = load_store(str(tmp_path / "x.vspec.arrow")), load_ms2(str(tmp_path / "x.vspec.arrow"))
    assert a.run == b.run
    for k in ("scan", "index", "rt", "precursor_mz", "charge", "mz", "intensity", "offsets"):
        x, y = getattr(a, k), getattr(b, k)
        assert x.dtype == y.dtype and np.array_equal(x, y, equal_nan=x.dtype.kind == "f"), k
