"""Store schema/writer/reader invariants on synthetic data (fast, no external files)."""
import json
import os

import numpy as np
import pyarrow as pa
import pytest

from vulcan.io import schema as S
from vulcan.io import store as St


def _synthetic(n=10, seed=0):
    rng = np.random.default_rng(seed)
    lvl = np.array([1, 2, 2, 2] * ((n + 3) // 4))[:n].astype(np.uint8)
    lens = rng.integers(0, 50, n)
    off = np.zeros(n + 1, np.int64)
    np.cumsum(lens, out=off[1:])
    mz = np.sort(rng.uniform(100, 2000, off[-1])).astype(np.float32)
    it = rng.uniform(1, 1e6, off[-1]).astype(np.float32)
    ms2 = lvl == 2
    prec = np.where(ms2, rng.uniform(400, 1200, n), np.nan)
    ch = np.where(ms2, 2, 0).astype(np.int16)
    filt = [f"FTMS + c NSI Full ms{'2' if l == 2 else ''} [{i}]" for i, l in enumerate(lvl)]
    return dict(scan=np.arange(1, n + 1), ms_level=lvl, rt_min=np.linspace(0, 10, n), injection_time_ms=np.full(n, 20.0),
                precursor_mz=prec, isolation_target_mz=prec, mono_mz=np.full(n, np.nan),
                isolation_lower_offset=np.where(ms2, 0.7, np.nan), isolation_upper_offset=np.where(ms2, 0.7, np.nan),
                precursor_charge=ch, precursor_charge_mask=ch <= 0, master_scan=np.where(ms2, 1, -1),
                master_scan_mask=~ms2, activation=np.where(ms2, 2, 0).astype(np.uint8), activation_mask=~ms2,
                nce=np.where(ms2, 28.0, np.nan), analyzer=np.ones(n, np.uint8), polarity=np.ones(n, np.int8),
                peak_source=np.zeros(n, np.uint8), filter_string=filt, peak_offsets=off, mz=mz, intensity=it)


def _write(tmp_path, batches, compression=None, run_id=7):
    p = str(tmp_path / "x.vspec.arrow")
    w = St.StoreWriter(p, dict(run_name="x", run_id=run_id, vendor="test", instrument_model="m", source_path="s",
                                source_size=0, converter="test", converter_version="0", acquisition="DDA",
                                peak_order="mz_ascending"), compression=compression)
    for b in batches:
        w.write(b)
    return p, w.close()


@pytest.mark.parametrize("compression", [None, "lz4"])
def test_roundtrip_and_invariants(tmp_path, compression):
    b1, b2 = _synthetic(10, 1), _synthetic(7, 2)
    p, meta = _write(tmp_path, [b1, b2], compression)
    t = St.read_table(p)
    assert t.schema.equals(S.SPECTRUM_SCHEMA, check_metadata=False)
    n = 17
    assert t.num_rows == n and meta["n_spectra"] == n
    idx = t.column("index").to_numpy()
    assert np.array_equal(idx, np.arange(n))
    sid = t.column("spectrum_id").to_numpy()
    assert np.array_equal(sid, (np.uint64(7) << np.uint64(32)) | idx.astype(np.uint64))
    P = St.peaks_numpy(t)
    exp_mz = np.concatenate([b1["mz"], b2["mz"]])
    assert np.array_equal(P.mz, exp_mz)
    assert np.array_equal(np.diff(P.offsets), np.concatenate([np.diff(b1["peak_offsets"]), np.diff(b2["peak_offsets"])]))
    assert meta["n_peaks"] == len(exp_mz)
    lvl = t.column("ms_level").to_numpy()
    assert meta["n_ms1"] == int((lvl == 1).sum()) and meta["n_ms2"] == int((lvl == 2).sum())
    # MS1 rows: precursor & activation null; MS2 rows non-null
    for f in ("precursor_mz", "precursor_charge", "activation", "master_scan", "nce"):
        nulls = np.asarray(t.column(f).is_null())
        assert np.array_equal(nulls, lvl == 1), f
    assert t.column("im_index").null_count == n
    assert t.column("filter_string")[3].as_py() == b1["filter_string"][3]
    side = json.load(open(p + ".json"))
    assert side["n_spectra"] == n
    rm = St.run_meta(p)
    assert rm["run_id"] == 7 and rm["store_version"] == S.STORE_VERSION
    assert not os.path.exists(p + ".partial")


def test_offsets_are_chunked_int32(tmp_path):
    p, _ = _write(tmp_path, [_synthetic(5, 3), _synthetic(5, 4)])
    t = St.read_table(p)
    assert t.column("mz").num_chunks == 2
    assert t.column("mz").type == pa.list_(pa.float32())


def test_default_run_id_deterministic():
    assert S.default_run_id("abc") == S.default_run_id("abc")
    assert 0 <= S.default_run_id("abc") < 2 ** 32


def test_mmap_policy(tmp_path):
    # never mmap outside /tmp (Lustre); test only the decision, path under tmp_path may be anywhere
    p, _ = _write(tmp_path, [_synthetic(3, 5)])
    r = St.open_store(p, mmap=False)
    assert r.num_record_batches == 1
