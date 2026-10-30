"""Write / read the Vulcan spectral store (Arrow IPC file, see schema.py and STORE_SPEC.md).

Writers hand over *columnar* batches (numpy arrays + offsets); nothing here loops per spectrum.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
from dataclasses import dataclass

import numpy as np
import pyarrow as pa

from . import schema as S

# float sentinel -> null for these columns (NaN means missing)
_FLOAT_NULL_NAN = {
    "injection_time_ms", "precursor_mz", "isolation_target_mz", "mono_mz", "isolation_lower_offset",
    "isolation_upper_offset", "precursor_intensity", "nce", "supplemental_energy", "resolution",
    "faims_cv", "scan_window_lower", "scan_window_upper", "tic", "im",
}


def configure_threads(n: int | None = None) -> int:
    """Pin Arrow CPU/IO pools (os.cpu_count() lies under affinity/cgroup limits; see vulcan.hw)."""
    n = int(n or os.environ.get("VULCAN_THREADS", os.environ.get("OMP_NUM_THREADS", 1)))
    pa.set_cpu_count(max(1, n))
    pa.set_io_thread_count(max(1, n))
    return n


_NP = {pa.uint8(): np.uint8, pa.uint16(): np.uint16, pa.uint32(): np.uint32, pa.uint64(): np.uint64,
       pa.int8(): np.int8, pa.int16(): np.int16, pa.int32(): np.int32, pa.int64(): np.int64,
       pa.float32(): np.float32, pa.float64(): np.float64}


def _arr(name: str, v, n: int, mask=None) -> pa.Array:
    typ = S.SPECTRUM_SCHEMA.field(name).type
    if v is None:
        return pa.nulls(n, type=typ)
    if isinstance(v, pa.Array):
        return v
    v = np.asarray(v)
    if mask is None and name in _FLOAT_NULL_NAN and v.dtype.kind == "f":
        mask = np.isnan(v)
        if not mask.any():
            mask = None
    if mask is not None and not np.any(mask):
        mask = None
    return pa.array(v.astype(_NP[typ], copy=False), type=typ, mask=mask)


def _list(values: np.ndarray, offsets: np.ndarray, typ: pa.DataType, mask=None) -> pa.Array:
    off = np.asarray(offsets)
    off = (off - off[0]).astype(np.int32)  # chunked int32 offsets
    vals = pa.array(np.ascontiguousarray(values[int(offsets[0]):int(offsets[-1])]), type=typ.value_type)
    return pa.ListArray.from_arrays(pa.array(off, pa.int32()), vals, type=typ,
                                    mask=None if mask is None else pa.array(mask, pa.bool_()))


class StoreWriter:
    """Streams record batches into one run's .vspec.arrow file.

    Static run metadata goes into the schema metadata (written first); final counts are
    merged into the footer metadata at close().
    """

    def __init__(self, path: str, run_meta: dict, compression: str | None = None):
        self.path = path
        self.meta = dict(run_meta)
        self.meta.setdefault("store_version", S.STORE_VERSION)
        self.meta.setdefault("created_utc", _dt.datetime.now(_dt.timezone.utc).isoformat())
        self.meta.setdefault("run_id", S.default_run_id(self.meta["run_name"]))
        for k in ("n_spectra", "n_ms1", "n_ms2", "n_peaks"):
            self.meta.setdefault(k, None)
        self.meta["compression"] = compression or "none"
        self.run_id = int(self.meta["run_id"])
        self.schema = S.make_schema(self.meta)
        opts = pa.ipc.IpcWriteOptions(compression=compression, use_threads=False)
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        self._tmp = path + ".partial"
        self._sink = pa.OSFile(self._tmp, "wb")
        self._w = pa.ipc.new_file(self._sink, self.schema, options=opts)
        self.n = 0
        self.n_peaks = 0
        self.level_counts: dict[int, int] = {}

    def write(self, cols: dict) -> None:
        """cols: numpy arrays keyed by schema field names (scalar fields, NaN = null for floats),
        plus peaks: 'peak_offsets' (n+1, any int), 'mz', 'intensity', optional 'im_index';
        optional 'sps_offsets'+'sps'; optional 'filter_offsets'+'filter_bytes' or 'filter_string' (list/pa.Array).
        Optional '<name>_mask' (bool, True = null) for integer nullable fields."""
        n = len(cols["scan"])
        if n == 0:
            return
        idx = np.arange(self.n, self.n + n, dtype=np.uint64)
        cols = dict(cols)
        cols["index"] = idx.astype(np.uint32)
        cols["spectrum_id"] = (np.uint64(self.run_id) << np.uint64(32)) | idx
        cols["run_id"] = np.full(n, self.run_id, np.uint32)
        arrays = []
        po = np.asarray(cols["peak_offsets"], dtype=np.int64)
        for f in S.SPECTRUM_FIELDS:
            name = f.name
            if name in ("mz", "intensity"):
                arrays.append(_list(cols[name], po, f.type))
            elif name == "im_index":
                arrays.append(_list(cols[name], po, f.type) if cols.get(name) is not None else pa.nulls(n, f.type))
            elif name == "sps_mz":
                if cols.get("sps_offsets") is not None:
                    so = np.asarray(cols["sps_offsets"], dtype=np.int64)
                    arrays.append(_list(np.asarray(cols["sps"], np.float64), so, f.type, mask=np.diff(so) == 0))
                else:
                    arrays.append(pa.nulls(n, f.type))
            elif name == "filter_string":
                if cols.get("filter_offsets") is not None:
                    fo = np.asarray(cols["filter_offsets"], dtype=np.int64)
                    data = cols["filter_bytes"][: int(fo[-1])]
                    arrays.append(pa.StringArray.from_buffers(n, pa.py_buffer(np.ascontiguousarray((fo - fo[0]).astype(np.int32))),
                                                              pa.py_buffer(np.ascontiguousarray(data)), None, 0))
                elif cols.get("filter_string") is not None:
                    arrays.append(pa.array(cols["filter_string"], pa.string()))
                else:
                    arrays.append(pa.nulls(n, pa.string()))
            else:
                arrays.append(_arr(name, cols.get(name), n, cols.get(name + "_mask")))
        batch = pa.RecordBatch.from_arrays(arrays, schema=self.schema)
        self._w.write_batch(batch)
        lv, cnt = np.unique(np.asarray(cols["ms_level"]), return_counts=True)
        for a, b in zip(lv.tolist(), cnt.tolist()):
            self.level_counts[a] = self.level_counts.get(a, 0) + b
        self.n += n
        self.n_peaks += int(po[-1] - po[0])

    def close(self, extra_meta: dict | None = None) -> dict:
        final = dict(self.meta)
        final.update(extra_meta or {})
        final.update(n_spectra=self.n, n_peaks=self.n_peaks, n_ms1=self.level_counts.get(1, 0),
                     n_ms2=self.level_counts.get(2, 0),
                     n_by_level={str(k): v for k, v in sorted(self.level_counts.items())})
        self._w.close()
        self._sink.close()
        # footer metadata cannot be added after the writer is created; rewrite footer cheaply is not
        # supported, so final counts go to a sidecar JSON *and* are recomputable from the file.
        with open(self.path + ".json", "w") as fh:
            json.dump(final, fh, indent=1, default=str)
        os.replace(self._tmp, self.path)
        return final


# ---------------------------------------------------------------- reading
@dataclass
class Peaks:
    mz: np.ndarray          # float32, flat
    intensity: np.ndarray   # float32, flat
    offsets: np.ndarray     # int64, n+1
    im_index: np.ndarray | None = None


def open_store(path: str, mmap: bool | None = None) -> pa.ipc.RecordBatchFileReader:
    """mmap only on local disk (DESIGN: never mmap on Lustre). Default: mmap iff path is under /tmp."""
    if mmap is None:
        mmap = os.path.abspath(path).startswith("/tmp/")
    src = pa.memory_map(path, "r") if mmap else pa.OSFile(path, "rb")
    return pa.ipc.open_file(src)


def read_table(path: str, columns=None, mmap: bool | None = None) -> pa.Table:
    r = open_store(path, mmap)
    t = r.read_all()
    return t.select(columns) if columns else t


def run_meta(path: str) -> dict:
    meta = S.run_meta_of(open_store(path).schema)
    side = path + ".json"
    if os.path.exists(side):
        with open(side) as fh:
            meta.update(json.load(fh))
    return meta


def peaks_numpy(table: pa.Table, with_im: bool = False) -> Peaks:
    """Flatten the chunked list columns into contiguous numpy (one memcpy per chunk)."""
    mz_chunks, it_chunks, im_chunks, lens = [], [], [], []
    for mzc, itc, imc in zip(table.column("mz").chunks, table.column("intensity").chunks,
                             table.column("im_index").chunks if with_im else [None] * table.column("mz").num_chunks):
        off = mzc.offsets.to_numpy()
        mz_chunks.append(mzc.values.to_numpy()[off[0]:off[-1]] if len(mzc.values) else np.empty(0, np.float32))
        it_chunks.append(itc.values.to_numpy()[off[0]:off[-1]] if len(itc.values) else np.empty(0, np.float32))
        if with_im and imc is not None and imc.null_count < len(imc):
            im_chunks.append(imc.values.to_numpy()[off[0]:off[-1]])
        lens.append(np.diff(off))
    lens = np.concatenate(lens) if lens else np.empty(0, np.int64)
    offsets = np.zeros(len(lens) + 1, np.int64)
    np.cumsum(lens, out=offsets[1:])
    return Peaks(np.concatenate(mz_chunks) if mz_chunks else np.empty(0, np.float32),
                 np.concatenate(it_chunks) if it_chunks else np.empty(0, np.float32),
                 offsets, np.concatenate(im_chunks) if im_chunks else None)


def column_numpy(table: pa.Table, name: str) -> np.ndarray:
    """Zero-copy when there are no nulls; otherwise a float64 copy with NaN for nulls."""
    col = table.column(name).combine_chunks()
    if col.null_count:
        return col.to_numpy(zero_copy_only=False)
    return col.to_numpy()
