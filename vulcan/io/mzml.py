"""Fast streaming mzML -> store reader.

Strategy: read the mzML in ~32 MB blocks cut at </spectrum>; within a block every field is
found with ONE compiled-regex scan over the whole block (C speed) and assigned to its
spectrum with np.searchsorted on the spectrum start offsets. The only per-spectrum work is
base64 + zlib decoding of the binary arrays (C calls; no XML tree, no per-spectrum dicts).
Handles TRFP and msconvert attribute orders (accession before value). numpress is rejected.
"""
from __future__ import annotations

import binascii
import os
import re
import time
import zlib

import numpy as np

from . import schema as S
from .store import StoreWriter

BLOCK = 32 << 20

_RE_SPEC = re.compile(rb'<spectrum\s[^>]*?id="([^"]*)"[^>]*>')
_RE_CV = re.compile(rb'<cvParam\s[^>]*?accession="(MS:\d+)"[^>]*?\svalue="([^"]*)"')
_RE_SCANREF = re.compile(rb'<scan\s[^>]*?instrumentConfigurationRef="([^"]+)"')
_RE_PRECREF = re.compile(rb'<precursor\s[^>]*?spectrumRef="[^"]*?scan=(\d+)"')
_RE_MONO = re.compile(rb'Monoisotopic M/Z:?"[^>]*?value="([^"]+)"')
_RE_BDA = re.compile(rb'<binaryDataArray\b[^>]*>(.*?)<binary>([^<]*)</binary>', re.S)
_RE_SCANNUM = re.compile(rb'scan=(\d+)')

# accession -> store column for numeric scalar cvParams
_NUM = {
    b"MS:1000511": "ms_level", b"MS:1000016": "rt", b"MS:1000927": "injection_time_ms",
    b"MS:1000827": "isolation_target_mz", b"MS:1000828": "isolation_lower_offset",
    b"MS:1000829": "isolation_upper_offset", b"MS:1000744": "precursor_mz", b"MS:1000041": "precursor_charge",
    b"MS:1000042": "precursor_intensity", b"MS:1000045": "nce", b"MS:1002680": "supplemental_energy",
    b"MS:1000285": "tic", b"MS:1000501": "scan_window_lower", b"MS:1000500": "scan_window_upper",
    b"MS:1001581": "faims_cv",
}
_ACT = {  # primary activation accessions
    b"MS:1000422": 2, b"MS:1000133": 1, b"MS:1000598": 3, b"MS:1000250": 4, b"MS:1000262": 5,
    b"MS:1000599": 6, b"MS:1003246": 7, b"MS:1003247": 8, b"MS:1002631": 3,
}
_SUPP = {b"MS:1002678": 1, b"MS:1002679": 2, b"MS:1003181": 2}  # supplemental CID / HCD
_ETHCD = b"MS:1002631"
_ANALYZER_ACC = {b"MS:1000484": 1, b"MS:1000264": 2, b"MS:1000082": 2, b"MS:1000083": 2, b"MS:1000291": 2,
                 b"MS:1000084": 3, b"MS:1000081": 5, b"MS:1003379": 4}


def _header_info(head: bytes) -> dict:
    ic = {}
    for m in re.finditer(rb'<instrumentConfiguration\s[^>]*id="([^"]+)"(.*?)</instrumentConfiguration>', head, re.S):
        code = 0
        for a in re.finditer(rb'<analyzer\b.*?</analyzer>', m.group(2), re.S):
            for acc in re.findall(rb'accession="(MS:\d+)"', a.group(0)):
                code = _ANALYZER_ACC.get(acc, code)
        ic[m.group(1)] = code
    model = re.search(rb'<referenceableParamGroup\b.*?<cvParam\s[^>]*?name="([^"]+)"', head, re.S)
    default_ic = re.search(rb'<run\s[^>]*?defaultInstrumentConfigurationRef="([^"]+)"', head)
    run_id = re.search(rb'<run\s[^>]*?id="([^"]+)"', head)
    sha1 = re.search(rb'accession="MS:1000569" value="([0-9a-f]+)"', head)
    return dict(ic=ic, model=model.group(1).decode() if model else "unknown",
                default_ic=default_ic.group(1) if default_ic else (next(iter(ic)) if ic else b""),
                run=run_id.group(1).decode() if run_id else "", source_sha1=sha1.group(1).decode() if sha1 else None,
                trfp=b"ThermoRawFileParser" in head)


def _decode(header: bytes, payload: bytes) -> tuple[int, np.ndarray | None]:
    """Returns (kind, values) where kind 0 = m/z, 1 = intensity, -1 = other."""
    if b"MS:1000514" in header:
        kind = 0
    elif b"MS:1000515" in header:
        kind = 1
    else:
        return -1, None
    if b"MS:1002312" in header or b"MS:1002313" in header or b"MS:1002314" in header:
        raise NotImplementedError("numpress-compressed mzML is not supported")
    raw = binascii.a2b_base64(payload) if payload else b""
    if raw and b"MS:1000574" in header:
        raw = zlib.decompress(raw)
    dt = np.float64 if b"MS:1000523" in header else np.float32
    return kind, np.frombuffer(raw, dtype=dt)


def _parse_block(buf: bytes, hdr: dict) -> dict:
    starts, scans, nat = [], [], []
    for m in _RE_SPEC.finditer(buf):
        starts.append(m.start())
        sm = _RE_SCANNUM.search(m.group(1))
        scans.append(int(sm.group(1)) if sm else -1)
    n = len(starts)
    st = np.asarray(starts, np.int64)
    nan64 = np.full(n, np.nan)
    c = {k: nan64.copy() for k in set(_NUM.values())}
    act = np.zeros(n, np.uint8)
    supp = np.zeros(n, np.uint8)
    pol = np.ones(n, np.int8)
    # one scan over all cvParams
    pos, acc, val = [], [], []
    for m in _RE_CV.finditer(buf):
        pos.append(m.start()); acc.append(m.group(1)); val.append(m.group(2))
    pos = np.asarray(pos, np.int64)
    si = np.searchsorted(st, pos, "right") - 1
    acc = [a if i >= 0 else b"" for a, i in zip(acc, si)] if (si < 0).any() else acc  # params before 1st spectrum
    acc_a = np.asarray(acc, dtype="S12")
    val_a = np.asarray(val, dtype=object)
    for a, col in _NUM.items():
        sel = np.nonzero(acc_a == a)[0]
        if len(sel):
            c[col][si[sel]] = np.asarray(val_a[sel].astype("S32"), dtype="S32").astype(np.float64)
    for a, code in _ACT.items():
        sel = np.nonzero(acc_a == a)[0]
        if len(sel):
            act[si[sel]] = code
            if a == _ETHCD:
                supp[si[sel]] = 2
    for a, code in _SUPP.items():
        sel = np.nonzero(acc_a == a)[0]
        if len(sel):
            supp[si[sel]] = code
    sel = np.nonzero(acc_a == b"MS:1000129")[0]
    pol[si[sel]] = -1
    # analyzer via scan instrumentConfigurationRef
    analyzer = np.full(n, hdr["ic"].get(hdr["default_ic"], 0), np.uint8)
    for m in _RE_SCANREF.finditer(buf):
        analyzer[np.searchsorted(st, m.start(), "right") - 1] = hdr["ic"].get(m.group(1), 0)
    master = np.full(n, -1, np.int32)
    for m in _RE_PRECREF.finditer(buf):
        master[np.searchsorted(st, m.start(), "right") - 1] = int(m.group(1))
    mono = nan64.copy()
    for m in _RE_MONO.finditer(buf):
        mono[np.searchsorted(st, m.start(), "right") - 1] = float(m.group(1))
    # filter strings (MS:1000512, value is text)
    filt = [None] * n
    for m in re.finditer(rb'accession="MS:1000512"[^>]*?value="([^"]*)"', buf):
        filt[np.searchsorted(st, m.start(), "right") - 1] = m.group(1).decode()
    # binary arrays
    mz_parts = [None] * n
    it_parts = [None] * n
    for m in _RE_BDA.finditer(buf):
        k, arr = _decode(m.group(1), m.group(2))
        if k < 0:
            continue
        i = np.searchsorted(st, m.start(), "right") - 1
        if k == 0:
            mz_parts[i] = arr
        else:
            it_parts[i] = arr
    empty = np.empty(0, np.float32)
    lens = np.fromiter((0 if a is None else len(a) for a in mz_parts), np.int64, n)
    off = np.zeros(n + 1, np.int64)
    np.cumsum(lens, out=off[1:])
    mz = np.concatenate([empty] + [a for a in mz_parts if a is not None]).astype(np.float32, copy=False)
    it = np.concatenate([empty] + [a for a in it_parts if a is not None]).astype(np.float32, copy=False)
    if len(mz) != len(it):
        raise ValueError("m/z and intensity arrays differ in length")
    lvl = c["ms_level"]
    lvl = np.where(np.isnan(lvl), 1, lvl).astype(np.uint8)
    ms1 = lvl < 2
    charge = np.nan_to_num(c["precursor_charge"], nan=0).astype(np.int16)
    out = dict(
        scan=np.asarray(scans, np.int64), ms_level=lvl, rt_min=c["rt"], injection_time_ms=c["injection_time_ms"],
        precursor_mz=c["precursor_mz"], isolation_target_mz=c["isolation_target_mz"], mono_mz=mono,
        isolation_lower_offset=c["isolation_lower_offset"], isolation_upper_offset=c["isolation_upper_offset"],
        precursor_charge=charge, precursor_charge_mask=charge <= 0, precursor_intensity=c["precursor_intensity"],
        master_scan=master, master_scan_mask=master < 0, activation=act, activation_mask=ms1,
        nce=c["nce"], supplemental_activation=supp, supplemental_activation_mask=supp == 0,
        supplemental_energy=c["supplemental_energy"], analyzer=analyzer, resolution=None, faims_cv=c["faims_cv"],
        polarity=pol, filter_string=filt if any(f is not None for f in filt) else None,
        scan_window_lower=c["scan_window_lower"], scan_window_upper=c["scan_window_upper"], tic=c["tic"],
        peak_source=np.full(n, 3, np.uint8), peak_offsets=off, mz=mz, intensity=it,
    )
    return out


def iter_blocks(path: str, block: int = BLOCK):
    """Yields (header_bytes_or_None, spectra_bytes) cut at </spectrum> boundaries."""
    with open(path, "rb", buffering=0) as fh:
        try:
            os.posix_fadvise(fh.fileno(), 0, 0, os.POSIX_FADV_SEQUENTIAL)
        except (AttributeError, OSError):
            pass
        rest = b""
        header = None
        done = False
        while not done:
            chunk = fh.read(block)
            if not chunk:
                done = True
            buf = rest + chunk
            if header is None:
                k = buf.find(b"<spectrumList")
                if k < 0:
                    if done:
                        raise ValueError("no <spectrumList> in mzML")
                    rest = buf
                    continue
                header, buf = buf[:k], buf[k:]
                yield header, b""
            end = buf.find(b"</spectrumList>")
            if end >= 0:
                yield None, buf[:end]
                return
            cut = buf.rfind(b"</spectrum>")
            if cut < 0:
                rest = buf
                continue
            cut += len(b"</spectrum>")
            yield None, buf[:cut]
            rest = buf[cut:]
            if done and rest.strip():
                raise ValueError("truncated mzML")


def convert(mzml_path: str, out_path: str, run_name: str | None = None, run_id: int | None = None,
            compression: str | None = None, source_raw: str | None = None, block: int = BLOCK) -> dict:
    t0 = time.perf_counter()
    hdr = None
    w = None
    t_parse = t_write = c_parse = c_write = 0.0
    for head, buf in iter_blocks(mzml_path, block):
        if head is not None:
            hdr = _header_info(head)
            run_name = run_name or hdr["run"] or os.path.basename(mzml_path).rsplit(".", 1)[0]
            src = source_raw or mzml_path
            meta = dict(run_name=run_name, vendor="Thermo" if hdr["trfp"] else "unknown",
                        instrument_model=hdr["model"], source_path=os.path.abspath(src),
                        source_size=os.path.getsize(src), source_sha1=hdr["source_sha1"],
                        intermediate_path=os.path.abspath(mzml_path), intermediate_size=os.path.getsize(mzml_path),
                        converter="vulcan.io.mzml (regex/numpy streaming)", converter_version="1",
                        acquisition="unknown", peak_order="as_given")
            if run_id is not None:
                meta["run_id"] = run_id
            w = StoreWriter(out_path, meta, compression=compression)
            continue
        if not buf:
            continue
        a = time.perf_counter(); ac = time.process_time()
        cols = _parse_block(buf, hdr)
        t_parse += time.perf_counter() - a; c_parse += time.process_time() - ac
        a = time.perf_counter(); ac = time.process_time()
        w.write(cols)
        t_write += time.perf_counter() - a; c_write += time.process_time() - ac
    return w.close({"timing_s": {"parse": t_parse, "arrow_write": t_write, "total": time.perf_counter() - t0},
                    "cpu_s": {"parse_incl_read": c_parse, "arrow_write": c_write}})
