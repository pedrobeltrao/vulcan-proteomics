"""Bruker timsTOF (.d, TDF) -> store, via timsrust_pyo3.

Two modes (STORE_SPEC.md, "Bruker"):
  frames  (default, lossless): one row per MS1 frame and one row per (MS2 frame, quadrupole
          window) -- the diaPASEF window or ddaPASEF precursor scan range -- holding *every* TDF
          event of that row: m/z (f32), intensity (f32) and the per-peak TIMS scan index
          (`im_index`, u16; 1/K0 via the run-level `im_lut`). Peaks are ordered by (scan, m/z).
  spectra (derived, lossy): timsrust's SpectrumReader output -- IM-collapsed, centroided MS2
          spectra (DDA: one per precursor; DIA: one per frame x window). MS1 is not emitted.
m/z and 1/K0 use timsrust's converters (sqrt-linear TOF->m/z from MzAcqRange and
DigitizerNumSamples; linear scan->1/K0 from OneOverK0AcqRange); parameters are recorded in
run metadata so TOF indices are recoverable. Recalibration is the L0.5 pilot's job.
"""
from __future__ import annotations

import os
import sqlite3
import time

import numpy as np

from . import schema as S
from .store import StoreWriter


def _tdf_meta(d: str) -> dict:
    con = sqlite3.connect(f"file:{os.path.join(d, 'analysis.tdf')}?mode=ro", uri=True)
    g = dict(con.execute("select Key, Value from GlobalMetadata").fetchall())
    n_frames, max_scans = con.execute("select count(*), max(NumScans) from Frames").fetchone()
    msms = dict(con.execute("select MsMsType, count(*) from Frames group by MsMsType").fetchall())
    tables = {r[0] for r in con.execute("select name from sqlite_master where type='table'")}
    windows = []
    if "DiaFrameMsMsWindows" in tables:
        windows = con.execute("select WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy "
                              "from DiaFrameMsMsWindows order by WindowGroup, ScanNumBegin").fetchall()
    con.close()
    acq = "DIA-PASEF" if msms.get(9) else ("DDA-PASEF" if msms.get(8) else "MS1")
    return dict(global_meta=g, n_frames=n_frames, max_scans=max_scans, msms_counts=msms, windows=windows, acquisition=acq)


def _converters(g: dict, max_scans: int):
    mz_min, mz_max = float(g["MzAcqRangeLower"]), float(g["MzAcqRangeUpper"])
    tof_max = int(g["DigitizerNumSamples"])
    icpt = np.sqrt(mz_min)
    slope = (np.sqrt(mz_max) - icpt) / tof_max
    im_lo, im_hi = float(g["OneOverK0AcqRangeLower"]), float(g["OneOverK0AcqRangeUpper"])
    im_lut = (im_hi + (im_lo - im_hi) / max_scans * np.arange(max_scans + 1)).astype(np.float32)
    return icpt, slope, im_lut, dict(tof_intercept=icpt, tof_slope=slope, tof_max_index=tof_max,
                                     mz_formula="mz = (tof_intercept + tof_slope*tof)^2",
                                     im_formula="1/K0 = im_hi + (im_lo-im_hi)/max_scans*scan", im_lo=im_lo, im_hi=im_hi,
                                     max_scans=max_scans)


class _Acc:
    """Accumulates rows + peaks and flushes columnar batches."""

    def __init__(self, w: StoreWriter, with_im: bool):
        self.w, self.with_im = w, with_im
        self.reset()

    def reset(self):
        self.rows = {k: [] for k in ("scan", "ms_level", "rt_min", "isolation_target_mz", "isolation_lower_offset",
                                     "isolation_upper_offset", "nce", "window_index", "im_scan_begin", "im_scan_end",
                                     "im", "precursor_mz", "precursor_charge", "precursor_intensity", "tic")}
        self.mz, self.it, self.imi, self.lens = [], [], [], []

    def add(self, row: dict, mz, it, imi=None):
        for k, v in self.rows.items():
            v.append(row.get(k, np.nan))
        self.mz.append(mz); self.it.append(it); self.lens.append(len(mz))
        if self.with_im:
            self.imi.append(imi)
        if len(self.lens) >= S.BATCH_SPECTRA:
            self.flush()

    def flush(self):
        n = len(self.lens)
        if not n:
            return
        r = {k: np.asarray(v, dtype=np.float64) for k, v in self.rows.items()}
        off = np.zeros(n + 1, np.int64)
        np.cumsum(self.lens, out=off[1:])
        lvl = r["ms_level"].astype(np.uint8)
        ch = np.nan_to_num(r["precursor_charge"], nan=0).astype(np.int16)
        wi = np.nan_to_num(r["window_index"], nan=-1).astype(np.int16)
        sb = np.nan_to_num(r["im_scan_begin"], nan=0).astype(np.uint16)
        se = np.nan_to_num(r["im_scan_end"], nan=0).astype(np.uint16)
        cols = dict(
            scan=r["scan"].astype(np.uint32), ms_level=lvl, rt_min=r["rt_min"], injection_time_ms=None,
            precursor_mz=r["precursor_mz"], isolation_target_mz=r["isolation_target_mz"], mono_mz=None,
            isolation_lower_offset=r["isolation_lower_offset"], isolation_upper_offset=r["isolation_upper_offset"],
            precursor_charge=ch, precursor_charge_mask=ch <= 0, precursor_intensity=r["precursor_intensity"],
            master_scan=None, activation=np.where(lvl >= 2, 1, 0).astype(np.uint8), activation_mask=lvl < 2,
            nce=r["nce"], supplemental_activation=None, analyzer=np.full(n, 3, np.uint8), polarity=np.ones(n, np.int8),
            tic=r["tic"], peak_source=np.full(n, 4 if self.with_im else 5, np.uint8),
            window_index=wi, window_index_mask=wi < 0, im_scan_begin=sb, im_scan_begin_mask=se == 0,
            im_scan_end=se, im_scan_end_mask=se == 0, im=r["im"],
            peak_offsets=off, mz=np.concatenate(self.mz).astype(np.float32, copy=False),
            intensity=np.concatenate(self.it).astype(np.float32, copy=False),
            im_index=np.concatenate(self.imi).astype(np.uint16, copy=False) if self.with_im else None,
        )
        self.w.write(cols)
        self.reset()


NATIVE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "tools", "vulcan_tims")


def _native():
    """vulcan_tims (vulcan/io/rust/vulcan_tims, built by scripts/ingest/build_vulcan_tims.sh) or None."""
    import sys
    p = os.path.abspath(NATIVE_DIR)
    if os.path.isdir(p) and p not in sys.path:
        sys.path.insert(0, p)
    try:
        import vulcan_tims
        return vulcan_tims
    except ImportError:
        return None


def _block_cols(x: dict, icpt: float, slope: float, im_lut: np.ndarray, group_keys: np.ndarray,
                group_of_frame: np.ndarray, win_keys: np.ndarray, win_ids: np.ndarray) -> tuple[dict, int]:
    """Vectorised rows for one block of decoded frames (vulcan_tims.read output)."""
    fidx = x["frame_index"].astype(np.int64)
    nf = len(fidx)
    lvl = x["ms_level"]
    ev = x["event_offsets"].astype(np.int64)
    sp = x["scan_ptr"].astype(np.int64)
    so = x["scan_offsets"].astype(np.int64)
    qp = x["quad_ptr"].astype(np.int64)
    nq = np.diff(qp)
    nscan = np.diff(sp) - 1                       # scans per frame
    # MS1 (or MS2 without quad info): whole frame
    whole = (lvl != 2) | (nq == 0)
    f_whole = np.nonzero(whole)[0]
    # MS2 windows
    f_win = np.repeat(np.arange(nf), np.where(whole, 0, nq))
    q_sel = np.concatenate([np.arange(qp[i], qp[i + 1]) for i in np.nonzero(~whole)[0]]) if (~whole).any() else np.empty(0, np.int64)
    qs = x["quad_scan_start"].astype(np.int64)[q_sel]
    qe = x["quad_scan_end"].astype(np.int64)[q_sel]
    lo_w = ev[f_win] + so[sp[f_win] + np.minimum(qs, nscan[f_win])]
    hi_w = ev[f_win] + so[sp[f_win] + np.minimum(qe, nscan[f_win])]
    # merge rows in (frame, window) order
    r_frame = np.concatenate([f_whole, f_win])
    r_lo = np.concatenate([ev[f_whole], lo_w])
    r_hi = np.concatenate([ev[f_whole + 1], hi_w])
    r_q = np.concatenate([np.full(len(f_whole), -1), np.arange(len(f_win))])
    order = np.lexsort((r_lo, r_frame))
    r_frame, r_lo, r_hi, r_q = r_frame[order], r_lo[order], r_hi[order], r_q[order]
    n = len(r_frame)
    lens = r_hi - r_lo
    off = np.zeros(n + 1, np.int64)
    np.cumsum(lens, out=off[1:])
    gather = np.repeat(r_lo - off[:-1], lens) + np.arange(off[-1])
    # per-event TIMS scan index (scan_offsets are per-frame, n_scans+1 entries)
    counts = np.diff(so)
    counts[sp[1:-1] - 1] = 0                      # drop frame-boundary diffs
    keep = np.ones(len(counts), bool)
    keep[sp[1:-1] - 1] = False
    scan_id = (np.arange(len(so)) - np.repeat(sp[:-1], np.diff(sp)))[:-1][keep]
    ev_scan = np.repeat(scan_id, counts[keep]).astype(np.uint16)
    tof = x["tof"]
    inten = x["intensity"]
    tof_g = tof[gather].astype(np.float64)
    mz = ((icpt + slope * tof_g) ** 2).astype(np.float32)
    it = inten[gather].astype(np.float32)
    im_index = ev_scan[gather]
    is_w = r_q >= 0
    qq = np.where(is_w, r_q, 0)
    def qv(name, dtype=np.float64):
        v = np.full(n, np.nan)
        if len(f_win):
            v[is_w] = x[name][q_sel][qq[is_w]]
        return v
    iso = qv("quad_mz")
    wdt = qv("quad_width")
    ce = qv("quad_ce")
    sb = np.zeros(n, np.int64); se = np.zeros(n, np.int64)
    if len(f_win):
        sb[is_w] = qs[qq[is_w]]; se[is_w] = qe[qq[is_w]]
    grp = group_of_frame[np.searchsorted(group_keys, fidx[r_frame]).clip(0, max(len(group_keys) - 1, 0))] if len(group_keys) else np.full(n, -1)
    if len(group_keys):
        grp = np.where(group_keys[np.searchsorted(group_keys, fidx[r_frame]).clip(0, len(group_keys) - 1)] == fidx[r_frame], grp, -1)
    wk = grp * 100000 + sb
    wpos = np.searchsorted(win_keys, wk).clip(0, max(len(win_keys) - 1, 0))
    widx = np.where(is_w & (grp >= 0) & (len(win_keys) > 0) & (win_keys[wpos] == wk if len(win_keys) else False), win_ids[wpos] if len(win_ids) else -1, -1).astype(np.int16)
    lv = np.where(is_w, 2, np.where(lvl[r_frame] == 0, 1, lvl[r_frame])).astype(np.uint8)
    tic = np.add.reduceat(it, off[:-1]) if off[-1] else np.zeros(n)
    tic = np.where(lens > 0, tic, 0.0)
    mid = ((sb + se) // 2).clip(0, len(im_lut) - 1)
    cols = dict(
        scan=fidx[r_frame].astype(np.uint32), ms_level=lv, rt_min=x["rt_s"][r_frame] / 60.0, injection_time_ms=None,
        precursor_mz=None, isolation_target_mz=iso, mono_mz=None, isolation_lower_offset=wdt / 2,
        isolation_upper_offset=wdt / 2, precursor_charge=None, precursor_intensity=None, master_scan=None,
        activation=np.where(lv >= 2, 1, 0).astype(np.uint8), activation_mask=lv < 2, nce=ce,
        analyzer=np.full(n, 3, np.uint8), polarity=np.ones(n, np.int8), tic=tic, peak_source=np.full(n, 4, np.uint8),
        window_index=widx, window_index_mask=widx < 0, im_scan_begin=sb.astype(np.uint16), im_scan_begin_mask=~is_w,
        im_scan_end=se.astype(np.uint16), im_scan_end_mask=~is_w,
        im=np.where(is_w, im_lut[mid], np.nan), peak_offsets=off, mz=mz, intensity=it, im_index=im_index,
    )
    dropped = int(ev[-1] - ev[0] - off[-1])
    return cols, dropped


def convert(d: str, out_path: str, mode: str = "frames", run_name: str | None = None, run_id: int | None = None,
            compression: str | None = None, max_frames: int | None = None) -> dict:
    import timsrust_pyo3 as tr
    t0 = time.perf_counter(); c0 = time.process_time()
    d = d.rstrip("/")
    tm = _tdf_meta(d)
    g = tm["global_meta"]
    icpt, slope, im_lut, conv = _converters(g, tm["max_scans"])
    # window table: global window index by (group, scan_begin)
    win_index = {(int(w[0]), int(w[1])): i for i, w in enumerate(tm["windows"])}
    run_name = run_name or os.path.basename(d).rsplit(".", 1)[0]
    src_size = sum(os.path.getsize(os.path.join(d, f)) for f in ("analysis.tdf", "analysis.tdf_bin"))
    meta = dict(run_name=run_name, vendor="Bruker", instrument_model=g.get("InstrumentName", "timsTOF"),
                instrument_serial=g.get("InstrumentSerialNumber"), instrument_software=g.get("AcquisitionSoftwareVersion"),
                acquisition_date=g.get("AcquisitionDateTime"), source_path=os.path.abspath(d), source_size=src_size,
                converter=f"vulcan.io.bruker[{mode}] (timsrust_pyo3)", converter_version="1",
                acquisition=tm["acquisition"], peak_order="im_scan_then_mz" if mode.startswith("frames") else "mz_ascending",
                bruker_mode=mode, tims_converters=conv, im_lut=im_lut.tolist(),
                dia_windows=[dict(index=i, group=int(w[0]), scan_begin=int(w[1]), scan_end=int(w[2]), isolation_mz=w[3],
                                  isolation_width=w[4], collision_energy=w[5]) for i, w in enumerate(tm["windows"])])
    if run_id is not None:
        meta["run_id"] = run_id
    w = StoreWriter(out_path, meta, compression=compression)
    t_dec = 0.0
    c_dec = 0.0
    dropped = 0
    native = _native() if mode == "frames" else None
    if mode == "frames_pyo3":
        mode = "frames_list"
    if mode == "frames" and native is None:
        mode = "frames_list"
    if mode == "frames":
        con = sqlite3.connect(f"file:{os.path.join(d, 'analysis.tdf')}?mode=ro", uri=True)
        tabs = {r[0] for r in con.execute("select name from sqlite_master where type='table'")}
        gf = np.array(con.execute("select Frame, WindowGroup from DiaFrameMsMsInfo order by Frame").fetchall()
                      if "DiaFrameMsMsInfo" in tabs else np.empty((0, 2)), dtype=np.int64).reshape(-1, 2)
        con.close()
        wins = np.array([(int(w[0]) * 100000 + int(w[1]), i) for i, w in enumerate(tm["windows"])], dtype=np.int64).reshape(-1, 2)
        wo = np.argsort(wins[:, 0]) if len(wins) else np.empty(0, np.int64)
        rdr = native.TdfFrames(d)
        n_frames = len(rdr) if max_frames is None else min(max_frames, len(rdr))
        BLOCK_FRAMES = 1024
        for s0 in range(0, n_frames, BLOCK_FRAMES):
            a = time.perf_counter(); ac = time.process_time()
            x = rdr.read(s0, min(s0 + BLOCK_FRAMES, n_frames))
            t_dec += time.perf_counter() - a; c_dec += time.process_time() - ac
            cols, dr = _block_cols(x, icpt, slope, im_lut, gf[:, 0], gf[:, 1], wins[wo, 0] if len(wins) else np.empty(0, np.int64),
                                   wins[wo, 1] if len(wins) else np.empty(0, np.int64))
            dropped += dr
            w.write(cols)
        w.meta["converter"] = "vulcan.io.bruker[frames] (vulcan_tims native: timsrust 0.4.2 + numpy)"
    elif mode == "frames_list":
        acc = _Acc(w, with_im=True)
        fr = tr.FrameReader(d)
        n_frames = tm["n_frames"] if max_frames is None else min(max_frames, tm["n_frames"])
        group_of_frame = {}
        con = sqlite3.connect(f"file:{os.path.join(d, 'analysis.tdf')}?mode=ro", uri=True)
        tabs = {r[0] for r in con.execute("select name from sqlite_master where type='table'")}
        if "DiaFrameMsMsInfo" in tabs:
            group_of_frame = dict(con.execute("select Frame, WindowGroup from DiaFrameMsMsInfo").fetchall())
        con.close()
        for i in range(n_frames):
            a = time.perf_counter(); ac = time.process_time()
            f = fr.read_frame(i)
            tof = np.asarray(f.tof_indices, dtype=np.float64)
            inten = np.asarray(f.intensities, dtype=np.float32)
            so = np.asarray(f.scan_offsets, dtype=np.int64)
            t_dec += time.perf_counter() - a; c_dec += time.process_time() - ac
            mz = ((icpt + slope * tof) ** 2).astype(np.float32)
            scan_of = np.repeat(np.arange(len(so) - 1, dtype=np.uint16), np.diff(so))
            fid = int(f.index)
            rt = float(f.rt) / 60.0
            if f.ms_level == tr.MSLevel.MS1 or f.quadrupole_settings is None or not len(f.quadrupole_settings.scan_starts):
                acc.add(dict(scan=fid, ms_level=1, rt_min=rt, tic=float(inten.sum())), mz, inten, scan_of)
                continue
            q = f.quadrupole_settings
            grp = group_of_frame.get(fid)
            dropped += len(tof) - sum(int(so[min(int(e), len(so) - 1)] - so[min(int(b), len(so) - 1)])
                                      for b, e in zip(q.scan_starts, q.scan_ends))
            for k in range(len(q.scan_starts)):
                b0, b1 = int(q.scan_starts[k]), int(q.scan_ends[k])
                lo, hi = so[min(b0, len(so) - 1)], so[min(b1, len(so) - 1)]  # [b0, b1) scans
                wdt = float(q.isolation_width[k])
                it_w = inten[lo:hi]
                acc.add(dict(scan=fid, ms_level=2, rt_min=rt, isolation_target_mz=float(q.isolation_mz[k]),
                             isolation_lower_offset=wdt / 2, isolation_upper_offset=wdt / 2,
                             nce=float(q.collision_energy[k]),
                             window_index=win_index.get((grp, b0), np.nan) if grp is not None else np.nan,
                             im_scan_begin=b0, im_scan_end=b1, im=float(im_lut[min((b0 + b1) // 2, len(im_lut) - 1)]),
                             tic=float(it_w.sum())),
                        mz[lo:hi], it_w, scan_of[lo:hi])
        acc.flush()
    elif mode == "spectra":
        acc = _Acc(w, with_im=False)
        sr = tr.SpectrumReader(d)
        n = len(sr)
        for i in range(n):
            a = time.perf_counter(); ac = time.process_time()
            s = sr.get(i)
            mz = np.asarray(s.mz_values, dtype=np.float32)
            it = np.asarray(s.intensities, dtype=np.float32)
            t_dec += time.perf_counter() - a; c_dec += time.process_time() - ac
            p = s.precursor
            wdt = float(s.isolation_width)
            acc.add(dict(scan=int(p.frame_index), ms_level=2, rt_min=float(p.rt) / 60.0,
                         precursor_mz=float(p.mz) if tm["acquisition"] != "DIA-PASEF" else np.nan,
                         isolation_target_mz=float(s.isolation_mz), isolation_lower_offset=wdt / 2,
                         isolation_upper_offset=wdt / 2, nce=float(s.collision_energy), im=float(p.im),
                         precursor_charge=float(p.charge) if p.charge else np.nan,
                         precursor_intensity=float(p.intensity) if p.intensity else np.nan, tic=float(it.sum())),
                    mz, it)
        acc.flush()
    else:
        raise ValueError(mode)
    return w.close({"timing_s": {"timsrust_decode": t_dec, "total": time.perf_counter() - t0},
                    "cpu_s": {"timsrust_decode": c_dec, "rows_and_arrow_write": time.process_time() - c0 - c_dec},
                    "events_outside_windows": dropped, "converter": w.meta["converter"]})
