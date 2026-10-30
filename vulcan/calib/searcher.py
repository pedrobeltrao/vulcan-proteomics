"""Pilot-pass searcher interface (DESIGN §4 L0.5) and the interim Sage implementation.

The pilot needs one thing from a search engine: rank-1 PSMs with q-values for a batch of spectra taken from a store.
`Searcher.search(batch, cfg, workdir)` returns a pyarrow Table with PSM_SCHEMA. `SageSearcher` (tools/bin/sage via an
MGF) is the interim implementation; the GPU search (vulcan.gpu) replaces it by implementing the same method.

Spectra cross the interface as a `SpectraBatch` (columnar numpy, no per-spectrum Python objects in the hot path) so
calibrated spectra can be searched without touching the store.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from dataclasses import asdict, dataclass, field, replace
from typing import Protocol

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pcsv

from ..io import store

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SAGE_BIN = os.environ.get("VULCAN_SAGE", os.path.join(ROOT, "tools", "bin", "sage"))

PSM_SCHEMA = pa.schema([
    pa.field("spectrum_id", pa.uint64(), nullable=False),
    pa.field("peptide", pa.string(), nullable=False),          # modified sequence, Sage notation (X[+mass])
    pa.field("charge", pa.int8(), nullable=False),
    pa.field("isotope_error", pa.int8(), nullable=False),      # observed = M + isotope_error neutrons
    pa.field("exp_mass", pa.float64(), nullable=False),        # (precursor_mz - PROTON) * z as searched
    pa.field("calc_mass", pa.float64(), nullable=False),
    pa.field("score", pa.float32(), nullable=False),           # engine primary score (Sage: hyperscore)
    pa.field("q_value", pa.float32(), nullable=False),         # spectrum-level q (engine TDC on its own decoys)
    pa.field("peptide_q", pa.float32(), nullable=False),
    pa.field("is_decoy", pa.bool_(), nullable=False),          # pilot-internal only; never leaves vulcan.calib
    pa.field("rt", pa.float32(), nullable=False),
    pa.field("matched_peaks", pa.int16(), nullable=False),
    pa.field("proteins", pa.string()),
])


@dataclass
class SpectraBatch:
    spectrum_id: np.ndarray      # u64
    rt: np.ndarray               # minutes
    precursor_mz: np.ndarray     # f64 (NaN -> spectrum skipped)
    charge: np.ndarray           # i64, 0 = unknown
    mz: np.ndarray               # f32/f64 flat peaks
    intensity: np.ndarray
    offsets: np.ndarray          # i64 (n+1)
    ms_level: int = 2

    def __len__(self):
        return len(self.spectrum_id)

    def take(self, idx: np.ndarray) -> "SpectraBatch":
        idx = np.asarray(idx, np.int64)
        lens = np.diff(self.offsets)[idx]
        off = np.zeros(len(idx) + 1, np.int64)
        np.cumsum(lens, out=off[1:])
        src = np.repeat(self.offsets[idx], lens) + (np.arange(off[-1]) - np.repeat(off[:-1], lens))
        return replace(self, spectrum_id=self.spectrum_id[idx], rt=self.rt[idx], precursor_mz=self.precursor_mz[idx],
                       charge=self.charge[idx], mz=self.mz[src], intensity=self.intensity[src], offsets=off)


def load_batch(path: str, ms_level: int = 2, precursor: str = "precursor_mz") -> SpectraBatch:
    """All spectra of one level from a store. `precursor`: store column used as the precursor m/z
    (precursor_mz = instrument mono if reported, else the isolation target; see STORE_SPEC §5.3)."""
    t = store.read_table(path, columns=["spectrum_id", "ms_level", "rt_min", "precursor_mz", "isolation_target_mz",
                                        "precursor_charge", "mz", "intensity"])
    t = t.filter(pc.equal(t.column("ms_level"), ms_level))
    pk = store.peaks_numpy(t)
    z = np.asarray(store.column_numpy(t, "precursor_charge"), np.float64)
    pmz = np.asarray(store.column_numpy(t, precursor), np.float64)
    meta = store.run_meta(path)
    mz, it = pk.mz, pk.intensity
    if meta.get("peak_order") != "mz_ascending":
        lens = np.diff(pk.offsets)
        o = np.lexsort((mz, np.repeat(np.arange(len(lens)), lens)))
        mz, it = mz[o], it[o]
    return SpectraBatch(np.asarray(t.column("spectrum_id").to_numpy(), np.uint64),
                        np.asarray(store.column_numpy(t, "rt_min"), np.float64), pmz,
                        np.where(np.isnan(z), 0, z).astype(np.int64), mz, it, pk.offsets, ms_level)


@dataclass
class SearchConfig:
    fasta: str
    precursor_ppm: tuple = (-50.0, 50.0)
    fragment_ppm: tuple | None = (-50.0, 50.0)
    fragment_da: tuple | None = None             # used instead of ppm for low-resolution MS2
    precursor_da: tuple | None = None            # open search (overrides precursor_ppm); error window, Da
    isotope_errors: tuple = (-1, 3)
    precursor_charge: tuple = (2, 4)
    static_mods: dict = field(default_factory=lambda: {"C": 57.021464})
    variable_mods: dict = field(default_factory=lambda: {"M": [15.9949], "^": [42.010565]})
    max_variable_mods: int = 2
    missed_cleavages: int = 2
    min_len: int = 7
    max_len: int = 50
    peptide_min_mass: float = 500.0
    peptide_max_mass: float = 5000.0
    min_peaks: int = 15
    max_peaks: int = 150
    min_matched_peaks: int = 4
    max_fragment_charge: int = 1
    report_psms: int = 1
    predict_rt: bool = True
    threads: int | None = None

    def to_dict(self):
        return asdict(self)


class Searcher(Protocol):
    name: str

    def search(self, batch: SpectraBatch, cfg: SearchConfig, workdir: str) -> tuple[pa.Table, dict]:
        """Return (PSM table with PSM_SCHEMA, accounting dict with wall_s / core_s / gpu_s)."""
        ...


# ---------------------------------------------------------------------------------------------------- Sage
def write_mgf(batch: SpectraBatch, path: str) -> int:
    """MGF with TITLE=<spectrum_id>; CHARGE only if the instrument gave one (Sage then tries cfg.precursor_charge)."""
    ok = np.flatnonzero(np.isfinite(batch.precursor_mz) & (np.diff(batch.offsets) > 0))
    mz = np.asarray(batch.mz, np.float64)
    it = np.asarray(batch.intensity, np.float64)
    with open(path, "w", buffering=1 << 22) as fh:
        for i in ok.tolist():          # text serialisation only; one formatted block per spectrum
            a, b = batch.offsets[i], batch.offsets[i + 1]
            ch = f"CHARGE={int(batch.charge[i])}+\n" if batch.charge[i] > 0 else ""
            peaks = "\n".join(f"{m:.6f} {x:.2f}" for m, x in zip(mz[a:b].tolist(), it[a:b].tolist()))
            fh.write(f"BEGIN IONS\nTITLE={int(batch.spectrum_id[i])}\nPEPMASS={batch.precursor_mz[i]:.7f}\n{ch}"
                     f"RTINSECONDS={batch.rt[i] * 60.0:.4f}\n{peaks}\nEND IONS\n")
    return len(ok)


def sage_config(cfg: SearchConfig, mgf: str, outdir: str) -> dict:
    """SearchConfig windows are error windows, error = experimental - theoretical. Sage 0.14.7 centres the precursor
    window on the experimental mass and accepts peptide masses in [exp + lo, exp + hi], i.e. error in [-hi, -lo]:
    the precursor window is negated and swapped here. Sage's fragment window is centred on the experimental peak in
    the fragment index (database.rs) but on the theoretical ion in scoring (scoring.rs), so asymmetric fragment
    windows are not well defined; the fragment window is symmetrised to +-max(|lo|, |hi|)."""
    if cfg.precursor_da:
        prec = {"da": [-cfg.precursor_da[1], -cfg.precursor_da[0]]}
    else:
        prec = {"ppm": [-cfg.precursor_ppm[1], -cfg.precursor_ppm[0]]}
    fw = cfg.fragment_da if cfg.fragment_da else cfg.fragment_ppm
    fw = max(abs(fw[0]), abs(fw[1]))
    frag = {"da": [-fw, fw]} if cfg.fragment_da else {"ppm": [-fw, fw]}
    return {
        "database": {
            "bucket_size": 8192 if not cfg.precursor_da else 32768,
            "enzyme": {"missed_cleavages": cfg.missed_cleavages, "min_len": cfg.min_len, "max_len": cfg.max_len,
                       "cleave_at": "KR", "restrict": "P", "c_terminal": True},
            "peptide_min_mass": cfg.peptide_min_mass, "peptide_max_mass": cfg.peptide_max_mass,
            "ion_kinds": ["b", "y"], "min_ion_index": 2,
            "static_mods": cfg.static_mods, "variable_mods": cfg.variable_mods,
            "max_variable_mods": cfg.max_variable_mods, "decoy_tag": "rev_", "generate_decoys": True,
            "fasta": cfg.fasta,
        },
        "precursor_tol": prec, "fragment_tol": frag,
        "precursor_charge": list(cfg.precursor_charge), "isotope_errors": list(cfg.isotope_errors),
        "deisotope": False, "chimera": False, "wide_window": False, "predict_rt": cfg.predict_rt,
        "min_peaks": cfg.min_peaks, "max_peaks": cfg.max_peaks, "min_matched_peaks": cfg.min_matched_peaks,
        "max_fragment_charge": cfg.max_fragment_charge, "report_psms": cfg.report_psms,
        "output_directory": outdir, "mzml_paths": [mgf],
    }


def _time_v(path: str) -> dict:
    txt = open(path).read()
    g = lambda pat: float(re.search(pat, txt).group(1)) if re.search(pat, txt) else 0.0
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", txt)
    w = sum(float(p) * 60 ** i for i, p in enumerate(reversed(wall.group(1).split(":")))) if wall else None
    u, s = g(r"User time \(seconds\): ([\d.]+)"), g(r"System time \(seconds\): ([\d.]+)")
    return {"wall_s": w, "user_s": u, "sys_s": s, "core_s": u + s,
            "max_rss_mb": g(r"Maximum resident set size \(kbytes\): (\d+)") / 1024}


class SageSearcher:
    """Interim pilot searcher: Sage 0.14.7 on an MGF written from the batch. Thread count from cfg.threads or
    vulcan.hw.n_threads(); timing from /usr/bin/time -v (core-s = user + sys)."""
    name = "sage"

    def __init__(self, binary: str = SAGE_BIN, annotate: bool = False):
        self.binary = binary
        self.annotate = annotate    # also write matched_fragments.sage.tsv (per-fragment errors for calibration)

    def search(self, batch: SpectraBatch, cfg: SearchConfig, workdir: str) -> tuple[pa.Table, dict]:
        from .. import hw
        os.makedirs(workdir, exist_ok=True)
        mgf = os.path.join(workdir, "spectra.mgf")
        t0 = time.time()
        n = write_mgf(batch, mgf)
        t_mgf = time.time() - t0
        conf = sage_config(cfg, mgf, workdir)
        cpath = os.path.join(workdir, "sage.json")
        with open(cpath, "w") as fh:
            json.dump(conf, fh, indent=1)
        nthr = cfg.threads or hw.n_threads()
        env = dict(os.environ, **hw.thread_env(nthr))
        tv = os.path.join(workdir, "time_v.txt")
        with open(os.path.join(workdir, "sage.log"), "w") as log:
            cmd = ["/usr/bin/time", "-v", "-o", tv, self.binary, cpath, "--batch-size", "1",
                   "--disable-telemetry-i-dont-want-to-improve-sage"]
            if self.annotate:
                cmd.append("--annotate-matches")
            r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0:
            raise RuntimeError(f"sage failed ({r.returncode}); see {workdir}/sage.log")
        acct = _time_v(tv)
        acct.update(n_spectra=n, mgf_write_s=t_mgf, threads=nthr, searcher="sage")
        psms = read_sage(os.path.join(workdir, "results.sage.tsv"))
        frag = os.path.join(workdir, "matched_fragments.sage.tsv")
        if self.annotate and os.path.exists(frag):
            acct["fragments"] = frag
        return psms, acct


def read_sage(path: str) -> pa.Table:
    t = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"),
                      convert_options=pcsv.ConvertOptions(column_types={"scannr": pa.string(), "peptide": pa.string()}))
    sid = np.array([int(s.split("=")[-1]) if "=" in s else int(s) for s in t.column("scannr").to_pylist()],
                   dtype=np.uint64)
    f = lambda c, d: np.asarray(t.column(c).to_numpy(zero_copy_only=False), d)
    return pa.table({
        "spectrum_id": pa.array(sid, pa.uint64()),
        "peptide": t.column("peptide").cast(pa.string()),
        "charge": pa.array(f("charge", np.int8)),
        "isotope_error": pa.array(np.rint(f("isotope_error", np.float64) / 1.0033548).astype(np.int8)),  # Sage: Da
        "exp_mass": pa.array(f("expmass", np.float64)),
        "calc_mass": pa.array(f("calcmass", np.float64)),
        "score": pa.array(f("hyperscore", np.float32)),
        "q_value": pa.array(f("spectrum_q", np.float32)),
        "peptide_q": pa.array(f("peptide_q", np.float32)),
        "is_decoy": pa.array(f("label", np.int64) < 0),
        "rt": pa.array(f("rt", np.float32)),
        "matched_peaks": pa.array(f("matched_peaks", np.int16)),
        "proteins": t.column("proteins").cast(pa.string()),
    }, schema=PSM_SCHEMA).append_column("psm_id", pa.array(f("psm_id", np.int64)))
