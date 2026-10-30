"""Vulcan L0 spectral store schema (store format version 1).

One Arrow IPC *file* per run. Each row is one spectrum; peaks live in list columns
(`mz`, `intensity`, optional `im_index`) whose int32 offsets are chunked per record
batch (<= BATCH_SPECTRA spectra per batch, far below 2^31 peaks), which is the
"chunked offsets" option of DESIGN.md section 4 L0. Run-level metadata is stored as
JSON in the Arrow schema metadata under the key ``vulcan``.

See vulcan/io/STORE_SPEC.md for field semantics, invariants and vendor mappings.
"""
from __future__ import annotations

import hashlib
import json

import pyarrow as pa

STORE_VERSION = 1
FILE_SUFFIX = ".vspec.arrow"
BATCH_SPECTRA = 4096  # spectra per record batch (offset chunk)

# ---- enumerations (u8 codes; keep in sync with vulcan/io/dotnet/VulcanRaw/Extractor.cs) ----
ACTIVATION = {
    0: "none", 1: "CID", 2: "HCD", 3: "ETD", 4: "ECD", 5: "MPD", 6: "PQD", 7: "UVPD",
    8: "NETD", 9: "PTR", 255: "other",
}
ACTIVATION_CODE = {v: k for k, v in ACTIVATION.items()}
ANALYZER = {0: "unknown", 1: "OT", 2: "IT", 3: "TOF", 4: "Astral", 5: "Quad", 6: "Sector"}
ANALYZER_CODE = {v: k for k, v in ANALYZER.items()}
PEAK_SOURCE = {
    0: "vendor_centroid_stream",      # Thermo centroid (label) stream, incl. reference/exception peaks
    1: "vendor_centroid_scan",        # Thermo scan acquired in centroid mode (e.g. ITMS centroid)
    2: "vendor_centroider_profile",   # Thermo profile scan centroided by RawFileReader
    3: "mzml_as_given",               # peaks copied from an mzML (centroided upstream)
    4: "tims_raw_events",             # Bruker TDF events (tof index -> m/z), lossless, one row per frame/window
    5: "tims_centroided",             # Bruker spectra centroided by timsrust (IM collapsed)
}

# ---- the spectrum table ----
SPECTRUM_FIELDS = [
    # identity
    pa.field("spectrum_id", pa.uint64(), nullable=False),     # (run_id << 32) | index
    pa.field("run_id", pa.uint32(), nullable=False),
    pa.field("index", pa.uint32(), nullable=False),           # 0-based row number within the run
    pa.field("scan", pa.uint32(), nullable=False),            # native scan number (Thermo) / frame id (Bruker)
    pa.field("ms_level", pa.uint8(), nullable=False),
    # timing
    pa.field("rt_min", pa.float64(), nullable=False),         # scan start time, minutes
    pa.field("injection_time_ms", pa.float32()),
    # precursor (null for MS1)
    pa.field("precursor_mz", pa.float64()),                   # best selected-ion m/z (mono if instrument-determined, else isolation target)
    pa.field("isolation_target_mz", pa.float64()),
    pa.field("mono_mz", pa.float64()),                        # instrument monoisotopic m/z (Thermo trailer), null if absent/0
    pa.field("isolation_lower_offset", pa.float32()),         # mzML semantics: window = [target - lower, target + upper]
    pa.field("isolation_upper_offset", pa.float32()),
    pa.field("precursor_charge", pa.int8()),                  # instrument charge, null if unknown/0
    pa.field("precursor_intensity", pa.float32()),
    pa.field("master_scan", pa.int32()),                      # scan of the precursor (MS1 for MS2), null if unknown
    pa.field("sps_mz", pa.list_(pa.float64())),               # SPS ions (MS3), null otherwise
    # acquisition
    pa.field("activation", pa.uint8()),                       # ACTIVATION code; null for MS1
    pa.field("supplemental_activation", pa.uint8()),          # ACTIVATION code (EThcD: activation=ETD, supplemental=HCD)
    pa.field("nce", pa.float32()),                            # collision energy as reported (NCE for Thermo HCD, eV for Bruker)
    pa.field("supplemental_energy", pa.float32()),
    pa.field("analyzer", pa.uint8(), nullable=False),         # ANALYZER code
    pa.field("resolution", pa.float32()),                     # nominal resolving power setting (Thermo trailer)
    pa.field("faims_cv", pa.float32()),                       # FAIMS compensation voltage, null if FAIMS off/absent
    pa.field("polarity", pa.int8(), nullable=False),          # +1 / -1
    pa.field("filter_string", pa.string()),                   # Thermo scan filter, null for other vendors
    pa.field("scan_window_lower", pa.float32()),
    pa.field("scan_window_upper", pa.float32()),
    pa.field("tic", pa.float32()),
    pa.field("peak_source", pa.uint8(), nullable=False),      # PEAK_SOURCE code
    # Bruker / ion mobility (null for non-IM data)
    pa.field("window_index", pa.int16()),                     # DIA: index into the run's window table (metadata); null otherwise
    pa.field("im_scan_begin", pa.uint16()),                   # TIMS scan range of this row [begin, end)
    pa.field("im_scan_end", pa.uint16()),
    pa.field("im", pa.float32()),                             # precursor/window-centre 1/K0 (V s cm^-2)
    # peaks
    pa.field("mz", pa.list_(pa.float32()), nullable=False),
    pa.field("intensity", pa.list_(pa.float32()), nullable=False),
    pa.field("im_index", pa.list_(pa.uint16())),              # per-peak TIMS scan index (lookup 1/K0 via run metadata im_lut); null for non-IM
]

SPECTRUM_SCHEMA = pa.schema(SPECTRUM_FIELDS)
FIELD_NAMES = [f.name for f in SPECTRUM_FIELDS]
PEAK_LIST_FIELDS = ("mz", "intensity", "im_index")
NULLABLE = {f.name for f in SPECTRUM_FIELDS if f.nullable}

# run-level metadata keys (JSON under schema.metadata[b"vulcan"])
RUN_META_REQUIRED = (
    "store_version", "run_name", "run_id", "vendor", "instrument_model", "source_path",
    "source_size", "converter", "converter_version", "created_utc", "n_spectra",
    "n_ms1", "n_ms2", "n_peaks", "peak_order", "acquisition",
)


def default_run_id(run_name: str) -> int:
    """Deterministic 32-bit run id from the run name.

    Only a default for single-run work: the repository catalog (L6) must assign
    collision-free run ids (32-bit hashes collide with p~1% at 10^4 runs).
    """
    return int.from_bytes(hashlib.sha1(run_name.encode()).digest()[:4], "little")


def spectrum_id(run_id: int, index):
    """Global u64 spectrum id; ``index`` may be an int or a numpy array."""
    return (int(run_id) << 32) | index


def make_schema(run_meta: dict) -> pa.Schema:
    missing = [k for k in RUN_META_REQUIRED if k not in run_meta]
    if missing:
        raise ValueError(f"run metadata missing keys: {missing}")
    meta = dict(run_meta)
    meta.setdefault("activation_codes", ACTIVATION)
    meta.setdefault("analyzer_codes", ANALYZER)
    meta.setdefault("peak_source_codes", PEAK_SOURCE)
    return SPECTRUM_SCHEMA.with_metadata({b"vulcan": json.dumps(meta, default=str).encode()})


def run_meta_of(schema: pa.Schema) -> dict:
    return json.loads(schema.metadata[b"vulcan"])
