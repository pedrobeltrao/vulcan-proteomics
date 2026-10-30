# Vulcan spectral store, format version 1 (L0)

Status: implemented in `vulcan/io/` (schema in `schema.py`, writer and reader in `store.py`). DESIGN.md §4 L0 is the requirement. E1 measurements are in `results/2026-10-26_E1_ingest/`.

## 1. Container

- **One file per run:** `<run>.vspec.arrow`. It is an Arrow IPC **file**, not a stream, so it has a footer and can be read randomly by record batch.
- **Compression:** uncompressed by default. `--compression lz4|zstd` is optional; IPC compresses per buffer. E1 measured only about 4% saving for lz4 on f32 peaks.
- **Rows:** one row per spectrum, in acquisition order (`index` 0..n-1).
- **Peaks:** stored as list columns (`mz`, `intensity`, optional `im_index`).
- **Offsets:** int32, one set per record batch ("chunked offsets", DESIGN §4). A batch holds a converter-chosen 10³–10⁴ spectra and always fewer than 2³¹ peaks:
  - Thermo direct: 4096 scans per batch;
  - mzML: one 32 MB block;
  - Bruker: 1024 frames.
- **Flattened read:** `store.peaks_numpy()` returns contiguous f32 `mz` and `intensity` plus int64 `offsets` (n+1). It costs one memcpy per batch.
- **Static run metadata:** JSON under the schema-metadata key `vulcan`, written first.
- **Final counts:** `n_spectra`, `n_ms1`, `n_ms2`, `n_peaks`, `n_by_level` and converter timings go to a sidecar `<file>.json`, because pyarrow fixes the footer metadata when the writer opens. Counts can always be recomputed from the file.
- **Atomic write:** the file is written as `<file>.partial` and renamed on close, so a crashed conversion never leaves a valid-looking store.
- **I/O policy (DESIGN):**
  - Stores are staged to local `/tmp` and read sequentially. `open_store()` memory-maps only paths under `/tmp/`; on Lustre it uses a buffered sequential `OSFile`.
  - **Never mmap on Lustre.** The same rule applies to the RAW input (§5.1).

## 2. Spectrum table

Columns may be nullable (N); null means "not applicable or not reported". Codes are u8 enums, and their dictionaries are also stored in the run metadata.

| Group | Field | Type | Semantics |
|---|---|---|---|
| Identity | `spectrum_id` | u64 | `(run_id << 32) \| index`. It is globally unique if `run_id` is. |
| | `run_id` | u32 | Assigned by the repository catalog (L6). The default is the first 4 bytes of sha1(run_name) little-endian. **Collision risk is about 1% at 10⁴ runs, so the catalog must assign it.** |
| | `index` | u32 | 0-based row number in the run |
| | `scan` | u32 | Thermo: native scan number. Bruker: TDF frame id. mzML: `scan=` from the nativeID. |
| | `ms_level` | u8 | 1, 2, 3, ... |
| Timing | `rt_min` | f64 | Scan start time, in minutes |
| | `injection_time_ms` | f32 N | Thermo trailer "Ion Injection Time (ms)" / mzML MS:1000927 |
| Precursor (null for MS1) | `precursor_mz` | f64 N | Best selected-ion m/z: `mono_mz` if the instrument reported one, else `isolation_target_mz`. On the mzML path it is the mzML selected-ion m/z (MS:1000744). See §5.3. |
| | `isolation_target_mz` | f64 N | Isolation window target (Thermo `Reaction.PrecursorMass` + `IsolationWidthOffset`) |
| | `mono_mz` | f64 N | Thermo trailer "Monoisotopic M/Z" (null if 0 or absent) |
| | `isolation_lower_offset`, `isolation_upper_offset` | f32 N | mzML semantics: the window is [target − lower, target + upper] |
| | `precursor_charge` | i8 N | Instrument-assigned charge; null if 0 or unknown |
| | `precursor_intensity` | f32 N | mzML MS:1000042 if present (null on the Thermo direct path) |
| | `master_scan` | i32 N | Precursor scan, from trailer "Master Scan Number". If that is absent (e.g. QE), the most recent scan of level ms_level−1 in acquisition order. |
| | `sps_mz` | list\<f64\> N | SPS ions (trailer "SPS Masses" + "SPS Masses Continued") |
| Acquisition | `activation` | u8 N | 1 CID, 2 HCD (beam-type CID), 3 ETD, 4 ECD, 5 MPD, 6 PQD, 7 UVPD, 8 NETD, 9 PTR, 255 other. Bruker TIMS CID → 1. |
| | `supplemental_activation` | u8 N | Same codes (EThcD = ETD + supplemental HCD) |
| | `nce` | f32 N | Collision energy as the instrument reports it: normalized CE for Thermo HCD/CID, eV for Bruker |
| | `supplemental_energy` | f32 N | |
| | `analyzer` | u8 | 0 unknown, 1 Orbitrap, 2 ion trap, 3 TOF, 4 Astral, 5 quadrupole, 6 sector |
| | `resolution` | f32 N | Thermo trailer "FT Resolution" / "Orbitrap Resolution" (the setting) |
| | `faims_cv` | f32 N | Scan-filter CV (`cv=`) or trailer "FAIMS CV" when FAIMS is on |
| | `polarity` | i8 | +1 / −1 |
| | `filter_string` | utf8 N | Thermo scan filter (direct and TRFP-mzML paths) |
| | `scan_window_lower`, `scan_window_upper` | f32 N | Scan range |
| | `tic` | f32 N | Total ion current (Bruker: sum of the row's intensities) |
| | `peak_source` | u8 | 0 Thermo centroid stream; 1 Thermo scan acquired centroided; 2 RawFileReader centroider on a profile scan; 3 mzML as given; 4 Bruker raw TDF events; 5 timsrust-centroided |
| Bruker / IM | `window_index` | i16 N | diaPASEF: index into the run's `dia_windows` table |
| | `im_scan_begin`, `im_scan_end` | u16 N | TIMS scan range [begin, end) of this row |
| | `im` | f32 N | 1/K0 at the window centre (DIA) or precursor (DDA spectra mode) |
| Peaks | `mz` | list\<f32\> | Peak m/z |
| | `intensity` | list\<f32\> | Peak intensity |
| | `im_index` | list\<u16\> N | Per-peak TIMS scan index; 1/K0 = `im_lut[im_index]` (run metadata). Null for non-IM data. |

**Peak order** is given by the run-metadata key `peak_order`:
- `mz_ascending`: Thermo direct and timsrust spectra;
- `as_given`: mzML;
- `im_scan_then_mz`: Bruker frames. Peaks are grouped by TIMS scan, and the TOF index (hence m/z) ascends within a scan.

The GPU view (K0) does its own binning, so no consumer may assume m/z order without checking `peak_order`.

**Run metadata (JSON):**
- `store_version`, `run_name`, `run_id`, `vendor`, `instrument_model`, `instrument_serial`, `instrument_software`;
- `acquisition_date`, `source_path`, `source_size`, `source_sha1` (mzML path, from TRFP), `converter`, `converter_version`;
- `acquisition` (DDA/DIA-PASEF/unknown), `peak_order`, `compression`, `include_reference_peaks`;
- code dictionaries;
- Bruker only: `tims_converters` (TOF→m/z and scan→1/K0 formulas and parameters), `im_lut`, `dia_windows`, `events_outside_windows`;
- sidecar only: counts, timings (`timing_s` wall, `cpu_s` process CPU).

**What is deliberately not here:**
- No peak preprocessing. Top-k, deisotoping and precursor removal are GPU search-time views (DESIGN K0).
- No separate TMT reporter block. Every centroid is kept, so the reporter region is intact and the reporter block becomes a K0 view taken before top-k.
- No noise, baseline or per-peak charge arrays. They are a possible v2 addition for Thermo, from the centroid stream.

## 3. Invariants (tested in `tests/io/`)

1. `spectrum_id == (run_id << 32) | index`, and `index == arange(n)`.
2. `len(mz) == len(intensity)` per row. If `im_index` is non-null, it has the same length.
3. MS1 rows have null precursor and activation fields. MS2 rows on Thermo have non-null `isolation_target_mz` and `activation`.
4. Counts in the sidecar equal the counts recomputed from the file.
5. Thermo direct vs TRFP-mzML+parser for the same RAW are identical in:
   - every peak (m/z bit-equal after f32 cast, intensity bit-equal);
   - `scan`, `ms_level`, `rt`, `charge`, `master_scan`, activation, NCE, isolation, injection time, analyzer, polarity and filter string.

   The only field that differs is `precursor_mz` semantics (§5.3).

## 4. Converters

`python -m vulcan.io.convert <in> <out> [--via direct|mzml] [--compression lz4] [--bruker-mode frames|frames_pyo3|spectra] [--stage-dir DIR] [--run-id N]`

| Input | Path | Implementation |
|---|---|---|
| Thermo .raw | **direct (default)** | `thermo_direct.py`: pythonnet hosts .NET 8 in-process and loads `VulcanRaw.dll` (C#, `vulcan/io/dotnet/VulcanRaw`). That wraps Thermo RawFileReader 8.0.6, the DLLs shipped with TRFP 2.0 in `conda_envs/trfp/bin` (Thermo licence; not redistributed). One .NET call returns a columnar batch of 4096 scans; Python copies each array with one memmove. Peaks come from the centroid stream, including reference/exception peaks (the TRFP default). |
| Thermo .raw | `--via mzml` | TRFP 2.0 → indexed mzML (zlib, 64-bit) on local disk → `mzml.py` |
| .mzML | `mzml.py` | Streaming 32 MB blocks; each field is found with one regex scan per block and assigned with `searchsorted`; base64/zlib decoding runs in C. TRFP and msconvert attribute orders are handled; numpress is rejected. |
| Bruker .d | `bruker.py` frames (default) | `vulcan_tims` (Rust, `vulcan/io/rust/vulcan_tims`: timsrust 0.4.2 + numpy) decodes 1024 frames per call into flat numpy arrays; rows are built vectorised. `frames_pyo3` is the same output via `timsrust_pyo3` (Python lists, slower). |
| Bruker .d | `--bruker-mode spectra` | timsrust `SpectrumReader`: IM-collapsed centroided MS2 (DDA: per precursor; DIA: per frame × window). MS1 is not emitted. |

**Build helpers.**
- `scripts/ingest/build_vulcanraw.sh` builds the C# helper into `tools/vulcanraw/`. It needs the .NET 8 SDK; `dotnet-install.sh` puts it in `tools/dotnet-sdk/`.
- `scripts/ingest/build_vulcan_tims.sh` builds the Rust extension into `tools/vulcan_tims/`.
- Both outputs are gitignored under `tools/`.

**.NET hosting.**
- `thermo_direct.load_runtime()` uses the SDK's shared runtime if present. Otherwise it synthesises a `DOTNET_ROOT` from TRFP's self-contained runtime: it symlinks the runtime files and copies TRFP's `deps.json` as `Microsoft.NETCore.App.deps.json`, because hostfxr ignores a framework directory without one.
- Globalization is invariant with `PredefinedCulturesOnly=false`. The node has no libicu, and RawFileReader's filter `ToString()` builds `CultureInfo("en-US")`.

## 5. Vendor notes and known deviations

### 5.1 Thermo RAW on Lustre: prefetch is mandatory

RawFileReader memory-maps the RAW and touches pages in a non-sequential order. With a cold page cache on Lustre, every 4 KiB fault is a round trip.

| Case | Decode wall | Total wall | CPU |
|---|---|---|---|
| Velos 335 MB, cold | 53 s | 65 s | about 5 core-s |
| QE HF-X 960 MB, cold | — | about 330 s | — |
| Velos 335 MB, after sequential prefetch | 1.2 s | 6.9 s | — |

`thermo_direct.convert` therefore runs a sequential `posix_fadvise(SEQUENTIAL)` + `readv` pass first (`prefetch`). Alternatively, stage the RAW to local disk (`--stage-dir`). Both cost about 1 sys core-s per GB.

### 5.2 TRFP Parquet is not a store substitute

TRFP 2.0 `-f 3` writes "mzparquet", one row per peak, with columns `scan`, `level`, `rt`, `mz`, `intensity`, `ion_mobility`, `isolation_lower`/`upper`, `precursor_scan`, `precursor_mz` and `precursor_charge`.

- **Missing:** activation, NCE, analyzer, filter string, injection time, FAIMS CV and resolution. It therefore fails the L0 metadata contract.
- **Cost:** it is cheaper than mzML for TRFP but far more expensive than the direct path (E1).

### 5.3 Precursor m/z semantics

TRFP's mzML selected-ion m/z follows its own rule.

- On PXD001468 01A it equals the isolation target for 86.5% of MS2 and the trailer monoisotopic m/z for 54.9%; 45% of MS2 have no trailer mono.
- The direct path's `precursor_mz` (mono if present) differs from TRFP's value for 32% of MS2, by up to 6 Th. These are isotope-peak picks that the instrument corrected.

Both raw ingredients (`isolation_target_mz`, `mono_mz`) are stored, so L0.5 or the search can choose. For Sage parity on mzML, use the mzML path or reproduce TRFP's rule from those fields.

### 5.4 Bruker diaPASEF: events outside windows are not stored

In MS2 frames, the store keeps one row per quadrupole window (the scan range [begin, end) from `DiaFrameMsMsWindows` / timsrust `QuadrupoleSettings`).

- Events in scans covered by no window have no defined isolation. They are counted in `events_outside_windows` but not stored. On A9_R1 that was **35.2M of 335.6M events (10.5% of all; 18.6% of MS2 events)**, the only deviation from "nothing discarded".
- The deviation is reversible by re-ingest. An alternative, if they are needed (e.g. for background modelling), is to store whole MS2 frames as one row each; `im_index` would then still allow window assignment.

### 5.5 Other Bruker notes

- **Conversions:** TOF→m/z and scan→1/K0 use timsrust's boundary-based converters, not Bruker's full calibration polynomials. Recalibration belongs to the L0.5 pilot. TOF indices are recoverable from m/z with the stored parameters.
- **ddaPASEF** frames are handled by the same code path (timsrust's quadrupole settings give per-precursor scan ranges) but are **untested**: no DDA .d file was available. The precursor id per row is not yet linked; it is a TODO for DDA, via `PasefFrameMsMsInfo`.

### 5.6 Thermo direct, not yet exercised

- **Peak-source fallbacks.** If a scan has no centroid stream (ion-trap analyzer), the scan's own centroids are used (`peak_source=1`), or RawFileReader's centroider for profile scans (2). The PXD001077 "Velos" file turned out to be HCD-OT, so ITMS has not been exercised yet; it is a TODO with an ion-trap MS2 dataset.
- **FAIMS, SPS and EThcD** are implemented from the RawFileReader API but have not been exercised on data. No FAIMS, SPS or EThcD files are in Tier S yet.
