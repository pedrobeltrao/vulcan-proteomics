# vulcan.ms1 output schemas (features v0.2, hypotheses v1.0)

These tables are what the L2 search consumes (DESIGN §4 L2). The Python schemas are the source of truth:
`FEATURE_SCHEMA` in `features.py` and `HYPOTHESIS_SCHEMA` in `hypotheses.py`. `tests/ms1/test_hypotheses.py`
pins the hypothesis schema, so any change has to update this file and bump `HYPOTHESIS_VERSION`.

Both tables are one Arrow IPC **file** per run, written atomically (`.partial`, then rename). Producer:
`python -m vulcan.ms1 STORE OUTDIR [--calib CALIB.json]`, which writes:
- `<run>.features.arrow`
- `<run>.hypotheses.arrow`
- `<run>.ms1.json` (timings and statistics)

Neither table holds a decoy, pair or entrapment label, and nothing that produces them reads one (DESIGN §3).

## 1. Conventions

| Name | Value |
|---|---|
| PROTON | 1.007276466812 Da |
| NEUTRON (isotope spacing) | 1.0033548378 Da (13C − 12C) |
| neutral mass | (m/z − PROTON) · z |
| isotope error e | The observed precursor m/z is the M+e peak. The searched peptide mass is (m/z_cal − PROTON)·z − e·NEUTRON. This is the same convention as Sage and `vulcan.gpu.engine.make_hypotheses`. |
| RT | minutes (store `rt_min`) |
| MS1 ordinal | 0-based index over the run's MS1 spectra in acquisition order |

## 2. Feature table (`FEATURE_SCHEMA`, one row per isotope-envelope feature)

| Field | Type | Meaning |
|---|---|---|
| `feature_id` | u32 | 0..n−1 in (apex MS1 ordinal, mono m/z) order. Assigned from the data, so it is stable for identical input and params. |
| `run_id` | u32 | store `run_id` |
| `mono_mz` | f64 | Intensity-weighted mean m/z of the monoisotopic hill. Uncalibrated. |
| `charge` | i8 | 1..6 |
| `neutral_mass` | f64 | (mono_mz − PROTON) · charge |
| `rt_apex` | f32 | RT of the apex scan of the most intense isotope hill |
| `rt_start`, `rt_end` | f32 | RT bounds of the mono hill |
| `ms1_apex`, `ms1_start`, `ms1_end` | u32 | The same three positions as MS1 ordinals |
| `intensity_apex` | f32 | Sum of the isotope hills' apex intensities |
| `intensity_mono_apex` | f32 | Apex intensity of the mono hill. Equals biosaur2 `intensityApex` at `-iuse 0`. |
| `intensity_area` | f32 | Sum of the isotope hills' trapezoid areas (intensity × min) |
| `n_isotopes` | u8 | Accepted isotope peaks, mono included (≥ 2) |
| `iso_cos` | f32 | Cosine between the isotope apex intensities and the averagine envelope |
| `trace_cos` | f32 | Mean cosine of each isotope's elution trace against the mono trace |
| `n_scans` | u16 | Length of the mono hill in MS1 scans |

The schema metadata key `vulcan_ms1` holds the JSON parameters (`FeatureParams`), the version, the hill count and
the candidate count. `FeatureResult.peak_feature` and `peak_iso` map every flat MS1 centroid to its feature and
isotope index (−1 if none). The hypothesis builder uses that map.

## 3. Precursor-hypothesis table (`HYPOTHESIS_SCHEMA`)

There is one row per (MS2 spectrum, precursor group, isotope error). A *precursor group* is one
(source, feature, charge) explanation of the spectrum's precursor. Each group is expanded over
`isotope_errors` (default −1, 0, +1, +2, +3).

Rows are sorted by:
1. `spectrum_id`;
2. group priority: source, then chimeric rank (feature groups) or charge (other sources);
3. isotope error.

All hypotheses of a spectrum compete together in L4.

| Field | Type | Null | Meaning |
|---|---|---|---|
| `spectrum_id` | u64 | no | Store `spectrum_id` of the MS2 |
| `hypothesis_id` | u16 | no | 0..n−1 within the spectrum, in the row order above |
| `precursor_group` | u16 | no | 0..g−1 within the spectrum. Rows that differ only in isotope error share it. |
| `source` | u8 | no | 0 = `ms1_feature`, 1 = `instrument`, 2 = `fallback` (see §3.1) |
| `charge` | i8 | no | Precursor charge of the group |
| `isotope_error` | i8 | no | e, as defined in §1 |
| `precursor_mz` | f64 | no | Observed mono m/z of the group, uncalibrated. Feature groups use the feature's `mono_mz`; other groups use the store `precursor_mz` (instrument mono if reported, else the isolation target). |
| `precursor_mz_cal` | f64 | no | `precursor_mz / (1 + ppm(rt, mz)·1e−6)` with the run's precursor offset curve (`vulcan.calib`). Equal to `precursor_mz` if no calibration was given. |
| `neutral_mass` | f64 | no | (precursor_mz_cal − PROTON)·charge − isotope_error·NEUTRON. **This is the mass to search.** |
| `isotope_error_prior` | f32 | no | P(e \| source). Taken from the pilot pass when available (table metadata `priors`), otherwise from `DEFAULT_PRIORS`. Sums to 1 over a group. |
| `feature_id` | u32 | yes | MS1 feature (feature table) for `ms1_feature` groups |
| `feature_intensity` | f32 | yes | Feature `intensity_apex` |
| `feature_window_intensity` | f32 | yes | The feature's centroid intensity inside the isolation window at the master MS1 scan |
| `isolation_fraction` | f32 | yes | `feature_window_intensity` / total centroid intensity in the window (precursor purity) |
| `chimeric_rank` | u16 | yes | Label-free rank of the feature among **all** features (any charge) with peaks in the window, by `feature_window_intensity`. 1 is the dominant precursor. Ties are broken by `feature_id`. |
| `feature_iso_cos` | f32 | yes | Feature `iso_cos` |
| `feature_rt_apex` | f32 | yes | Feature `rt_apex` |
| `instrument_agrees` | bool | no | True for the one feature group whose charge equals the instrument charge and whose mono equals the instrument m/z up to an isotope offset (≤ 3) within 20 ppm |
| `instrument_isotope_offset` | i8 | yes | (instrument m/z − feature mono)·z / NEUTRON, rounded. Set only when `instrument_agrees`. +1 means the instrument reported the M+1 peak. |
| `n_features_in_window` | u16 | no | Features of any charge with ≥ 1 centroid in the window. Repeated on every row of the spectrum. |

The schema metadata key `vulcan_hypotheses` holds JSON with:
- `version` and `run_id`;
- `source_codes` and `isotope_errors`;
- `params` (`HypothesisParams`);
- `priors` per source;
- `calibration`: the `Calibration.summary()` that was applied, or null;
- `feature_params`.

### 3.1 Sources and rules (`HypothesisParams` defaults)

The isolation window is [target − lower, target + upper], using the store's isolation offsets (±1 Th if absent).
It is sliced from the MS2's master MS1 scan, or from the preceding MS1 if the store has no master scan.

1. **`ms1_feature`**: every feature with centroids in the window that meets all of:
   - charge 2–6;
   - `chimeric_rank` ≤ 5;
   - `isolation_fraction` ≥ 0.01.

   A charge-1 feature qualifies only if the instrument charge is also 1. These groups are the co-isolated
   precursors for the chimeric search: restrict to `chimeric_rank == 1` for a non-chimeric search.
2. **`instrument`**: the instrument charge is known and no feature group agrees with it. The group uses the
   instrument m/z and charge.
3. **`fallback`**: there is no instrument charge and no feature group. Charges 2, 3 and 4 are tried at the
   instrument m/z.

When the instrument agrees with a feature, only the feature group is emitted. Its mono (from MS1) replaces the
instrument's value, which corrects instrument mono picks: `instrument_isotope_offset ≠ 0`.

### 3.2 Consumer notes

- The search reads `spectrum_id`, `charge`, `neutral_mass` and `isotope_error` and needs nothing else. Precursor
  tolerance windows are applied to `neutral_mass` by the searcher, using the run's empirical tolerance
  (calibration JSON `tolerances.precursor_ppm_after`).
- For a non-chimeric search that is equivalent to using the instrument precursor (but with the MS1-corrected mono),
  keep groups with `source != 0 or instrument_agrees`, then filter `isotope_error` to the engine's range.
- Calibration is never baked into an index. The hypothesis table carries both the raw and the calibrated m/z, and
  the offsets live in `<run>.calib.json` (`vulcan.calib.calibration.Calibration`).
