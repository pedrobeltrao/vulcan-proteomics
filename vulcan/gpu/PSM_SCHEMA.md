# vulcan.gpu output schemas (PSM table v1, family top-2 table v1)

Producers: `vulcan.gpu.engine.psm_table` (used by `engine.search`, `python -m vulcan.gpu.run`, `vulcan.gpu.batch`,
`python -m vulcan.cohort`) and `vulcan.gpu.family.family_table` (`vulcan.cohort --family-top2`). One Parquet file per run:
`<run>.psms.parquet`, `<run>.family_top2.parquet`. Neither table carries a target/decoy/entrapment label; join
`pf_id` to the index's `labels.arrow` in L4 only (DESIGN §3).

M2b (2026-10-27) did **not** change the PSM table: same columns, same meaning, same order. Additions go at the end
and are listed in the change log below.

## 1. PSM table (`<run>.psms.parquet`)

One row per (spectrum, rank ≤ `report_psms`). Candidates with `matched_b + matched_y < min_matched_peaks` are not
ranked. Ranking is per spectrum over all its precursor hypotheses by hyperscore, descending, ties to the lower
`pf_row`. It is label-blind.

| Column | Type | Meaning |
|---|---|---|
| `run` | string | Run name (store `run_name`) |
| `spec_idx` | i64 | 0-based index over the run's MS2 spectra in store order |
| `scan` | u32 | Native scan number |
| `rt` | f64 | Retention time, minutes |
| `rank` | i32 | 1 = best |
| `pf_row` | i64 | Row in the L1 index (`peptidoforms.arrow`). This is index-specific. Use `pf_id` across indexes. |
| `pf_id` | u64 | Peptidoform id: xxhash64 of the canonical peptidoform. The join key to `labels.arrow`. |
| `charge` | i8 | Precursor charge of the hypothesis |
| `isotope_error` | i8 | Isotope error of the hypothesis |
| `expmass` | f32 | (precursor m/z − proton)·z, Sage f32 arithmetic |
| `calcmass` | f32 | Peptidoform mass (Sage f32) |
| `delta_mass_ppm` | f32 | Sage's symmetric ppm mass error (isotope-corrected) |
| `hyperscore` | f64 | Sage hyperscore (K2) |
| `delta_next` | f64 | hyperscore − next candidate's hyperscore in the same spectrum, regardless of label (equals hyperscore if there is no other candidate) |
| `delta_best` | f64 | best hyperscore in the spectrum − this one |
| `matched_b`, `matched_y`, `matched_peaks` | i16 | Matched fragment ions (b, y, sum) |
| `summed_intensity` | f32 | Matched intensity (b + y) |
| `matched_intensity_pct` | f32 | 100 · summed_intensity / TIC of the processed peaks (the K0 view) |
| `average_ppm` | f32 | Intensity-weighted mean fragment ppm error |
| `longest_b`, `longest_y` | i16 | Longest consecutive b / y ladder |
| `longest_y_pct` | f32 | longest_y / peptide length |
| `peptide_len` | i16 | Peptide length |
| `missed_cleavages` | u8 | Missed cleavages of the forward target peptide (a pair attribute) |
| `n_var` | u8 | Variable modifications |
| `k1_count` | i32 | K1 preliminary count (−1 if the row entered only through pair closure) |
| `scored_candidates` | i32 | Candidates with K1 count ≥ 1 over all hypotheses of the spectrum (Sage semantics) |
| `poisson` | f64 | log10 Poisson probability (Sage) |

`scored_candidates` and `poisson` depend on the whole searched database, so an entrapment-augmented search gives
different values than a production search.

### 1b. M4 additions (mass-offset / open modes, `vulcan.ptm.search.search_modes`)

Closed-mode PSM tables (`engine.search`) are unchanged. Tables written by `search_modes` carry all §1 columns with
the same meaning, plus, at the end:

| Column | Type | Meaning |
|---|---|---|
| `mode` | u8 | Hypothesis kind of the PSM: 0 closed (spectrum x charge x isotope error), 1 mass-offset, 2 open |
| `hyp_offset` | f32 | Nominal offset (Da) of a mass-offset hypothesis; 0 for closed; NaN for open |
| `mass_shift` | f32 | expmass − isotope_error·NEUTRON − calcmass (Da), i.e. the unexplained delta |

For mode 1/2 PSMs `hyperscore`, `matched_*`, `summed_intensity` come from K2's localization-aware shifted scoring
(the delta is placed on the residue that maximises the hyperscore); `average_ppm`, `longest_b/y` use the union of
shifted and unshifted matches. `delta_mass_ppm` keeps its §1 definition (so it is large for shifted PSMs).

**Variable-mod indexes (M4).** `peptidoforms.arrow:mask` bit 63 (`vulcan_core.NTERM_BIT`) marks a variable
peptide-N-term modification (Sage `^`, e.g. acetyl). Residue bits stay 0..len−1. Consumers that decode masks must
ignore/handle bit 63; it only occurs in indexes built with a `^` variable mod.

## 2. Family top-2 table (`<run>.family_top2.parquet`, `vulcan.gpu.family`)

This table is written for a multi-replicate index (`vulcan.index.build.build_multi_decoy`: T, E_s, dec_j(T),
dec_j(E_s)). Each row is one (spectrum, retention group, `grank` ∈ {1, 2}).

| Column | Type | Meaning |
|---|---|---|
| `run`, `spec_idx`, `scan` | | As in §1 |
| `fgroup` | i16 | Opaque retention group ((f·J) + j)·2 + side. f: 0 = T family, 1 + i = entrapment seed i (sorted); j: decoy context; side 0 = forward, 1 = decoy. Use `FamilyIndex.replicate_groups(j, i)`. |
| `grank` | i8 | 1 or 2 within (spectrum, fgroup), by hyperscore desc then `pf_row` asc |
| `pf_row`, `pf_id` | | As in §1 |
| `hyperscore` | f64 | As in §1 |
| `charge`, `isotope_error` | i8 | Canonical hypothesis: the lowest-index hypothesis whose window contains the row |
| `matched_b`, `matched_y` | i16 | As in §1 |

**Reconstructing replicate (j, s).** Use `family.reconstruct(tbl, fi.replicate_groups(j, s))`: per spectrum, take the
top-2 by (hyperscore desc, `pf_row` asc) over the four groups. This equals a separate search of that replicate under
the per-family top-K rule (`scripts/m2b/family_validate.py`, exact). `delta_next` per replicate is top1 − top2.

## Change log
- v1 (E2, 2026-10-27): PSM table as above.
- M2b (2026-10-27): no PSM column changed. Added the family top-2 table (new file, new producer).
- M4 (2026-10-27): no existing column changed. `search_modes` tables add `mode`, `hyp_offset`, `mass_shift` (§1b);
  mask bit 63 = peptide-N-term var mod. Closed search is bit-identical to E2 (checked on D0 01A).
