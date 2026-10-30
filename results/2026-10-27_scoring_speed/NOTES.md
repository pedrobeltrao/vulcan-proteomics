# Scoring speed: the frozen static_gbm_pred production path, made faster with identical outputs (dev only)

**Scope.**
- The scorer is frozen for T1 (DESIGN §6.0 A1.1). The model, the features, the NCE grid and rule, the competition and the picked-protein estimators are unchanged, and no model was retrained.
- Data: D0 (2 runs), D1 (1 run) and the D8 holdout (10 runs; M3e stores and RAWs). Nothing sealed was touched.
- Files:
  - `metrics.json`: summary.
  - `profile_before.json`: frozen path, per-sub-stage timers.
  - `bench_after_library.json`: new path, GPU-resident library. `bench_after_library_host.json`: host-resident library. `bench_after_firstcontact.json`: new code, M3e cache state.
  - `gpuprof_{frozen,fast}.json`: CUDA kernel seconds.
  - `equality.json`, plus `equality_mmap_placement.json` and `equality_reference_cache.json`.
  - `predictor_precision.json`, `library_build.json`.
  - `cohort_*.json`, `cohort_summary.json`.
  - `config.yaml`, `git_sha`.
- Code:
  - `vulcan/predict/library.py`: consolidated, memory-mapped / host / GPU-resident predicted library with the PredCache lookup semantics.
  - `vulcan/rescore/fast.py`: GPU peaks, the once-per-run NCE-independent block, and `GpuHGB`.
  - `vulcan/rescore/post.py`: torch-free competition and proteins.
  - `vulcan/rescore/prod.py`, `vulcan/rescore/cohort_scoring.py`.
  - `vulcan/cohort.py`: `--score-lib`.
  - `scripts/scoring_speed/`.
  - `tests/rescore/test_fast_equal.py`: 5 tests.

## Bottom line
- **Per run, scoring is 7.2× faster in wall time and 8.7× in core-s.** This is the production configuration in a persistent process, under the lock, median of 12 measured runs (10 D8 + D0 run 2 + D1).
  - Wall: 3.89 s → 0.54 s. Core-s: 4.60 → 0.53. Max: 11.3 s → 1.7 s. Sum over the 12 runs: 48.0 s → 7.95 s (6.0×).
  - GPU kernel time: 1.04 → 0.087 GPU-s per run.
- **Outputs are bit-identical, so no float tolerance was needed.** On all 13 D0/D1/D8 runs:
  - all 10 feature columns are bitwise equal to the stored reference features (M3/M3e), and to the frozen FeatureBuilder on the same library;
  - the GPU GBM equals sklearn bitwise, and the scores are bitwise equal;
  - the NCE choice, the competition winners and q, and the q≤1% (and 0.1%, 5%) ID sets are identical;
  - the picked-protein counts are identical.

  In the cohort, RAW→PSM tables are identical to the reference, the `score` column is bitwise equal, and the q≤1% peptide sets and protein counts are identical for every base run.
- **End to end inside the cohort (RAW → PSMs → scored peptides → picked protein groups):** 192 runs (D8 + D0 cycled), 24 decode workers, node loadavg 60–70.

  | config | runs/h | steady state |
  |---|---|---|
  | **fast scoring** | **5,658** | 10,491 |
  | frozen scoring | 642 | |
  | search only (no scoring) | 7,485 | 11,989 |

  - That is 8.8× the frozen-scoring throughput.
  - GPU-thread capacity with scoring is 7,913 runs/h, against 789 runs/h for the frozen components.
  - Cost per run: 6.77 core-s with fast scoring, 6.18 without scoring, 11.57 with frozen scoring.
- **The ≥5× depends on the precomputed library.** The new code on the M3e "first-contact" cache (only NCE 29 shipped, other NCEs predicted lazily) gives 1.57 s median, i.e. 2.5×. In that state, prediction of misses dominates, and it can only be made faster by predicting ahead.

## Per sub-stage (median wall per run, s; before = frozen path, M3e state; after = new path, GPU-resident library)

| sub-stage | before | after | × |
|---|---|---|---|
| predictor inference RT+MS2 (incl. NCE grid misses) | 1.67 | 0 (all hits; z≥5 lazy) | — |
| NCE selection (6 NCEs × ≤2000 PSMs) | 1.34 | 0.051 | 26 |
| fragment matching + spectral features, chosen NCE (lookup + match) | 0.57 | 0.065 | 8.7 |
| RT lookup + calibration | 0.14 | 0.010 | 14 |
| peak table | 0.17 | 0.016 | 11 |
| **features total** | 2.48 | 0.158 | 15.7 |
| GBM inference (incl. feature matrix) | 0.58 | 0.035 | 16.7 |
| double competition (estimator, unchanged) | 0.047 | 0.048 | 1 |
| picked protein (estimator, unchanged) | 0.30 | 0.30 | 1 |
| **scoring total** | **3.89** (4.60 core-s) | **0.54** (0.53 core-s) | **7.2** (8.7 core) |
| I/O: store load (not scoring) | 0.63 | 0.59 | page-cache state |

- Picked protein (`picked_group_fdr`, an estimator) is now 55% of inline scoring. In the cohort it runs in the post-worker processes, off the GPU thread: 0.50 s wall and 0.47 core-s per run.
- With a host-resident library instead of GPU, scoring is 0.60 s (features 0.19 s).
- `end_to_end_total` single-process (1.23 s after) is inflated by the store load (0.59 s), which is page-cache eviction on the node and not part of scoring.

## What changed (all outputs identical by construction and by test)
1. **Predicted library, content-hashed per (model, instrument, NCE, mods), shared across runs and processes.**
   - `scripts/scoring_speed/build_library.py`:
     - copies the M3 cache root, which is the one the reference features were computed from; it is append-only and first-copy-wins, so every reference value is preserved;
     - runs the frozen `PredCache.precompute` (misses only, fp32, batch 16384) for every row of the human and *P. furiosus* prod indexes, for RT and MS2 at z = 2, 3, 4 and all 6 grid NCEs;
     - consolidates each table (stable merge, first copy wins = `PredCache` lookup semantics).
   - Build cost: 2,555 s one-time, run unlocked on the shared GPU. Disk: 50 GB (24.9 GB consolidated + source shards).
   - At run time the library is mirrored in GPU memory (24.9 GB; lookups by `searchsorted` + gather on the GPU). It loads in 5.7 s from the page cache and 45 s cold.
   - z ≥ 5 (~1% of PSMs) stays lazy. All charges' misses go in one predictor call, and misses are appended as ordinary shards.
2. **NCE choice:** the same per-run label-blind rule.
   - The NCE-independent half (fragment m/z, observed intensities, slot mask) is computed once per run, not 6 times.
   - Each grid NCE's spectral-angle block keeps the frozen shapes (≤2000 rows), so the medians are bitwise equal.
   - Per-instrument or per-project caching of the NCE was **not** applied: the frozen pipeline chooses per run, and caching would change semantics.
3. **Features:**
   - The peak table is built on the GPU (float64 key, stable sort). TIC and m/z min/max are computed with a numba sequential sum, equal to `np.bincount`.
   - Predictions are scattered on the GPU, and the features come back in one D2H copy.
   - The frozen chunking (200k rows) is kept, so every reduction sees the same shapes.
4. **GBM:** `GpuHGB` evaluates sklearn's 150 HistGradientBoosting trees on the GPU in float64. It uses sklearn's split rule (NaN → missing_go_to_left, else `x <= thr`) and sklearn's summation order (zeros + baseline, then tree by tree), so the output is bitwise identical; this is tested on random data, on values at the thresholds, and on the frozen model.
5. **Cohort integration** (`python -m vulcan.cohort ... --params '{"report_psms": 5}' --score-lib DIR`):
   - After each search, the GPU-stage thread runs features + GBM and adds a `score` column to the PSM table.
   - Double competition + picked proteins run in torch-free post-worker processes (default 2; 3 here), which write `<run>.peptides.parquet` and `<run>.proteins.parquet`.
   - The scorer setup overlaps the first decodes.
   - `report_psms` must be 5 (the frozen scorer's PSM table); the driver refuses anything else.
   - `--score-impl frozen` runs the old components for comparison.

## Rejected or not applied
- **Reduced-precision predictor (autocast fp16 / bf16):** 19% / 24% of stored MS2 values change (max 0.26 / 0.65) and RT changes. Rejected; the predictor stays fp32.
- **Larger predictor batches (65k, 262k):** outputs identical but no faster (49–51k vs 57k MS2/s). Batch 16384 is kept.
- **The predictor is not batch-composition invariant:** in a side check (30k peptides predicted whole vs in 7 random sub-batches), 1,058 of 2.44M stored fp16 MS2 values differed by 1 ulp and RT by up to 4.9e-4. The frozen cache has the same property. Values therefore depend on cache history. This is why the library is built on top of the reference cache and never rebuilt from scratch.
- **RT keyed independently of NCE:** would change values, not applied. **Per-project NCE:** see above.
- **Two GPU-stage threads:** 6,196 vs 5,658 runs/h, within noise; each run's GPU-thread time doubles (GIL). Kept at 1.

## Measurement notes (caveats)
- **Setup:**
  - Timed sections ran under `flock /tmp/vulcan_bench.lock`: before profile, two after benches, two gpuprofs, each cohort.
  - Node loadavg was 33–90 on 32 affinity cores (other agents: sage, ray jobs, E9).
  - The job cgroup sat at about 245/256 GB, with heavy page-cache refaults.
  - Absolute cohort runs/h are therefore load-dependent: the same search-only config gave 2,513 and 7,485 runs/h 30 min apart. Compare configs measured back to back (`fast` then `noscore2`).
- **Host-resident library (24.9 GB anon):** under this cgroup pressure, decode was slower in the cohort (ingest wall 3× core). The default is therefore GPU placement. `--score-placement host|mmap` remain.
- **GPU memory:** the cohort process peaks at 31–34 GB (library 24.9 GB) on a shared A100. Use `--score-placement host` when GPU memory is scarce.
- **Cohort sample:** 192 runs is 16 copies of 12 RAWs. From the second copy on, the frozen config's lazy cache is warm, so its 642 runs/h is a best case. The GBM on the host also competes with the decode workers for CPU: 1.4 s per run vs 0.58 s alone.
- **First cohort attempts** (superseded, kept in `cohort_*_hostlib.json` / `/tmp`): one used the M2b default `report_psms 1`, which is not the frozen scorer's table, so its IDs differed. Another built the scorer before starting the decode feeder. Both were fixed before the reported runs.
- **GPU-s:** CUDA kernel time from torch.profiler in a separate pass. The per-run GPU-s from cohort events covers search only.
- **hwinfo:** in `metrics.json` (`hw`; AMD EPYC 9354, 32 affinity cores, A100 80GB).
