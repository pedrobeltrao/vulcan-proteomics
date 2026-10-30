# E11: GPU MS1 feature finder (opt-in `ms1_backend=gpu`; dev only: D0, D1, D8h)

## Bottom line
- **Bitwise exact.** The GPU backend (hill detection + isotope grouping) gives tables identical to the CPU path on all
  13 dev runs (D0 2, D1 1, D8h 10). That covers every Hills array, the feature table, `peak_feature`/`peak_iso`, params
  and the hypothesis table, and the stored E7 tables. No stage needed a "deterministic but inexact" fallback.
  - So the ID-set Jaccard is 1.0 and the v2_sep ID ratio change is 0 by construction.
  - In the 96-run cohort, the 12 base runs are bitwise equal to E7b's CPU cohort on all 5 output tables. The 84 renamed
    copies differ only by run identity.
- **Per run (E7b protocol, 8 threads, full chimeric v2_sep path, 13 runs):**
  - MS1: 90.8 → 13.8 core-s (6.6×) and 49.4 → 14.6 s wall; the MS1 kernels use 0.80 GPU-s in total.
  - Whole pipeline: 116.7 → 40.0 core-s (2.9×).
  - Per D0 run: MS1 12.4 → 1.8 core-s; total 14.3 → 4.4 core-s and 8.3 → 4.5 s wall.
  - Counts (features, hypotheses, PSMs, IDs, protein groups) are identical.
- **Cohort (96 runs, RAW → protein groups, same session, back to back):**
  - Throughput: 3,307 → **3,647 runs/h** including setup (1.10×); steady state 5,718 → **11,892 runs/h** (2.08×).
  - Cost: 16.8 → 11.5 core-s/run; MS1 7.9 → 1.85 core-s/run.
  - The cohort's own GPU memory peaked at 48 GB of 80.
- **Opt-in only.** The frozen/default path is unchanged: `backend="cpu"` is the default everywhere, `--ms1-backend`
  defaults to `cpu`, and FeatureParams and the table metadata are unchanged. The CPU path stays the reference.

Nothing sealed was touched: no T1–T5 data, no `bench/sealed_run.py`, no frozen model or its training. Only dev stores
were used (`/tmp/vulcan_t1_<uid>/{D0,D1,D8h}`, the E7b cohort RAWs). Node: gpu-node (A100 80GB, 32 affinity cores).
It was shared with other jobs throughout (load 8–42), and every record carries `hw` and the load averages.

## What was built
| file | role |
|---|---|
| `vulcan/ms1/ms1_kernels.cu` | CUDA kernels + host orchestration (CUB scans/radix sorts, thrust sort). Session API over ctypes; static cudart, **no torch**, so CPU decode workers can load it (LESSONS 4). |
| `scripts/e11/build_ms1_kernels.sh` | nvcc `-fmad=false`, sm_80 + sm_86 + PTX. It is auto-run by `vulcan/ms1/gpu.py` when the source hash changes. Output goes to `tools/vulcan_ms1_kernels/` (gitignored). |
| `vulcan/ms1/gpu.py` | Binding: `Session.hills/isotopes/isotopes_raw/peak_map`, `finish_isotopes` (CPU greedy), the memory check, and a cross-process gate. |
| `vulcan/ms1/features.py` | `find_features(..., backend="cpu"\|"gpu", hill_peaks=True)`. `gpu_supported()` covers gap 0, full-norm and max_iso ≤ 16 (the production config); anything else runs the CPU path. |
| `vulcan/cohort.py` | `--ms1-backend {cpu,gpu}` (default cpu) and `--ms1-gpu-slots N` (default 2). Each decode worker gets its own CUDA context (0.4 GiB). |
| `tests/ms1/test_gpu_exact.py` | 15 test cases: synthetic runs (4 seeds), raw isotope arrays, `SCAN_BIN=2`, 7 parameter variants (m/z window, min_intensity, no ban, max_iso 4, tolerances, min_overlap 3, charges), gap-1 fallback, empty run, hypotheses windows, `hill_peaks=False`, gated early release, memory fallback, and that the CPU default is untouched. |
| `scripts/e11/*` | Benchmark and equality scripts: `equal_gpu.py`, `bench_run.py`, `run_bench.py`, `cohort_bench.py`, `equal_cohort.py`, `gpu_mem.py`, `make_metrics.py`. |

### Where each step runs and why it is exact
**On the GPU:**
- Linking.
- Labels by list ranking.
- Counting-sort placement.
- Valley split.
- Length filter and compaction.
- Hill summaries.
- m/z argsort.
- Per-hill norms.
- Scan-block bins.
- First-isotope pre-screen.
- Isotope chains with the averagine prefix.
- `nn ≥ 2` filter and charge ban.
- Priority order.
- Per-peak feature map (in the single-process mode).

**On the CPU:**
- The sequential greedy acceptance, using the same numba kernel.
- The per-feature table.
- Hypotheses.
- The store load.

**Why it is exact:**
- **Floating-point reductions.** Each one runs sequentially in one GPU thread, in the CPU's loop order. Numba's typing is replicated: f32·f32 and f32+f32 round to f32 before widening, e.g. the trace dot products, the norms and the trapezoid area. FMA contraction is off (`-fmad=false`), and sqrt/div are IEEE on both sides.
- **Link claims.** The one-to-one link claim is (min distance, then min index), done with two `atomicMin` passes. This is the CPU's ascending strict-< scan.
- **Order-dependent integer steps** are order-exact:
  - Peak placement uses chain depth, which equals index order.
  - The argsort and the block bucketing use CUB's stable radix sort.
  - Candidates come out in row-major `nonzero` order via a scan.
  - The priority order is a total order (ties end at (mono, charge)), so any correct sort reproduces `np.lexsort`.
- **Tie rule.** The explicit isotope tie rule (best cosine, then lowest global m/z rank) is kept.

## Equivalence evidence
- **`equality_gpu.json`, final build: `all_ok` true, 13/13 runs.** For each run, `find_features(backend="gpu")` was compared with the CPU default. Equal:
  - all 12 Hills arrays (bytes + dtype);
  - the feature table (every column + schema metadata), `peak_feature`, `peak_iso` and `params`;
  - the `build_hypotheses` table (every column + metadata);
  - both tables against E7's stored files (`<work>/<ds>/ms1/`), 13/13. E9_0001 has no MS1 peaks, and its empty tables are equal.
- **`equality_cohort_gpu.json`: 96-run GPU cohort vs E7b's CPU cohort `base12`.**
  - The **12 base-named runs are bitwise equal on all 5 tables** (features, hypotheses, chimeric PSMs incl. score, peptides incl. q, proteins) and on all metadata.
  - In the **84 renamed copies (`_rNN`)**, the only differences are run identity, and that is verified per run (`identity_only`). Each copy gets its own `run_id` from store meta, so `spectrum_id = (run_id << 32) | index` changes in its high bits. Everything else is equal: the low 32 bits, every other column, and the hypothesis metadata once `run_id` is removed.
  - This is the same situation as E7b's copies check, which excluded the `run` column.
- **Bench counts.** `counts_equal` is true in both benches: n_features, n_hypotheses, n_psms, NCE, `ids_q_sep_1pct` and protein groups per run are identical CPU vs GPU.
- **Tests.** `pytest tests/ms1 tests/gpu/test_cohort_load.py`: 35 passed.

## Per-run cost (`bench_{cpu,gpu}_t8.json`, `bench_{cpu,gpu}_t1_ms1.json`)
**Setup:** the protocol of `scripts/e7b/bench_run.py` ("after"). One warm process per backend, 1 warm-up + 3 repeats per run, each run's repeats under `flock /tmp/vulcan_bench.lock`, run back to back. The GPU backend runs in production mode (`hill_peaks=False`).

**Sums over the 13 runs, 8 threads, full chimeric path:**

| stage | wall s CPU → GPU | core-s CPU → GPU | × core-s |
|---|---|---|---|
| **MS1 features + hypotheses** | 49.4 → 14.6 | 90.8 → 13.8 (GPU-s 0.80) | **6.6** |
| search | 4.5 → 4.7 | 4.4 → 4.7 | – |
| features | 10.3 → 8.8 | 10.1 → 8.7 | – |
| model | 1.0 → 1.3 | 1.0 → 1.3 | – |
| compete | 3.1 → 3.1 | 3.1 → 3.1 | – |
| proteins | 6.4 → 7.9 | 6.4 → 7.9 | – |
| **total** | **75.5 → 41.5** | **116.7 → 40.0** | **2.9** |

- E7b's `bench_after` (same protocol) had MS1 92.4 and total 116.8 core-s, which matches this CPU run.
- The non-MS1 stages are the same code in both arms. Their ±20% swings are node load: load average was 8→42 during the CPU run and 41→30 during the GPU run.

**MS1 sub-stages, summed wall over 13 runs (1 thread, MS1 only, the decode-worker setting):**

| sub-stage | CPU (numba, 1 thr) | GPU |
|---|---|---|
| store load | 3.97 | 3.95 |
| hills: link + split/filter + summaries | 10.87 + 9.74 + 1.27 = 21.9 | 0.84 (incl. H2D) + 0.28 D2H |
| isotopes | 51.3 | 2.77 = 0.43 GPU + D2H, then **2.0 greedy** (CPU) |
| feature table | 5.69 | 3.30 (GPU peak map) |
| **MS1 core-s total** | **83.3** | **12.7** (6.6×) |

- **Per D0 run (mean of 2):**
  - 1 thread: MS1 13.1 s / 12.6 core-s → 2.0 s / 1.8 core-s.
  - GPU kernel time per D0 run is about 0.13 s (hills 85–105 ms, isotopes 30–40 ms, peak map 5 ms).
- **Largest run (E9_0000, 79.9M MS1 peaks):** 0.92 GPU-s.
- **CPU waits use blocking-sync events**, so host waits do not spin a core.

## Cohort throughput (`cohort_chim_{gpu,cpu}.json`, `cohort_runs.jsonl`; `scripts/e11/cohort_bench.py`)
**Setup:** the E7b design: 10 D8 + 2 D0 RAWs cycled to 96 runs, human_contam_prod, the 3i predicted library (GPU-resident), 24 decode workers, 3 post workers, `--chimeric`. Each cohort ran under `flock /tmp/vulcan_bench.lock`, GPU first, then CPU.

**The /tmp free-space guard:**
- The driver refuses to start below 50 GB free and kills the process group below 20 GB.
- It also samples GPU memory per process, split into this cohort vs others.

**One difference from E7b: RAWs are decoded in place (`--no-stage-raw`) in both arms.**
- The RAWs already sit on local NVMe.
- With staging, the first decode wave copies ~24 RAWs (~55 GB). On this /tmp (58–75 GB free, shared), that tripped the guard at 13 s.

| config | runs/h (incl. setup) | steady runs/h | core-s/run | MS1 core-s/run | MS1 wall/run (median) | GPU-s/run (search+score) | MS1 GPU-s/run | cohort GPU mem max |
|---|---|---|---|---|---|---|---|---|
| chimeric, **MS1 GPU** (2 slots) | **3,647** | **11,892** | **11.5** | **1.85** | 1.64 s | 0.22 | 0.14 | 48.4 GB |
| chimeric, MS1 CPU (E7b default) | 3,307 | 5,718 | 16.8 | 7.88 | 5.69 s | 0.29 | – | 38.3 GB |

- **Steady-state throughput doubles** (2.08×) and core-s/run falls 1.46×.
- **The incl.-setup figure gains only 10%.** At 96 runs the ~95–105 s makespan is dominated by setup and tail:
  - Scorer-ready comes at 32 s with the GPU backend vs 21 s with the CPU one, because 24 worker CUDA contexts initialise while the 25 GB library warms up.
  - The decode pool is only 35% busy.
  - A 192-run cohort, which would show this better, was attempted for both arms. The /tmp guard stopped both at 89 s and 110 s, because other jobs had left only 58 GB free. That is recorded in `cohort_runs.jsonl`.
- **Earlier GPU cohort variants (`prelim/`, same equality result):**
  - Holding the GPU slot through the CPU greedy pass and feature table made workers wait 592 s in total on 2 slots. Throughput was 3,167 runs/h (steady 6,734).
  - With 3 slots it was 3,213 runs/h (steady 10,032).
  - The final code frees the device right after the isotope arrays are copied out (only `p_src` is fetched) whenever a cross-process gate is active. Slot waits fell to 36 s in total.
  - In single-process use (bench) the peak map stays on the GPU.
- **GPU memory** (`gpu_mem_worker.json`):
  - Per worker: 0.42 GiB idle context, 4.4 GiB peak on an 80M-peak run, 3.4 GiB on D0.
  - A run starts on the GPU only if `cudaMemGetInfo` shows ≥ 96 B/peak + 256 MiB free. Otherwise, or on a CUDA allocation error, it takes the bitwise-equal CPU path, recorded as `cpu_fallback`. There were 0 fallbacks in the cohorts.
  - Node-wide used memory once read 79.4 GB, but that included other jobs: 24 GB was in use before our cohort allocated anything. The per-process split shows our share peaking at 48–51 GB.

## Caveats
- **The node was shared and loaded** (load 8–42), so absolute times are noisy. Core-s comparisons hold up better than wall.
- **The per-run benches used the build before the memory check and early-release changes.** Neither affects the single-process path, apart from one `cudaMemGetInfo` call. The equality record and all tests were rerun on the final build.
- **The cohort comparison is not E7b's configuration** (`--no-stage-raw`). Compare the two arms here with each other, not with E7b's 2,987 / 4,643 runs/h.
- **Worktree runs used E7b's Python env** (`/tmp/vulcan_env_<uid>/pred`) plus symlinks into the main repo's gitignored built artifacts (`tools/vulcan_core`, `vulcan_kernels`, `vulcanraw`, `dotnet*`, `conda_envs/trfp`, `data/cache/peptdeep`). These were read-only use, and no source in them changed.

## What remains
1. **The MS1 stage is no longer the bottleneck.** About 1.85 core-s per run is left, mostly on the CPU:
   - the sequential greedy acceptance (≈0.15–0.25 s per D0 run);
   - the per-feature table (≈0.25–0.45 s);
   - the store load (≈0.3–0.6 s);
   - hypotheses (≈0.23 s).
   A GPU greedy would need conflict-free rounds with an exactness proof; it is not attempted here.
2. **Per D0 run the chimeric pipeline is now ≈4.4 core-s.** The main costs are now decode (5.0–5.6 core-s/run in the cohort), proteins and features.
3. **Cohort startup.** Initialising the worker CUDA contexts slows library warm-up by about 11 s. Staggering context creation, or creating it at first use after the scorer is ready, should recover most of the incl.-setup gap.
4. **Steady-state throughput at ≥192 runs** still needs a node with more free /tmp, or backpressure on decode-ahead (stores pile up while the GPU stage waits).
5. **Landing.** This is a post-freeze speed path. Before any production default switch, it needs review and a decision entry. The CPU path stays the reference, and nothing here changes A3-frozen behaviour.
