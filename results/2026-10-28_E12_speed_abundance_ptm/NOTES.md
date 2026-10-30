# E12: speed for protein abundance (LFQ) and PTMs (phospho), exact speedups only (dev data only)

## Bottom line
- **Every change is either byte/bitwise-identical in output, or opt-in with the default off.** All new paths are
  opt-in flags; no frozen default changed. No frozen model was touched, and no sealed set (T1–T5) or
  `bench/sealed_run.py` was read.
- **Equivalence evidence.** All checks are on dev data.
  - 96-run cohorts, fast vs default: all output tables are bitwise equal in closed, LFQ front end and chimeric.
    - Tables compared: PSMs incl. score, peptides incl. q, proteins, MS1 features, hypotheses. That is 288 or 480
      files per comparison, in 4 rounds.
    - The outputs are also identical across rounds/sessions.
  - Index build: sha256-identical files at full scale (184.5M peptidoforms).
  - directLFQ: sha256-identical protein and ion tables on D9 and D11.
  - D9 quant metrics JSON: identical, MBR off and on.
  - Localization peak prep and summaries: bitwise identical.
- **Headline (96-run cohorts, round 4, back to back, same node session, load 33–40; `metrics.json` `headline`):**

| workflow | runs/h before → after | core-s/run before → after | vs Sage (core-s/run) | equivalence |
|---|---|---|---|---|
| closed (search + static_gbm_pred + picked proteins) | 4,570 → **5,753** (+26%) | 9.9 → 9.1 | Sage 70 (4 thr) – 120 (32 thr), + 250 for RAW→mzML (M1e, D0) | 288/288 tables bitwise |
| LFQ, MBR off (cohort + D9 quant) | 561 → **1,280** (2.3×)* | 19.7 → **15.3** | Sage search+LFQ 266 (D9, 32 thr, mzML given); FlashLFQ alone 153 (M6) | cohort 480/480 bitwise; directLFQ/MaxLFQ outputs identical |
| LFQ, MBR on | 572 → **1,194** (2.1×)* | 19.6 → **15.5** | (Sage has no MBR FDR) | as above |
| chimeric (v2_sep) | 2,567 → **3,867** (+51%) | 16.2 → 14.9 | Sage chimera not timed (≥ closed) | 480/480 bitwise |
| phospho + L3 localization + Ala-decoy FLR | 2,362 → **3,172** (+34%)* | 10.7 → 9.8 | Sage phospho on D0: 3,087 core-s / 2 runs incl. 188 s DB build (§6; no localization) | peaks/summaries bitwise; scoring: see §4 |

\* LFQ and phospho runs/h are a conservative estimate. They assume the experiment-level step runs **serially after**
the cohort: makespan + 96 × the per-run quant (or post-search) wall. The cohort parts are measured.
MSFragger: there is no measured number in `results/`. It is unlicensed, and M1e/M4 cite literature only.

- **Biggest wins:**
  1. **directLFQ**: 22–28 s → 4.7–5.0 s per D9 experiment (whole run 21.1 → 2.9 s), byte-identical.
     - It was 2/3 of LFQ quant, so the D9 quant step goes from 32.9 → 12.5 s per experiment (5.5 → 2.1 s per run).
  2. **Phospho index build**: 444 s → **81 s** (locked), byte-identical.
     - The old code took 549 s in this session. Arrays fell 416 → 15.6 s and labels 58 → 6.9 s.
  3. **Cohort start and GPU-stage parallelism** (`--fast-setup --search-threads 2`): +26% / +28% / +51% runs/h
     (closed / LFQ / chimeric).
     - Scorer ready at 18–23 s instead of 37–52 s.
  4. **Localization host prep**: 2.5 → 0.19 s wall per D2 batch (prepare_peaks), summarize 0.63 → 0.21 s.
- **New finding (affects T2 reproducibility; for the PI):**
  - The frozen L3 localization scoring is **not bit-reproducible run to run**: its float64 GPU `index_add_`
    atomics vary the summation order. On PXD009449, 313 of 45,845 PSMs change their best configuration between
    two identical reference runs (every one an exact score tie broken by last-bit noise).
  - An opt-in `deterministic=True` mode is reproducible bit for bit. It also implements the documented tie rule
    (lowest configuration index) exactly.
  - It differs from the reference by about as much as the reference differs from itself in best configurations
    (342–351 PSMs). But it accepts **1.7% fewer** mono sites at FLR̂ ≤ 1% on PXD009449 (40,747 vs 41,433–41,460):
    with exact ties, tied sites get delta = 0 exactly instead of random ±1e-14.
  - PXD000138 and PXD007058 counts are unchanged.
  - Not adopted. It needs a decision before any T2 claim, see §4.

## 1. Setup
- **Node and data:**
  - Node: gpu-node (32 Zen 4 cores, A100 80GB). It was shared, at load average 18–55 throughout.
  - Every record carries the hw fingerprint and the load.
  - Worktree branch `worktree-agent-a0df6277cef40ec85`, with gitignored tools, envs and caches symlinked
    read-only from the main repo.
  - Cohort python: `/tmp/vulcan_env_<uid>/pred` (as E11).
- **Cohorts (`scripts/e12/cohort_bench.py`):**
  - Design is E11's: 10 D8 + 2 D0 QE RAWs cycled to 96 runs, human_contam_prod, the 3i predicted library, 24 decode
    workers, 3 post workers, RAWs decoded in place.
  - Each config runs under `flock /tmp/vulcan_bench.lock`. The page cache is prewarmed (untimed) before each config.
  - **Stores and outputs live in /dev/shm.** Decode run-ahead therefore cannot fill the shared /tmp.
  - /tmp still hit 100% at 12:55 from other jobs (A5/T2 look dry runs, A6). My /tmp footprint was < 10 MB. The guard
    killed one of my configs (round 3 `chim`, recorded with rc −15).
- **Configs:**
  - `closed`: production closed path.
  - `lfq`: closed + `--ms1 --ms1-backend gpu`, the LFQ front end (MS1 features + hypotheses).
  - `chim`: `--chimeric --ms1-backend gpu`.
  - `_fast` adds `--fast-setup`; `_st2` adds `--search-threads 2`.
- **Quant (`scripts/e12/quant_bench.sh`, `quant_bench2.sh`):**
  - D9, 6 HYE Astral runs, `scripts/m6d/lfq_eval.py d9 all m6d [--mbr]`.
  - 2 locked repeats per configuration, best of repeats per stage.
  - Production stages only: link, align, MBR, XIC, siblings, MaxLFQ, directLFQ.
- **Phospho (`scripts/e12/phospho_bench.py`):**
  - M4 phospho PSMs (D2: PXD000138 dev pools, PXD009449, PXD007058; 9 runs).
  - Stages: index open → rank-1 → labels/decode → peptidoform FDR → spectra → localization (S/T/Y + Ala) → FLR curve.
  - One process, ref and fast alternately, 2 repeats; the second (warm) repeat is reported. Locked.

## 2. Per-stage profile, before → after (core-s per run unless noted)

**Cohort (round 4, per run; `round4/cohort_*.json`):**

| stage | closed | LFQ front end | chimeric |
|---|---|---|---|
| decode (.NET RawFileReader + store write + search prep) | 7.65 → 6.96 | 6.84 → 6.41 | 6.73 → 5.71 |
| MS1 features + hypotheses (GPU) | – | 3.12 → 3.06 | 4.05 → 4.58 |
| decode-worker CPU outside jobs | 1.36 → 1.24 | 3.30 → 2.79 | 3.72 → 2.90 |
| GPU process (search + scoring + loaders) | 0.91 → 0.90 | 1.00 → 0.97 | 1.71 → 1.71 |
| GPU-s (kernels) | 0.11 → 0.20 | 0.14 → 0.24 | 0.27 → 0.43 |
| scorer ready (s from start) | 36.9 → 22.8 | 39.1 → 18.4 | 52.4 → 21.6 |
| makespan, 96 runs (s) | 75.6 → 60.1 | 89.3 → 69.9 | 134.6 → 89.4 |

How to read it:
- The per-run core-s changes come from node noise, not code. The setup changes do not alter per-run work.
- GPU-s/run rises with 2 GPU-stage threads because the event timers overlap, not because there is more work.

**Setup (unchanged code paths, now overlapped or deferred):**
- Predictor (AlphaPeptDeep) load: 11–22 s on the critical path → 0.05–0.3 s.
  - It now loads in a background thread. A fully covered library never waits for it.
- Worker CUDA contexts: their creation inflated main `cuda_init` to 4.5–6.7 s → 0.2 s (created at the first GPU MS1
  run, inside the slot gate).
- Library warm: 25 GB H2D in 8–29 s, mostly unchanged.
  - The 24-way thread pool did not help.
  - A pinned double-buffered staging upload (round 2) was slower (13 s) and was **reverted**.
  - Page-cache state dominates.

**What limits a cohort now:**
- After setup, the critical path is the **single GPU-stage thread's host time**: search 0.09 + scoring 0.21 s per
  run closed; 0.29 + 0.32 s chimeric.
- It degrades under CPU contention: all decodes finished at ~42 s in a 89 s chimeric run.
- A second GPU-stage thread (`--search-threads 2`, an existing knob) recovers this. It is exact (bitwise tables) and
  gives +43% on chimeric in round 3 (99.6 → 69.5 s).
- The cohort's own "bottleneck: decode_pool" label is misleading at 96 runs.

**LFQ quant, D9 experiment (6 runs), s (single-threaded; wall ≈ core-s):**

| stage | MBR off ref → fast | MBR on ref → fast |
|---|---|---|
| link PSMs→features | 1.88 → 2.56 | 2.03 → 2.02 |
| RT alignment | 0.20 → 0.21 | 0.21 → 0.21 |
| MBR | – | 1.15 → 2.28 (noise; same code) |
| XIC re-quant | 1.36 → 1.77 | 1.68 → 1.26 |
| charge siblings | 1.01 → 1.18 | 1.42 → 1.01 |
| MaxLFQ (`impl="fast"`) | 2.77 → 1.73 | 3.12 → 1.99 |
| **directLFQ (`impl="fast"`)** | **25.67 → 5.05** | **22.56 → 4.95** |
| **total** | **32.9 → 12.5** | **32.2 → 13.7** |

Notes on the quant numbers:
- The first `quant_*` run was cold. XIC took 11.5 s because the MS1 peaks were read from a cold page cache.
- Warm numbers are reported throughout.

**Phospho post-search, 9 D2 runs, warm repeat (wall / core-s, s):**

| stage | ref | fast |
|---|---|---|
| L1 index open (184.5M rows) | 4.13 / 4.12 | 0.00 / 0.00 (mmap; pages read on use) |
| labels + decode rank-1 | 0.67 / 0.67 | 0.74 / 0.74 |
| peptidoform FDR (paired TDC) | 0.46 / 0.46 | 0.47 / 0.47 |
| spectra | 0.33 / 0.33 | 0.34 / 0.34 |
| localization: prepare_peaks | 2.51 / 2.51 | **0.19 / 1.66** (numba, 32 thr) |
| localization: score_configs (GPU) | 1.58 / 1.58 | 1.84 / 1.84 (deterministic sum) |
| localization: summarize | 0.63 / 0.63 | 0.21 / 0.21 |
| FLR curve | 0.23 / 0.23 | 0.56 / 0.56 (noise; same code) |
| **total** | **10.8 / 10.9** | **4.8 / 6.3** |

- The search itself is 0.2 s per run on the phospho index (M4; as fast as closed).
- The index build is a one-off: 444 s → 81 s.
- prepare_peaks scales with threads: 1 thread 2.38 s, 32 threads 0.32 s, on 50k synthetic spectra; ref 3.33 s.
- The first call compiles the kernel (about 8 s, cached on disk).

## 3. What was built (all opt-in; defaults unchanged)

| change | file(s) | flag | exactness |
|---|---|---|---|
| fast arrays/labels step of the L1 index build | `vulcan/index/build.py`, new `vulcan/index/_fastbuild.py` | `build(impl="fast")`, `VULCAN_INDEX_BUILD_IMPL=fast`, CLI `--impl fast` (default stays `reference`) | 6 files sha256-identical at full scale; `tests/index/test_build_fast_exact.py` (10 tests) |
| memory-mapped L1 index load | `vulcan/index/build.py` | `L1Index(path, mmap=True)` | arrays equal (test) |
| per-protein directLFQ on numpy | new `vulcan/quant/dlfq_fast.py`, `lfq.py` | `directlfq_rollup(impl="fast", num_cores=1)` | protein + ion TSV sha256-identical (D9, D11); `tests/quant/test_dlfq_fast.py` |
| MaxLFQ positional rows | `vulcan/quant/lfq.py` | `maxlfq(impl="fast")` | bitwise (test + D9 metrics) |
| prepare_peaks as a parallel numba kernel | `vulcan/ptm/localize.py` | `prepare_peaks_fast` | bitwise, both weight modes, ties, NaN, float32, empty (`tests/ptm/test_localize_fast.py`) |
| segmented summarize for k = 1 | `vulcan/ptm/localize.py` | `summarize_fast` | bitwise |
| deterministic localization score sum | `vulcan/ptm/localize.py` | `score_configs(deterministic=True)` | reproducible; **not** equal to the reference (which is not reproducible), §4 |
| cohort fast setup | `vulcan/cohort.py`, `vulcan/rescore/prod.py`, `cohort_scoring.py` | `--fast-setup` | 96-run tables bitwise, 4 rounds |
| bench / equality / metrics scripts | `scripts/e12/*`, `scripts/m6d/lfq_eval.py` (`--dlfq-impl`, `--maxlfq-impl`) | – | – |

What `--fast-setup` changes:
- A lazy AlphaPeptDeep proxy (`_LazyAPD`): the cache key needs only the model sha.
- A thread-pooled library warm.
- Worker CUDA contexts created at first use.

Tests: tests/index 14, tests/ptm 16, tests/quant 17, tests/gpu/test_cohort_load + tests/a3/test_a3_freeze 6 — all pass.

## 4. Localization reproducibility (finding)
- **The cause.** `score_configs` sums matched-peak weights per configuration with `torch.index_add_` on float64 CUDA
  tensors. Atomic accumulation order varies, so configurations that tie exactly (same matched peak set) get scores
  that differ in the last bits. "First maximum" then picks a random member of the tie, and delta_site of a tied
  best site becomes a random ±1e-14 instead of 0.
- **The reference is not reproducible** (`phospho_repro_unlocked.json`; ref run 2 vs ref run 1, same inputs):

  | dataset | PSMs | different best configuration | mono sites accepted at FLR̂ ≤ 1% |
  |---|---|---|---|
  | PXD009449 | 45,845 | 313 | 41,460 vs 41,433 |
  | PXD000138 | 3,858 | 16 | equal (2,949) |
  | PXD007058 | 506 | 1 | equal |

- **What `deterministic=True` changes.** Sums happen in a fixed order on the host (`np.bincount` over the sorted
  (row, peak) pairs). Tied configurations then get identical scores, and the documented tie rule (lowest
  configuration index, which depends on positions only) applies. The mode is residue-blind as before.
- **Fast vs reference:** 342–351 PSMs differ in best configuration on PXD009449, about the same as reference vs
  reference.
  - The accepted count at FLR̂ ≤ 1% is 40,747 vs 41,433–41,460 (−1.7%), from the exact-zero deltas of tied sites
    near the threshold.
  - PXD000138 and PXD007058 are unchanged.
- **Recommendation (not done here):** the PI should decide before any T2 claim whether the frozen T2 path should
  use the deterministic sum. Then G0/FLR calibration should be re-checked on dev. The T2 look as frozen would
  otherwise give non-reproducible site lists (≈0.7% of PSMs on Ymod).

## 5. Exactness details
- **Index build** (`index_build.json`):
  - The fast path replaces np.unique (252 s), lexsort (103 s), gathers and Python string columns with numba kernels:
    - a bucketed parallel sort equal to `np.lexsort((pf_id, mass32))`, since pf_id is unique, so the order is total;
    - a parallel uniqueness check;
    - parallel gathers;
    - Arrow string columns built from buffers, which reproduce `pa.array`'s 2³¹−2-byte chunking.
  - The new files are byte-identical to both the M4 index and a same-session old-code rebuild.
- **directLFQ** (`dlfq_fast.json`):
  - The per-protein estimation is re-run on numpy slices with directLFQ's own kernels: get_normfacts, merge_distribs,
    numba nanmedian/nanvar.
  - It is installed only inside a context manager.
  - Ion and protein tables are identical. Per protein, 5,135 (D9) and 5,300 (D11) proteins agree bitwise.
  - With 8 cores: wall 1.1 s, 4.6 core-s.
- **Cohort:** setup changes cannot alter any per-run computation. That was verified per table, including NaN-aware
  float bit compare (`scripts/e12/equal_dirs.py`), on every config in rounds 1–4. Round 3 vs round 4 default closed
  outputs are identical too.

## 6. Comparators
- **Sage 0.14.7** (from earlier locked records):
  - closed D0: 46 s / 140 core-s per 2 runs at 4 threads; 34 s / 239 at 32 threads. RAW→mzML (TRFP) adds 501 core-s
    per 2 runs (M1e).
  - search + LFQ on D9: 17.6 s / 266 core-s per Astral run at 32 threads (M6, mzML given).
  - FlashLFQ: 72 s / 153 core-s per run (M6).
- **Sage phospho** (`scripts/e12/sage_phospho.py`): same search space as d2_phospho_prod, on the 2 D0 runs,
  32 threads, unlocked, load recorded. Result (`timings.jsonl` `E12_sage_phospho_D0_t32`, load ≈ 38):
  - 386 s wall / **3,087 core-s for 2 runs**, peak RSS 123 GB.
  - Its own log: fragment-index build 188 s (5.18G fragments, 184.6M peptides, the same size as our index), total
    260 s.
  - Per run that is ~1,540 core-s with the build, or ~36 s wall of search per run without it.
  - Sage does no site localization.
  - Ours: 9.8 core-s per run end to end (decode + search + scoring + FDR + localization + FLR), plus a one-off
    81 s index build (444 s before).
- **On the comparison itself:**
  - Our per-run numbers come from QE runs (D8/D0 mix) for the cohort part and D9 Astral for the quant part. Sage's
    LFQ number is Astral (123k MS2 per run), so it is not like-for-like on data. The core-s gap (15 vs 266) is far
    larger than that difference.
  - Accuracy: E12 changes no output, so the accuracy comparisons of M6/M6d (LFQ), M4 (localization), E7b (chimeric)
    and D17 (IDs) carry over unchanged. For localization this holds for the reference path; the deterministic mode
    is opt-in, see §4.

## 7. Caveats and next targets
- **The node was loaded** (load 18–55; other agents were running gates and looks). Absolute runs/h vary ±25% between
  rounds. Compare before/after within a round:
  - round 1: closed 6,381 → 7,439, LFQ 5,913 → 6,752, chim 3,876 → 3,689 (fast-setup only);
  - round 2: staged-warm variant, reverted, mixed;
  - round 3: closed 5,197 → 6,237 (fast) → 6,872 (fast+st2), chim_fast 3,469 → 4,976 (st2);
  - round 4: the headline.
  - `fast-setup` alone does not always win: when the GPU stage starts earlier it competes with decode, and the
    GPU-stage thread is the critical path. Combine it with `--search-threads 2`.
- **Recommended production config** (opt-in until reviewed): `--fast-setup --search-threads 2`, plus the fast
  quant/index/localization opt-ins. All are bitwise-identical except `deterministic`.
- **Next bottlenecks:**
  1. Decode, 6–7 core-s/run (vendor .NET reader).
  2. GPU-stage host time, 0.3–0.6 s/run of Python scoring/search glue. CUDA graphs, or moving features to the
     loader thread.
  3. Decode-worker CPU outside jobs, 1.2–3 core-s/run, not yet attributed. A thread sample shows .NET tiered-JIT
     and threadpool threads ≈ 0.7 core-s/run (`diag/`).
  4. GPU MS1 host remainder, 3 core-s/run.
  5. Link/XIC/siblings quant, about 1 s per run, Python.
- **Not done (time):** a GPU localization kernel; a quiet-node headline session; LFQ cohort on D9 Astral RAWs
  end-to-end (the predicted library covers only human_contam_prod).
