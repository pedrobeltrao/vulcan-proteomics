# M2b: cohort driver, startup, family top-2, wire format, ingest (dev data only)

- **Data:** dev data only. The cohort list is `scripts/m2b/dev_dda_raws.txt`: all 35 dev Thermo DDA RAWs (D0 ×2, D1, D2 ×12, D3, D3b, D4 ×12, D7 ×6; 25.2 GB, mean 0.72 GB and 28.7 k MS2 per run). **To reach 50 or 200 runs the list is cycled.** Repeats get unique run names, so 200 runs means each file about 5.7 times. No sealed data was touched.
- **Search space:** every run is searched against the human prod index (`/tmp/vulcan_e2/index/human_contam_prod`, 7.06 M peptidoforms) with the E2 parameters. This emulates throughput only. It is not the right database for the Pfu, pool, XL or TMT files, so these runs have no ID meaning.
- **Node:** gpu-node, 32 Zen 4 affinity cores, A100 80GB (`metrics.json:hw`). The batch-scheduler job has a 256 GB cgroup memory limit.
- **Contention:** the node was **shared and noisy**. Other agents and users ran Comet, training jobs and others, and host busy cores averaged 11–34 during our runs. Every timed config held `flock /tmp/vulcan_bench.lock`, one hold per config of 30–125 s.
- **Page cache:** **RAW reads were page-cache warm** (repeats, earlier reads). Uncached Lustre was measured separately (§3).

Files:
- `metrics.json`: curated numbers, from `scripts/m2b/summarize.py`.
- `runs/*.json`: full per-run stats and timelines of the final code.
- `runs/v1_*`, `runs/v2_*`, `runs/v3_*`, `runs/discarded_lustre_mmap/`: earlier rounds. They are kept for the record, and every one is superseded.
- `config.yaml`.

## 1. Cohort driver (`python -m vulcan.cohort`, `vulcan/cohort.py`)

One long-lived process runs the pipeline:

- **Decode workers** are a supervised pool of N spawn processes, 1 thread each, kept alive for the whole cohort. Each worker loads the .NET runtime and the RawFileReader once. Per run, a worker:
  - copies the RAW from Lustre to local NVMe;
  - decodes it to the Arrow store on `/tmp`;
  - builds search-ready MS2 arrays and the precursor-hypothesis table (`engine.pre_hypotheses`) as `.npy` files.
- **The GPU process** builds the L1/K1 index once and keeps it resident. A loader thread memory-maps the `.npy` files. The GPU stage then runs K0 → K1 → top-K → K2 → features (`engine.search(pre_hyp=…)`).
- **Writer threads** write `<run>.psms.parquet` and delete the staged files.
- **Overlap:** decoding of run i+k overlaps the search of run i. Jobs are dispatched largest file first.
- **Supervision:** a job with no CPU progress for 30 s, or one past its timeout, gets its worker killed and respawned, and the job is retried (at most 2 retries).
- **Optional `--ms1`:** the vulcan.ms1 features and hypotheses run in the worker. Their output is written but not yet consumed by the search.
- **Optional `--family-top2`:** see §5.

**Identity (`ingest_identity.json`, `determinism_threads.json`):**
- Stores written by the cohort path are column-for-column identical to the E2 one-shot `vulcan.io.convert` stores for Velos, 01A and 01B. Only timing metadata differs.
- The PSM tables are identical to the E2 PSM tables.
- All 35 runs give identical PSM tables across three configurations: (a) E2-style, (b) 3 GPU threads, (c) worker-side prep.

### Throughput: 200-run cohorts, final code, makespan including setup

| decode workers | makespan | **runs/h** | runs/h, 10–90% steady* | core-s per run (all processes) | core-s per 1e5 MS2 | decode core-s per run | GPU busy (NVML mean) |
|---|---|---|---|---|---|---|---|
| 4 | 122.7 s | **5,870** | 7,380 | 2.52 | 8.8 | 2.33 | 5% |
| 8 | 59.5 s | **12,100** | 16,000 | 2.60 | 9.1 | 2.38 | 10% |
| 16 | 39.3 s | **18,340** | 27,800 | 3.04 | 10.6 | 2.73 | 16% |
| 24 (3 reps: c200, abB1, abB2) | 34.9 / 31.5 / 28.6 s | **20,600 / 22,900 / 25,100** | 37–42 k | 3.5–3.6 | 12.1–12.6 | 3.1–3.2 | 20–23% |
| 28 | 27.5 s | **26,230** | 46,200 | 3.78 | 13.2 | 3.33 | 25% |

\*The 10–90% steady rate is biased upward, because the largest files are dispatched first.

- **50-run cohorts** (`c50_*`, the requested size): 5.5 k / 9.8 k / 15.1 k / 12.8 k runs/h at 4 / 8 / 16 / 24 workers. At this size the makespan is latency-bound: 10–14 s set by the slowest single D0 decode (~8 s) plus about 2–4 s of setup. Time to first PSM table is 6–12 s.
- **Scaling:** close to linear up to 8 workers, about 0.8× per added worker up to 16, then flattening. Decode core-s per run grows from 2.3 to 3.3 between 4 and 28 workers. This reflects shared-node contention (memory bandwidth, other jobs), not idle time: the decode pool stays 70–99% busy.
- **Compared with E2 one-shot** (7.5 s and 4.9 core-s per D0 file from the store, plus about 8 core-s for a one-shot ingest): the cohort does RAW → PSM in **2.5–3.6 core-s per run** (dev mix) and **~0.06 GPU-s per run** (CUDA events, an upper bound).

### Where the bottleneck is: per-stage capacity at 24 workers (abB2; `metrics.json:configs.*.stage_capacity`)

| Stage | Parallelism | Mean s per run | Capacity, runs/h |
|---|---|---|---|
| **Decode workers** (RAW copy + decode + Arrow + prep) | 24 processes | 3.05 | **~28 k** (24–30 k across reps) |
| GPU stage (K0..features, serial) | 1 thread | 0.070 | 52 k |
| GPU busy (CUDA events; NVML says 20–25% utilized) | — | 0.069 | 52 k |
| Writer (parquet + cleanup) | 2 threads | 0.084 | 85 k |
| Loader (mmap of the worker-prepared arrays) | 2 threads | 0.007 | ≫ 100 k |
| *Lustre RAW read, uncached (measured, §3)* | 16 streams | 0.72 GB per run at 1.27 GB/s | **~6.5 k** |

- **With page-cache-warm RAWs, the CPU decode stage binds.** The GPU stage has roughly 2× headroom, and the GPU itself is about 75% idle.
- **For a real cohort streamed cold from this Lustre, the RAW read binds at about 6.5 k runs/h** for this file mix. That is 4× below decode capacity. To run thousands of PRIDE runs, the next lever is I/O: more read streams or nodes, staging during download, or prefetching RAWs to NVMe ahead of decode. More cores or GPU work will not help.
- **Worker-side prep matters for headroom.** Interleaved A/B at 24 workers (A = arrays and hypotheses built in the GPU process): GPU-stage capacity is 30–31 k runs/h for A vs 47–52 k for B. Makespan is about the same (A 34.0 / 34.4 s, B 31.5 / 28.6 s), because decode binds in both.
- Two GPU-stage threads in one process (round v2) made every stage slower through GIL contention, so it was reverted. The default is `--search-threads 1`.

### MS1 stage (`--ms1`, `ms1_w24`; 50 runs, 24 workers)

| | Without MS1 (c50_w24) | With MS1 |
|---|---|---|
| core-s per run | 4.6 | 10.2 |
| MS1 core-s per run | — | 4.0 mean (D0 ≈ 30, Velos ≈ 4) |
| Share of decode-worker busy time | — | **53%** |
| Worker fixed cost (numba load/compile) | — | about +4 core-s per worker |
| runs/h, 50-run makespan | 12.8 k | **4.5 k** (tail: one D0 MS1 job is ~30 s on 1 core) |
| Decode-stage capacity | 26 k | **11.4 k** |

**The MS1 stage roughly halves the CPU-bound runs/hour budget.** It is still above the cold-Lustre ceiling.
- In this setting, MS1 should use more than one numba thread for the large runs, or split D0-sized runs across workers.
- The GPU hill linker (0.26 s vs 0.55 s) would only pay if MS1 ran in the GPU process. Here it runs in CPU workers that do not import torch, by design (§2).

## 2. Startup (`startup.json`; fresh process: imports, CUDA init, kernels .so, L1 load, GPU index build)

| | Lustre env (`conda_envs/pe`) | **Node-local NVMe env** (`scripts/m2b/stage_env.sh`) |
|---|---|---|
| Process wall | 13.2 s (cold: 42–183 s seen today) | **2.6 s** |
| numpy + pyarrow | 5.8 s | 0.12 s |
| torch import | 4.1 s | 1.1 s |
| CUDA init | 0.3 s | 0.3 s |
| L1 load + GPU index | 1.1 + 0.8 s | 0.3 + 0.3 s |

**Fixes:**
1. `scripts/m2b/stage_env.sh` stages the env and the .NET runtime, RawFileReader and VulcanRaw onto local NVMe. It runs as parallel tar pipes and is idempotent; the one-time stage took about 7 min on loaded Lustre. `scripts/m2b/vpy.sh` runs the staged env and falls back to the Lustre env.
2. **Persistent processes.** The cohort pays this setup once, overlapped with the first decodes (2.2–3.8 s), not once per run.
3. **Decode workers never import torch.** The torch-free host code moved to `vulcan/gpu/host.py`; engine re-exports every name. Importing torch in each worker had cost about 7 core-s per worker (fixed cost 65 → 14 core-s for 8 workers).

## 3. Ingest: direct Thermo reader at cohort scale (`ingest_profile.json`, `lustre_read.json`)

**Is it the bottleneck?** With warm caches, yes: the decode stage binds, at 2.3–3.3 core-s per run. With cold Lustre, no: the Lustre read binds first.

### Profile: one 1-thread process like a decode worker, local RAW copy; Velos, D4 SPS (Eclipse), D0 01A (QE), D7 LM2480 (Lumos)

| | core-s per 1e5 MS2 |
|---|---|
| First conversion in a fresh process | 17.3 (= the old "19 core-s" figure: .NET JIT, runtime load, lazy Lustre imports) |
| **Steady state in a persistent worker** | **6.5** |
| RawFileReader peak fetch (centroid decode) | 3.12 (48%) |
| Scan event + MS2 reactions/trailer | 0.68 |
| double → float peak copy (C#) | 0.62 |
| Python batch copy (pythonnet → numpy) | 0.49 |
| Arrow write | 0.47 |
| Filter string | 0.27 |
| Filter object | 0.13 |
| Scan stats and trailer | 0.09 |

- **Implemented:**
  - persistent workers (2.7× fewer core-s);
  - RAW copied to local NVMe before decode. In round 1, RawFileReader mmap on Lustre plus page-cache reclaim under the 256 GB cgroup limit turned 5 s decodes into 77 s;
  - supervision against hangs;
  - .NET section profiling knobs (`Profile`, `FastTrailer`) with the default behaviour unchanged. `FastTrailer` (one `GetTrailerExtraValues` call) was slower (0.10 vs 0.015 s), so it is off; its output was identical.
- **Hangs found:** RawFileReader hangs sporadically, in a futex wait at 0% CPU. We saw 4 hangs in about 6 k runs: 2 before supervision existed, and 2 in the 4.1 k supervised runs (`v3_torch_in_worker/abB1_w24`, `c50_w8`).
  - It keeps a per-RAW named mutex and shared memory under `$HOME/.dotnet/corefx/Runtime/SharedMemory`. On this cluster that is NFS, and neither HOME nor TMPDIR moves it. This is the suspected cause, not proven.
  - The supervisor recovered both supervised hangs, each by one kill, respawn and retry. There were none in the final round.
  - A hang now costs at most about 30 s of one worker (the stall check).
- A race in `thermo_direct.load_runtime` was fixed: all workers wrote one shared runtimeconfig path, causing hostfxr `InvalidConfigFile`. The path is now per-process and written atomically.
- **Proposed, not done:**
  1. Write Arrow from C# (Apache.Arrow) to drop the Python copy and Arrow write (−15%).
  2. Emit the filter string only for MS1/MSn events that need it (−4%).
  3. A RawFileReader ThreadAccessor per file, to cut per-run latency (the tail of finite cohorts), not core-s.

  The remaining ~50% is RawFileReader's own decode, which we cannot change.
- **Uncached Lustre read** (O_DIRECT, distinct RAWs): 0.30 GB/s with 1 stream, 0.80 with 4, 1.15 with 8, **1.27 GB/s with 16**. This is a shared filesystem and was measured once.

## 4. Compact wire format (`h2d.json`): measured, not worth it

- Pageable H2D of one D0 run's MS2 peaks (71 MB f32) takes **5.1 ms, 14 GB/s, which is 4.5% of the run's 115 ms search**. Pinned memory gives 4.6 ms.
- The u32 bin + f16 format would move 53 MB in 2.9 ms. It saves about 2 ms per run but costs 0.22 core-s to pack in numpy, and the GPU stage is not the binding stage.
- **Not implemented.** Revisit only if the GPU stage binds, for example with entrapment or open search on a pre-staged cohort.

## 5. Family top-2 (`vulcan/gpu/family.py`; `family_validation_D1.json`, `family_validation_D0.json`, `tests/gpu/test_family.py`)

**Index:** `vulcan.index.build.build_multi_decoy` builds T ∪ E_s ∪ dec_j(T) ∪ dec_j(E_s) for entrapment seeds 2001 and 2002 and decoy seeds 1001 and 1002. The decoys are C-term-fixed hash shuffles, drawn against a shared used-set so that no two sets collide. The index has one `partner_j` per decoy set, and `labels.arrow` has `decoy_set` and `decoy_seed`. The original `build()` is unchanged.

**Search:**
- K1 top-K is taken **per family** (T, E_s, dec_j(T), dec_j(E_s)).
- Pair closure runs over all decoy sets.
- **Admissibility per decoy context j:** a row counts for j only if it, or its set-j partner, was shortlisted in its own family. This is what a separate replicate search would score.
- **Canonical hypothesis:** the lowest-index hypothesis whose window contains the row.
- Top-2 is kept per (spectrum, family × context).
- A replicate (j, s) is reconstructed by merging its 4 groups (`family.reconstruct`).

**Validation:** separate searches of each replicate's row-subset index (same peptidoforms) under the same per-family top-K rule.

| | Replicates | Exact equality (spec, rank, pf_id, hyperscore bits) | Combined search | 4 separate searches (incl. 4 GPU index builds) |
|---|---|---|---|---|
| D1 Velos (2.30 M rows) | 4 (2 decoy × 2 entrapment seeds) | **4/4 exact** | 0.16 s | 0.41 s |
| D0 01A (31.7 M rows, 789 M postings) | 4 | **4/4 exact** | 0.35 s | 4.8 s |

- **Against the default engine mode** (one global top-K = 500 in a separate search), the reconstructed top-1 agrees on **99.97%** of D0 spectra and 100% of D1 spectra. The rule therefore barely changes results.
- **Not reconstructable:** `scored_candidates` and `poisson` depend on the whole database, so the family table omits them (`PSM_SCHEMA.md` §2).
- **Cost:** building the D0 multi-replicate index takes 5 min at 4 threads. Its GPU index uses 3.1 GB on the A100.
- **Label use:** retention needs the family/side codes from `labels.arrow`. This is the family-tag exception sanctioned in DESIGN §4 L1, and the rule is symmetric in the target and decoy roles.

## 6. Correctness

- **Tests:** `tests/index` and `tests/gpu` (K1 == brute force, G0 role swap and permutation, family exactness, multi-decoy structure, loader equality) and `tests/io` all pass. The slow mzML-equality test was run separately with `VULCAN_SLOW=1`.
- **Two determinism fixes in `engine.py`, found by running several CUDA streams concurrently:**
  1. The per-spectrum TIC came from an f64 global prefix-sum difference. It was not bit-stable from run to run (scan association), which showed up as last-bit changes in `matched_intensity_pct` on 5 of 35 runs. It is now an f64 row sum over a padded matrix.
  2. When a (spectrum, row) came from several hypotheses, its hypothesis was picked arbitrarily (`index_put` with duplicates). It now takes the lowest index. There were 0 such duplicates on D0 and D1.

  E2's D0 PSM tables are reproduced bit-identically after both fixes (`psm_vs_E2_after_tic_fix.json`).
- **PSM schema:** `vulcan/gpu/PSM_SCHEMA.md` documents it. No column changed. The family table is a new file.

## Caveats

- Throughput uses repeated dev files with page-cache-warm RAWs, on a shared, noisy node. Reps at 24 workers spread from 20.6 k to 25.1 k runs/h.
- Runs/h depends on the file mix (mean 28.7 k MS2 and 0.72 GB per run here).
- "GPU-s" are CUDA-event phase times. They include host work inside phases, so they are an upper bound. NVML utilization is the cross-check.
- `/tmp` state (staged env, stores, indexes) is rebuilt by `scripts/m2b/stage_env.sh`, `scripts/m2/e2_prepare.py` and `scripts/m2b/family_validate.py`.

## Next

1. **I/O first.** Measure a truly cold cohort with a download/stage pipeline, plus multi-stream RAW prefetch to NVMe ahead of decode. It binds at about 6.5 k runs/h on this Lustre.
2. Arrow writing in C# and lighter MS1 handling for large runs (multi-threaded or split).
3. Wire the MS1 hypothesis table into the search (it replaces `pre_hypotheses`; the precursor-hypothesis sharing is ready).
4. E3 proper: micro-cohort mass-sorted batching only pays once the GPU stage binds, which is not yet the case for closed search. Run it with the entrapment/multi-replicate index (4.5× the rows) or open search.
5. Use family top-2 for the M1e/M3 entrapment replicates (5 decoy × 3 entrapment at milestones).
