# A5.11.3: bench/t5_armM.py (arm M on T5; shared dev/sealed module)

- **Date:** 2026-10-28. **Agent:** A5 implementer, fork "T5". T1-blind and T5-blind: no file under data/sealed and no T5 sample sheet was read. Nothing sealed was run.
- **Code:** `bench/t5_armM.py`. Tests: `tests/bench/test_t5_armM.py` (17 pass).

## What it does (A5.6, unchanged design)
- **Databases.**
  - Targets: yeast_UP000002311 + contaminants (fixed in sealed mode; no spiked protein; no sample sheet).
  - Entrapment: human Swiss-Prot peptides, neighbour-filtered by `scripts/entrapment/build_dbs.py --arms C` (unchanged), restricted to the effective search space, and subsampled to r = 1 with `default_rng(seed)`. This is the `armM_prepare.py` rule.
  - Sealed: recipe `shuffle` (forced). Per entrapment seed 2004–2006, one `build_multi_decoy` index with the 5 DECOY_SEEDS C-term-fixed shuffles, materialised by `bench/t1/dbs.materialize` into replicates M_<j>_<s>. That gives **15 replicates**.
  - Dev: recipe `dev_foreign` (`build_foreign`, decoy seed 1001), which is the dev design.
- **v1.** Closed search (report_psms 5) through `bench/t1/ours.search_replicates`, with stores ingested by `vulcan.cohort` (the T1 path; checked bit-for-bit against the engine on the first replicate). Scoring goes through **`bench/t1/ours.py Scorer`** (the T1 path). Paired double competition (seed 11).
- **v2_sep.** `vulcan.ms1`, then the chimeric search, then `ScorerV2` (frozen v2), then `chimeric_competition` with q_sep.
  - Merged strata are keyed by block label when every replicate has the same merge map, otherwise by "contains_k". They decide voids.
  - Raw strata {1, 2, >=3} are reported.
- **Hyperscore** is run through each pipeline as a reference.
- **Estimator and verdicts.**
  - Combined FDP_M with r = peptide-level effective r (1.0).
  - Units: peptidoform (deciding) and sequence.
  - α ∈ {0.1%, 1%, 5%}.
  - CI convention (a) decides; (b) is reported.
- **Prediction caches** are fresh and live under --work (`predcache_t5_v1`, `predcache_t5_v2`). They are never shared with T1.
- **Sealed-mode refusals:**
  - it must run under sealed_run (SEALED_DATASET = T5), with raws equal to $SEALED_FILES;
  - seeds are fixed;
  - `--stores`, `--targets` and the dev recipe are refused;
  - `--out` must be under $SEALED_LOOK_ROOT;
  - `--work` must be under $SEALED_LOOK_ROOT or $SEALED_LOOK_SCRATCH/T5;
  - both must be absent or empty (no resume).
- **Sealed console.** Only data-free stage lines are printed (`CONSOLE_RE`). Everything else goes to the work log, which is deleted.
- **Completeness.** Every replicate × run must be present in every list. Otherwise `complete=false` and the module exits 3.
- **Sealed clean-up.** `work/` is deleted except `db/` and `db_*/` (FASTA-derived indexes). `<out>/manifest.sha256` is written.

## Dev equivalence (A5.6, A5.11.3.3): PASS, exact
`t5_dev_equivalence.json`. Setup: D7 (LM2486/2487/2490), seeds 2001–2003, decoy seed 1001, the dev target FASTA (yeast + contaminants + P37840, sha256 15716a2c…), and the dev stores.
- **Rebuilt databases:**
  - the foreign lists of every seed and the arm-C peptide TSV are byte-identical to the dev ones;
  - the seed-2001 index (`peptidoforms.arrow`, `labels.arrow`) is byte-identical.
- **v2_sep** (overall and raw strata 1, 2, >=3, α = 0.1/1/5%) and **hyperscore_sep** are **identical** to `results/2026-10-28_A3_preconditions/armM.json`: every verdict, mean, CI, LB, LB CI, mean N_T and per-seed count.
- **v1 through the T1 path** (`ours.py Scorer`) and the hyperscore reference are **identical** to 3e's `results/2026-10-27_M3e_hardening/armM.json`. So `ours.py` and `scripts/m3/features.py` give the same result on D7, and no difference needs explaining.
- **Engineering change (no analysis change).** Each (list, unit, α) now makes one `evaluate_arm` call (convention a), and reads convention (b)'s CI from the same bootstrap (`ci95_b`, `lb_ci95_b`, `ENT.verdict`). The calls run in a fork pool.
  - Both versions give identical output on D7 (`single_vs_two_bootstrap_calls_identical`), and a unit test checks this.
  - Evaluation time dropped from about 10.5 min to about 1.5 min on D7.
- **Wall time on D7:** 18.4 min for the first run (3 replicates; the db build took 4 min). See the 15-replicate timing below.

## Deviations (recorded for config/a5.yaml)
- **A5.11.3.1.** `scripts/m3e/armM.py`, `scripts/m3e/armM_prepare.py` and `scripts/a3/dev_armM.py` are **not** refactored to call the module. This keeps their committed outputs' provenance. Equivalence is shown by the exact dev reproduction above instead.
- **Arm-C neighbour filter.** It rejects collisions with the shuffled decoys of DECOY_SEEDS 1001–1005 and the pseudo-reverse set (`build_dbs.py`, unchanged). The same list serves every decoy seed.

## Sealed-like dev run (timing for look.time_limit)
`t5_dev_D7_15rep_armM.json`. Setup: `--mode dev` on D7 with seeds 2001–2003 × decoys 1001–1005, giving **15 replicates**. It used the `shuffle` recipe and default targets (yeast + contaminants), with a fresh `vulcan.cohort` ingest from RAW. That is the sealed code path in every respect except the seeds and the data.
- **Result:**
  - complete, nothing missing;
  - the cohort output is bit-identical to the engine search on the first replicate;
  - merge maps are identical across replicates.
- **Merged-block keys** are the block labels "1", "2" and "3". "3" is the pooled rank ≥ 3 block, and a merged block such as "2_3" would be keyed that way.
- **Wall time: 26.2 min for 3 runs**:
  - arm-C list 5 min;
  - 3 × 15-replicate indexes 7 min;
  - ingest 0.3 min;
  - v1 search+score 4.5 min;
  - MS1 + chimeric search + v2 score 9 min;
  - evaluate < 1 min.

  The node was loaded (load average ≈ 40 on 32 cores; the GPU was shared), and work was in /dev/shm.
- **/tmp was full during this run** (100%, other projects' /tmp/vulcan_* dirs). The first attempt died with ENOSPC. The 15-replicate work dir needs about 20 GB. The look's node-local scratch (/tmp/vulcan_look_<id>) needs free space checked before the look.
