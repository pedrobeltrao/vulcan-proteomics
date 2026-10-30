# E7b: fair chimeric comparator + v2_sep speed (dev only: D0, D1, D8h)

## Bottom line
- **Fair comparator.** Against Sage with chimera mode plus the same AlphaPeptDeep features and seeded mokapot, v2_sep still gives **+14–15% peptidoforms at matched arm-P FDP 1%**. With CIs: D0 1.149 [1.147, 1.152] (R = 5), D1 1.137 [1.120, 1.140], **D8h holdout 1.150 [1.069, 1.203]**. In E7, without chimera, the advantage was 1.17–1.22.
  - **There is no protein-group gain on D8h: 0.999 [0.951, 1.039].** D0 is 1.016.
  - At 0.1% the advantage shrinks: D8h 1.096 [1.024, 1.150] at matched FDP. On D1 at own q it is 0.964 [0.959, 0.971].
- **Comparator calibration.** Sage-chimera is calibrated overall: arm-P FDP at 1% is 0.91–1.09%, nothing violated. Its rank ≥ 2 IDs are not: 4–13% FDP at nominal 1%, and on D8h the arm is "violated" (5.3%, 2,631 IDs). So the ratios above favour the comparator.
- **Speed (all exact).**
  - Per-run v2_sep cost: 2.07× fewer core-s (summed over 13 runs, 241.5 → 116.8). Per D0 run, 16.9 s / 32.5 core-s → 9.3 s / 15.5 core-s.
  - By stage: MS1 1.9×, features 3.6×, competition 5.6×.
  - Chimeric cohort RAW → protein groups: **2,987 runs/h (steady 4,643)** at 19.0 core-s/run, against 2,468 / 3,217 with the old MS1 kernels. Closed v1 in the same session runs at 5,433 / 13,501.
  - ID sets at q ≤ 0.1/1/5% are identical on D0/D1/D8h. Features, scores and MS1 tables are bitwise equal, except 260 of 26M D8h PSM scores that trace to E7's prediction-cache history; the old code path reproduces the new values bitwise.

Pre-registration: `config.yaml`, commit 7a2d6ee. It was committed before any Sage chimera search, any comparator evaluation and any speed change. Nothing sealed was touched (no T1–T5 data, no `bench/sealed_run.py`). No DN2 data and nothing in `/tmp/vulcan_a2` was read. The frozen v1 and v2_sep models, their training, the metrics, controls, labels, estimators and tests are unchanged. Seeds follow DESIGN §6.1.

## 1. Fair comparator: Sage chimera + AlphaPeptDeep features + seeded mokapot

### Configuration (pre-registered, `config.yaml`)

The arm is the existing strongest comparator with `chimera: true`. Everything else is identical: replicate FASTAs, tolerances, predict_rt, Sage 0.14.6, the same 9 AlphaPeptDeep features from the same PredCache, and the seeded mokapot wrapper.

Two adaptations were needed so that Sage's chimeric PSMs reach the output at all. Both were fixed before any run.

1. **`report_psms` = R ≥ 2.** Sage's `score_chimera_fast` iterates "best PSM, then remove its matched peaks" until `report_psms` PSMs exist, so with R = 1 chimera mode is a no-op. R = 2 is primary; R = 5 is a sensitivity arm on D0/D1. The precursor window is unchanged (±10 ppm, isotope 0..2, no wide_window).
2. **mokapot spectrum key = ScanNr·16 + rank.** mokapot keeps one PSM per (FileName, ScanNr, ExpMass), which would delete every rank ≥ 2 PSM. The prediction-feature calibrations use Sage rank-1 PSMs, the analogue of v2's primary-only calibration.

Checks:
- On D1 prod, Sage's rank-1 PSMs are unchanged by chimera mode (11,489 = 11,489).
- Rank ≥ 2 PSMs make up 11–14% of mokapot's target PSMs.

**R rule** (pre-registered): R = 5 also runs on D8h only if it beats R = 2 at matched FDP 1% on both D0 and D1.
- D0: R = 5 wins (26,470 vs 26,074).
- D1: R = 5 does not win (6,238.3 vs 6,238.9).

So D8h ran R = 2 only. The headline comparator per dataset is the R with more IDs: R = 5 on D0, R = 2 on D1 and D8h. In the driver, the methods are `sage_chim2`, `sage_chim5` (Sage LDA, all ranks) and `sage_chim{2,5}+pred+mokapot`.

### Ratios: v2_sep / Sage-chimera+pred+mokapot

Driver estimator, unchanged: run-level IDs summed over runs, mean over 5 decoy shuffles, run × shuffle bootstrap 95% CI. The matched 0.1% level was enabled through the opt-in `VULCAN_T1_FMATCH_0001`.

| dataset (R) | level | matched FDP 0.1% | **matched FDP 1%** | matched FDP 5% | own q 1% (comparator verdict) | own q 0.1% |
|---|---|---|---|---|---|---|
| D0 (R=5) | peptidoform | 1.062 [0.996, 1.085] | **1.149 [1.147, 1.152]** | 1.171 [1.165, 1.175] | 1.150 [1.147, 1.154] (inconclusive) | 1.065 [1.034, 1.096] |
| D0 (R=5) | protein groups | 1.002 [0.989, 1.024] | **1.016 [1.009, 1.021]** | 1.015 | 1.015 (inconclusive) | 0.999 |
| D0 (R=2) | peptidoform / protein | 1.084 / 1.042 | 1.166 [1.164, 1.171] / 1.060 [1.058, 1.064] | 1.195 / 1.059 | 1.168 / 1.060 | 1.078 / 1.039 |
| D1 (R=2) | peptidoform | 1.013 [0.977, 1.098] | **1.137 [1.120, 1.140]** | 1.145 [1.142, 1.148] | 1.124 [1.109, 1.138] (supported, counts) | **0.964 [0.959, 0.971]** |
| D1 (R=2) | protein groups | 1.379 [1.103, 2.448] (few groups) | **1.100 [1.094, 1.106]** | 1.093 | 1.092 (not-evaluable) | – |
| D1 (R=5) | peptidoform | 1.003 | 1.138 [1.121, 1.141] | 1.144 | 1.123 (supported) | 0.961 |
| **D8h (R=2, holdout)** | peptidoform | 1.096 [1.024, 1.150] | **1.150 [1.069, 1.203]** | 1.173 [1.085, 1.231] | 1.149 [1.074, 1.196] (inconclusive) | 1.096 [1.043, 1.133] |
| **D8h (R=2)** | sequence | – | 1.150 [1.070, 1.203] | – | – | – |
| **D8h (R=2)** | **protein groups** | 0.975 [0.932, 1.018] | **0.999 [0.951, 1.039]** | 0.999 | 0.997 (supported, counts) | 0.972 |

**Reproduction.** The non-chimeric reference arm reproduces E7 exactly: D0 1.211 [1.208, 1.214] and D1 1.171 [1.155, 1.175].

**What the fair comparator changes:**

- **Peptidoforms.**
  - The v2_sep advantage at matched FDP 1% falls from +17 to +22% (vs Sage + predicted features + mokapot without chimera) to **+14 to +15%** (D0 1.149, D1 1.137, D8h 1.150).
  - Every lower bound is ≥ 1.069.
  - Chimera mode is worth only +3 to +6% to the comparator: 1.211/1.149 on D0, 1.171/1.137 on D1, and about 1.222/1.150 on D8h (the last uses E7's evaluation of the non-chimeric arm).
- **Protein groups.**
  - The v2_sep advantage largely disappears: D0 1.016 (R = 5) or 1.060 (R = 2), D1 1.100 (few groups), and **D8h 0.999 [0.951, 1.039]**.
  - Compared with E7's 1.19 against the non-chimeric comparator on D8h, chimeric peptides mostly add peptidoforms to proteins that are already identified.
- **0.1%.**
  - On D1 at own q 0.1%, v2_sep is **0.964 [0.959, 0.971]**× the chimeric comparator (lower bound below 0.97; the comparator is "not-evaluable" there).
  - At matched FDP 0.1%, the ratio is 1.013 [0.977, 1.098] on D1, 1.062 [0.996, 1.085] on D0 and 1.096 [1.024, 1.150] on D8h.
  - So the advantage is concentrated at 1–5%.

### Comparator entrapment verdicts and calibration (arm P / arm C mean FDP at own q; convention (a) decides)

| dataset | method | 0.1% P / C | 1% P / C (verdict) | 5% P / C (verdict) |
|---|---|---|---|---|
| D0 | sage_chim2 (Sage LDA) | 0.11 / 0.10% | 1.09 / 0.87% (inconclusive) | 5.29 / 4.68% (inconclusive) |
| D0 | sage_chim2+pred+mokapot | 0.11 / 0.10% | 1.02 / 0.86% (inconclusive) | 5.03 / 4.18% (inconclusive) |
| D0 | sage_chim5 | 0.11 / 0.08% | 1.04 / 0.86% (inconclusive) | 5.25 / 4.59% (inconclusive) |
| D0 | sage_chim5+pred+mokapot | 0.10 / 0.11% | 1.07 / 0.84% (inconclusive) | 4.96 / 4.06% (supported) |
| D1 | sage_chim2 | 0.05 / 0.07% | 1.02 / 0.88% (inconclusive) | 4.89 / 5.38% (inconclusive) |
| D1 | sage_chim2+pred+mokapot | 0.12 / 0.03% | 0.92 / 0.85% (supported) | 4.49 / 4.41% (supported) |
| D1 | sage_chim5 | 0.06 / 0.06% | 1.07 / 0.84% (inconclusive) | 5.06 / 5.60% (inconclusive; sequence level, convention b: violated) |
| D1 | sage_chim5+pred+mokapot | 0.12 / 0.04% | 0.91 / 0.76% (supported) | 4.51 / 4.22% (supported) |
| D8h | sage_chim2 | 0.09 / 0.14% | 0.99 / 0.85% (supported) | 4.89 / 4.63% (supported) |
| D8h | sage_chim2+pred+mokapot | 0.09 / 0.19% | 0.97 / 1.02% (inconclusive) | 4.80 / 4.52% (supported) |
| ref. D0 / D1 / D8h | vulcan_v2_sep | 0.10 / 0.07 / 0.10% | 1.07 / 0.90 / 1.06% | 5.09 / 5.12 / 5.07% |

**Overall, Sage's chimeric output stays calibrated.** At 1% the arm-P FDP is 0.91–1.09%, no verdict is violated under convention (a), and "inconclusive" on D0 is the known upward bias of the paired estimator (as for every method).

**By Sage rank it does not.** Arm-P FDP at own q, by the Sage rank of each peptidoform's best PSM (`sage_rank_strata_<ds>.json`; descriptive, pre-registered as such):

| dataset | method | rank 1: IDs / P FDP at 1% (verdict) | **rank ≥ 2: IDs / P FDP at 1% (verdict)** | rank ≥ 2 at 5% |
|---|---|---|---|---|
| D0 | sage_chim2 | 23,535 / 1.03% (inconclusive) | 409 / 4.3% (not evaluable) | 810 / 34.8% (violated) |
| D0 | sage_chim2+pred+mokapot | 24,406 / 0.77% (supported) | 887 / 8.3% (not evaluable) | 1,444 / 40.1% (violated) |
| D0 | sage_chim5 | 23,369 / 0.97% (supported) | 453 / 5.3% (not evaluable) | 963 / 37.8% (violated) |
| D0 | sage_chim5+pred+mokapot | 24,336 / 0.70% (supported) | 1,122 / 9.4% (not evaluable) | 1,925 / 38.4% (violated) |
| D1 | sage_chim2 / +pred+mokapot | 6,014 / 0.98% / 6,011 / 0.56% | 92 / 3.7% / 214 / 11.5% (not evaluable) | 31.6% / 45.2% (not evaluable) |
| D1 | sage_chim5 / +pred+mokapot | 5,990 / 0.86% / 6,011 / 0.47% | 130 / 10.2% / 220 / 13.0% (not evaluable) | 46.4% / 44.7% (not evaluable) |
| **D8h** | sage_chim2 | 61,929 / 0.93% (supported) | 828 / 5.6% (not evaluable) | 1,894 / 24.0% (violated) |
| **D8h** | **sage_chim2+pred+mokapot** | 69,492 / 0.80% (supported) | **2,631 / 5.3% (violated)** | 4,088 / 26.0% (violated) |

How to read this:
- **The pattern.** The rank-1 list is conservative. The 3–4% of IDs that come from rank ≥ 2 PSMs carry 4–13% FDP at nominal 1% and 24–46% at nominal 5%. This is the same failure mode E7 found for Group-walk over strata: the pooled list is controlled, but the error moves into the secondary stratum.
- **v2_sep has no such failure.** Its per-stratum TDC had no violated stratum on D0, D1 or D8h (E7 §3, §5).
- **What this means for the ratios.** Part of the comparator's chimeric gain is false discoveries in its rank ≥ 2 output. The matched-FDP ratios above (overall arm P) therefore *favour* the comparator. A per-rank-controlled Sage chimera arm (e.g. per-rank mokapot) would gain less.
- **Status.** This is descriptive only, not a verdict on the comparator, as pre-registered.

**Answer to "does Sage-chimera stay calibrated?"** Overall yes; by rank no. Sage's own chimeric secondary PSMs are not FDR-controlled as a list on D0 or D8h.

## 2. Speed: v2_sep production path

### What changed (every change is exact; the old path is kept as the reference)

| stage | change | equality evidence |
|---|---|---|
| MS1 isotope grouping (`vulcan/ms1/isotopes.py`) | RT-binned isotope lookups: hills are bucketed by 32-MS1-scan blocks, m/z-sorted within each block, and visited in (block, m/z) order with a forward pointer. On D0 the original m/z-window sweep visited 98 window hills per lookup, of which 0.16% overlapped in scans. Full-norm trace cosines use per-hill sums of squares precomputed with the same loop. The existence pre-screen stops at the first passing isotope. The tie rule (best cosine, first in m/z order) is kept explicitly. Order-only argsorts are replaced by counting passes. `exact_reference=True` runs the original kernels. | `tests/ms1/test_fast_exact.py`; `equality_ms1.json` |
| MS1 hills (`vulcan/ms1/hills.py`) | Fused numba passes: hill ids, stable counting sort plus gathers, and length filter plus compaction. `reference=True` keeps the original path. | as above |
| v2 features (`vulcan/rescore/fast.py`) | `FastFeatureBuilder.run(primary=group 0)`: the 3i GPU fast path (GpuPeaks, NCE-independent blocks, predicted library mirrored on the GPU). The NCE choice and RT calibration use primary PSMs only, exactly as `FeatureBuilderV2.run_v2` computes the v1 columns. The 19 unused v2 extras are no longer computed. | `equality_scoring.json`, `check_history_D8h.json` |
| chimeric competition (`vulcan/chimeric/fast_compete.py`) | Ion-bin sets in numba with numpy's exact arithmetic (sequential cumsum, numpy pairwise sum, round-half-even), built only for spectra with ≥2 group winners. A numba CONGA filter in the same visiting order. Group-walk q is skipped in production (q_sep does not depend on it). | `tests/rescore/test_fast_compete.py`; competition columns in `equality_scoring.json` |
| cohort (`vulcan/cohort.py --chimeric`, `vulcan/rescore/cohort_scoring.py`) | MS1 features and hypotheses come from the decode workers (existing `--ms1` stage, now the fast exact code). Precursor groups are built in the loader threads. The secondary search, v2 features and model run in the GPU-stage thread. The chimeric competition (q_sep) and picked proteins run in the post workers. | `equality_cohort_reference.json`, `equality_cohort_copies_*.json` |

### Equivalence

- **MS1** (`equality_ms1.json`): on all 13 dev runs (D0 2, D1 1, D8h 10), every column of the feature and hypothesis tables is bitwise identical to the tables E7 wrote with the old code. Schema metadata is equal, and fast and reference Hills/isotope arrays are bitwise equal. One run (E9_0001) has no MS1 hills and gives empty, identical tables.
- **Scoring + competition** (`equality_scoring.json`), on the stored E7 chimeric PSM tables for all 25 replicates of each dataset:
  - D0 and D1: every v1 feature column and the v2 score are bitwise equal to E7, the NCE is equal, and there are 0 library misses. Winners, strata, q_sep and q_pooled are identical. The v2_sep ID sets at q_sep ≤ 0.1%, 1% and 5% are identical on all 5 prod replicates.
  - D8h: the ID sets at q ≤ 0.1%, 1% and 5% are identical on all 5 prod replicates. 62 of 250 (replicate, run) pairs are not bitwise equal to E7, in 260 of 26.0M PSM scores (|Δscore| ≤ 8e-6), almost all in the RT features. The cause is cache history, not the port:
    - `check_history_D8h.json` reruns the OLD slow path (FeatureBuilderV2 + PredCache) on the same cache for all 62 pairs. It is bitwise equal to the fast path everywhere.
    - E7 scored D8h in parallel processes. These appended predictions of the same peptides with batch-dependent last bits (3i: the predictor is not batch-composition invariant), so the first stored copy no longer matches what some processes used in memory.
    - Because of the 260 changed scores, q_sep differs slightly on some prod replicates (`all_prod_competition_equal` false), while every thresholded ID set is identical. Jaccard is 1.0, and the ratio change is 0.
- **Cohort** (`equality_cohort_reference.json`): a 12-run chimeric cohort (the 10 D8 + 2 D0 RAWs, with stores kept) was checked against the per-run reference pipeline built only from pre-E7b components, on the same stores and library root. That pipeline is: original MS1 kernels, closed primary search without worker-prepared hypotheses, `FeatureBuilderV2` + PredCache, the frozen v2 model, `chimeric_competition`, and picked proteins. For 12/12 runs the following are identical: the MS1 hypothesis table, the PSM table, `score` (bitwise), the peptide table (pf_id, spec_idx, score, is_decoy, q), the ID sets at 0.1%, 1% and 5%, and the protein groups at 1%.
- **Cohort copies** (`equality_cohort_copies_{chim,chim_ms1ref}.json`): in both timed 96-run cohorts, every run's PSM, peptide and protein tables equal its base run in the 12-run reference cohort, for all columns except the run name (`<base>_rNN` by construction). For the `chim` check the name column was excluded post hoc, after the first comparison flagged only that column.

### Per-run cost (`bench_before.json`, `bench_after.json`; `scripts/e7b/bench_run.py`)

Setup: one persistent process per implementation, warm, 8 numba/torch threads, the dataset's T1 prod_1001 index, and the dataset's consolidated prediction-cache copy (`data/cache/e7b_predlib/<ds>`). Each run gets 1 warm-up plus 3 measured repeats, taken while holding `flock /tmp/vulcan_bench.lock` per run. The two implementations ran back to back. Load average was 46 at the start of "before" (it fell to 19) and 21 at the start of "after".

Sums over the 13 runs (medians per run):

| stage | wall s before → after | core-s before → after | × core-s |
|---|---|---|---|
| MS1 features + hypotheses | 70.3 → 49.0 | 178.0 → 92.4 | 1.9 |
| chimeric search | 4.6 → 4.5 (GPU-s 3.28 → 3.13) | 4.6 → 4.5 | 1.0 |
| v2 features | 32.1 → 8.9 (GPU stream incl. model 33.0 → 9.7 s) | 31.7 → 8.9 | 3.6 |
| v2 model | 1.1 → 1.1 | 1.1 → 1.1 | 1 |
| competition (q_sep) | 16.8 → 3.0 | 16.7 → 3.0 | 5.6 |
| picked proteins (estimator, unchanged) | 8.0 → 7.0 | 8.0 → 7.0 | 1.1 |
| **total** | **133.8 → 74.2** | **241.5 → 116.8** | **2.07** |

- **Per D0 run (mean of the 2 runs):** 16.9 s → 9.3 s wall, and 32.5 → 15.5 core-s.
  - MS1: 26.8 → 13.1 core-s.
  - Features: 3.0 → 0.8 s.
  - Competition: 1.4 → 0.3 s.
- **For comparison:** closed v1 was 2.5 s / 3 core-s per D0 run in E7. The E7 runtime table (13 s / 28 core-s) left out competition and proteins (≈2.3 s per D0 run before).
- **MS1 is now 79% of the remaining per-run core-s.**

### End-to-end cohort, RAW → protein groups, chimeric mode (`cohort_{chim,chim_ms1ref,closed}.json`)

Setup: the 3i design (10 D8 + 2 D0 RAWs on local NVMe, cycled), human_contam_prod index, 3i predicted library (24.9 GB GPU-resident), 24 decode workers, 3 post workers. Each config was one cohort process under the lock, run back to back.

The runs were cut from 192 to 96: the first 192-run attempt filled `/tmp` (see Caveats).

| config | runs/h (incl. ~50 s setup) | steady runs/h | core-s/run | of which MS1 / decode in workers | GPU-s/run (search) | bottleneck |
|---|---|---|---|---|---|---|
| **chimeric v2_sep (after)** | **2,987** | **4,643** | **19.0** | 8.4 / 6.4 | 0.32 | decode pool |
| chimeric, MS1 reference kernels | 2,468 | 3,217 | 28.2 | 16.9 / 6.6 | 0.33 | decode pool |
| closed v1 (3i config, same session) | 5,433 | 13,501 | 9.1 | 0 / 7.1 | 0.11 | decode pool |

- In the cohort, MS1 runs single-threaded in the decode workers: 16.9 → 8.4 core-s per run, which raises steady throughput 1.44×.
- Scoring was already fast in both chimeric configs. Before E7b, scoring in the cohort was not implemented for chimeric mode; the per-run bench above covers the scoring "before".
- **Chimeric mode costs ≈2.1× the closed core-s per run.** Peptide IDs at 1% summed over the 96 runs: 1.30M vs 1.08M closed (different scorers, v2 vs v1).
- Load average at the start: 12 / 29 / 29 (chim / chim_ms1ref / closed). During the `chim_ms1ref` and `closed` cohorts, a light equality checker and the D8h history check ran on the node. So the absolute runs/h are load-dependent, as in 3i.

## Caveats

- **Disk incident.** The first 192-run chimeric cohorts filled node-local `/tmp` to 100% for a few minutes.
  - Decode ran ahead while the GPU stage waited ~50 s for the library warm-up: 63 GB of stores. The setup overlap had not been budgeted.
  - It killed both comparator jobs; they resumed idempotently, and temporary shards are dot-prefixed and were removed. It may also have affected other agents' jobs on the node in that window.
  - The cohorts were re-run at 96 runs with a free-space guard, and LESSONS.md has a new rule.
- **Before/after load differed** (46 → 21 at start); the core-s comparison is less load-sensitive than wall time.
- **The "MS1 reference" cohort is not a full "before" cohort:** its scoring and competition are already the fast ones.
- **Dev only.** These are speed changes to a candidate scorer that is not cleared for sealed use (D21).
