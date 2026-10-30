# E7 chimeric search + scorer v2 (dev only)

## 0. Pre-registration and amendment
- **Plan (commit 8be8f45, before any chimeric search, v2 training or dev evaluation).** `config.yaml` fixes the chimeric-search parameters, the competition and strata, the v2 feature list, the training-set rule and the **model-choice rule**.
- **Training set frozen (bd6e215, before training):** `train_v2_folds.tsv`, 43 runs from 43 projects (M3 13 + M7 30, raw sha256). Three runs were excluded label-independently (`train_prep_log.json`): two RAWs RawFileReader cannot decode, plus MN-WT1-Frac8, which is not OT-HCD (as in M3). That leaves **40 runs**.
- **Amendment E7-A1 (c7bdaae).** It is dated after the D0/D1 `vulcan_chim_v1` stratum results were seen, but before any v2 model existed and before any chimeric result on D8h. It records the outcome of the pre-registered E7 rule and adds the per-stratum-TDC variant `sep`. D8h is its confirmation set.
- **v2 frozen by the rule (6607ecd) before any dev evaluation of v2.**
- Untouched: the frozen v1 `static_gbm_pred` (sha 31bfe685), §6, the metrics, the labels and the estimators. Not used: sealed data, DN2.

## 1. What was built
- **`vulcan/chimeric/search.py`: E7 chimeric search.**
  - Group 0 is exactly the closed search: instrument precursor, isotope 0..2. It is chimeric rank 1.
  - Every co-isolated MS1 feature group that is *not* the instrument's precursor becomes its own query unit. This uses the vulcan.ms1 hypothesis table: charge 2–6, MS1 chimeric_rank ≤ 5, isolation fraction ≥ 1%, feature mono m/z, isotope 0..1. The chimeric rank is 1 + its order by MS1 chimeric_rank, which is label-free.
  - Each group has its own K1 top-K (K=500), pair closure and K2 ranking. This is implemented as *virtual spectra*: copies of the real peaks carrying the group's precursor, run through the unchanged K0/K1/K2 kernels.
  - **No kernel changed. With chimeric mode off, the output is bit-identical by construction and by test** (`tests/gpu/test_chimeric.py`). In the T1 driver, the reused primary equals a fresh closed search on prod_1001.
  - On D0 there are 1.02–1.04 secondary groups per spectrum, and 59% of spectra have at least one. On D8h it is 0–2.1 per run; one run (E9_0001) has no MS1 hills in its store, so for that run the chimeric search equals the closed search. Hypotheses go up 1.7×.
- **`vulcan/chimeric/compete.py`: competition.** It runs in this order:
  1. per-group winner;
  2. the CONGA neighbour filter: drop a winner that shares >5% of its b/y ions with a higher-scoring kept winner of the same spectrum. Ions are b/y with index ≥ 3, binned at 0.02 Da. b1/b2/y1/y2 are excluded because y1/y2 coincide between all K- or R-terminal peptides; this choice was fixed in config.yaml before any run. On D0 prod the filter drops 11.8% of group winners;
  3. paired double competition;
  4. strata = chimeric rank {1, 2, ≥3}, with ≥3 pooled and merged label-blind (≥200 pairs; no merge was needed);
  5. q-values: Group-walk (pre-registered), pooled TDC, and `sep` = TDC within each stratum (A1).
- **`vulcan/rescore/features_v2.py`.** These are the v1 features, with the run-level calibrations (NCE, RT, hs_runz) computed on primary PSMs only. On top come 19 label-blind additions:
  - per-ion-series predicted-vs-observed SA/Pearson;
  - cleavage coverage and the longest b/y ladder;
  - precursor ppm against a label-blind calibration curve;
  - MS1 feature quality (iso/trace cosine, n isotopes, intensity run-z, isolation fraction, RT offset in the feature);
  - chimeric context;
  - instrument-charge agreement.
- **T1 driver** (dev only):
  - new methods `vulcan_chim_v1`, `vulcan_chim_v1_sep`, `vulcan_v2` and `vulcan_v2_sep` (`bench/t1/ours_chimeric.py`);
  - a `chim` stage and dev dataset `D8h` (the first 10 D8 runs; it reuses D0's identical DB);
  - `--ours` with several numerators.

  T1 behaviour is unchanged: default methods, numerator "vulcan".

## 2. v2 training (pre-registered rule; `models/training_log.json`)
Training used 40 runs, 6.06 M chimeric PSMs, 3 project folds, leave-one-fold-out. The metric is held-out peptidoforms at 1% after the full chimeric pipeline. Hyperscore alone gives 286,288.

| features | linear C=0.01 / 0.1 / 1 | GBM 50 / 150 / 300 |
|---|---|---|
| v1 set (31) | **407,941** / 405,974 / 404,444 | 394,786 / 400,482 / 398,703 |
| v2 set (50) | 407,393 / 404,591 / 403,018 | 394,720 / 398,909 / 397,053 |

**Chosen: linear, C = 0.01, on the v1 feature set** (sha256 3f2c69cc). It was the best, and the rule prefers it anyway: the v1 feature set first, then linear, within 0.5%.
- **The 19 richer label-blind features gave no held-out gain** (−0.13% for linear, −0.4% for GBM).
- The v2 gain therefore comes from the chimeric search, plus retraining on 3× more projects with chimeric PSMs.

## 3. E7 gain (frozen v1 model; closed `vulcan` vs chimeric). Dev, D0/D1 = selection set
Peptidoform ratio at matched arm-P FDP 1%, with the run × decoy-shuffle bootstrap CI (T1-driver estimator):

| | D0 | D1 |
|---|---|---|
| chim_v1 (Group-walk) / closed v1 | **1.208 [1.206, 1.211]** | 1.171 [1.170, 1.172] |
| chim_v1_sep / closed v1 | **1.204 [1.203, 1.208]** | 1.155 [1.150, 1.160] |
| protein groups, chim / closed | 1.140 [1.127, 1.156] | 1.114 [1.111, 1.118] |

Overall entrapment at 1% on D0 is unchanged:
- closed v1: P 1.01%, C 0.83%;
- chim_v1: P 1.05%, C 0.87%;
- chim_v1_sep: P 1.00%, C 0.88%.

**By rank stratum** (D0, at matched FDP 1%; IDs = production peptidoforms by the stratum of the pair winner; FDP̂ restricted to the stratum; verdicts under convention (a)):

| stratum | IDs GW / sep | GW: P / C / verdict / decoy D/T | sep: P / C / verdict / decoy D/T |
|---|---|---|---|
| 1 | 22,152 / 22,347 | 0.73% / 0.60% / supported / 0.56% | 0.99% / 0.85% / supported / 0.99% |
| 2 | 4,595 / 4,527 | 1.50% / 1.20% / inconclusive / 1.54% | 0.94% / 0.84% / supported / 0.97% |
| ≥3 | 2,378 / 2,162 | 3.09% / 2.75% / **violated** (LB CI [1.23, 1.94]%) / 3.94% | 1.30% / 1.21% / inconclusive / 0.93% |

The closed baseline is 24,324.

On D1 (`strata_D1_*.json`):
- stratum 1 is supported for both variants;
- strata 2 and ≥3 are not evaluable at 1% (α·N_T < 20);
- at 5%, GW stratum 2 arm C is "violated" (16.8%), while every stratum under sep is supported or not evaluable.

**E7 verdicts.**
- **The pre-registered rule fails for Group-walk:** stratum ≥3 is violated at 1% on D0. This is recorded as final, so chimeric + Group-walk is not the default.
- Inside each stratum the entrapment FDP tracks the decoy estimate (3.1% vs 3.9%). So TDC is calibrated within strata. Group-walk controls only the *overall* list, and it moves the error into the high-rank strata.
- **The `sep` variant passes the D0/D1 part of A1:** ratio lower bound 1.203 > 1, and no stratum violated on D0 or D1. It is confirmed on D8h (§5).
- v2 (Group-walk) repeats the pattern: stratum ≥3 is violated on D0 (P 2.54%).

## 4. v2 vs comparators (D0/D1 = selection set)
At matched FDP 1%, peptidoform level, `vulcan_v2_sep`:

| vs | D0 | D1 |
|---|---|---|
| **Sage+pred+mokapot (strongest)** | **1.211 [1.208, 1.214]** | **1.171 [1.155, 1.175]** |
| vulcan v1 (closed, frozen) | 1.222 [1.220, 1.226] | 1.173 [1.158, 1.176] |
| Sage (RT on) | 1.257 [1.253, 1.261] | 1.178 [1.161, 1.181] |
| Sage+mokapot | 1.254 | – |
| Comet+mokapot | 1.264 | – |
| protein groups vs strongest | 1.144 [1.129, 1.161] | – |

D0 IDs at 1%:
- `vulcan_v2_sep` 29,484 peptidoforms / 27,746 sequences / 6,979 protein groups;
- Sage+pred+mokapot 24,576 / 23,180 / 6,465.

The overall verdicts are "inconclusive" for every method on D0, Sage included (the paired estimator is upward-biased, see 1c). For `vulcan_v2_sep` at 1%: P 1.07% [1.01, 1.14], C 0.77%. **v2_sep strata on D0 at 1%:** P 1.08 / 1.02 / 1.13% with C 0.74 / 0.86 / 0.83%, all inconclusive or supported; none violated. v2 (GW) gives 1.217 [1.215, 1.221] against the strongest comparator.

v2 vs chim_v1 isolates the retraining effect: +1.5% on D0 (29,484 vs 29,056 sep).

## 5. D8h confirmation (first 10 D8 runs; not used for any selection)
D8h = the A1.1 holdout, 10 runs from 10 projects. It went through the same driver, 25 replicate DBs and estimators, and was looked at only after v2 and A1 were frozen. Ratios are run-level (summed over runs) at matched arm-P FDP 1%, with the run × shuffle bootstrap CI.

| ratio | peptidoform | sequence | protein groups |
|---|---|---|---|
| **v2_sep / Sage+pred+mokapot (strongest)** | **1.222 [1.138, 1.279]** | 1.226 [1.142, 1.281] | **1.192 [1.133, 1.252]** |
| v2_sep / vulcan v1 closed (frozen) | 1.252 [1.149, 1.321] | 1.256 [1.152, 1.324] | 1.221 [1.156, 1.282] |
| chim_v1_sep / vulcan v1 closed (E7 gain, frozen model) | 1.191 [1.113, 1.244] | 1.194 [1.115, 1.246] | 1.198 [1.138, 1.261] |

- At own q 1%, v2_sep / strongest = 1.220 [1.140, 1.272]. That comparator is "supported" at 1% on D8h, so this own-q comparison counts under §1.
- Dataset-level IDs at 1%: v2_sep 79,307 peptidoforms / 6,375 protein groups; Sage+pred+mokapot 70,194 / 6,149; v1 closed 67,392 / 5,984.
- Derived (not bootstrapped): v1 closed / strongest ≈ 0.976 on D8h, matching the D0 dry run (0.990).

**Overall entrapment at 1%:**
- v2_sep: P 1.06% [1.01, 1.10], C 0.87% [0.82, 0.91], "inconclusive";
- v1 closed: P 1.02%, C 1.04%;
- Sage+pred+mokapot: P 0.97%, C 0.99%, "supported".

**By stratum, at 1%** (IDs at the matched-FDP-1% threshold; FDP̂ restricted to the stratum):

| stratum | v2_sep IDs | P | C | verdict (a) | chim_v1_sep: P / C / verdict |
|---|---|---|---|---|---|
| 1 | 55,780 | 1.06% | 0.89% | inconclusive (C supported) | 1.03% / 1.04% / inconclusive |
| 2 | 18,217 | 1.02% | 0.80% | inconclusive (C supported) | 0.90% / 0.96% / supported |
| ≥3 | 6,041 | 1.15% | 0.93% | inconclusive (C supported) | 1.07% / 1.03% / inconclusive |

At 5% every stratum is supported or inconclusive. The in-stratum decoy estimates are 0.97–1.00%.

**E7-A1 verdict: confirmed.** On D8h no stratum and no overall verdict is "violated" at 1%, for v2_sep or chim_v1_sep. Per-stratum TDC (`sep`) is therefore the v2 default competition, and chimeric search + sep is the E7 outcome.

The Group-walk variants were not evaluated on D8h. They were already rejected on D0 by the pre-registered rule, and a first parallel D8h evaluation of all six methods was lost to the job memory limit, so the work was trimmed.

Note: one D8h run (E9_0001) has no MS1 hills in its store, so its chimeric search equals the closed search. It is kept, which makes the ratio conservative.

## 6. Gates on v2 (`gates_v2.json`)
- **G0: PASS.** D1 + E2001 against the E2 role-swapped, row-shuffled index. The chimeric PSM tables (42,833 PSMs, 9,704 secondary) are bit-identical in every column, as are all 50 v2 features (fresh caches) and the v1 and v2 scores.
  - The first run printed FAIL only because `pyarrow.Table.equals` treats NaN ≠ NaN (`ms1_feature_rt_apex` is NaN when there is no MS1 feature).
  - Re-run with a bitwise column comparison, which was recorded.
- **G2-original: PASS.** This is the M3 mixture design, reproduced from its seed and checked against the saved mixtures, run through the chimeric search. Offset copies shift every group.
  - 208,150 null competitions (73,496 secondary units); fixed p(target win) 0.494.
  - v2 vs hyperscore: diff −0.0010, UCB 0.0005 (≤ 0.01); top-decile UCB 0.0053 (≤ 0.02).
  - Diagnostics:
    - the offset-null subset alone (n = 3,763) has diff +1.06%, UCB 2.1%, top 2.8%;
    - the secondary units pass (UCB 0.003).
  - The frozen v1 on the chimeric search also passes (UCB 0.0003).
- **G5: PASS.** AUC 0.504 (≤ 0.55); the fixed-score reference is 0.50.
- **Pending:** tail-G2, G1 and G3 wait for the A2 DN2 nulls. No DN2 data were read. So under D12/A1, FDR claims for v2 stay blocked, exactly as for v1.

## 7. Runtime (`runtime.json`)
Setup: D0, one persistent process, warm prediction cache, under `flock`, median of 3; loadavg 17–25 on 32 cores; slow-path PredCache features for both columns.

| per run | closed v1 | chimeric v2 |
|---|---|---|
| MS1 features + hypotheses | – | 9.8–10.7 s (≈ 25 core-s, 8 threads) |
| search (GPU-s) | 0.14 s (0.14) | 0.35 s (0.25) |
| hypotheses | 125k–140k | 213k–236k (1.7×) |
| PSMs scored | 147k–152k | 287k–307k (2.0×) |
| features + model | 2.2–2.5 s | 2.2 s |
| **total wall / core-s** | **2.4–2.6 s / 2.9–3.4** | **12.5–13.3 s / 28–29** |

**The MS1 feature finder dominates:** about 25 of the 28 core-s, i.e. about 9× the closed per-run core-s. Search GPU-s rise 1.8×.
- v2 uses only v1 features, so the 7.2× fast path (GPU library, `vulcan.rescore.fast`) applies directly. The only change it needs is primary-only calibration, which has not yet been ported.
- At cohort scale, MS1 cost (open item 2e) is the lever: GPU hill-linking, or computing MS1 only for chimeric runs.

## 8. Caveats
- **Dev only.** D0 and D1 are the selection set; D8h confirms.
- **Verdicts at 1% are "inconclusive" for every method on D0.** This is the known ~20% upward bias of the paired estimator. Only matched-FDP ratios compare methods.
- **The D0 CIs are narrow:** 2 runs, so the bootstrap is mostly over decoy shuffles.
- **Arm C on D1 uses human as the foreign proteome.**
- **The comparators did not search chimerically.** Sage has a chimera option, but it is not in the pre-registered comparator set. A "Sage chimera + pred + mokapot" arm would be the fair next comparator.
- **Two training runs are not tryptic** (LysC, chymotrypsin; few PSMs). They were kept because the pre-registered rule did not exclude them.
- **Engineering notes.**
  - The D8h evaluation of all methods in one process ran into the job's 256 GiB limit while other jobs were running. It was re-run as separate processes under `VULCAN_T1_NO_LABEL_CACHE=1` (labels loaded per replicate; same numbers, less memory).
  - The neighbour-filter ion cache is keyed by object identity, not `id()`, so a recycled id can never return stale ion sets. Earlier runs held every label object alive, so they were not affected.
- **v2 is not cleared for sealed use.** It is a new scorer version. A1.1 requires D8 confirmation (done here), and tail-G2/G1/G3 are pending.

## Files
| File | Contents |
|---|---|
| `metrics.json` | collected numbers (`scripts/e7/make_metrics.py`) |
| `config.yaml` | pre-registration + A1 |
| `t1/{D0,D1,D8h}/` | driver outputs: metrics.json, tables.md, logs |
| `t1/{D0,D1}_e7`, `t1/{D0,D1}_e7sep` | the E7 runs vs closed v1 |
| `strata_<ds>_<method>.json` | per-stratum read-outs |
| `gates_v2.json` | gates |
| `runtime.json` | runtime |
| `models/` | v2 model + training log |
| `train_v2_folds.tsv` | fold map |
| `train_prep_log.json` | training-run preparation log |

Code: `vulcan/chimeric/`, `vulcan/rescore/features_v2.py`, `bench/t1/ours_chimeric.py`, `scripts/e7/`, `tests/gpu/test_chimeric.py`.

Tests: tests/index, tests/gpu (incl. role-swap G0 and test_chimeric), tests/bench and tests/rescore: 32/32 pass.
