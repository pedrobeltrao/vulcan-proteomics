# Independent review of M3a (static scorer, predicted features, RESET, gates), E2 and M1e, before any sealed look

Reviewer: independent agent (statistics and computational proteomics). Date: 2026-10-27.

Scope:
- PROMPT.md; DESIGN.md §1, §3, §4 L2–L4 and §6 (prereg-v0.4); RESEARCH_TREE.md.
- results/2026-10-27_{E2_gpu_search, M3_scoring, M1e_baselines}; vulcan/{index, gpu, predict, rescore, stats}; scripts/m3.

Only dev data were used; nothing under data/sealed was touched. The check scripts and their outputs are in `reviews/checks/`:

| Script | What it checks |
|---|---|
| `rt_null_asymmetry.py` | RT range and features of the G2 null sources vs the host run; in-situ null proxy |
| `g2_null_tail.py` | null winners above the production q threshold, by label |
| `g2_null_tail_breakdown.py` | the same, split by provenance and by near-homology to the foreign proteome |
| `estimator_checks.py` | rank-1-only ablation, arm-C r units and per-decile r, bootstrap-only CIs |
| `train_independence.py` | accession, lab, sha256 and spectrum-hash overlap; fold map; per-feature null AUCs |

One untimed Sage rerun with `predict_rt: true` (D0/D1, E2 oracle configuration) was made under the bench lock. Its outputs are outside the repo.

## Verdict

**A sealed T1 look is not justified yet.**
- The pre-registered T1 comparison (§1 non-inferiority, §6.3 milestone replication) cannot be computed with the code that exists today.
- The gate set is incomplete, and the gates that were run test less than they claim:
  - G1, G3 and G4 were not run (M2).
  - G2 is nearly vacuous for the default GBM, and its nulls neither exercise the RT channel (M1) nor are clean nulls (M0).
  - Arm M is required by §6.3 for predicted-feature scorers and was skipped (M3).
- My own tail analysis found **no** evidence of label leakage. The target excess among learned-score null winners is explained by near-homologous yeast spectra (M0).
- None of the issues is fraud or leakage through labels: L0–L3 are label-blind by construction, and G0 is sound. Training-set independence holds. The problems are the protocol not yet being complete, gates that test less than they claim, and an unfair Sage-LDA reference.

---

## BLOCKERS

### B1. The pre-registered T1 protocol is not implemented
DESIGN §1 says: "Both methods are averaged over the same 5 shuffled-decoy FASTAs; lower 95% bootstrap CI of the ratio (runs × decoy shuffles) ≥ 0.97." §6.3 says: "At milestones: the full 5 decoy × 3 entrapment".

What exists:
- M3 used one deterministic pseudo-reverse decoy set. `DECOY_SEEDS` (`vulcan/stats/_util.py:8`) is only used for collision repair.
- The M1e sensitivity runs use Sage's internal decoys: `results/2026-10-27_M1e_baselines/configs/sage_D0_prod_t32.json` has `generate_decoys: true`, and Comet uses `decoy_search=1`.
- No code computes the ID-ratio bootstrap CI (grep for the 0.97 threshold and the ratio finds nothing).
- No T1 driver exists in scripts/ or bench/.
- M1e is unfinished. The D0 `pipelines` entry is `{}` in `metrics.json`, D0/D4 rows are empty in `tables/main_tables.md`, and D4 is reduced to 2 files in the uncommitted `scripts/baselines/run_searches.py`.

Why it blocks: every scorer, Sage included, is "inconclusive" at 1%. By §1 only the **matched-entrapment-FDP** comparison then counts, and that needs comparator entrapment runs on the same DBs, which do not exist on D0.

**Fix:** before tagging a milestone, write and dry-run on D0/D1 a single T1 driver that:
- builds the 5 shuffled-decoy DBs plus 3 sealed-seed entrapment DBs (family-max makes this one search);
- runs ours plus Sage+LDA (predict_rt on), Sage+mokapot and Comet+mokapot on the identical external-decoy DBs;
- computes the ratio CI and the matched-FDP curves.

Only then tag and open T1.

### B2. (Procedural) §6.4 preconditions for any FDR claim are unmet
§6.4 says the gates "all must pass before any FDR claim", and §6.3 lists the arms required for this model class. A T1 entrapment verdict is an FDR claim, so majors **M0–M3** below are preconditions for T1, not follow-ups:
- the G1/G3/G4 gates;
- the tail-G2 on clean, RT-mapped nulls;
- dev arm M.

No individual one of them is evidence of a violation.

---

## MAJOR

### M0. G2 nulls are contaminated by partly-true matches, and G2 lacks the tail test (was a candidate blocker; downgraded after the breakdown)
G2 compares PSM-level target-win rates over *all* null spectra, but FDR is decided in the tail. On the G2 mixture runs (production search, prod index, homolog filter as in `gates.py`), I applied the production estimator: paired double competition, q ≤ 1%. I then counted null-provenance winners whose score is at or above the lowest accepted target (`reviews/checks/g2_null_tail.json`):

| score | null T ≥ thr | null D ≥ thr | binomial p | of which "far" (no shared 6-mer, Hamming > 2): T / D |
|---|---|---|---|---|
| hyperscore | 106 | 91 | 0.32 | 12 / 14 |
| static_linear_base | 141 | 91 | 0.0012 | 15 / 17 |
| static_gbm_base | 179 | 132 | 0.009 | 20 / 18 |
| static_linear_pred | 136 | 96 | 0.010 | 17 / 21 |
| **static_gbm_pred** | **159** | **115** | **0.009** | **19 / 24** (p = 0.54) |

Source: `g2_null_tail_breakdown.json`. Homology is judged on the pair's target sequence, so it is label-symmetric.

What the breakdown shows:
- The whole target excess sits in **yeast nulls whose candidate shares a 6-mer with the yeast/contaminant proteome**. For gbm_pred that is 121 T vs 64 D; even hyperscore shows 89 vs 68.
- Those are partly-true matches: conserved sub-sequences that explain real fragment ladders. A better scorer legitimately lifts them.
- On "far" nulls every scorer is balanced. **So I find no evidence of label leakage in the tail.** However, only ~40 far nulls clear the 1% threshold, so the power to detect a tail leak is low.

The consequences are real nonetheless:
1. G2's null set is contaminated. The §6.4 homolog filter (I/L, single substitution, isobaric) is too weak for yeast vs human, so p_fixed and every G2 difference are computed partly on non-null spectra. That inflates the "learned − fixed" difference for good scorers and can mask a leak.
2. Nothing in the gates looks at the tail where the 1% cut is made.

**Fix (dated §6 note before re-running):**
- Strengthen the null homolog filter: shared k-mer (k = 6) and Hamming ≤ 2, label-symmetric on the pair target.
- Add a **tail-G2**: null target-minus-decoy winners above the production 1%/5% thresholds, with a one-sided bound and a margin fixed now.
- Add enough null runs from several projects that the far-null tail has ≥ ~200 winners at 1%. 3 yeast runs carry 93% of the nulls today, and the G2 SE treats 142k spectra as independent.

### M1. G2 is nearly vacuous for the default scorer, and its nulls do not exercise the RT channel
- **Vacuity.** G2 tests only whether the learned score changes *which candidate wins* a null spectrum. static_gbm_pred changes the winner on **0.19%** of nulls (`gates.json:G2_G5.scorers.static_gbm_pred.frac_winner_changed = 0.0019`), so its G2 diff of +0.0001 is guaranteed whether or not the score leaks. The rank-1-only ablation confirms the effect: 24,518 vs 24,512 IDs (`reviews/checks/estimator_checks.out`).
  - Leakage in a static re-scorer acts on the **between-spectrum ranking**. That is G5 (an AUC over all nulls, diluted) and the tail test in M0.
- **The RT channel is untested.** The species nulls keep their source-run retention times:

  | Spectra | RT q1 / q50 / q99 (min) |
  |---|---|
  | Pfu nulls | 6.6 / 40 / 86 |
  | yeast nulls | 14.6 / 72 / 125 |
  | D0 host run | 30–34 / 108–123 / 198–201 |

  (`reviews/checks/rt_null_asymmetry.json:G2_mixture`). Their RT residuals against the D0 calibration are therefore ~45–57 min for **both** labels, i.e. saturated.
  - In production, null spectra elute inside the host gradient at times consistent with natural (human) peptides, and `rt_pred` separates targets from pseudo-reverse decoys (AUC 0.56–0.61 on all candidates here; 0.71 per M3 on D1).
  - So `rt_abs_resid`/`rt_resid_z` could favour null targets in production, and G2 as built cannot see it. On yeast nulls the RT residual even slightly favours *decoys* (AUC 0.506), which would mask an opposite production effect.
- **Fix:** put each null source's RT on the host time scale before mixing (quantile-map source RT to host RT, or calibrate each source run with its own label-blind calibration set and map through iRT). Report G2/G5/tail-G2 with and without the RT features.

### M2. Required gates G1, G3 and G4 were not run, yet M3 reports "gates pass"
§6.4: "all must pass before any FDR claim".
- `gates.json` contains only G0 and G2_G5. The `gates.py` docstring mentions G1's diagnostic, but no G1 code exists.
- **G1** gate (|Spearman ρ| between scorer output and a sequence-only target/decoy classifier logit, top decile, nulls): not done. It is directly relevant because the predictor is a sequence model and `rt_pred` is excluded for exactly this reason.
- **G3** (retrieval symmetry on nulls, binned by score): K1 top-K (K = 500) is a pruning stage. E2 measured recall on all spectra (target 99.90 / decoy 99.78% on D0). §6.4 asks for null spectra binned by score.
- **G4**: not done (it is needed before any strata are used).
- **Fix:** run G1, G3 and G4 on the G2 mixture, then re-state the gate status.

### M3. The memorization arm (M) is required, not optional, for predicted-feature scorers
§6.3 requires arms P, C **and M** for "models with sequence capacity (… any peptide-token model)". static_gbm_pred's strongest inputs (pred_sa, pred_entropy; D1 AUC 0.97) are outputs of AlphaPeptDeep, a peptide-token transformer trained on large human corpora. The composite scorer has sequence capacity through its features. NOTES ("not required for static low-capacity scorers") contradicts §6.3.

Arm P cannot detect this failure. In D0_E2001 low-score spectra, the shuffled entrapment E and both decoy families look alike on pred_sa (0.229–0.237) and RT residual (50.5–56.3 min; `rt_null_asymmetry.json:D0_E2001`). A "natural-sequence bonus" for real human peptides would favour T over E and decoys together, and arm P would read it as control.

**Fix:** run arm M on dev before T1. D7 yeast (PXD031616, never in training) searched against yeast + human entrapment, r ≈ 1, with the frozen gbm_pred. Use T5 only at the milestone.

### M4. The Sage-LDA reference is handicapped, so "+3.4% / +10.5%" is overstated
- `vulcan/gpu/sage_bridge.py:88` sets `"predict_rt": False`. Sage's default and M1e (`scripts/baselines/run_searches.py:107`) use `True`.
- Rerun with only that change, through our estimator:

  | Set | Sage-LDA, RT off | Sage-LDA, RT on | gbm_pred vs RT-on Sage |
  |---|---|---|---|
  | D0 @1% | 23,703 | **23,836** | **+2.8%** |
  | D0 @0.1% | 21,324 | **21,717** | **+8.5%** |
  | D1 @1% | 6,044 | 6,057 | parity |

- The no-cleavage FASTA bridge also makes Sage's `missed_cleavages` feature constant, and `report_psms=1` makes `rank` constant.
- The remaining gain is almost surely the predicted MS2 features. The fair state-of-the-art comparator for that claim is a predicted-feature rescorer (MS2Rescore/Oktoberfest + mokapot) on the identical DB, not bare Sage-LDA.
- There is no file-level CI for the difference (n = 2 files).
- **Fix:** re-issue sage_lda_ours with predict_rt on. Add Sage+mokapot and an MS2-predicted rescorer on the same index. Quote the difference with a CI, in `uniq seq` units for cross-tool comparisons (the M1e rows are sequences, ours peptidoforms).

### M5. D0/D1 have become a selection set; quote the default's numbers as post-selection
D7 picked gbm_pred over linear_pred and RESET by looking at D0/D1 IDs and entrapment (4 static variants + 2 RESET). Model files predate D0/D1 evaluation (`train_independence` check), so no retraining leaked. Variant selection did, mildly.

The margin used to rank them is noise:
- gbm_pred's "cwc" rests on a P-arm lower CI of 0.0097.
- The CI adds bootstrap variance plus between-seed variance (`vulcan/stats/entrapment.py:251-253`), which double-counts the same binomial noise.
- With the bootstrap-only CI the lower bound is 0.01004, so cwc is lost (`estimator_checks.out`).
- At 5% on D0 all four static scorers fail cwc on arm P (gbm_pred mean 0.0536, lower CI 0.0504; hyperscore passes). NOTES says "P arm cwc on both datasets" without that qualifier.

**Fix:**
- Freeze exactly one scorer (`static_gbm_pred.pkl`, sha256 31bfe685…) for T1, and label the others exploratory.
- Report both CI variants.
- Do not use cwc to rank scorers.
- Correct the NOTES wording.
- Pre-register (dated) whether the between-seed term is added or replaces part of the bootstrap. Simulation (§6.7) should decide this; the result it produces should not.

### M6. Speed numbers in M3 and E2 are not the production configuration, and the headline ratios are favourable
- M3's "0.43 s warm features / ~2.1 s per run":
  - `scripts/m3/bench_scoring.py` fixes the library NCE (29). Production chose NCE 23 on D0, which gave 3,583 / 147,151 MS2 cache hits and 0 RT hits, so the real feature stage was ~27 s.
  - The RT cache is keyed by NCE even though RT does not depend on NCE (`vulcan/predict/cache.py:22-25`).
  - Library load (6.9 s/process), precompute (328 s) and start-up are excluded.
- E2's "7.5 s / 4.9 core-s, 4.0× wall, 48× core-s" is computed: process minus repeated warm passes (`scripts/m2/e2_summarize.py:70-75`). The measured single-process numbers are 8.59 s / 6.0 core-s (3.5× and 40×).
  - Sage's wall includes its index build (~7.8 s) and MGF parsing. Ours excludes the L1 build (40 s prod, 2 min with entrapment) and RAW→store ingest.
  - Sage-32 core-s on a node at loadavg 14–22 is inflated by rayon spinning (Sage-4 used 114–160 core-s).
- The two locked `io1` pipelined runs both stalled (217/227 s) and were replaced by an unlocked rerun. Two out of two is a reproducible behaviour, not "environmental".
- **Fix:**
  - Benchmark the exact production path (per-run NCE). Key the RT library independently of NCE, and precompute the NCE grid or snap to the library NCE.
  - Report measured (not derived) process numbers, with a "including index build" row.
  - Use Sage at 4–8 threads as the core-s reference, and give ranges.
  - Rerun under the lock on a quiet node.

### M7. The T4 selector does not mechanically exclude the training projects
`scripts/data/select_t4.py` never reads `data/manifests/train.tsv`. The exclusion exists only as prose in DATA.md, and a failing project is replaced "from the same cell", which could draw a training project.

**Fix:** add train.tsv to the selector's exclusions, and make `bench/sealed_run.py` assert T4∩train = ∅.

---

## MINOR

1. **Arm units.** Arm P pairs by target *sequence* (`entrap_pair`, `vulcan/index/build.py:~203`), but IDs are peptidoforms. Arm C uses peptide-level r (D0 1.022 vs 1.096 for peptidoforms; D1 12.94 vs 13.54). The §6.3 per-mass-decile r is not used, although entrapment hits concentrate in the low-mass decile (D0 decile-1 r = 0.77). The per-decile estimator gives gbm_pred C = 0.82% vs 0.78% reported. State the units and report the per-decile estimator.
2. **G5.** Every learned scorer sits ~2σ above 0.5 (0.503–0.504, SE ≈ 0.0015; hyperscore 0.4996). This is below the alarm but detectable. Report its CI and a top-decile G5.
3. **Verdict precedence.** `entrapment.py` checks "not evaluable" before "violated", so a short list can hide an LB violation. Check "violated" first.
4. **Offset nulls are near-empty by design** in closed search (595 of 177k): an n+0.5 Da shift lands in the forbidden mass-defect zone at ±10 ppm. This was predictable, and is not a result-driven discovery.
5. **Lazy vs library predictions.** G0 used a fresh cache on each side, so it never compares a library hit against a lazy prediction of the same peptidoform. In production, T/dec(T) can come from the bulk library while E/dec(E) are lazy (different batch composition). Add a test that lazy and library predictions agree to fp16 ulp.
6. **hs_runz** is computed per search, so augmented (P/C) searches use a slightly different normalization than production. It is label-symmetric and so valid, but note it under §6.3 "Reporting".
7. **Training data** (9/12 human, HEK-like lines) is close to D0. Gains on D0 may be optimistic for other species and instruments. The training PSM/feature tables live only in /tmp; commit their digests.
8. The `scripts/baselines` changes (D4 cut to 2 files) and `vulcan/gpu` edits are uncommitted. Commit them with a note before any milestone tag (sealed_run requires a clean tree anyway).

## What checks out (no action)
- **Label blindness of L0–L3 by construction.**
  - Labels live only in `labels.arrow`.
  - ids are xxhash64 of the canonical peptidoform; rows are sorted by (f32 mass, id).
  - Top-K ties go to the higher row, never to a label.
  - `missed_cleavages`, length and n_var are pair attributes copied from the forward target (`vulcan/index/build.py:159-164`), so pseudo-reverse decoys cannot differ on them.
  - The prediction cache is keyed by content (pf_id) plus (model, instrument, NCE, mods).
  - NCE choice and RT calibration use rank-1 PSMs of both labels.
- **G0** is a real test of positional and flag leaks: shuffled FASTA, permuted pairs, swapped roles, fresh caches, bit-identical outputs. By design it cannot detect *content* leaks (reversed sequences look different); that is G1/G2/G5's job, hence M0/M1/M2.
- **Estimators.**
  - `double_competition` follows §4 L4 (hypothesis → PSM → own-pair; (D+1)/T, monotone, decoys first at ties, hash tie-break).
  - The paired FDP formula matches FDRBench.
  - RESET uses c_train 3 / c_final 2 with +1 inside.
  - The verdict table matches §6.3.
  - `tests/stats`: 60 passed.
- **Training independence.**
  - The 13 training accessions overlap none of D*/T*, the T4 seed-11 selection, the predictor manifests, D7 yeast or the T1 source. No yeast, no shared labs, no sha256 overlap; spectrum-hash overlap is 4 of 457k.
  - The fold map is project-disjoint and actually used for hyper-parameters.
  - Static models predate any D0/D1 static evaluation.
- **Re-ranking top-5 vs Sage's top-1 / K=50** does not explain the gain (rank-1-only ablation: same IDs).

## M3's proposed §6 changes

| # | Proposal | Ruling |
|---|---|---|
| 1 | Make wrong-species nulls primary; offset nulls widened or dropped | **Accept, in modified form.** The near-emptiness is a structural property of closed search, not a result. Keep offset nulls as specified (they do no harm). Add, by a dated note before re-running, (a) RT mapping of species nulls to the host gradient (M1), (b) the stronger homolog filter, and (c) more null projects. Do **not** drop a source after seeing that the gate passes on the other. |
| 2 | Peptide/pair-level G2 for RESET | **Accept, and extend it to the static scorers** as the tail-G2 in M0, with margins fixed before re-running. It makes the gate stricter, so it is not result-driven. |
| 3 | Report expected LB α·r/(1+r) next to arm C for large r | **Accept as a diagnostic only.** It must not change any verdict. It was motivated by the uncomfortable D1 linear_pred number (LB 1.17%), so any use to soften a verdict would be result-driven and is rejected. |
| 4 | State that a static scorer may re-rank the top-K if G2 passes | **Accept the clarification**, conditioned on the tail-G2/G1/G3 gates, not the PSM-level G2 alone (which is vacuous for re-ranking that seldom changes winners, M1). |

**Not proposed by M3, and to be rejected if raised:**
- dropping the P arm or downgrading the cwc rule because "everyone is inconclusive";
- treating "consistent with control" as "supported";
- choosing the CI variant that keeps cwc for the default scorer after seeing which one does (M5).

## Minimum path to a justified T1 look
1. B1: T1 driver (5 decoy × 3 sealed-entrapment family-max, comparators on identical external-decoy DBs, ratio CI, matched-FDP), dry-run on D0/D1; finish M1e D0.
2. Dated §6 note (before re-running): tail-G2 with margins, RT-mapped and more diverse species nulls, Hamming-2/6-mer homolog filter, CI convention, arm-C diagnostic.
3. Re-run G1/G2(tail, M0)/G3/G4/G5 and dev arm M (D7 yeast + human entrapment) on the frozen `static_gbm_pred`.
4. Correct M3 NOTES and RESEARCH_TREE wording (Sage RT, cwc at 5%, speed configuration, units).
5. Clean tree; milestone tag; T1 via `bench/sealed_run.py`.
