# M3a / E4: rescoring with predicted features (static-first scorers + RESET arm), dev only

Data: D0 = PXD001468 (2 HEK files, QE HCD), D1 = PXD001077 (1 Pfu file, Velos HCD/OT). These are development curves, not claims (DESIGN §6.6). No sealed file was touched. §6 was not changed; proposals are at the end.
Files: `metrics.json` (everything), `eval.json` (full entrapment output), `gates.json`, `bench_scoring.json`, `predictor_library.json`, `sage_lda_ours.json`, `models/` (frozen scorers + `training_log.json`), `config.yaml`, `git_sha`.

## What was built
- **Predictor** (`vulcan/predict/apd.py`): AlphaPeptDeep generic MS2 + RT (peptdeep 1.5.1, `pretrained_models_v3.zip`, sha256 345cb241…). It runs in `conda_envs/pred`, a venv over `pe` with peptdeep installed `--no-deps`; peptdeep pins numpy<2/pandas<3, and its own `predict()` crashes under pandas 3 copy-on-write.
  - We feed the torch models directly. Inputs are built on the GPU from the L1 index layout (residues u8 + var-mod bitmask), length-bucketed exactly like peptdeep.
  - Output vs peptdeep's own featurization: RT is identical, MS2 differs by ≤2.4e-4 (fp16 storage). Tests: `tests/predict`.
  - The pandas path manages 3.4k peptides/s (RT) and 11k precursors/s (MS2); the fast path is about 20× faster.
- **Decision D4** is enforced: `check_allowed` refuses D2/T2 and PXD009449/PXD000138.
- **Leakage check:** no D0/D1/training-run or D7 accession appears in `predictor_training_accessions.tsv`.
- **Cache** (`vulcan/predict/cache.py`):
  - Content-hashed by (model sha, instrument, NCE, mods) and keyed by peptidoform id. Because the key is content-defined, one library serves prod, P and C indexes and decoys.
  - Append-only shards; misses are predicted lazily. Hit rates are in `metrics.json:feature_runs`.
  - Human tryptic library: 7.06M peptidoforms × z = 2, 3 at NCE 29 (QE class) → 2.6 GB.
- **L3 features** (`vulcan/rescore/features.py`, GPU, label-blind):
  - Fragment features: predicted b/y (z1, plus z2 when the precursor is ≥3+) vs all observed centroids (±20 ppm). Spectral angle (raw and sqrt), Pearson r, entropy similarity, predicted-intensity fraction observed, top-6 observed, observed/TIC.
  - RT: |residual| and its robust z after a per-run monotone calibration. The calibration is fitted on the top-decile-hyperscore rank-1 PSMs of **both** labels.
  - NCE: chosen per run from {20,…,35} by the median spectral angle of the same label-blind set. Results: D0 23, D1 29, training runs 26–32.
  - `rt_pred` itself is **excluded** from the models. It is a pure sequence property with target/decoy AUC 0.71 on D1, so it is a G1-type leak risk.
  - Diagnostic AUCs on D1 rank-1 PSMs: hyperscore 0.92; pred_sa_sqrt 0.97, pred_entropy 0.97, |RT resid| 0.80.
- **Search:** the E2 GPU search with `report_psms=5`, so static scorers can re-rank a spectrum's top-5.
- **Arm-C index:** `vulcan/index/foreign.py` is a copy of `build()` (commit 166e6de) plus M1e's neighbour-filtered foreign peptides as family E. It is kept separate because M2b is concurrently editing `build.py`; fold it in later.
- **Static-first scorers** (`vulcan/rescore/static.py`, `scripts/m3/train_static.py`): frozen, per instrument class (OT-HCD only).
  - Variants: {standardized logistic regression, HistGBM with depth 3 / 8 leaves / min leaf 200} × {base = K2 label-blind features, pred = base + predicted}.
  - **Training data:** 13 runs from 13 PRIDE projects (9 human, 2 mouse, 2 E. coli), selected with seed 11 from the T4 frame. Excluded: T4 selection, D*/T*, predictor training sets, yeast, Arabidopsis. Details in `DATA.md` (role `train`), `data/manifests/train.tsv`, `train_folds.tsv` (3 project-disjoint folds) and `train.sha256`.
  - One run (PXD019362 `MN-WT1-Frac8`) turned out to be IT-CID and was excluded automatically. That leaves 12 runs and 829k PSMs.
  - **Labels:** Percolator-style target-decoy on the *training* runs only.
    - Positives: rank-1 targets at PSM q ≤ 1%.
    - Negatives: rank-1 decoys plus rank≥2 candidates of positive spectra.
    - 3 iterations.
  - **Hyper-parameters:** leave-one-project-fold-out. Then the model is refit on all folds and frozen; it is never fit on D0/D1.
  - Training-set PSMs at 1%: hyperscore 137,119 → linear_pred 174,949, gbm_pred 174,149.
- **RESET arm:** `vulcan.stats.reset` exactly per §6.4: single-decoy, seed 11, LDA, init = hyperscore, run on the hyperscore double-competition winners.
  - Pooled peptide level only; no protein FDR.
  - Fallback: 0/20 RESET fits (2 datasets × 2 feature sets × 5 searches) needed the α_train = 5% retry or the fixed-score fallback.

## IDs (peptidoform units, paired double competition; q ≤ 0.1 / 1 / 5 %) and dev entrapment verdicts

Sage reference row: Sage 0.14.7's LDA discriminant from E2's run on the identical index and spectra, put through **our** estimator (`sage_lda_ours.json`). Sage's own peptide_q gives 23,637 (D0) and 6,037 (D1).

Verdict notation: the overall §6.3 verdict at α = 1% / 5%, followed by "cwc" (paired lower CI ≤ α, the weaker "consistent with control" claim). A scorer marked "cwc only on C" means the P-arm lower CI exceeds α: no claim.

| D0 | 0.1% | 1% | 5% | uniq seq @1% | prot @1% | verdict 1% (P / C) | verdict 5% |
|---|---|---|---|---|---|---|---|
| hyperscore | 12,635 | 19,600 | 23,640 | 18,615 | 5,560 | inconclusive (P .0109 cwc / C .0100) | supported |
| Sage LDA (our TDC) | 21,324 | 23,703 | 25,423 | 22,410 | – | M1e D0 pending | – |
| static_linear_base | 18,686 | 23,369 | 25,447 | 22,142 | 6,245 | inconclusive (P .0113 cwc / C .0085 sup) | inconclusive |
| static_gbm_base | 21,157 | 23,695 | 25,173 | 22,419 | 6,342 | inconclusive (P .0119 no-cwc / C .0091 sup) | inconclusive |
| **static_linear_pred** | 22,600 | **24,736** | **26,127** | 23,387 | **6,477** | inconclusive (P .0121 no-cwc / C .0073 sup) | inconclusive |
| **static_gbm_pred** | **23,566** | 24,512 | 25,700 | 23,162 | 6,455 | inconclusive (P .0112 cwc / C .0078 sup) | inconclusive |
| reset_lda_base | 17,028 | 22,922 | 25,022 | 21,718 | n/a (§6.4) | inconclusive (P .0124 no-cwc / C .0078) | inconclusive |
| reset_lda_pred | 23,235 | 24,535 | 25,747 | 23,172 | n/a | inconclusive (P .0122 cwc / C .0074) | inconclusive |

| D1 | 0.1% | 1% | 5% | prot @1% | verdict 1% (P / C) |
|---|---|---|---|---|---|
| hyperscore | 5,303 | 5,768 | 6,297 | 800 | **supported** (P .0102 / C .0062) |
| Sage LDA (our TDC) | 5,883 | 6,044 | 6,299 | 834 (M1e) | M1e: inconclusive (P inconcl / C sup) |
| Sage+mokapot (M1e, seq units) | 5,680 | 5,899 | 6,157 | 837 | M1e: inconclusive |
| static_linear_base | 5,860 | 6,061 | 6,299 | 842 | inconclusive (P .0122 / C .0094 sup) |
| static_gbm_base | 5,806 | 5,960 | 6,298 | 828 | inconclusive (P .0122 / C .0093 sup) |
| static_linear_pred | **5,992** | **6,089** | 6,309 | 843 | inconclusive (P .0109 / **C .0126 inconclusive**, LB .0117) |
| static_gbm_pred | 5,917 | 6,052 | 6,297 | 838 | inconclusive (P .0107 / C .0095 sup) |
| reset_lda_base | 5,705 | 5,956 | 6,297 | n/a | inconclusive (P .0123 / C .0111) |
| reset_lda_pred | 5,868 | 6,062 | 6,297 | n/a | inconclusive (P .0114 / C .0125) |

How to read the verdicts:
- **No scorer is "violated".** For every scorer and level, the lower 95% CI of the mean LB is ≤ α.
- **At 1%, every learned scorer is "inconclusive".** This includes Sage+LDA and Sage+mokapot in M1e, and the hyperscore on D0. The paired shuffled arm sits at about 1.1–1.2% for everyone.
  - That is the known ~20% overstatement of the paired estimator (§6.3).
  - The natural-proteome arm C is well inside α on D0 for all predicted-feature scorers (0.73–0.78%).
- **D0, lower CI of P > 1%** for static_linear_pred, static_gbm_base and reset_lda_base. For these, "consistent with control" cannot be claimed. static_gbm_pred and reset_lda_pred keep cwc.
- **D1 arm C (r = 12.9 human entrapment in a Pfu sample):**
  - static_linear_pred and reset_lda_pred accept 71 random human peptides among about 6,080 IDs: LB 1.17%, CI [0.90, 1.43]%. Hyperscore accepts 31 (0.58%) and Sage-LDA 0.70%.
  - The accepted human peptides are singletons from unrelated genes, i.e. genuine false positives, not a contaminating human protein.
  - This is within about 2σ of what a tight 1% list would give (expected N_E ≈ 0.93·α·N), so it is not a violation. Still, it is the least comfortable number in this report.
  - static_gbm_pred stays at 0.88% (supported on C).
- **Protein level (picked groups, arm P entrapment proteins):** D0 is inconclusive at 1% for all scorers (paired 1.1–1.4%) and supported at 5% for all. D1 is not evaluable at 1% (α·N < 20) and supported at 5%.
  - Protein IDs at 1%: D0 5,560 (hyperscore) → 6,455 (gbm_pred, +16%) / 6,477 (linear_pred); D1 800 → 838 / 843 (Sage 834).
- **E4 (value of predicted features), same estimator, D0 at 1%:**
  - base → pred adds +1,367 (linear), +817 (GBM) and +1,613 (RESET).
  - At 0.1% the gain is large: +3,914 for linear and +6,207 for RESET.
  - Versus Sage+LDA on identical PSMs: +4.4% (linear_pred) and +3.4% (gbm_pred) at 1%; +10.5% (gbm_pred) at 0.1%.
  - On D1 the gain over Sage is small (+0.1 to +0.7% at 1%; +1.9% for linear_pred at 0.1%).
- **D4: not feasible this pass.** It is TMT11 SPS-MS3 with **IT-CID MS2** (M1e). That needs an IT-CID static model and TMT mods in the index; neither exists, and we have no IT-CID training runs yet.

## Gates (`gates.json`)
- **G0 PASS.** Setup: D1 searched against the E2001 index and against E2's swap index (FASTA order shuffled, pair list permuted, target/decoy roles swapped).
  - Bit-identical: the top-5 PSM tables (33,129 rows), the chosen NCE, all 10 L3 features, and all 4 frozen static scores.
  - Each side used its own fresh prediction cache.
- **G2 PASS** for all 4 static scorers (class OT-HCD, one-sided 95%).
  - **Null set:** n = 142,501 null competitions, above the `g2_sample_size` budget of 107,570.
    - Sources: 595 offset, 9,541 Pfu and 132,365 yeast (D7, PXD031616 LM2480/81/86, Lumos OT-HCD; dev data, not training).
    - Mixture design: added into the two D0 runs.
    - Homolog filter removed 1,082 Pfu and 20,702 yeast spectra.
  - **Fixed target-win rate:** p_fixed = 0.4944.
  - **Learned − fixed** (diff / UCB; budget 0.01 overall, 0.02 top decile):

    | Scorer | Overall diff / UCB | Top-decile diff / UCB |
    |---|---|---|
    | linear_base | +0.0011 / 0.0023 | +0.0013 / 0.0038 |
    | gbm_base | +0.0003 / 0.0005 | +0.0011 / 0.0017 |
    | linear_pred | −0.0024 / −0.0002 | +0.0063 / 0.0121 |
    | gbm_pred | +0.0001 / 0.0003 | +0.0010 / 0.0017 |

  - linear_pred re-ranks the winner on 24% of null spectra; gbm_pred on 0.2%.
  - **RESET:** its PSM and peptide competitions are never redone (§6.4), so its G2 difference is 0 by construction. The test is vacuous for RESET (see proposals).
  - **The offset source is nearly empty.** n + 0.5 Da shifts leave almost no candidate within ±10 ppm: 595 competitions from 177k spectra. G2 therefore rests on wrong-species nulls.
- **G5 PASS.** Null winner AUC by learned score:

  | Scorer | AUC |
  |---|---|
  | linear_base | 0.504 |
  | gbm_base | 0.503 |
  | linear_pred | 0.502 |
  | gbm_pred | 0.503 |
  | RESET-pred | 0.501 |
  | hyperscore (reference) | 0.4996 |

## Runtime (A100 80GB + 32 Zen4 cores; `bench_scoring.json` under the bench lock, `predictor_library.json`)
Per D0 run: second file, 46,816 MS2 and 152,335 top-5 PSMs, steady state (warm-up excluded). Wall / CPU core-s:

| Stage | Wall | Core-s | Notes |
|---|---|---|---|
| Predictor inference, cold (no cache: RT + MS2 for 128.9k precursors, lazy) | 2.30 s | ~2.3 | ≈ 56k precursors/s incl. RT |
| Feature assembly, cold total (peaks + fragment matching + RT calibration + inference) | 2.78 s | 2.83 | |
| Feature assembly with the **resident** tryptic library (100% hits) | **0.43 s** | 0.43 | fragment features 0.16 s, RT 0.03 s, peak upload 0.25 s |
| Library load into a fresh process (2.6 GB, 14.1M precursors) | 6.9 s | 6.9 | once per process; a persistent worker amortises it |
| NCE grid choice (6 NCEs × 2,000 PSMs) | 1.05 s | 1.05 | |
| Feature matrix | 0.04 s | 0.04 | |
| Static linear apply, GPU / CPU | 0.006 s / 0.016 s | – | |
| Static GBM apply, CPU (sklearn) | 0.64 s | 1.29 | |
| Paired double competition | 0.08 s | 0.16 | |
| RESET-LDA (5 iterations) | 0.07 s | 0.14 | |

- **Comparison:** for the default scorer, scoring with the resident library plus NCE choice plus GBM apply costs about 2.1 s wall / 2.8 core-s per run. With the linear model on GPU it is about 1.5 s. Ingest is ~8.3 core-s per run (E1: 17.8 core-s / 10⁵ MS2), so scoring is well below ingest. It is still ~10× the warm GPU search itself (0.12 s per run, E2).
  - The next wins are a GPU GBM (e.g. a treelite/FIL-style kernel), caching the NCE choice per instrument, and fusing peak upload with K0, which uploads the same peaks.
- **Predictor throughput, A100:**

  | Condition | RT (peptides/s) | MS2 (precursors/s) |
  |---|---|---|
  | Clean, single process, 500k random human peptidoforms | 211k | 65k |
  | Under GPU contention (`predictor_library.json`) | 122k | 48k |

  - Full human library precompute (7.06M peptidoforms × z = 2, 3): 328 s wall, under contention.
  - The first call per peptide length costs ~4–50 s of kernel set-up; it is excluded above.
- **Process start-up is environmental:** Python import 179 s and peptdeep import + model load 272 s in this locked run. Lustre-hosted env, node loadavg 10–16. A persistent worker is required (as found for E2).

## Recommendation
**Default scorer: `static_gbm_pred`** (frozen OT-HCD GBM over label-blind and AlphaPeptDeep features).
- **IDs:**
  - D0: 24,512 at 1% (+25% vs hyperscore, +3.4% vs Sage-LDA under the same estimator) and the best at 0.1% (23,566, +10.5% vs Sage).
  - D1: 6,052 (≈ Sage).
  - Proteins: +16% on D0.
- **Statistics are simple:** it is static, so it needs no RESET. Plain paired TDC, picked protein and global FDR share one score.
- **Gates and entrapment:** best behaviour among the learned scorers.
  - G2 diff ≈ 0 and G5 0.503.
  - Arm C supported on both datasets.
  - P arm consistent-with-control on both datasets.
- **linear_pred** has the most IDs at 1%, but it fails the cwc rule on D0 arm P and is inconclusive on D1 arm C. Keep it as an alternative, not the default.
- **RESET-LDA-pred** gives similar IDs but only peptide-level FDR. It stays an arm (E4): it does not beat the static scorer at verified FDP.
- **Before any claim:**
  - Sealed T1 at the M3 milestone.
  - The M arm (T5 yeast) is not required for static low-capacity scorers, but cheap to add.
  - An IT-CID model is needed for D4/Lumos-IT data.

## Caveats
- The node was shared: loadavg 10–16 and other agents' GPU jobs (M2b cohort). The library precompute ran **unlocked-concurrent** with feature jobs, so its throughput numbers are under contention. Clean single-process numbers are 211k RT peptides/s and 65k MS2 precursors/s.
- The first `python import` from the Lustre env took 200–300 s several times (E2 saw 4 s). Startup is environmental; a persistent worker is needed (M2 2b).
- D0 RT calibration MAD is about 9–11 min on a 180-min gradient (Spearman 0.93 on the calibration set). It is not fine-tuned; RT transfer learning on label-blind calibration PSMs is a cheap next step.
- Units: ours and Sage-LDA (our TDC) are peptidoforms; M1e reports peptide sequences. The `uniq seq` column is the comparable one (D0 gbm_pred 23,162 vs Sage-LDA 22,410).
- The training runs are in the T4 *frame*, though not in the current selection. Any future T4 redraw must exclude them (recorded in DATA.md).

## Proposed changes (NOT applied; §6 frozen at prereg-v0.4)
1. **G2 nulls.** With ±10 ppm precursor tolerance, n+0.5 Da offsets yield ~0.3% usable competitions. Propose making wrong-species nulls (with the homolog filter) the primary source, and either:
   - running offset nulls with a wide precursor window (open-search style); or
   - dropping them.
2. **G2 for RESET is vacuous**, because competition is never redone. Propose a peptide-level G2 variant for RESET: the null target-win rate among reported pair winners at q ≤ α under the learned score vs the fixed score.
3. **Arm C for small proteomes** (D1, r = 12.9). With many foreign peptides per target, LB ≈ FDP. Propose reporting the expected LB under exact control (α·r/(1+r)) next to the verdict.
4. **Static-scorer re-ranking.** The frozen model scores all top-5 candidates, so the PSM-level winner can change (this is what G2 tests). State explicitly in §4 L3 that a static scorer may re-rank the top-K, provided G2 passes.

## Erratum (2026-10-27, added by M3e after the independent review; the original text above is left unchanged)
Source: reviews/m3_independent_review.md (M1, M3, M4, M5, M6, minor 1). New measurements: results/2026-10-27_M3e_hardening/.
None of the numbers above has been rewritten. Where a number is wrong or mislabelled, the correction is given here.

1. **The Sage-LDA reference ran with Sage's RT prediction OFF.** `vulcan/gpu/sage_bridge.py:88` sets `"predict_rt": False`; Sage's default and M1e use `True`. Every "Sage LDA (our TDC)" row and every "vs Sage-LDA" gain above therefore refers to a handicapped comparator. With `predict_rt: true` (review rerun, same index, spectra and estimator):

   | Set | RT off (above) | RT on | static_gbm_pred vs RT-on |
   |---|---|---|---|
   | D0 @1% | 23,703 | 23,836 | +2.8% (not +3.4%) |
   | D0 @0.1% | 21,324 | 21,717 | +8.5% (not +10.5%) |
   | D1 @1% | 6,044 | 6,057 | parity |

   - The no-cleavage FASTA bridge makes Sage's `missed_cleavages` feature constant, and `report_psms=1` makes `rank` constant.
   - The fair comparator for a predicted-feature claim is a predicted-feature rescorer on the comparator's own PSMs (A1.9; tree item 3f). It has not been run, so no gain over the strongest comparator is established yet.
   - n = 2 files, so there is no CI on the difference.
2. **"Consistent with control" (cwc) at 5% and the CI convention.**
   - The Recommendation's "P arm consistent-with-control on both datasets" holds only at 1%. At 5% on D0, all four static scorers fail cwc on arm P (static_gbm_pred mean 0.0536, lower CI 0.0504); hyperscore passes.
   - At 1%, static_gbm_pred's cwc on D0 holds only under the CI convention that adds between-seed variance: lower CI 0.0097 under (a), 0.01002 under bootstrap-only (b).
   - Per A1.6, cwc is never used to rank scorers, so "best entrapment behaviour among the learned scorers" is not supported as a ranking criterion.
   - The convention chosen by the pre-registered simulation, and the D0/D1 verdicts recomputed under it, are in results/2026-10-27_M3e_hardening (NOTES, `d0d1_ci_conventions.json`).
3. **Arm M is required, not optional.** The Recommendation says "The M arm (T5 yeast) is not required for static low-capacity scorers". That contradicts §6.3. static_gbm_pred's strongest inputs come from AlphaPeptDeep, a peptide-token model, so the scorer has sequence capacity through its features; arms P, C and M are required (A1.5). Dev arm M (D7 yeast) has now been run in M3e.
4. **The speed numbers are not the production configuration.**
   - `scripts/m3/bench_scoring.py` fixes the library NCE at 29. Production chooses the NCE per run (D0: NCE 23), which gave 3,583 / 147,151 MS2 cache hits and 0 RT hits, so the real feature stage on D0 was ~27 s, not 0.43 s.
   - The RT cache is keyed by NCE, although RT does not depend on NCE.
   - Library load (6.9 s/process), library precompute (328 s) and process start-up were excluded.
   - "Scoring ... about 2.1 s wall / 2.8 core-s per run" is therefore a best case for a cache-resident NCE, not the production cost.
   - The measured production-configuration cost (per-run NCE choice, first-contact cache, persistent process) is in results/2026-10-27_M3e_hardening (`bench_prod.json`, NOTES).
   - The E2 figures quoted in the tree ("7.5 s / 4.9 core-s, 4.0x / 48x") were derived (process minus warm passes). The measured single-process figures are 8.59 s / 6.0 core-s (3.5x wall, 40x core-s), and they exclude the L1 build and RAW ingest.
5. **Units.**
   - Ours and "Sage LDA (our TDC)" count peptidoforms; M1e rows (e.g. "Sage+mokapot (M1e, seq units)") count peptide sequences. For cross-tool comparisons, use the `uniq seq` column.
   - Arm P pairs by target sequence while IDs are peptidoforms.
   - Arm C used peptide-level r (D0 1.022; peptidoform-level 1.096; D1 12.94 vs 13.54), and the §6.3 per-mass-decile r was not used (D0 decile-1 r = 0.77). The per-decile estimator gives static_gbm_pred C = 0.82% at 1% on D0, not 0.78%.
6. **"Gates pass" overstated what was tested.**
   - Only G0, G2 and G5 were run; G1, G3, G4 and arm M were not.
   - G2 is nearly vacuous for static_gbm_pred, which changes the winner on 0.19% of null spectra.
   - The G2 nulls were contaminated by near-homologous yeast spectra and kept their source-run RTs (RT channel untested).
   - G5 for every learned scorer was ~2σ above 0.5 (detectable, below the alarm).
   - The A1 re-run of all gates is in results/2026-10-27_M3e_hardening.
