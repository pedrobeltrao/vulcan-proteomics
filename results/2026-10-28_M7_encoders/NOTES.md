# M7: split-then-joint encoders: E2b (size), E6 (cross-encoder reranker), E5 (bi-encoder retrieval). Dev data only

Only dev data were used: training-role runs, D0, D1, D7 (arm M), the D8 holdout (first 10 runs) and the existing dev nulls
(M3e G2 mixtures). Nothing under data/sealed was read, and DN2 data were not touched. §6, A1, A2, the metrics, labels
and estimators are unchanged, and the frozen T1 scorer (static_gbm_pred.pkl, sha256 31bfe685…) is only read.

Sections 0–2 were written and committed **before any model was trained or any dev set was scored**. The commit
history is the record.

## 0. Hardware and conditions
- Hardware: gpu-node, AMD EPYC 9354 (32 affinity cores) and an A100 80GB PCIe. hwinfo is in `e2b.json:hwinfo` and
  `metrics.json:hw`.
- The node was shared throughout, with loadavg 30–70. Other agents' jobs held the GPU at 50–100% utilization at the
  start of most E2b sections (recorded per row as `gpu_before`), so all throughputs are lower bounds under
  contention.
- Each timed section held `/tmp/vulcan_bench.lock` for a few seconds only.

## 1. E2b: throughput vs parameter count; size decision (fixed before training)
`scripts/m7/e2b_bench.py`, `e2b.json`. The architecture is `vulcan/encoders/model.py`:
- **Split stage:** a spectrum encoder over 63 peaks plus a charge token, and a peptide encoder over residue tokens
  with Fourier b/y mass features. Pooled 128-d outputs are the E5 embeddings.
- **Joint stage:** cross-attention of peptide tokens to spectrum tokens, with a learned fragment-match bias
  (b1/y1/b2/y2 vs peaks, 20 ppm).

The forward pass is the full split + joint stage in inference mode, on synthetic but shape-realistic batches of
8,192 candidates.

| size | params | dtype | K=5 spectra/s | K=10 | K=20 | K=50 | candidates/s | split: spectra/s, peptides/s |
|---|---|---|---|---|---|---|---|---|
| xs | 0.29M | bf16 | 71,225 | 36,525 | 18,584 | 7,488 | 3.6–3.7e5 | 1.0e6, 1.4e6 |
| **s** | **1.66M** | bf16 | **41,859** | 21,917 | 11,256 | 4,539 | **2.1–2.3e5** | 4.3e5, 6.3e5 |
| m | 5.5M | bf16 | 19,566 | 10,226 | 5,216 | 2,108 | 1.0e5 | 2.3e5, 3.3e5 |
| l | 22.5M | bf16 | 10,201 | 5,419 | 2,798 | 1,138 | 5.1–5.7e4 | 8.7e4, 1.3e5 |
| xl | 50.4M | bf16 | 6,495 | 3,457 | 1,789 | 731 | 3.2–3.7e4 | 5.1e4, 7.9e4 |

What the table shows:
- FP16 is within ±3% of BF16 at every size.
- INT8 was not measured: this torch build has no GPU int8 path for these layers (that would need TensorRT).
- Cost is linear in K. The spectrum encoder is negligible, and the joint stage dominates.
- One implementation fix was needed: the SDPA fast path requires the key length to be a multiple of 8, so we use
  63 peaks plus 1 charge token, and peptides are padded. Without it, the masked cross-attention fell back to a kernel
  2.2× slower.
- `torch.compile` is unavailable on this node (Triton cannot build its helper: gcc header path).

**Pipeline rate.** From results/2026-10-27_scoring_speed, persistent process, library resident:
- search takes 0.13 s per run (median) and scoring 0.60 s, so 0.80 s end to end;
- the page-cached cohort runs at about 10.5k runs/h, decode-bound, with the GPU otherwise about 90% idle;
- the median run has 48.8k MS2 spectra, and the production top-5 table has about 3.5 candidate rows per spectrum
  (D0: 147,151 / 41,820). That is about 1.7e5 candidates per run.
- PRIDE streaming runs at about 580 runs/h (download-bound, E9).

The page-cached cohort therefore feeds about 1.4e5 spectra/s and streaming about 8e3 spectra/s.

The DESIGN's assumption (about 2M parameters at 5×10⁴ spectra/s with K=20) is about 4.5× optimistic for this unfused
PyTorch implementation on a contended GPU. The measured s model does 1.1×10⁴ spectra/s at K=20.

**Decision rule, fixed now:** use the largest size whose reranker GPU time at the production candidate set (top-5,
about 1.7e5 candidates per median run) is **≤ 1.0 GPU-s per run**. That is at most about one extra copy of the
current search + scoring cost. Under that budget the GPU stage stays above about 3.5k runs/h, which is more than 6×
the streaming rate, and the GPU fraction of a streamed run stays below 20%.
- xs: 0.46 s/run;
- **s: 0.78 s/run, chosen;**
- m: 1.7 s;
- l: 3.1 s;
- xl: 5.3 s.

**Model of record: size s (1.66M parameters), bf16 inference, K = 5 (the same top-5 table the frozen static scorer
re-ranks).** K = 20 costs about 3 s/run and is reserved for streaming-only deployments. The 20–40M class (l, 22.5M)
would fit only the streaming rate, and with about 40 training runs it is not expected to be data-supported. It is not
trained in this pass.

## 2. Pre-stated go/no-go criteria (written before training and before any dev scoring)

### Training (task 2)
- **Runs:** the 12 OT-HCD M3 training runs plus 30 new training-role runs. These were drawn by
  `scripts/m7/select_train_m7.py` (seed 23, same frame and filters as M3) and exclude every accession in any
  manifest: T4, D*/T*, DN, DN2 walk, D8, IT-CID sets and predictor training.
- **Project-disjoint 5-fold map:** `data/manifests/train_m7_folds.tsv`. Fold 4 is the only validation fold, used for
  early stopping and LR choice. No dev set (D0/D1/D7/D8/nulls) is used for any choice.
- **Labels:**
  - Positives are PSM winners that are targets, whose peptidoform is accepted at peptide-level q ≤ 0.1% under the
    frozen static_gbm_pred with paired double competition per run, and whose PSM score is at least that list's
    threshold.
  - Negatives are **incorrect targets only**:
    - (1) other target rows of the spectrum's top-20, I/L-different;
    - (2) shifted-precursor: confident target peptides of other spectra at ±5–60 Da, so every confident peptide also
      appears as a negative (prior balancing);
    - (3) wrong-species: E. coli-only index hits for human/mouse spectra, minus any I/L sequence in the human or
      mouse target space.
  - Decoys are never candidates. This is asserted against labels.arrow when the dataset is built.
- **Dedup:** raw sha256 against all manifests, plus an exact spectrum hash within training and against every
  evaluation store.
- **Loss:** BCE over candidates (absolute, comparable across spectra, so that one global TDC threshold is meaningful)
  plus 0.5 × InfoNCE on the pooled split embeddings (E5).

### E6: frozen cross-encoder reranker
- **Arm.** `ce_s` is the frozen size-s cross-encoder logit, used as the score in exactly the production estimator: the
  top-5 PSM table, paired double competition, peptidoform units and picked protein. It is the same code path as
  static_gbm_pred in scripts/m3/evaluate.py and scripts/m3e/{d8,armM}.py.
- **Comparators:** frozen static_gbm_pred (A1.1), RESET-LDA-pred (§6.4, M3 recipe, pooled peptide level only) and
  hyperscore.
- **Sets:**
  - D0 and D1 (arms P with seeds 2001–2003, and C);
  - the D8 holdout, first 10 runs (per-run experiments, pooled verdict; arms P and C);
  - D7 arm M (yeast + human entrapment, r = 1, seeds 2001–2003).

  These are evaluated **once** with the single model frozen by training-fold validation. Anything after that is
  labelled exploratory.
- **Primary metric.** Peptidoform IDs at matched arm-P FDP 1%: the largest own-q threshold on a 200-point grid where
  the mean arm-P FDP̂ over the 3 seeds is ≤ 1%, counted on the production search. We also report IDs at own q
  0.1/1/5%, the matched arm-C FDP, and the §6.3 verdicts under CI convention (a).
- **GO as default candidate** requires all of the following:
  1. **Sensitivity.** CE / static_gbm_pred at matched arm-P FDP 1% is ≥ 1.00 on D0, on D1 and on pooled D8, **and**
     ≥ 1.02 on pooled D8 with a per-run bootstrap lower 95% bound > 1.00. CE / RESET-LDA-pred at matched FDP 1% on
     pooled D8 must be ≥ 1.00.
  2. **Entrapment.** No arm (P, C, M) is "violated" at α = 1% or 5% on any set. Arm M is a memorization screen: the
     CE's mean FDP̂_M at 1% must be ≤ static_gbm_pred's on the same searches + 0.25 percentage points. M3e had static
     at 1.05%, re-computed here on identical searches. We also report the enrichment of training-positive human
     peptides among accepted arm-M entrapment hits.
  3. **Gates.**
     - G0: bit-identical CE scores under label swap and DB permutation (D1 vs the E2 swap index).
     - G2-original on the M3e clean nulls: CE − hyperscore target-win rate, one-sided 95% UCB ≤ 0.01 overall and
       ≤ 0.02 in the top decile.
     - G5: null winner AUC ≤ 0.55.
     - Tail-G2, G1 and G3 need A2's fresh nulls (DN2, in progress). They are marked **pending**, so even a GO is
       conditional and no FDR claim is made.
  4. **Runtime.** Measured CE reranking ≤ 1.0 GPU-s per median run (the E2b budget).

  If (1) fails while (2)–(4) pass, the CE is a **NO-GO as default** and stays a research arm. If (2) or (3) fails, it
  is a **NO-GO on validity**.

### E5: bi-encoder retrieval for open/HLA-scale indexes
- **Reference.** D0 open mode as in M4: the human_contam_prod index, delta [−150, +500] Da, K1 top-K 500, report 5.
  We take the rank-1 **target** PSMs accepted at PSM-level TDC q ≤ 1% (hyperscore), split into narrow (|Δm| < 0.5 Da)
  and shifted.
- **Bi-encoder retrieval.** Cosine top-K over all peptidoforms (targets and decoys) whose mass lies in the same open
  window, by exact brute force on the GPU. This gives the best case for any ANN.
- **Metric.** Recall@K for K ∈ {10, 100, 500, 1000}.
- **GO for open/HLA-scale retrieval** requires both:
  1. recall@K ≥ 95% at some K ≤ 1000, on all accepted hits and separately on shifted hits (lit/01 R4);
  2. G3-style symmetry on the existing dev nulls, i.e. the M3e G2 mixture null spectra. For every decile of the
     exhaustive best score, |recall_T − recall_D| ≤ 1%, where recall_T and recall_D are the bi-encoder top-K
     (K = 10, 100, within the closed precursor window) recall of the exhaustive best target and best decoy
     (hyperscore).

  Otherwise it is a **NO-GO**. The retrieval cost is reported against K1 open-mode cost.

---
# Results (written after evaluation; sections 0–2 above are unchanged)

## 3. Training (task 2)
- **Runs.** 43 runs were selected: 13 from M3 and 30 new ones (seed 23). Of these:
  - 40 ingested as OT-HCD;
  - MN-WT1-Frac8 is IT-CID and was excluded;
  - 2 new RAW files are unreadable by RawFileReader (uzorak_07, Sh-271) and were excluded as recorded;
  - 7 runs gave 0 positives at q ≤ 0.1%. Examples are a chymotrypsin run (qExHF01_…_Chymo) and runs with few tryptic
    IDs. They contribute nothing.
- **Dataset.** 492,411 positive spectra and 3,062,800 candidates:
  - 1,375,046 hard top-K negatives;
  - 984,820 shifted-precursor negatives;
  - 210,523 wrong-species negatives.

  Spectra per fold are 151,665 / 95,457 / 114,259 / 39,328 / 91,702 (fold 4 is validation).
- **Dedup.** Raw sha256 overlaps none of the 1,228 hashes in the other manifests. The spectrum hash removed 0
  duplicates within training and 0 spectra that also occur in the 31 evaluation stores (1.48M spectra).
- **Decoys.** No decoy-recipe sequence entered training; this is asserted. Files: `dataset_meta.json`,
  `label_summary.json`, `train_spectrum_hashes.parquet`, `data/manifests/train_m7{.tsv,_folds.tsv,.sha256}`.
- **Model.** `models/ce_s.pt`, sha256 aa19c8f3…, best epoch 5 of 8 (early stop).
  - Fold-4 validation: BCE 0.021, top-1-among-candidates 99.66%, candidate AUC 0.9993.
  - Training took 504 s on the A100.
  - These validation numbers concern confident spectra only and say nothing about FDR behaviour.

## 4. E6: frozen cross-encoder vs static_gbm_pred and RESET (`e6.json`; dev curves, not claims)
Units are peptidoforms. "matched P1%" is the number of IDs where the mean arm-P FDP̂ reaches 1%. Verdicts are under CI
convention (a), and none is "violated" anywhere.

| set | scorer | IDs q 0.1% | IDs q 1% | IDs q 5% | matched P1% | matched P5% | P 1% (FDP̂, verdict) | C 1% |
|---|---|---|---|---|---|---|---|---|
| D0 | hyperscore | 12,635 | 19,600 | 23,640 | 19,441 | 23,561 | 1.09% inconcl. | 1.00% inconcl. |
| D0 | static_gbm_pred | 23,566 | 24,512 | 25,700 | 24,429 | 25,618 | 1.12% inconcl. | 0.78% supp. |
| D0 | **ce_s** | 20,762 | 24,285 | 25,912 | 24,234 | 25,844 | 1.13% inconcl. | 0.77% supp. |
| D0 | RESET-LDA-pred | 23,235 | 24,535 | 25,747 | 24,452 | 25,707 | 1.22% inconcl. | 0.74% supp. |
| D1 | static_gbm_pred | 5,917 | 6,052 | 6,297 | 6,046 | 6,294 | 1.07% inconcl. | 0.95% supp. |
| D1 | **ce_s** | 5,982 | 6,078 | 6,317 | 6,082 | 6,317 | 0.90% **supp.** | 0.75% supp. |
| D1 | RESET-LDA-pred | 5,868 | 6,062 | 6,297 | 6,060 | 6,288 | 1.14% inconcl. | 1.25% inconcl. |
| D8 (10 runs) | hyperscore | 41,517 | 68,506 | 91,533 | 68,233 | 90,749 | 1.02% inconcl. | 1.11% inconcl. |
| D8 | static_gbm_pred | 97,091 | 110,105 | 117,665 | 109,977 | 117,400 | 1.03% inconcl. | 0.87% supp. |
| D8 | **ce_s** | 73,495 | 107,930 | 119,632 | 108,067 | 119,679 | 0.98% **supp.** | 0.94% supp. |
| D8 | RESET-LDA-pred | 89,660 | 111,600 | 118,122 | 111,273 | 117,758 | 1.05% inconcl. | 0.83% supp. |

Protein groups at 1% (picked; not for RESET): D0 is 6,455 for static vs 6,386 for CE, and D1 is 838 vs 840.

**CE / comparator at matched FDP.** On D8 the bracket is a 2,000-draw run bootstrap; D0 and D1 are point estimates.

| | D0 P1% | D1 P1% | D8 P1% | D8 P5% | D8 C1% | D8 C5% |
|---|---|---|---|---|---|---|
| vs static_gbm_pred | 0.992 | 1.006 | **0.983 [0.968, 0.999]** | 1.019 [1.005, 1.030] | 0.980 [0.968, 0.997] | 1.015 [1.000, 1.025] |
| vs RESET-LDA-pred | 0.991 | 1.004 | 0.971 [0.955, 1.000] | 1.016 [1.001, 1.026] | 0.966 [0.949, 0.996] | 1.009 [0.995, 1.019] |
| vs hyperscore | 1.247 | 1.055 | 1.584 [1.386, 1.752] | 1.319 | 1.609 | 1.301 |

**Arm M** (D7 yeast with human entrapment, r = 1, 3 seeds, combined estimator; `e6.json:D7_armM`):

| scorer | FDP̂_M at 1% [95% CI] | verdict | LB | FDP̂_M at 5% | verdict |
|---|---|---|---|---|---|
| hyperscore | 1.046% [0.83, 1.26] | inconclusive | 0.52% | 5.27% | inconclusive |
| static_gbm_pred | 1.049% [0.80, 1.30] | inconclusive | 0.52% | 5.01% | inconclusive |
| **ce_s** | **1.045% [0.87, 1.22]** | inconclusive | 0.52% | 4.77% | **supported** |

The memorization screen (CE ≤ static + 0.25 pp) passes: 1.045% vs 1.049%.

Memorization diagnostic: the share of accepted human entrapment peptides that are training positives.

| scorer | share among accepted entrapment hits | share among all competing entrapment winners | enrichment |
|---|---|---|---|
| ce_s | 16–24% | 11.3–11.6% | about 1.8× |
| hyperscore | 13–25% | 8.4–8.7% | about 2.2× |

The enrichment is about the same as hyperscore's, so it is not a CE-specific memorization signature. Accepted
entrapment hits are easy-to-match human peptides whatever the scorer.

One CE-specific signal: among *competing* entrapment winners, the CE picks training-positive human peptides more
often (11.4% vs 8.6%). That is a weak peptide prior. It does not raise the accepted N_E (CE 92 vs static 95, mean
over 3 seeds).

**Gates** (`gates_ce.json`):
- **G0 PASS.** The D1 PSM tables are bit-identical between pfu_contam_E2001 and the swap/permuted index (33,129
  rows), and the CE scores are bit-identical. Repeat scoring on the shared GPU is also bit-identical.
- **G2-original PASS** on 41,918 clean null competitions (45 M3e mixtures).
  - CE − hyperscore null target-win rate: −0.0015 overall (UCB 0.0016 ≤ 0.01) and +0.0021 in the top decile (UCB
    0.0124 ≤ 0.02).
  - Unlike static_gbm_pred, which changes the null winner on 0.1% of spectra, the CE changes it on **15.2%**, so for
    the CE this G2 is a real test, not a vacuous one.
- **G5 PASS.** AUC 0.507 [0.501, 0.513], no alarm; static_gbm_pred gives 0.515.
- **Pending: tail-G2, G1 and G3.** They need A2's fresh nulls (DN2). As a diagnostic only, on the D12-contaminated
  nulls where static already fails:
  - CE tail-G2: n_tail 709, target fraction 0.772 (UCB 0.797), Δ vs hyperscore +0.044 (UCB 0.091);
  - static tail-G2: 0.766, Δ +0.038;
  - CE G1: ρ 0.072 [0.032, 0.111], against 0.055 for static. On all nulls it is 0.099 vs 0.061, so the CE's score
    correlates somewhat more with the sequence-only target classifier. This is the weak peptide prior seen in arm M.

  These would fail as A1 gates on these nulls, just as static does. They are not verdicts.

**Runtime** (`bench_ce.json`, under lock, persistent process, D0 run 2 + 10 D8 runs, production top-5 tables):
- median 0.56 GPU-s, 0.68 s wall and 0.68 core-s per run;
- 1.74e5 candidates per GPU-s and 7.7e4 spectra per GPU-s;
- model load 0.35 s, once.

That is within the 1.0 GPU-s budget and comparable to the current fast static scoring (0.60 s median,
scoring_speed). The first implementation selected peaks on the CPU, which cost 4.9 core-s per run (lexsort over all
centroids). It was replaced by a GPU version that is bit-identical; this was verified, and all reported scores are
unchanged (`/tmp/vulcan_m7/bench_ce_cpu_peaks.json` holds the old timing).

**E6 verdict under the pre-stated rules: NO-GO as default. The CE stays a research arm.**
- **Criterion 1 fails.** CE/static at matched P1% is 0.992 on D0 and 0.983 [0.968, 0.999] on D8 (needed ≥ 1.02 with
  LB > 1). CE/RESET on D8 is 0.971.
- **Criteria 2–4 pass:** no violated arm, arm M passes, G0/G2/G5 pass, runtime passes.

**Reading:**
- A 1.66M-parameter cross-encoder trained from scratch on 40 runs, with **no predicted intensities, no RT and no
  hand-crafted features**, comes within about 2% of the AlphaPeptDeep-feature GBM at 1%. It is ahead at 5% on D8:
  +1.9%, CI above 1, a reported number and not a criterion.
- It is clearly behind at 0.1%: −24% on D8 and −12% on D0. Its scores separate the confident core well but lack the
  extreme-tail resolution that predicted intensities give.
- Its entrapment behaviour is at least as good as static's. It is "supported" on arm P at 1% on D1 and D8 and on
  arm M at 5%, where static is inconclusive.
- The obvious next version is a hybrid (static features + CE logit, or CE with predicted-intensity tokens). That would
  be a new pre-registered arm, with CE logits out-of-fold on the training runs. It was not run here.

## 5. E5: bi-encoder retrieval (`e5.json`)
**Open-search recall on D0 (M4 open mode).** The reference is 18,389 rank-1 targets at PSM q ≤ 1%: 6,560 shifted and
11,829 narrow. Bi-encoder retrieval is exact cosine top-K within the same open window.

| K | 10 | 100 | 500 | 1000 |
|---|---|---|---|---|
| recall, all | 54.5% | 67.5% | 75.9% | **79.5%** |
| recall, narrow | 71.4% | 83.0% | 89.6% | 91.9% |
| recall, shifted | 24.2% | 39.5% | 51.2% | **57.2%** |

Median rank is 1 for narrow hits and 430 for shifted ones.

Cost:
- Brute-force query is 0.28 ms per spectrum on the GPU, against 0.24 ms per spectrum for the whole K1 open search
  (21.0 s for 88,636 spectra). An ANN would be cheaper but recalls less.
- Embedding the 7.06M peptidoforms takes 29 s for two charge banks.

**G3-style symmetry on dev nulls.** M3e mixtures 00/01 contain 26,750 null spectra with exhaustive closed-window
candidates (median window of 10 candidates).
- **K = 10 FAILS.**
  - All nulls: the top decile has recall_T 0.824 vs recall_D 0.746 (+7.8 pp), and deciles 8–9 are +2.1 and +1.3 pp.
  - Clean species subset (1,529 spectra): 8 of 10 bins are off by 1.9–6.6 pp, in both directions.
- **K = 100 passes only trivially**, because K is at least the window size.
- Post-hoc diagnostic, *added after seeing this result and not a criterion*: the target share of open-window
  top-10 hits on the same null spectra is 0.523 (SE 0.0004).
  - It is 0.528 on offset nulls, which are shifted human spectra and expected to find human peptides.
  - On species nulls it is 0.494–0.509; on yeast it is below 0.5.
  - At K = 100 it is 0.504, and at K = 1000 it is 0.499.

**E5 verdict: NO-GO** for bi-encoder retrieval on open/HLA-scale indexes, by both pre-stated rules:
- recall@1000 is 79.5% overall and 57.2% on shifted hits, against the 95% required;
- the G3-style symmetry fails at the only non-trivial K.

The split-stage embeddings are useful as the CE's internal representation. As a *retriever* they (i) miss 20–40% of
the K1 hits within K ≤ 1000, worst on exactly the shifted PSMs that open search exists for, and (ii) show a
target-favouring retrieval asymmetry at small K on null spectra. That asymmetry is the failure mode lit/03 §9.2(5)
warns about. This agrees with lit/01 R4: keep index retrieval followed by reranking, and use embeddings only as side
channels (clustering/dedup) after their own symmetry tests.

## 6. Runtime per stage (A100 + 32 Zen 4 cores, shared node, loadavg 20–70)

| stage | wall | notes |
|---|---|---|
| training-run selection (PRIDE API) | ~10 min | metadata only |
| download of 30 runs (26.7 GB) | ~30 min | proxy, SHA-1 verified where PRIDE has one (23/30) |
| ingest (vulcan.io.convert, NVMe staged) | 210 s total for 27 runs, 6 parallel | |
| search, top-20 (56 run-searches incl. E. coli-only) | 4.1 s GPU-search wall total | |
| features (AlphaPeptDeep, lazy) + labels | ~10 min | |
| dataset build (negatives, hashes over 1.48M eval spectra) | 252 s | |
| training (8 epochs, size s) | 504 s | |
| CE scoring of all 119 dev PSM tables (18.9M rows) | 115 GPU-s | wall 983 s with the old CPU peak path |
| E6 evaluation (D0/D1/armM; D8) | ~25 min; ~4 min | CPU, estimator and bootstrap |
| gates (G0 + G2/G5 over 45 mixtures) | 17 min | |
| E5 (open search 21 s; embeddings 29 s; queries 5 s; G3-style about 5 min) | | |
| **CE reranking in production** (median per run) | **0.56 GPU-s / 0.68 s wall / 0.68 core-s** | under lock |

## 7. Caveats
- Everything here is dev. D0/D1 are a selection set for the static scorer, not for the CE. No CE choice used any
  dev set.
- The matched-FDP CIs on D8 come from bootstrapping runs, and there are only 10 runs.
- The E2b throughputs and the bench were measured on a GPU shared with other agents.
- The tail-G2/G1/G3 diagnostics use the D12-contaminated nulls and are not verdicts.
- The E5 G3-style test uses closed windows, which hold only about 10 candidates, so K ≥ 100 is uninformative.
- Training labels come from the frozen static_gbm_pred. The CE therefore learns "what static calls confident", which
  may cap it at the static scorer's view of the confident core. This is a plausible reason it does not exceed static
  at 1%.
- E2b throughput holds for unfused PyTorch. TensorRT/int8 was not tried.
