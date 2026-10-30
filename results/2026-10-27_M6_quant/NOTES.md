# M6: quantification (LFQ with MS1 features + opt-in MBR with transfer q-values; TMT SPS-MS3), dev only

- **Data:** D9 (new dev set, see `DATA.md` "D9": ProteoBench DDA-Astral HYE A/B × 3, plus 3 human-only DDA-Astral runs as the MBR negative control) and D4 (PXD016766, TMT11 SPS-MS3). **No sealed file was touched; T1 was not opened.** All numbers are development curves, not claims (§6.6). §6 was not changed.
- **Node:** gpu-node, 32 Zen 4 cores, A100 80GB (`hwinfo.json`). It was shared and very noisy: load average 14 at the start, then 100–113 during the quant runs (other agents: E9 scoring, XL xiSEARCH).
- **Files:**
  - `metrics.json`: everything, assembled by `scripts/m6/summarize.py`.
  - `eval_*.json`: comparator evaluations.
  - `timings.jsonl`: timed runs under the bench lock, each with an hwinfo fingerprint.
  - `configs/`: Sage configs.
  - `config.yaml`: every parameter.
  - `git_sha`.
- **Code:**
  - `vulcan/quant/{lfq,tmt,metrics}.py`, with tests in `tests/quant` (7 tests, synthetic).
  - `scripts/m6/*`.
  - `scripts/data/fetch_d9.py` plus `data/manifests/{d9,d4_ref}.tsv/.sha256`; `verify_data.py` now also reads both manifests.

## 1. What was built

**LFQ (`vulcan.quant.lfq`).** Accepted PSMs are linked to MS1 features, then a fixed-window XIC re-integration, then cross-run RT alignment, then opt-in MBR, then roll-up.
- **Accepted PSMs:** `scripts/m6/ids.py`.
  - Paired double competition pooled over the experiment, in peptidoform units.
  - Run-level evidence is a PSM winner of an accepted peptidoform with a score at or above the 1% threshold.
- **Linking PSMs to MS1 features:**
  - First through the precursor-hypothesis table: the spectrum's `ms1_feature` groups with the PSM charge, within ±clip(4 SD, 3, 10) ppm of the run's empirical offset.
  - Otherwise by a targeted feature lookup.
  - On D9 this linked **87% of PSMs** (85% hypothesis + 3% targeted). 87% of precursors get a feature.
- **XIC re-integration:**
  - Integrates M, M+1 and M+2 over a window of the same width in every run (the precursor's median hill width), centred on the run's apex.
  - Feature hill areas stop at the noise floor, so they shrink faster than linearly at low abundance. With hill areas, E. coli log2FC was −2.28 (expected −2) and yeast +1.12 (expected +1). XIC fixes both (−2.04 / +1.01).
  - Hill areas are kept as an arm (§2 table).
- **RT alignment:** star alignment to a consensus scale, 2 rounds, monotone piecewise-linear. Residual SD on D9 is 0.005 min (0.3 s), using 18k anchors per run.
- **MBR (opt-in):**
  - Each candidate transfer is scored with a fixed a-priori score (mass, RT, isotope and trace cosine, intensity consistency). Nothing is fitted.
  - The transfer q-value is PIP-ECHO-style, from three decoy sources:
    - RT-shifted decoy windows;
    - m/z-shifted decoy windows (+11 × 1.0005 Da);
    - decoy-peptide donors: decoy peptidoforms that reach the ID score threshold are transferred exactly like targets.
  - FDR(s) = (max(D_rt, D_mz) + 1)/T + D_pep/T, computed **per acceptor run**.
  - Features already explained by an MS2-identified precursor of the acceptor run are excluded from **every** window, target and decoy alike.
  - Transfers carry `provenance = mbr`, `transfer_q` (plus `transfer_q_peak` without the decoy-peptide term) and `n_donor_runs`. **They are never counted as identifications**: every ID count below is MS2-only.
- **Roll-up:**
  - Top3.
  - MaxLFQ: our own implementation (pairwise medians, least squares per connected run component), unit-tested on synthetic data.
  - directLFQ 0.3.3: pip, installed `--no-deps` into `data/cache/pylibs_m6`, plus `multiprocess`.
- **Normalisation:** label-free per-run median, no species information. Every tool gets the same post-processing (`scripts/m6/evaluate_lfq.py`).

**TMT (`vulcan.quant.tmt`):**
- **SPS-MS3 linkage:** MS3 `master_scan` equals the MS2 scan, verified against the filter-string parent m/z. 100% of accepted PSMs linked.
- **Per-plex tolerance:** a first pass at ±20 ppm, then offset ± clip(3 SD, 3, 15) ppm. D4 offset is −3.5 to −4.6 ppm, SD 0.45–0.52 ppm, final window ±3 ppm.
- **Isolation purity** comes from the MS1 feature (hypothesis-table `isolation_fraction`).
- **Impurity correction:** NNLS against the OpenMS TMT11 default matrix. The lot sheet is not deposited.

## 2. D9 dev LFQ benchmark (6 HYE runs, A/B × 3; ProteoBench convention log2(A/B): H 0, Y +1, E −2)

Identifications (MS2 only, 1%):
- **Ours:** 32,955 peptidoforms. Frozen M3 static_gbm_pred, paired TDC.
  - The model is OT-HCD-class and is applied **out of class** to Astral HCD.
  - Plain hyperscore gives 16,191.
- **Sage 0.14.6:** 30,188 peptides by its own LDA q-values.

**Dev sanity check of the static model on Astral.** Yeast or E. coli-unique peptidoforms accepted in the 3 human-only runs (9-run q):
- 60 of 26,512 at 1%;
- natural-entrapment estimate N_YE(1+1/r)/N = **0.80%** (r = 0.39);
- 0.21% at 0.1% and 3.2% at 5%.

This is not a §6.3 verdict, but it does not suggest that the doubling over hyperscore is FDR inflation.

Quant metrics. n is analytes with ≥1 value in each condition. ε is the median |ε| at that filter. CV is the median per-condition CV over analytes with 3/3 values. DE uses proteins with ≥2 values per condition: the pAUC is standardised over FPR ≤ 5%, and the TPR is at empirical FDP ≤ 5% (directLFQ).

| tool (level) | n (≥1/cond) | n (3/3) | completeness | median abs ε | H / Y / E | Y / E log2FC | median CV | DE pAUC Top3 / MaxLFQ / dLFQ | TPR@5%FDP | proteins |
|---|---|---|---|---|---|---|---|---|---|---|
| **Vulcan, MBR off** (precursor) | 18,933 | 7,280 | 0.59 | **0.130** | 0.128 / 0.133 / 0.162 | +1.01 / −2.04 | **5.0%** | 0.38 / 0.57 / 0.62 | 0.45 | 5,135 |
| **Vulcan, MBR on** (transfer q ≤ 1%) (precursor) | 23,695 | 14,876 | 0.76 | 0.134 | 0.128 / 0.146 / 0.221 | +1.01 / −2.06 | 5.7% | 0.37 / 0.53 / **0.63** | **0.49** | 5,135 |
| Vulcan, MBR off (peptidoform) | 18,007 | 7,078 | 0.59 | 0.136 | 0.134 / 0.141 / 0.177 | +1.02 / −2.05 | 5.1% | | | |
| Vulcan, MBR on (peptidoform) | 22,473 | 14,288 | 0.76 | 0.136 | 0.130 / 0.148 / 0.222 | +1.01 / −2.08 | 5.7% | | | |
| Vulcan arm: feature hill area, MBR off | 18,933 | 7,280 | 0.59 | 0.166 | 0.157 / 0.190 / 0.375 | +1.12 / −2.28 | 5.9% | 0.41 / 0.56 / 0.65 | 0.51 | 5,135 |
| **Sage 0.14.6 LFQ** (peptide, charges combined) | 20,027 | 19,993 | 1.00 | 0.134 | 0.125 / 0.150 / 0.252 | +1.00 / −2.04 | 6.1% | 0.36 / 0.42 / 0.50 | 0.34 | 4,380 |
| **FlashLFQ 2.1.5 on our PSMs**, MBR off (precursor) | 20,130 | 7,387 | 0.56 | 0.164 | 0.163 / 0.170 / 0.169 | +0.92 / −1.92 | 8.0% | 0.45 / 0.53 / 0.60 | 0.45 | 6,163 |
| FlashLFQ on our PSMs, MBR on (PIP-ECHO q ≤ 1%) (precursor) | 26,839 | 14,343 | 0.75 | 0.165 | 0.159 / 0.184 / 0.227 | +0.88 / −1.86 | 8.7% | 0.41 / 0.48 / 0.58 | 0.43 | 6,163 |
| FlashLFQ on Sage PSMs, MBR off (precursor) | 20,431 | 7,631 | 0.60 | 0.168 | 0.167 / 0.175 / 0.169 | +0.92 / −1.91 | 8.2% | 0.44 / 0.53 / 0.61 | 0.42 | 5,942 |
| FlashLFQ on Sage PSMs, MBR on (PIP-ECHO q ≤ 1%) (precursor) | 25,179 | 13,833 | 0.76 | 0.168 | 0.163 / 0.185 / 0.206 | +0.88 / −1.87 | 8.7% | 0.44 / 0.46 / 0.58 | 0.45 | 5,942 |

**On matched analytes** (`eval_common.json`; peptidoforms quantified by every tool, ≥1 value per condition; I/L and CAM normalised): 13,506 peptidoforms.

| tool | median abs ε | H / Y / E | median CV | complete (6/6) |
|---|---|---|---|---|
| Vulcan, MBR off | 0.132 | 0.128 / 0.140 / 0.177 | **4.8%** | 5,947 |
| Vulcan, MBR on | 0.120 | 0.117 / 0.127 / 0.179 | 5.0% | 10,322 |
| Sage LFQ | **0.114** | 0.113 / 0.119 / 0.122 | 5.3% | 13,501 |
| FlashLFQ (our PSMs), MBR off / on | 0.154 / 0.145 | 0.153 / 0.156 / 0.161, then 0.145 / 0.148 / 0.136 | 7.6% / 8.0% | 6,544 / 10,476 |
| FlashLFQ (Sage PSMs), MBR off / on | 0.153 / 0.145 | | 7.8% / 8.0% | 6,835 / 10,569 |

How to read it:
- **Accuracy and precision.** With XIC re-integration, our full-depth ε equals Sage's (0.130 vs 0.134) and our CV is lower (5.0% vs 6.1%).
  - **On the analytes both tools quantify, Sage is more accurate** (ε 0.114 vs 0.132 MBR off, 0.120 MBR on; E. coli 0.122 vs 0.177). Its means use 3/3 replicates for almost every analyte, while ours with MBR off often use 1–2. Our CV stays lower (4.8% vs 5.3%).
  - The remaining E. coli gap (the low-abundance side: 5% in A) is the next target.
  - On the same PSMs, FlashLFQ gives ε 0.164 and CV 8.0%: it compresses ratios (Y +0.92, E −1.92).
- **Depth and completeness.** Sage LFQ integrates every 1%-peptide in every file, i.e. it always borrows identity, so it is 100% complete with no transfer FDR. With MBR off we are 59% complete (Astral DDA undersamples at 15 min). MBR raises this to 76% at an ε cost of +0.004.
- **DE pAUC.** Our DE pAUC is the best in the table (directLFQ 0.62–0.63 vs Sage 0.50, FlashLFQ 0.58–0.60), as is our TPR at 5% FDP (0.45–0.49 vs 0.34).
- **Roll-up ranking.** directLFQ > MaxLFQ >> Top3 for DE on every tool. MaxLFQ gives the lowest protein CV (4.1–4.4%).
- **Caveat on DE.** Protein sets differ between tools (protein-group strings), so DE compares tool pipelines, not identical proteins.

## 3. MBR transfer-FDR behaviour (truth: the 3 human-only runs as acceptors; 9-run experiment)

There are 22,590 yeast or E. coli precursors that are candidates for transfer into the H runs, and every one of them is truly absent. Columns:
- *false/accepted*: the lower bound on FDP in the H runs.
- *combined*: N_YE(1 + N_H,cand/N_YE,cand)/N_acc, which assumes equal false-acceptance per candidate.
- *Y/E transferred*: the share of yeast or E. coli precursors that were MS2-identified in A/B runs and transferred into ≥1 H run (Lim 2019 metric; MaxQuant in Lim: 44% at protein level).

| version (development order) | transfer q ≤ 1%: accepted in H / false | FDP LB / combined | Y/E transferred | q ≤ 5%: LB / combined |
|---|---|---|---|---|
| v1: pooled q, donors 1%, claimed features allowed (`v1_allow_claimed`) | 0 (pooled q never < 1%) | – | – | 4.1% / 9.0% |
| v1, donors 0.1% | 11,800 / 189 | 1.6% / 3.1% | 0.94% | 3.8% / 7.5% |
| v2: + claimed-feature exclusion (pooled q, donors 0.1%) | 12,187 / 130 | 1.1% / 2.1% | 0.87% | 2.7% / 5.2% |
| **v3 (default): + per-acceptor-run q, donors 0.1%** | 10,806 / 68 | **0.63% / 1.23%** | **0.51%** | 2.65% / 5.19% |
| v3, donors 1%, q_total | 729 / 7 | 0.96% / 2.1% | 0.06% | 2.9% / 6.3% |
| v3, donors 1%, q_peak only (no decoy-peptide term) | 15,895 / 179 | 1.1% / 2.5% | 1.1% | 2.9% / 6.3% |
| FlashLFQ 2.1.5 PIP-ECHO (our 9-run PSMs; its own MBR candidates) | 15,714 / 298 | **1.9%** / n.a. | n.a. | 7.7% / n.a. (no PIP filter: 25.7%) |

Findings:
- **Donor ID error is the floor.** With donors at ID q ≤ 1%, the decoy-peptide term alone is D_pep/T ≈ 1.0%. Decoy-peptide donors transfer as often as target donors (51% vs 52% beat both decoy peaks), so the propagated ID error is about the donor FDR. On top of the peak term, the total q-value barely reaches 1%: 729 transfers into H runs, and 4.5k on HYE. The default is therefore donors at ID q ≤ 0.1% (D_pep/T = 0.09%), giving 31.8k HYE transfers at q ≤ 1%. q_peak with 1% donors gives 42.9k, but it ignores the donor error and is anti-conservative (LB 1.1%).
- **The claimed-feature rule matters.** 58% of the false yeast/E. coli transfers (v1) landed on a feature that the acceptor run had already MS2-identified as a different (human) peptide, versus 2.5% of human transfers. Excluding such features in both target and decoy windows removed about a third of the false transfers. This rule was added *after* looking at the H-run diagnostic. It is label-free and a priori sensible, but the numbers above are dev, not validation.
- **Run-specific q-values are needed.** Pooled q underestimates FDP in acceptor runs where the null fraction is high: the H runs, where 46% of candidates are null.
- **Default calibration (v3).** At nominal 1% the false-transfer FDP in the H runs is between 0.63% (lower bound) and 1.23% (combined estimate). At 5% it is 2.65–5.2%. Roughly calibrated.
- **Decoy exchangeability in the pure-null subset** (Y/E into H, above the lowest accepted score): target wins 196, RT-decoy wins 205, m/z-decoy wins 82. **The RT-shifted decoy is exchangeable with null transfers. The m/z-shifted decoy under-counts about 2.4×**, so it never binds through max(D_rt, D_mz). Fewer features sit at +11 Da off the averagine grid.
- **FlashLFQ on the same truth.** FlashLFQ's PIP-ECHO at q ≤ 1% accepts more H-run transfers (15.7k vs 10.8k), with a false-transfer lower bound of 1.9% vs our 0.63%. At 5% the bound is 7.7% vs 2.65%. With no PIP filter, 25.7% of its H-run MBR peaks are yeast/E. coli. The candidate sets differ, so only lower bounds are comparable.
- **Accuracy cost on HYE.** MBR (v3, 1%) raises completeness from 0.59 to 0.76 for ε +0.004 and CV +0.7 points; E. coli ε goes from 0.16 to 0.22. Before XIC re-integration, the same MBR cost ε +0.015 (0.166 → 0.181, `v3_area_perrunq`). Low-abundance transferred values are the noisy part.

## 4. TMT on D4 (PXD016766: one TMT11 plex of 11 human cell lines, 12 fractions, IT-CID MS2 + SPS-MS3)

**Search and scoring:**
- The search uses CAM plus TMT-K statics in the index. The TMT N-term static is added at load time (mass32 + TMT, mask bit 63, kernel N-term delta; equivalent to Sage's `^` static) because `vulcan_core` supports residue statics only.
- Fragment tolerance is ±600 ppm (the kernels take ppm only).
- **The scorer is plain hyperscore paired TDC: there is no IT-CID model (M3).** 50,230 PSMs, 44,138 peptidoforms at 1%.
  - Sage (LDA): 65,654 PSMs.
  - Submitter (SEQUEST + LDA, `passThreshold`): 69,653 PSMs.

**Agreement:**
- **Our PSMs vs the submitter:** 48,037 share a scan with the submitter, and **99.4% of those (47,730) have the same sequence** (I = L). 2,193 PSMs are ours only.
- **Our PSMs vs Sage:** 49,281 on the same scan, 99.98% the same sequence.
- The gap to both is depth from rescoring, as on D0 before M3.

**Reporters:**
- **MS3 linkage** works for 100% of accepted PSMs.
- **Reporter completeness:** all 11 channels > 0 in **92.2%** of accepted PSMs, ≥1 channel in 99.8%; per-channel non-zero rate 96.9–98.7% (126 and 131C lowest).
- **Agreement with Sage TMT11 level-3 reporters on the same scans** (541,798 channel values): Pearson r(log2) = 0.9999996, and 99.997% within 1%. Channel-share r = 0.995.
  - 1,835 values are zero in ours but not in Sage. Our window is ±3 ppm around the per-plex offset, while Sage takes the nearest peak in a wider window.
  - The submitter files carry no reporter intensities, so quant agreement with the submitter cannot be checked.
- **Isolation purity** (MS1 feature) is available for 88% of PSMs: median 1.0, 87% ≥ 0.5, 79% ≥ 0.7. The MS2 window is 0.5 Th.
- **Ratio compression vs purity.** The median log2(max/min) channel spread rises with purity: 2.42 (purity < 0.5), 2.39, 2.54, 2.65 (purity ≥ 0.9). This is the expected co-isolation compression.

**CV.** The design has no replicate channels (one channel per cell line), so a biological-replicate CV is impossible. Instead we report the technical reproducibility of within-plex channel profiles (channel shares after impurity correction):
- median CV across PSMs of the same peptide: **10.9%** (3,190 peptides; 10.7% at purity ≥ 0.5);
- median CV across fractions for peptides seen in ≥2 fractions: **12.8%** (1,509 peptides; 12.3% at purity ≥ 0.5).

## 5. Runtime (A100 + 32 Zen 4 cores; noisy shared node, load 14 → 110)

Per Orbitrap Astral run (123.5k MS2, 3.7k MS1).
- **Locked runs** (`timings.jsonl`, `flock /tmp/vulcan_bench.lock`, hwinfo per record): Sage (D9 and D4), FlashLFQ (D9) and one-shot ingest.
- **The other stages are not locked.** The numbers below come from the data-generation runs on the noisy node. The lock was held by the XL agent's xiSEARCH runs for more than 2 h, and `scripts/m6/bench.sh` waited 58 min for its first stage; I stopped it after ingest. Re-run it on a quiet node for locked per-stage numbers.

| stage | per run | notes |
|---|---|---|
| ingest RAW→store, one-shot, 4 threads (**locked**) | 15.1 s / 13.4 core-s (incl. .NET start-up); in-process decode 8.9 core-s (= 7.2 core-s/1e5 MS2) | E1/M2b: 6.5 core-s/1e5 MS2 in a persistent worker |
| **MS1 features + hypotheses** (quant prerequisite) | 4.8 s / **4.5 core-s** | 3 threads; Astral MS1 is sparse (QE D0: ~29 core-s) |
| GPU search (K0–K2, top-5) | 1.8–2.1 s / 1.9 core-s, 0.5 GPU-s | HYE index 9.83M peptidoforms (T + decoys); setup 5 s once |
| L3 features + static_gbm_pred | 10.5–12.9 s (prediction misses, cold cache) + 1.5 s GBM | GPU shared with E9 scoring |
| pooled paired TDC (6 runs) | 11.5 s total (1.9 s/run) | |
| **quant: link + align + MBR + XIC re-quant** (6 runs) | 2.6 + 0.4 + 1.9 + 3.0 = **7.9 s total (1.3 s/run)**, single-threaded | |
| roll-up (6 runs) | Top3 0.14 s, MaxLFQ 4.1 s, directLFQ 38 s (num_cores=1) | |
| TMT reporter quant (D4, 12 fractions) | 5.9 s total (0.5 s/run): extract 3.7, link 0.8, impurity NNLS 0.7, purity 0.2 | D4 search: 0.10–0.27 s/run GPU wall |
| *Sage 0.14.6 search + LFQ, 6 HYE mzML (locked, 32 thr)* | *105.8 s / 1,595.5 core-s (17.6 s / 266 core-s per run)* | *mzML conversion not included* |
| *Sage 0.14.6 TMT11 MS3, 12 D4 mzML (locked, 32 thr)* | *48.6 s / 792 core-s* | |
| *FlashLFQ 2.1.5, 6 HYE mzML, our PSMs (locked, 32 thr)* | *433.6 s / 915.5 core-s (72 s / 153 core-s per run)* | *4-thread unlocked runs: ~35 min each* |

**Quant adds little to the cohort cost.**
- On top of RAW→PSM it adds about **5.8 core-s per Astral run**: MS1 4.5 plus link/align/MBR/XIC 1.3, plus 0.7 s for MaxLFQ.
- The MS1 stage was already an optional M2b worker stage. M2b measured it at 53% of decode-worker time on the dev mix (QE-heavy); on Astral it is cheap.
- directLFQ in single-core mode is the one expensive piece (6 s/run). Use MaxLFQ or a multi-core directLFQ at cohort scale.
- Everything in quant except MS1 runs once per experiment on small tables, so it scales with runs × precursors and does not touch spectra. XIC re-quant is the exception: it reads MS1 peaks once per run (0.5 s/run).


## 6. Issues / open
- **The static scorer on Astral is out of class.** OT-HCD model, Astral HCD data; the H-run check gives 0.8% at 1%, but no §6.3 entrapment run was done on D9. An Astral class (or a paired-entrapment run on D9) is needed before any Astral claim.
- **FlashLFQ cannot read .raw on this node**, so the comparators use mzML.
  - The 6 HYE files use ProteoBench's mzML (the same input as Sage).
  - The H runs were converted with TRFP 2.0.
  - So FlashLFQ's peak picking input differs from ours (vendor centroids via RawFileReader).
- **The T6 boundary.** D9's MBR control uses ProteoBench `DDA-astral-entrapment` human runs. §6.6 T6 says "the Astral entrapment module", which we read as the DIA one (`select_t4.py` maps T6 to PXD046444). If the T6 owner disagrees, swap the control, e.g. to PXD014415 (Lim 2019). Recorded in DATA.md.
- **Lock contention.** The bench lock was held by XL xiSEARCH for a long time. Stage timings marked "unlocked" come from the data runs under load 15–110.
- **D4 caveats.** TMT ID depth trails Sage by 24%, because there is no IT-CID rescoring model (3c). The impurity matrix is a generic lot, not the one used.
- **MBR scoring is fixed a priori.** A learned transfer score would need the same exchangeability checks (G2-style) as the PSM scorers. The m/z-decoy under-count should be understood before a wider rollout (non-averagine mass offset? feature density?).
- **Roll-up cost.** directLFQ runs single-core in our call (`num_cores=1`, 38 s for 6 runs under load). It is the most expensive quant step. MaxLFQ takes 4 s and Top3 0.1 s.
