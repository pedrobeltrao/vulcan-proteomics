> **Notice.** This work and this report were produced by an AI agent team (Claude Opus 5.5 / Sonnet 5.5). They have not been reviewed or verified by human experts.

# Vulcan: fast reprocessing of repository-scale proteomics data

*Research briefing, 2026-10-30. Sources: REPORT.md, RESEARCH_TREE.md, PROMPT.md, reports/figures/CAPTIONS.md and SOURCES.md.*

> **Read this first: there is no held-out test-set result.**
> Every number in this briefing comes from **development (dev) data**, and several of those sets were also used for tuning. The single pre-registered look at the sealed test sets (T5 + T1) started on 2026-10-28 at 14:23 CET. T5 finished. T1 was killed in its final scoring stage when the compute node was wiped (around 00:30 on 10-29). The look was closed as "failed (unterminated)", with all outputs deleted unread, so nobody has seen any T1 or T5 number. The protocol (amendment A5) allows at most two retries from the same frozen commit (4fe2e2a); none has been used, so up to two retries remain. A retry needs about 13 h of uninterrupted node time, including a check of the rebuilt mokapot environment. T2–T4 were never opened.
> Error-rate statements rest on **entrapment estimates on dev data**. The frozen scorer **FAILED** the first null-spectrum tests. A post-hoc contamination explanation exists but is unconfirmed.

---

## 1. The problem

Public repositories such as PRIDE hold many thousands of mass spectrometry (MS) proteomics runs. Reprocessing them with one pipeline could turn them into a shared reference, but this is slow and costly with current engines. A faster engine must find as many peptides as established tools (Sage, MSFragger, DIA-NN), localize modification sites as accurately, and quantify as precisely. It must also keep false discoveries under control when thousands of runs are pooled, because per-run error rates add up across studies.

The project goal (PROMPT.md) was to build such an engine by co-designing algorithms for CPUs and GPUs, and to report runtime alongside accuracy.

## 2. The idea

Vulcan places each task on the device that finishes the pipeline first and avoids repeating work across runs.

- **Resident GPU index.** The proteome index is built once on the CPU and kept on the GPU, where candidates are found and rescored for all search modes.
- **Persistent cohort process.** One long-lived process keeps the index and a predicted-spectrum library loaded; CPU workers decode RAW files while GPU threads search and score.
- **Predicted-feature scoring.** The frozen default scorer (v1) is a gradient-boosted model over label-blind match features and predicted fragment intensities and retention times. A test checks that swapping all target/decoy labels leaves search and scoring outputs bit-identical.
- **Global FDR.** For repository-scale lists, Vulcan keeps the best evidence per peptidoform across all runs, then applies one global target-decoy competition and a picked protein step, instead of merging per-run lists.

## 3. Key results (dev data, 32-core EPYC + A100, shared node)

Timings come from a shared node (load 14–110; runs/h varied about ±25%). Core-s means user + sys **CPU** seconds; it **does not include GPU time**, which Vulcan uses and Sage does not.

| Workflow | Accuracy vs comparator (dev) | Speed vs comparator (dev) |
|---|---|---|
| Standard DDA search (frozen v1) | 0.990× [0.989, 0.992] peptidoforms vs Sage + predicted features + mokapot at matched 1% FDP (D0); +2.8% vs plain Sage | search stage about 23–40× fewer CPU core-s than Sage at 4–32 threads (measured single-process 6.0 core-s) (CPU only; excludes GPU time); RAW → proteins about 5,700 runs/h at about 9 core-s/run vs Sage 70–120 core-s + about 250 for RAW→mzML |
| Chimeric search (v2_sep, candidate) | 1.150× [1.069, 1.203] peptidoforms vs Sage chimera + predicted features + mokapot on 10 dev runs held out from scorer selection (D8h; not a sealed test set); protein groups 0.999× | about 3,900 runs/h at about 15 core-s/run |
| Ion-trap CID / TMT | 1.034× (D4) and 1.093× (D10cf) vs Sage + Prosit features + mokapot | about 0.9 core-s per fraction vs Sage about 54 core-s per file (16 threads, unlocked, indicative) |
| Label-free protein abundance | parity with Sage on the HYE benchmark (matched-set error 0.118 vs 0.1145) | about 1,200–1,280 runs/h at about 15 core-s/run vs Sage search+LFQ 266 core-s/run (+ about 250 for mzML conversion) |
| Phospho + site localization | Ala-decoy FLR point estimates near true; stricter upper-bound rule not met on PXD009449 (n = 183); mixed vs pyAscore | about 3,200 runs/h |
| Cross-linking | DSSO: about 73% of xiSEARCH heteromeric pairs, at lower true FDP | DSS: about 460× less wall time than xiSEARCH for the whole pipeline (about 8,900× for the search stage) |
| Repository-scale FDR | global peptide FDP̂ 1.16–1.21% (arms P/C) and protein 1.00% at N = 1,000 (982 runs completed), vs 9.3% (peptide) and 19.1% (protein) for the naive union | streaming from PRIDE is download-bound, about 580 runs/h |

![Figure 2. Compute cost and throughput (dev data, shared node)](../figures/fig2_speed.png)

*Figure 2. (a) CPU core-s per unit of work, log scale. (b) Cohort throughput before and after the last-day exact speed-ups (opt-in, not merged). (c) The 1 TB LFQ estimate below. Dev data only.*

In plain terms, the frozen scorer finds about as many peptides as the strongest comparator we ran, at a small fraction of the CPU cost; its gain over plain Sage comes from the predicted features. On 982 dev runs (140 projects), global FDR kept entrapment-estimated FDP at 1.0–1.2%, while a naive union of per-run 1% lists reached 9.3% (peptide) and 19.1% (protein).

![Figure 3. Identification sensitivity relative to the strongest comparators (dev data)](../figures/fig3_sensitivity.png)

*Figure 3. Vulcan/comparator identification ratios at matched entrapment-estimated FDP of 1% (95% bootstrap CIs; dotted line: non-inferiority margin 0.97). D0/D1 are selection sets. FDR claims are not gate-validated.*

### An example: 1 TB of label-free data on one machine (estimate)

This is **an estimate extrapolated from per-run costs, not a measured end-to-end run**. It assumes about 1 GB per RAW file, so about 1,000 runs.

| | Vulcan (32 cores + 1 A100 GPU) | Sage + mzML conversion (same 32 cores, CPU only) |
|---|---|---|
| CPU cost per run | about 15 core-s (GPU time not counted) | about 266 core-s (search + LFQ) + about 250 core-s (RAW→mzML) |
| Wall time for about 1,000 runs | **about 0.8–1 h** (1,200–1,280 runs/h measured on 96-run cohorts from local disk) | **about 4.5–7 h** (4.5 h assumes ideal 32-core scaling; 2.3–5 h if mzML already exists) |
| Ratio | about 5–7× less wall time per machine; about 17–34× fewer CPU core-s | — |

The comparison favours a GPU node, since Sage runs on CPU only. The runtimes come from different datasets (Sage on Astral D9; Vulcan on a QE cohort with D9 quantification), and streaming from PRIDE (about 580 runs/h) would make I/O the bottleneck. Quantitative accuracy is at parity; error-rate claims are not gate-supported. MSFragger and DIA-NN were not run.

## 4. Statistical status

We pre-registered null-spectrum tests (see glossary). The frozen scorer v1 **FAILED** all three tests under amendment A1 (tail-G2, G1, G3). We later proposed that the null data were contaminated with human material, but this explanation came after the failure and is unconfirmed, so the failure stands.

Later attempts to obtain a decisive null test did not produce one (Figure 6):

- **A2** (fresh null projects): not evaluable; only 6 of at least 10 required projects were admitted.
- **Ion-trap nulls:** futile; 581 usable spectra against 120k needed.
- **A4:** drafted and reviewed, not run, because in simulation the chance of reaching its minimum sample was only 0.08–0.43 (and a clean scorer would pass with probability about 0.06–0.4).
- **A6** (human removed-peptide nulls, v1 only): "not tagged: infeasible", because the measured throughput gave P(v1 minimum) = 0.63 against 0.8 required. Its acceptance test also found that the bootstrap behind earlier G2 "pass" results undercovers with about 37 clusters.

The chimeric candidate v2_sep and the ion-trap scorer have never had a null-spectrum test. All error-rate statements in this briefing are therefore entrapment estimates on dev data. No arm was "violated" for v1 or v2_sep, but the paired estimator overstates by about 20%, so many verdicts are "inconclusive". Amendment A5 knowingly relaxed the claim conditions to allow the single held-out look, which then failed for technical reasons.

![Figure 6. History of the pre-registered statistical gates](../figures/fig6_gate_history.png)

*Figure 6. Each row is a gate, amendment or attempt. Red = failed; amber = not evaluable / not run / infeasible / futile; green = GO in simulation only. No null-spectrum gate was passed.*

## 5. What did not work

- **Learned encoders.** Bi-encoder retrieval was no faster than the inverted index (0.28 vs 0.24 ms per spectrum) with 79.5% open-search recall@1000: NO-GO. The cross-encoder scored 0.983× v1 on D8h, a dev subset held out from selection: NO-GO as default.
- **Null-set attempts.** As above: A1 failed, A2 not evaluable, A4 not run, ion-trap nulls futile, A6 infeasible.
- **Other negatives.** Group-walk FDR over chimeric strata violated its rank ≥3 stratum (3.1%); per-run score calibration cost 13–72% of identifications.

## 6. Limitations and next steps

Beyond the missing held-out comparison:

- Several dev sets were used for tuning (D0/D1 scorer selection, localization settings, quantification choices on D9).
- No DIA result exists.
- Cross-linking sensitivity is below xiSEARCH.
- Site localization is not run-to-run reproducible for about 0.7% of PSMs, because GPU atomic sums break score ties at random.
- All timings come from a shared, loaded node.
- The fastest throughputs are on opt-in or unmerged paths, not the frozen path.

The next session would (1) retry the T1 + T5 look and run the T2 look, (2) run A6 on a quiet node after replacing the undercovering bootstrap, (3) merge the last-day speed-ups, and (4) re-validate a deterministic localization tie rule under a new pre-registration.

---

### Glossary for non-specialists

| Term | Meaning |
|---|---|
| **FDR** (false discovery rate) | The expected fraction of reported identifications that are wrong. "1% FDR" is a target we ask for; whether it is achieved must be checked separately. |
| **Entrapment** | Adding sequences that cannot be in the sample (shuffled, or from another species such as *Arabidopsis*) to the search. How often they are reported gives an estimate of the true error rate. |
| **Null-spectrum test** | Searching spectra that should contain no correct match, to check that the scorer does not invent confident hits. |
| **DDA** | Data-dependent acquisition: the instrument picks individual peptide ions and fragments them one at a time. |
| **LFQ** | Label-free quantification: comparing protein amounts between runs from signal intensity, without chemical labels. |
| **Chimeric spectra** | Fragment spectra that contain two or more peptides fragmented together, because they were isolated at the same time. |
| **Dev vs held-out** | Development data are used to build and tune the method. Held-out (sealed) data are kept unseen until a single final test; none has been read here. |
