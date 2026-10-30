# Draft figure captions

These captions apply to the research briefing, the technical report and the manuscript. Shorten them as each venue needs, but keep the qualifiers. Every number comes from **development (dev) data**. No held-out test-set (T1–T5) result exists. The single sealed look was started on 2026-10-28 and was closed as "failed (unterminated)" with its outputs unread (REPORT.md §3).

Terms used below:
- **Nominal**: the FDR threshold that was requested.
- **Entrapment-estimated FDP (FDP̂)**: the false-discovery proportion estimated from entrapment sequences. Arm P uses paired shuffled entrapment; arm C uses natural foreign (Arabidopsis) entrapment.
- **Verified**: an FDP measured against a ground-truth design. In this project that is the phospho-library truth (Fig. 5a) and the XL library truth (not plotted).

No FDR claim has passed the pre-registered null-spectrum gates (Fig. 6).

Timings come from a shared node: gpu-node, 2× EPYC 9354 (32 affinity cores) + A100 80 GB, load average 14–110.
- Core-s = user + sys CPU seconds. Core-s exclude GPU time: Vulcan also uses the A100, whereas the comparators are CPU-only, so per-run core-s comparisons favour the GPU node.
- Vulcan timings were taken under a bench lock where stated.
- Several comparator timings were taken without the lock and are indicative only.

Each plotted value is listed with its source file and key in `SOURCES.md`, and `make_figures.py` regenerates all figures.

---

**Figure 1 | Vulcan architecture and CPU/GPU placement.**
**(a)** Processing layers L0–L6, coloured by the device that runs them in production:
- blue = GPU;
- grey = CPU.

The proteome index (L1) is built on the CPU and kept resident on the GPU. Candidate generation (L2) and exact rescoring with label-blind and predicted features (L3) run on the GPU. Statistics, cross-run quantification and repository-scale FDR (L4–L6) run on the CPU. GPU MS1 feature finding (L5) is an exact opt-in path. A label-swap/permutation test (G0) checks that L2/L3 outputs do not depend on target/decoy labels.

**(b)** A single long-lived `vulcan.cohort` process keeps the index and the predicted-spectrum library on the GPU. Work flows through three stages:
- 24 CPU decode workers: RAW decode, per-run pilot, MS1 features;
- one or two GPU-stage threads: search and scoring;
- 3 CPU post workers: competition, picked protein groups, quantification.

Per-run outputs feed the L6 global FDR. The schematic follows REPORT.md §1; the worker counts are the E7b/E12 cohort configuration.

---

**Figure 2 | Compute cost and throughput (dev data, shared node).**

**(a)** CPU core-seconds per unit of work for Vulcan (blue; 32 cores + A100) and the comparator (orange) on a log scale. The fold label is comparator core-s ÷ Vulcan core-s; Vulcan's GPU time is not included. What each row measures:

| Row | Unit | Vulcan | Comparator | Note |
|---|---|---|---|---|
| Closed DDA | per run | RAW → protein groups (E12 cohort) | Sage 0.14 at 32 threads, search only | hollow orange circle = Sage + RAW→mzML conversion (250 core-s) |
| LFQ | per run | RAW → protein abundance, MBR off | Sage search + LFQ on Astral D9, mzML given | the two inputs are not like-for-like data; hollow circle = + conversion |
| Phospho | per run | localization + Ala-decoy FLR | Sage phospho, including its database build, without localization | — |
| Open search, semi-specific (LiP), cross-linking | per run | search stage only | search stage only | — |
| Ion-trap CID/TMT | per fraction | search + scoring | Sage at 16 threads, unlocked | — |

Cross-linking (DSS, D3) compares against xiSEARCH 1.8.12 at 32 threads under the bench lock. Its 1.26 × 10⁵-fold ratio sets Vulcan's search stage (0.43 core-s) against the full xiSEARCH run (53,671 core-s). Including Vulcan's site localization the source gives about 50,000× fewer core-s, and for the whole pipeline the wall-time ratio is about 460× (M5_xl; REPORT §2.6).

**(b)** RAW-to-output throughput of 96-run cohorts before and after the E12 exact speed-ups. These speed-ups are opt-in, bitwise-identical in output, and not merged into the frozen path. Runs/h varied by about ±25% between rounds on the shared node, so compare within a row. The LFQ and phospho rows assume the experiment-level step runs serially after the cohort, which makes them conservative.

**(c)** **Estimate, not a measured run:** wall time for 1 TB of RAW files (about 1,000 runs at about 1 GB per run) for label-free quantification on one machine, extrapolated from per-run costs.
- Vulcan: 0.8–1 h, from 1,194–1,280 runs/h measured on 96-run cohorts read from local disk.
- Sage with RAW→mzML conversion: 4.5–7 h. The lower bound assumes ideal 32-core scaling of 516k core-s.

The estimate favours a GPU node. It mixes timings from different datasets (QE cohort and Astral D9). It ignores network or PRIDE I/O, and streaming from PRIDE is download-bound at about 580 runs/h. MSFragger and DIA-NN were not run.

---

**Figure 3 | Identification sensitivity relative to the strongest comparators (dev data).**
Each point is the ratio of Vulcan identifications to comparator identifications at matched entrapment-estimated FDP of 1% (arm P). The comparison works as follows:
- Each method's own q threshold is moved until its arm-P FDP̂ reaches 1%.
- Identifications are summed over runs and averaged over 5 decoy shuffles.
- Bars show the 95% percentile bootstrap CI over runs × decoy shuffles. For IT-CID the bootstrap is over runs only. D1 has a single run.
- The solid line marks parity; the dotted line marks the pre-registered non-inferiority margin (0.97).

| Group | Vulcan scorer | Comparator | Datasets and units |
|---|---|---|---|
| Standard DDA | frozen v1 (`static_gbm_pred`) | Sage + the same AlphaPeptDeep features + seeded mokapot (strongest comparator) | D0, D1 peptidoforms. D0/D1 are selection sets. |
| Chimeric | candidate v2_sep (per-stratum TDC) | fair chimeric comparator: Sage chimera mode + predicted features + seeded mokapot, narrow precursor window | D0 vs R = 5; D1 and D8h vs R = 2 (D8h = 10 dev runs held out from scorer selection, not a sealed test set). Peptidoforms; hollow = protein groups on D8h. |
| Ion-trap CID/TMT | `itcid_linear_pred` | Sage + Prosit features + seeded mokapot | D4, D10cf, sequence units |
| Cross-encoder (grey) | research arm `ce_s` | frozen v1, **not** Sage | D8h (dev holdout). Fails its pre-registered GO rule (NO-GO). |

v1 is non-inferior to the strongest comparator, and its gain over plain Sage comes from the predicted features. v2_sep gains about 15% in peptidoforms, with no protein-group gain on D8h.

Sage chimera's rank ≥ 2 identifications exceed nominal FDP (Fig. 4c), so the chimeric ratios favour the comparator. FDR claims for all scorers are not gate-validated.

---

**Figure 4 | Statistical error control as estimated by entrapment (dev data; no null-spectrum gate passed).**

**(a)** Per-run peptidoform lists at nominal 1%:
- filled circles: arm-P FDP̂ with 95% CI;
- hollow diamonds: arm-C FDP̂.

Shown for frozen v1 (blue), v2_sep (aqua) and Sage + predicted features + mokapot (orange) on D0, D1 and D8h (dev holdout). The paired estimator overstates by about 20%, so arm-P values just above 1% read as "inconclusive" under the pre-registered rule, not as "violated". No list is "violated".

**(b)** Repository-scale FDR on up to 1,000 dev runs (D8; 982 completed, 140 projects), frozen v1 score, arm P, nominal 1%.
- L6 global procedure (blue): best evidence per peptidoform, then global paired TDC, then picked protein.
- Naive union of per-run 1% lists (orange).
- Circles = peptide level; squares = protein level.

At N = 1,000 the global lists stay at 1.0–1.2%, while the naive union reaches 9.3% (peptide) and 19.1% (protein). One entrapment and one decoy seed were used.

**(c)** Arm-P FDP̂ within chimeric rank strata at nominal 1%, with 95% CI. Crosses mark strata with the pre-registered verdict "violated". The methods are:
- Group-walk over strata (chim_v1, D0; rejected under the pre-registered rule because stratum ≥3 was violated);
- v2_sep per-stratum TDC (D0, D8h; no stratum violated);
- the Sage chimera comparator + predicted features + mokapot, by Sage rank (R = 5 on D0, R = 2 on D8h).

Vulcan strata are evaluated at the matched-FDP-1% threshold; Sage ranks at Sage's own q ≤ 1%. Sage rank-≥2 identifications run at 5–9% FDP̂. They are "violated" on D8h and "not evaluable" on D0 because there are too few targets. These rank-level read-outs were pre-registered as descriptive only.

All values are entrapment estimates, not verified FDP. v1 FAILED the only completed null-spectrum test (A1). No null-spectrum test of v2_sep exists.

---

**Figure 5 | PTM localization and label-free quantification (dev data).**

**(a)** True (verified) false-localization rate of accepted mono-phospho STY sites on synthetic-library truth. The rate is measured at nominal FLR thresholds set by the Ala decoy (filled, primary) or the Leu decoy (hollow, secondary), on PXD000138 and PXD009449 (Ymod).
- The horizontal bars show the nominal level.
- The whiskers run up to the one-sided 95% upper confidence bound, which the stricter T2 rule requires to be ≤ nominal.
- The labels give the Ala verdict under that rule.

The Ala point estimates are calibrated or conservative. The upper-bound rule is met only on PXD000138. PXD009449 at 1% is inconclusive, because only 183 single-Y peptides carry that set.

**(b)** Difference in site recall, Vulcan − pyAscore (documented defaults), at matched true FLR, with the Ala decoy. 95% CI from a paired sequence-cluster bootstrap. n.d. = no difference shown.

Some caveats apply to (a) and (b):
- Localization settings were tuned on PXD000138 and PXD009449.
- Localization is not run-to-run reproducible for about 0.7% of PSMs (GPU atomic ties). Verdicts were stable over 11 repeats; values come from one T2 driver dry run.
- The sealed T2 set was not run.

**(c)** Label-free quantification accuracy on the D9 Astral HYE benchmark (6 runs, MBR off). The metric is the median absolute log2-ratio error against the expected ratios on the 13,674 peptidoforms quantified by both Vulcan (blue) and Sage LFQ (orange). It is shown overall and per species. Overall the two are at parity (0.118 vs 0.1145). Vulcan is worse on low-abundance E. coli. With MBR on (not shown), Vulcan's matched-set error is 0.111 vs Sage's 0.1145 (REPORT §2.5). The quantification choices (M6d) were made on D9.

---

**Figure 6 | History of the pre-registered statistical gates (honest summary).**
The x-axis lists decision-log entries in RESEARCH_TREE.md in order, with their dates. Spacing is ordinal, not to time scale. Each row is one gate, amendment or attempt. Markers show:
- process steps (grey squares: adopted, opened);
- outcomes, with status colour plus icon and label: red cross = failed; amber diamond = not evaluable / not run / infeasible / futile; green circle = GO in simulation only, not a gate.

The outcomes, in the decision log's words:
- A1 (tail-G2, G1, G3 on frozen v1): **FAILED**. A post-hoc contamination explanation is unconfirmed.
- A2/DN2 fresh nulls: **not evaluable** (6 of ≥ 10 projects admitted, after a re-walk to fix an implementation bug).
- IT-CID nulls: **futile**.
- A3 v2_sep tiers: expected not evaluable.
- A4/DN3: **not run** (simulated P(v1 minimum) 0.08–0.43; P(clean v1 passes) about 0.06–0.4).
- C9 pilot: GO by simulation.
- A5: knowingly relaxed the claim conditions to allow a single T1 + T5 look.
- A6/HN1: **not tagged: infeasible** (P(v1 minimum) = 0.63 < 0.8 at the measured throughput).
- Sealed look (opened 2026-10-28 14:23 CET; ledger 13:23:46Z): T5 completed. T1 was killed when the compute node was wiped, and the look was recorded as **failed (unterminated)** with outputs deleted unread. A5 allows at most two retries from the frozen commit; one attempt has failed, so up to two retries remain.

No null-spectrum gate was passed, and no held-out number has been seen.
