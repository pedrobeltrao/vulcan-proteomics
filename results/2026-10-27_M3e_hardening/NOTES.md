# M3e: hardening of the frozen default scorer under Amendment A1 (DESIGN §6.0, prereg-v0.4.1). Dev data only

**Scope.**
- Scorer: `static_gbm_pred.pkl` (sha256 in `metrics.json`), frozen. It was not retrained, tuned or varied.
- Protocol: §6, A1, metrics, labels and estimators are unchanged. The one code change to an estimator file is an *optional* `ci_convention` argument for `vulcan.stats.entrapment.evaluate_arm`. Its default is unchanged and it is used to implement A1.6.
- Data: no sealed file was touched, and the T4 selection was not used.
- Files: `metrics.json` (summary), `gates_a1.json`, `g3.json`, `g4.json`, `armM.json`, `ci_convention_sim.json`, `d0d1_ci_conventions.json`, `d8_confirm.json`, `bench_prod.json`, `diag_source_tail.json`, `config.yaml`, `git_sha`.
- Code: `vulcan/gates/` (homolog filter, RT map, gate statistics), `scripts/m3e/`, `tests/stats/test_a1_gates.py` (14 tests, including the tail-G2 MC fire and pass checks).

## Bottom line
- **The gates required before any FDR claim do not all pass.** Tail-G2, G1 and G3 FAIL under the pre-specified rules. Per A1.2/A1.4, FDR claims for the frozen scorer stay blocked, and no T1 look is justified on the gate record.
- **The failures trace to the null material, not to the scorer, but that is a post hoc diagnosis and changes no verdict.**
  - *S. aureus* nulls: one of the six null projects (*S. aureus* PXD023811) contains real human skin and sweat proteins (desmoplakin, desmoglein-1, plakoglobin, dermcidin, ...).
  - Offset nulls: the offset nulls that reach the tail are copies of real human spectra.
  - Both sources carry real target evidence. Hyperscore is also 96% target there; the frozen scorer recognises more of it.
  - Without those two sources, tail-G2 and G1 would pass. That split was chosen after seeing the result, so it is not a verdict.
- **Arm M (D7 yeast, human entrapment, r = 1): inconclusive, not violated.** The frozen scorer accepts no more human entrapment than hyperscore (FDP̂_C 1.05% vs 1.05% at 1%), despite +27% yeast IDs.
- **CI convention (A1.6): (a)**, bootstrap + between-seed variance. Coverage was 0.937–0.947 for (a) and 0.817–0.841 for (b), which is below 93%. The D0/D1 verdicts under (a) equal M3's.
- **D8 holdout (10 runs): the ID gain is confirmed and larger than on D0.**
  - At 1%: 110,105 vs 68,506 peptidoforms; per-run geometric-mean ratio 1.55 [1.41, 1.71].
  - Pooled entrapment: arm P inconclusive (1.03%), arm C supported (0.87%); hyperscore is also P-inconclusive. Nothing is violated.
- **Production-config speed (measured):** scoring median 4.3 s wall / 5.1 core-s per run (range 0.8–10.8 s); end-to-end store→proteins median 4.9 s / 5.4 core-s.

## Gate table (frozen static_gbm_pred; "without RT" = the same frozen model with rt_abs_resid / rt_resid_z set to the run median, an ablation)

| Gate | with RT | without RT | numbers (with RT) |
|---|---|---|---|
| G1 (seq-classifier correlation, top decile of nulls) | **FAIL** | **FAIL** | ρ = 0.048, 99% CI [0.019, 0.076] ⊄ (−0.05, 0.05); n = 8,139. Hyperscore reference ρ 0.028 [−0.001, 0.056] |
| G2 original (PSM-level differential) | pass | pass | diff −0.00002 (UCB 0.00014; budget 0.01), top-decile UCB 0.0005 (budget 0.02); winner changed on 0.08% of nulls, so still nearly vacuous |
| tail-G2 (A1.2, 1% tail) | **FAIL** | **FAIL** | n_tail 751 (≥ 200, evaluable); target fraction 0.715, UCB 0.742 > 0.60; diff vs hyperscore +0.085, UCB 0.124 > 0.10. Both clauses fail |
| G3 (top-K 500 vs exhaustive, nulls by score decile) | **FAIL** | **FAIL** (scorer-independent) | All nulls (n = 26,750): max \|ΔR\| 0.12 pp, pass. Clean species nulls (n = 1,529, ~100 per bin): one bin with recall_T 0.9897 vs recall_D 1.000 (1 missed target of 97; Δ 1.03 pp > 1 pp; Newcombe CI [−5.6, +2.6] pp). The pre-specified rule required both sets |
| G4 (decoy KS across charge strata, training decoys) | pass | pass | τ = 0.538 (95th pct of within-stratum between-run KS; median 0.22); stratum-vs-rest KS 0.08–0.18 |
| G5 (leak AUC on null winners) | pass | pass | 0.506 [0.502, 0.510] ≤ 0.55; top decile 0.535; hyperscore 0.501 |
| Arm M (D7, r = 1) | inconclusive | n/a | see below |

RT ablation: removing the RT channel changes nothing material (tail fraction 0.715 both ways; G1 ρ 0.047; G5 0.505). The failures are not RT-driven, and the RT-mapped nulls do not reveal an RT leak.

## Nulls (A1.3; `scripts/m3e/g2mix.py`, `gates_a1.py`)
- **Sources (OT-HCD, 6 projects; A1.3 requires ≥ 3):**
  - D1 *P. furiosus* (1 run);
  - D7 yeast (6 runs);
  - new dev set **DN** (DATA.md; 4 projects × 3 runs, seeded selection `scripts/m3e/select_nulls.py`): *T. brucei* PXD015648, *S. aureus* PXD023811, *Synechococcus* PCC 7002 PXD021787, *B. cereus* PXD055392. The 12 RAWs and 5 proteomes are sha256-recorded.
- **Mixtures:** 894,579 species-null spectra, each used once, in 45 mixtures of 20,000 with a D0 host.
- **RT:** quantile-mapped onto the host gradient.
- **Offset nulls:** kept as specified, in the first two mixtures. 579 have a candidate.
- **Homolog filter (6-mer or Hamming ≤ 2, I/L collapsed):**
  - Checked against the source proteome plus contaminants, on the matched peptide and its pair partner, for all compared pipelines (label-symmetric). Verified exactly against brute force (tests).
  - It removes most nulls: 586,201 nulls have a candidate and 81,390 are clean (yeast 31.9k, *B. cereus* 24.6k, *S. aureus* 12.6k, *T. brucei* 4.7k, *Synechococcus* 3.7k, Pfu 3.3k, offsets 579).
  - The removals are mostly chance 6-mer sharing with whole proteomes, so the filter is very conservative.
- **Tail count:** the tail reached n = 751 without gathering further nulls.

### Post hoc diagnostics (after seeing the tail result; they change no verdict)
- **Tail by source** (frozen scorer tail / targets; hyperscore tail / targets):

  | Source | frozen scorer | hyperscore |
  |---|---|---|
  | *S. aureus* | 375 / 362 | 131 / 126 |
  | offsets | 17 / 17 | 2 / 2 |
  | yeast | 147 / 67 | 277 / 155 |
  | *B. cereus* | 145 / 60 | 72 / 30 |
  | Pfu | 32 / 11 | 19 / 5 |
  | *T. brucei* | 21 / 8 | 6 / 1 |
  | *Synechococcus* | 14 / 12 | 1 / 1 |

- ***S. aureus* tail target winners** (`diag_source_tail.json`) map to human epidermis and sweat proteins: DSP 54, DSG1 29, JUP 13, TGM3 10, PKP1 9, DCD 8, DSC1 8, CASP14 7, PIP 6, ... This is human skin material in PXD023811. The A1.3 filter (source proteome + keratin-type contaminants) cannot remove it by construction.
- **Offset tail hits** (`scripts/m3e/diag_offset_tail.py`) are shifted copies of D0 spectra whose *unshifted* precursor had no candidate (likely mis-assigned precursors). The shifted copy lands within a few ppm of a plausible human peptide. They are real human fragment ladders, not label-blind nulls.
- **Excluding *S. aureus* and offsets:**
  - tail-G2: n 359, target fraction 0.44 (UCB 0.485), diff −0.072 (UCB −0.018); would pass.
  - G1: ρ 0.004 [−0.027, 0.035]; would pass.
  - G5: 0.4996.
  - G1 by source: offset 0.21, *S. aureus* 0.14; all other sources |ρ| ≤ 0.07 (*T. brucei* −0.07 on few nulls).

  The leakage signal therefore sits entirely in the two sources that contain real human evidence, and the remaining "far" nulls lean towards decoys. This split was chosen after seeing the result. Under §6/A1 it is a diagnosis, not a pass.

## Arm M (A1.5; `armM.json`)
- **Setup:**
  - Sample: D7 TrP runs LM2486/2487/2490 as one experiment.
  - Targets: yeast + contaminants + the spiked human α-synuclein.
  - Entrapment: human Swiss-Prot peptides neighbour-filtered vs targets, subsampled to r = 1.00 on the effective space with seeds 2001–2003 (peptidoform r 1.035; lowest mass decile r 0.83).
- **Estimator:** combined (upper bound; with r = 1 it equals 2·LB).

| scorer | α | N_T (mean) | N_E (mean) | FDP̂_C | 95% CI (a) | LB (CI a) | verdict |
|---|---|---|---|---|---|---|---|
| static_gbm_pred | 1% | 18,017 | 95 | 1.049% | [0.80, 1.30]% | 0.52% [0.40, 0.65] | inconclusive |
| static_gbm_pred | 5% | 18,639 | 479 | 5.008% | [4.67, 5.35]% | 2.50% | inconclusive |
| hyperscore | 1% | 14,197 | 75 | 1.046% | [0.83, 1.26]% | 0.52% | inconclusive |
| hyperscore | 5% | 16,417 | 445 | 5.274% | [4.92, 5.63]% | 2.64% | inconclusive |

- **Not violated** (lower CI of LB ≪ α).
- No memorization signature: the entrapment acceptance rate per target ID is the same as or lower than hyperscore's.
- "Inconclusive" because the r = 1 combined estimator is an upper bound of ~2× the LB.

## CI convention (A1.6; `ci_convention_sim.json`; `scripts/m3e/ci_convention_sim.py`)
- **Model and rule:** fixed in the script docstring before the run.
  - The target side is shared across J = 3 entrapment seeds, and the library `evaluate_arm` is used (n_boot 2000).
  - 10⁴ replicates per scenario; MC SE of a coverage ≈ 0.002.
- **Coverage of E[mean FDP̂_P]:**

  | scenario | (a) boot + between | (b) boot only |
  |---|---|---|
  | D0-like 1% | 0.944 | 0.839 |
  | D0-like 5% | 0.947 | 0.833 |
  | D1-like 1% | 0.944 | 0.841 |
  | biased (p_null 0.55) 1% | 0.937 | 0.817 |

  - (b) misses the seed-to-seed redraw of the entrapment null competitions, which the pair bootstrap holds fixed. The review's "double counting" concern does not hold in this model.
  - **Rule outcome:** (b) < 93%, so it is ineligible, and **(a) is adopted**. That is the existing default; no code default changes.
- **D0/D1 recomputed (`d0d1_ci_conventions.json`):** every verdict under (a) equals M3's.

  | dataset | scorer | α | arm P (lower CI) | arm C | overall |
  |---|---|---|---|---|---|
  | D0 | static_gbm_pred | 1% | inconclusive (0.00974) | supported | inconclusive |
  | D0 | static_gbm_pred | 5% | inconclusive (0.05037) | supported | inconclusive |
  | D0 | hyperscore | 1% | inconclusive | inconclusive | inconclusive |
  | D0 | hyperscore | 5% | inconclusive | supported | supported (fixed: P or C) |
  | D1 | static_gbm_pred | 1% / 5% | inconclusive | supported | inconclusive |

  - "Consistent with control" at 1% on D0 holds under (a) and would not under (b) (0.01002). It is a weaker claim, not used to rank.
  - With arm M now required for this scorer (A1.5) and inconclusive, the overall verdict stays inconclusive.

## D8 confirmation (A1.1; `d8_confirm.json`)
- **Data:** the first 10 runs of `data/manifests/e9_d8.tsv`, fetched independently (`scripts/m3e/fetch_d8.py`); all 10 sha256 match the E9 record. These runs were never used for selection.
- **Pipeline:** identical to D0. Each run is its own experiment.

| q ≤ | hyperscore peptidoforms | static_gbm_pred | ratio of sums | per-run geo-mean ratio [run-bootstrap 95%] | proteins (hs → gbm) |
|---|---|---|---|---|---|
| 0.1% | 41,517 | 97,091 | 2.34 | 2.25 [1.94, 2.66] (9 runs) | 5,273 → 13,593 |
| 1% | 68,506 | 110,105 | 1.61 | 1.55 [1.41, 1.71] (9 runs; gbm > hs in 10/10) | 12,702 → 18,467 |
| 5% | 91,533 | 117,665 | 1.29 | 1.27 [1.19, 1.35] (10 runs) | 12,803 → 18,477 |

- **Comparator:** hyperscore only. The strongest comparator (A1.9: Sage with RT prediction plus a predicted-feature rescorer) is tree item 3f and was not run. **These gains are not gains over the state of the art.**
- **Entrapment, pooled over runs (CI (a)):**

  | scorer | α | arm P | arm C | overall |
  |---|---|---|---|---|
  | static_gbm_pred | 1% | 1.032% [0.967, 1.098], inconclusive, cwc | 0.867%, supported | inconclusive (P, C, M required) |
  | static_gbm_pred | 5% | 5.21%, inconclusive | 4.68%, supported | inconclusive |
  | hyperscore | 1% | 1.02%, inconclusive | 1.11%, inconclusive | inconclusive |
  | hyperscore | 5% | 5.25%, inconclusive | 4.86%, supported | supported |

  - LB 0.53% at 1%: nothing violated.
  - Per run at 1% (frozen scorer): arm P 4 supported / 5 inconclusive / 1 not evaluable; arm C 8 / 1 / 1. Run E9_0007 is not evaluable: 10k MS2, 200 IDs.
- **Reading:** D8 confirms the direction and size of the dev ID gain over the fixed score on runs never used for selection. The gain is larger than on D0 (+25%), because hyperscore is weak on these heterogeneous runs. The entrapment pattern matches D0: P inconclusive at the known ~20% paired overstatement, C supported, nothing violated.

## Speed, production configuration (review M6; `bench_prod.json`)
- **Setup:**
  - Persistent process; staged NVMe env; under `flock /tmp/vulcan_bench.lock`.
  - Node load average ~50 (shared; other agents' jobs).
  - Per-run label-blind NCE choice as in production.
  - Prediction cache = "first contact": only the shipped NCE-29 library; every other NCE fills lazily.
  - D0 run 1 was the warm-up and is excluded. Measured: 10 D8 runs + D0 run 2.
- **Per run:**

  | stage | median wall | median core-s | range wall |
  |---|---|---|---|
  | scoring = features incl. NCE choice + GBM + double competition + picked protein | 4.3 s | 5.1 | 0.8–10.8 s |
  | per 10⁵ MS2 | 9.3 s | 9.6 | |
  | end-to-end (store load + GPU search + scoring) | 4.9 s | 5.4 | 1.9–13.2 s |

  - Components (median): features 2.8 s; GBM apply on CPU 0.7 s / 1.4 core-s; search 0.14 s; store load 0.8 s.
- **D0 run 2 (NCE 23):** scoring 10.8 s / 9.1 core-s.
  - Features 7.3 s, of which NCE grid 3.8 s and predictions 2.8 s.
  - 14% MS2 / 15% RT cache hits, despite following a warm-up run at the same NCE.
  - M3's "2.1 s / 2.8 core-s" was a cache-resident best case, **~4–5× below** the production cost for this run.
- **One-time costs per process:** imports 2.2 s (staged env; minutes on the Lustre env), model load 11.1 s, NCE-29 library load 7.9 s, index upload 2.1 s.
- **Excluded:** RAW ingest (M2b: ~6.5 core-s/10⁵ MS2 warm) and the L1 index build (once per proteome).
- **Next, not done (no new variants allowed here):** RT cache keyed independently of NCE, a precomputed NCE grid, and a GPU GBM.

## Deviations, choices and caveats (all made before the corresponding result unless marked)
- **Specification choices.** Where A1/§6.4 leave a detail open, the choice is written in `vulcan/gates/stats.py` (docstring) and `config.yaml` before the run:
  - tail membership = the scorer's own per-run 1% threshold, at spectrum level;
  - exact CP bound;
  - delta-method paired difference;
  - G1 Fisher-z 99% CI;
  - G3 point-estimate rule on all nulls *and* clean species nulls;
  - G4 τ from between-run KS.
- **"Without RT"** is an ablation of the frozen model (RT features neutralised), not a retrained model.
- **"Host proteome" in A1.3** is read as the null *source's* proteome. The human host proteome would remove every target, so it cannot be meant.
- **G4.** The partition was specified as charge {2, 3, 4+}. The training data contain 404 charge-1 decoys, which are kept as their own stratum (stricter). The tolerance τ is loose (0.54) because training runs differ a lot from one another.
- **G3** has little power on the clean subset (~100 targets per bin, where one miss = 1.03 pp). The pre-specified rule makes this a FAIL. It is reported as such and not re-run with more data.
- **Correction before fetching.** The DN proteome lookup gave the wrong *Synechococcus* strain; this was corrected before any DN file was fetched (DATA.md).
- **G1 classifier:** held-out pair AUC 0.70 (0.70 without rt_pred); D0 rank-1 AUC 0.77.
- **Operational.** `pkill -f`/`pgrep -f` matched my own shell twice (LESSONS); no data were affected. A CUDA OOM from sharing the GPU with E9 workers cost one feature rerun.
- **E9's code and files were not modified.** I re-downloaded the 10 D8 RAWs myself.

## Proposals for the lead (NOT applied; each would need a dated amendment, and both are informed by these results, so result-driven risk applies)
1. **Null-source admission check.** Before a null project enters G2/tail-G2, screen it label-blind for real human material, e.g. a human-only search at 1% with a cut on human protein groups such as skin proteins, or an extended contaminant list. Any such rule must be fixed before re-running and applied to all sources. It must not simply drop *S. aureus* after seeing that the gate passes without it; the review rejected exactly that move ("do not drop a source after seeing that the gate passes on the other"). The offset-null tail hits suggest the same issue for offsets with mis-assigned precursors.
2. **G3** needs a minimum per-bin n (or a CI-based rule) for the clean subset.

Until an amendment is made, the frozen scorer **has not passed the A1 gates**, so FDR claims and the T1 look remain blocked.
