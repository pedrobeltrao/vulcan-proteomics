# A3.1.8 dev preconditions for v2_sep, A3.1.9 power check (dev / synthetic only)

- **Date:** 2026-10-28. **Agent:** the A3 implementation agent. This agent is not the step-4b or the step-5 agent.
- **Not read:** sealed data, anything under `/tmp/vulcan_a2`, any DN2 record or count, any E7b output.
- **Code state.** Every run used a clean git worktree at a committed state (implementation note 32):
  - G4 ran at 8bfb431, after the partition was committed in that same commit and before any estimate existed.
  - Arm M ran at f258a9d.
  - The Monte Carlo ran at f293526.
  - The V12 truncation rate ran at d5291cc. The first version (31af0be) was superseded after review M1.

  The first three ran before E7b landed its speed paths (bd2cfb2). The truncation run uses no module that bd2cfb2 changed.
- **MS1 tables** were regenerated with the committed `vulcan.ms1`. They are bit-identical (sha256) to E7's tables for D0 (2 runs), D1 and LM2486.
- **Results are reported as is.** Under A3.3.8, inclusion of v2_sep in T1 does not depend on these results; v2_sep's claims do.

## Verdicts

| Precondition | Rule | Result |
|---|---|---|
| (i) G4, partition chimeric stratum × (OT-HCD, HCD, charge) | §6.4 `g4_gate`: every stratum's KS vs the rest ≤ τ | **PASS**: τ = 0.298, max stratum KS = 0.078 |
| (ii) arm M, D7, seeds 2001–2003, chimeric search | no "violated" at 1% or 5%, overall or in any stratum | **PASS**: no violated verdict, under either CI convention |
| (iii) §6.7 MC of the full v2_sep competition and pooled picked protein groups | each gated list: mean FDP − 2.326·SE ≤ α, ≥ 10⁴ replicates | **PASS**: all 4 scenarios × 10⁴ replicates |
| A3.1.9 power | report pass probabilities, n_min, joint per tier, "expected not evaluable" flags | Reported below. **Tier 1 is not evaluable at the A2 design point** (V6 class 1, V12). With B ≤ |A| it is "expected not evaluable" at every reachable size |

## (i) G4 (`g4.json`; `scripts/a3/dev_g4.py`; partition `g4_partition.json`, frozen in 8bfb431)

**Setup.**
- Data: the 40 frozen v2 training runs; training decoys only.
- Decoys: stage-1 unit winners that are decoys under v2.
- Strata: r{1, 2, ≥3} × z{1, 2, 3, 4+}, OT-HCD/HCD.
- Statistic: `vulcan.gates.stats.g4_gate`, unchanged.
- τ = 0.298 is the 95th percentile of 4,398 within-stratum between-run KS values (10/50/90th percentiles 0.05 / 0.12 / 0.26).

| stratum | n decoys | KS vs rest |
|---|---|---|
| r1 z1 / z2 / z3 / z4 | 1,596 / 97,851 / 71,232 / 34,621 | 0.054 / 0.063 / 0.053 / 0.078 |
| r2 z1 / z2 / z3 / z4 | 95 / 52,220 / 40,161 / 19,644 | 0.137 / 0.032 / 0.044 / 0.041 |
| r3 z1 / z2 / z3 / z4 | 9 / 49,088 / 45,917 / 22,087 | 0.266 / 0.049 / 0.070 / 0.051 |

**Reported only:**
- hyperscore reference: pass (τ 0.660);
- v2 on the A1 charge-only partition: pass;
- pair-winner decoys after the full v2_sep competition: fail, only in r3 z1 (n tiny, KS 0.327);
- strata split by instrument model: fail, only in two charge-1 strata of Q Exactive HF (KS 0.39, 0.50).

**Caveat.** τ is large because the 40 training runs are very heterogeneous. G4 is therefore a weak test here. This is the rule as pre-specified (A1/§6.4), applied unchanged.

## (ii) Arm M (`armM.json`; `scripts/a3/dev_armM.py`)

**Setup.**
- Data: D7 yeast LM2486/2487/2490, searched against yeast + contaminants + P37840, with human entrapment at r = 1. These are the A1 arm-M indexes, one per seed (2001–2003).
- Search: chimeric, report 5. Each run has about 35–45k secondary groups.
- Scoring: v2 model; competition `chimeric_competition` with q_sep.
- Neighbour filter: drops 6.5% of group winners. No merge was needed; the pair counts per stratum are 65,702 / 26,992 / 15,993 (seed 2001).
- Estimator: combined, r = 1.0. CI convention (a) decides; (b) is reported.

**v2_sep results** (mean FDP̂_C, 95% CI, verdict):

| list | 1% | 5% |
|---|---|---|
| overall | 1.04% [0.84, 1.24], inconclusive (N_T ≈ 20,970) | 5.39% [5.07, 5.70], inconclusive |
| stratum 1 | 1.00% [0.74, 1.25], supported | 5.49% [5.09, 5.89], inconclusive |
| stratum 2 | 1.26% [0.77, 1.75], inconclusive (N_T ≈ 2,950) | 5.29% [4.28, 6.29], inconclusive |
| stratum ≥3 | not evaluable (α·N_T = 11.4 < 20) | 4.30% [2.88, 5.71], supported |

- No "violated" verdict appears under (a) or (b). The largest LB lower CI is 0.42% at 1% and 2.54% at 5%, both well below α.
- The hyperscore reference through the identical pipeline gives the same picture: overall 1.04% / 5.20%, inconclusive.
- Memorisation has no visible signature: v2_sep's FDP̂_C equals hyperscore's at 1%.

## (iii) §6.7 Monte Carlo (`mc_v2sep.json`; `scripts/a3/dev_mc_v2sep.py`)

**What is tested.**
- The committed `vulcan.chimeric.compete.chimeric_competition` is called unchanged. Its steps are:
  - stage-1 unit winners;
  - the CONGA neighbour filter over all group winners, in descending score order (rank-1 winners can be dropped);
  - paired double competition;
  - strata from the pair winner;
  - label-blind merge at < 200 pairs;
  - q_sep.
- Pooled picked protein groups use `picked_group_fdr` (subset grouping, rescue, seed 11) over the pair-winner scores pooled across strata, as `bench/t1/evaluate`.
- Synthetic data include:
  - stratum-heterogeneous null laws;
  - shadows (the rank-1 true peptide re-seen by a secondary group);
  - label-symmetric ion families, giving a neighbour-filter drop rate of 11.7–12.5%, the D0 level;
  - many null pairs whose best target and best decoy come from units of different strata (224–767 per replicate).
- Seed 131; 10⁴ replicates per scenario.

| scenario | merges | worst gated list: mean FDP − 2.326·SE (α) | protein 1% / 5% (mean FDP) |
|---|---|---|---|
| S1 dev-like (6,000 spectra) | none | union@1%: 0.00965 (0.01) | 0.86% / 0.88% |
| S2 small (700 spectra; strata 2+3 merged in 79% of replicates) | 1 \| 2_3 in 7,853 | union@1%: 0.00954 | 0.71% / 0.94% |
| S3 strong heterogeneity | none | union@1%: 0.00966 | 0.02% / 4.15% |
| S4 cross-rank pairs (small universe) | none | merged:1@1%: 0.00873 | 0.69% / 1.74% |

- Every merged-stratum list, every union and every protein list passes at 1% and 5% in every scenario.
- **Reported finding.** When strata 2 and ≥3 are merged (S2), the *raw* stratum ≥3 under the merged q is not controlled: 1.18% at 1%, 8.0% at 5%. q_sep controls the merged block, which is what A3.3.6 evaluates ("each merged stratum"). Any statement about a raw stratum inside a merged block would be uncontrolled, so T1 claims must name the merged block.

## A3.1.9 power check (`reviews/checks/a3_power.json`; `reviews/checks/a3_power.py`, synthetic, seed 131; revised after the independent review)

**n_min.** Values are per label per bin, by the V6 rule (a symmetric process at the dev rate passes one 90% TOST bin
with probability ≥ 0.90; 2·10⁴ MC replicates per n).

| row | n_min | basis |
|---|---|---|
| V6 | **22,700** | dev drop rate 0.118 |
| V12 | **19,200** | dev miss rate 0.0981 = the max-variance decile of `truncation_dev.json` |
| V14 | 7,200 | A2.6 rule (learned-CI pass probability 0.906 at 7,200) |
| V8 | 107,570 | `g2_sample_size()` |
| V11 | — | unused |

V12's dev source is a D0-host mixture with 12,000 D7 yeast nulls and production calibration (review M1). There the pooled miss rate is 0.0889 (T 0.087, D 0.091), and 64% of nulls have a candidate per label.

**Reachable DN2 sizes** are bounded by A2.9.3's stopping rule (review M3). The walk stops at the first count ≥ 10 projects at which |A| ≥ 120k, and only reaches 16 projects at low yield. The resulting reachable points are (J, |A|) = (16, 96k), (10, 120k), (10, 240k), (10, 450k), (10, 750k), for per-run A densities of 2k–25k. B/|A| ∈ {1, 3} is an assumption.

| \|A\| (J) | not evaluable (B/A 1 / 3) | tier 1 joint: independence [Fréchet lower, min upper] | tier 2 (independence) |
|---|---|---|---|
| 96k (16) | V2, V6 c1+c2, V12, V8 | 0 [0, 0] | 0 |
| 120k (10), A2 design point | V6 c1+c2, V12 (+V8 at B/A 1) | 0 [0, 0] | 0 |
| 240k (10) | same | 0 | 0 |
| 450k (10) | V6 c1+c2 | 0 | 0 |
| 750k (10), B/A 1 | V6 c1+c2 | 0 | 0 |
| 750k (10), B/A 3 | none | **0.891 [0.887, 0.948]** | 0.838 |

**Clean-scorer pass probabilities where evaluable:**
- V1: 0.75 (96k), 0.86 (120k), 0.95–0.96 above;
- V2: 0.88 at 120k, 0.99 above;
- V12: 0.83 at 450k, 0.995 at 750k;
- V13: 0.92–0.97. The sensitivity run at half the A-syn tail fraction gives 0.76–0.96;
- V14: 0.98–0.99;
- the G3, G2 and G5 rows: ≥ 0.96.

**Pass rate at the margin (MC SE about 0.003–0.01):**
- tail-G2 at fraction 0.60: 0.01;
- G1 at ρ 0.05: ≤ 0.007;
- TOST rows at a 1-pp gap: 0.048–0.056, i.e. at the nominal 0.05 within MC error;
- G2 at +0.01: 0.05;
- G5 at AUC 0.55: 0.5, because G5 is a point rule with no controlled level.

**Joint per tier.** Rows are simulated separately; no single synthetic DN2 per replicate was built. The json therefore gives:
- bounds valid under any dependence between rows (Fréchet lower, min upper);
- the independence product.

This is a disclosed departure from "simulated joint" (review M4).

**Flags (A3.1.9, read strictly as "can never pass": clean pass probability < 0.05 at every reachable point):**
- **At the A2 design point (|A| = 120k), tier 1 and tier 2 are not evaluable** (V6 class 1 and 2, V12).
- **With B ≤ |A|, V6 can never pass at any reachable DN2 size, so tier 1 and tier 2 are "expected not evaluable"** (json `B_over_A_1.0`).
- Only at the largest reachable yield with B ≥ 3|A| (750k A spectra at 10 projects) does tier 1 become passable (≈ 0.89), and tier 2 (≈ 0.84).
- The binding row is **V6 class 1**: about 1,400 top-decile winners per label at the design point, against 22,700 needed.
- Conclusion: the A3.5.2 risk is realistic. The most likely DN2 outcome for v2_sep is tier 1 "not evaluable". That blocks v2_sep FDR claims but is not final (A3.2.6).

## Files
| File | Contents |
|---|---|
| `g4_partition.json` | frozen partition (8bfb431, before estimation) |
| `g4.json` | G4 decision and reported variants |
| `armM.json` | arm M for v2_sep and hyperscore_sep, all lists, both conventions, merge maps |
| `mc_v2sep.json` | MC: 4 scenarios × 10⁴ replicates, every list |
| `truncation_dev.json` | dev top-5 truncation rate of the exhaustive v2-best candidate (V12 n_min input) |
| `metrics.json` | headline numbers |
| `logs/` | run logs |
