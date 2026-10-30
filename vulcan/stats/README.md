# vulcan.stats: FDR estimation and verification

This package implements DESIGN.md §4 L4 and the pre-registered protocol in §6 (v0.3). It is pure numpy/scipy. scikit-learn is optional and only needed by `LinearSVMLearner`.

Conventions used throughout:
- Higher scores are better, and NaN counts as -inf.
- Competition ties are broken by a seeded RNG, or by an explicit `tiebreak` key. The recommended key is `_util.hash_tiebreak(sequence, seed)`: unlike an RNG stream, it does not depend on row order, so permuting rows never changes a winner.
- FDR estimates count every item with score ≥ s, so all decoys tied with a target are counted when that target is accepted ("decoys first").
- q-values are monotone (the minimum over lower thresholds) and capped at 1. "Report q ≤ α" is identical to the SeqStep+ rule (the largest list with FDR̂ ≤ α).

Pre-registered seeds (§6.1) are in `_util.py`: `SEEDS = (11, 23, 47, 83, 131)`, `DECOY_SEEDS`, `ENTRAP_SEEDS_DEV/SEALED`, `CALIB_DECOY_SEED`. All Monte-Carlo simulations split their replicates evenly over `SEEDS`.

| Module | What it does |
|---|---|
| `tdc.py` | `compete`, `psm_competition`, `tdc_qvalues`, `double_competition`, `multi_decoy_competition` |
| `groupwalk.py` | `group_walk`, `group_walk_batch`, the merge rules, `label_blind_pair_counts`, `check_pair_invariant` |
| `reset.py` | `reset` (single/two-decoy), `LDALearner`, `LinearSVMLearner`, `two_decoy_competition` |
| `entrapment.py` | `lower_bound`, `combined`, `paired`, `paired_counts`, protein-level rule, `evaluate_arm`, `verdict`, `combine_arms` |
| `protein.py` | `picked_group_fdr` (with rescued subsets), `classic_protein_fdr`, `picked_protein_fdr_simple` |
| `flr.py` | `decoy_aa_flr`, `flr_curve`, `flr_strata` |
| `xl.py` | `xl_qvalues`, `merge_level`, `xl_level_fdr` |
| `g2.py` | `g2_differential`, `g2_absolute` (diagnostic), `ks_winners`, `leak_auc` (G5), `g2_power`, `g2_sample_size`, `simulate_g2` |
| `familymax.py` | `top_k_per_family`, `reconstruct_competition`, `replicate_specs` |
| `global_fdr.py` | L6: `psm_winners`, `global_competition`, `run_competitions`, `naive_union`, `score_comparability`, `runner_up_stats`/`calibrate`, `global_outputs` |
| `sim.py` | Monte-Carlo simulators (used by `tests/stats` and `scripts/run_stats_sims.py`) |

---

## 1. Target-decoy competition (`tdc.py`)

**Estimator.** FDR̂(s) = factor·(D(s)+1)/max(T(s),1), with factor = 1 for TDC.

**Assumptions.**
- *Equal chance:* a null competition is won by the target with probability ½, independently of the winning score.
- *Independence:* these coin flips are independent across competitions.

Under these, the "+1" gives finite-sample FDR control (Barber & Candès 2015, knockoff+; He et al. 2015, arXiv:1501.00537).

At PSM level, multiplicity breaks independence. Many spectra of one absent analyte share one coin flip. PSM lists are therefore labelled "not FDR-controlled".

**Paired double competition** (Lin, Short, Noble, Keich 2022, JPR 21:2023):
1. All hypotheses and candidates of a spectrum compete.
2. The best PSM per unit is kept: per peptide for DDA, per precursor for DIA (pass `unit_id` = (peptide, charge)).
3. Each target unit competes with its own decoy, using units that share a `pair_id`.
4. (D+1)/T is applied to the pair winners.

Each stage has its own seeded tie-breaking stream.

**Multiple decoys** (Emery, Hasam, Noble, Keich, RECOMB 2020; arXiv:1907.01458), optional.
- With d exchangeable decoys, the winning score is the maximum of the d+1 scores.
- Let r be the rank of the original. The competition is a target win if r ≤ k_t, a decoy win if r > d+1−k_d, and ignored otherwise.
- FDR̂ = (k_t/k_d)·(D+1)/T. With d = 1 this is plain TDC.

## 2. Group-walk (`groupwalk.py`)

Freestone, Short, Noble, Keich, Bioinformatics 38:i82 (2022). This follows the CRAN `groupwalk` reference code.

**The walk.**
- Each group is sorted by winning score, ascending.
- The walk starts with every hypothesis included. At each step it peels the lowest remaining hypothesis of one group.
- Before every peel it computes FDR̂ = factor·(D+1)/max(T,1) over the remaining hypotheses. The running minimum is the q-value of the hypothesis being peeled.

**Choosing the group to peel.**
- Phase 1: while some group has had fewer than K = 40 peels, peel the group with the fewest peels.
- Phase 2: peel the group with the most decoys among its last K peels.
- Ties between groups are broken at random (seeded), or by lowest group index (`tie_break="first"`, which `group_walk_batch` uses).
- Within a group, at equal scores, targets are peeled before decoys.

**Validity requirements.**
- The partition must be fixed before estimation, using label-blind information only.
- Membership must be identical for a target and its decoy (`check_pair_invariant`).
- The TDC null assumptions must hold within each group.

**Merge rules** (DESIGN v0.3 §4 L4). The minimum size of 200 pair competitions is a power rule. Counts are label-blind: `label_blind_pair_counts` counts distinct pairs, not winners.
- `merge_adjacent(order, counts, 200)` handles ordered factors, such as mass-offset bins sorted by mass.
  - Repeatedly take the smallest block (ties go to the earliest) and merge it into its smaller adjacent neighbour (ties go to the preceding one).
  - Stop when every block has ≥ 200, or when one block remains (the pooled-TDC fallback).
  - The rule is deterministic.
- `merge_rare_to_other` handles unordered PTM classes. Classes below 200 go to "other", which is kept even if it is still below 200.
- `pool_chimeric_rank` pools chimeric ranks ≥ 3.

`factor=2` gives the "RESET + Group-walk" adaptation (c = 2/3). It is **unproven**, and our simulation of it is exploratory.

## 3. RESET (`reset.py`)

Freestone, Käll, Noble, Keich. bioRxiv 10.1101/2023.10.26.564068; RECOMB 2024; JPR 2025 24(5):2266. In v0.3 RESET is arm E4. The default static-first scorer does not need it.

**Input.** Winners after PSM-level and paired peptide-level competition on the search-engine score. Competition is never redone with the learned score.

**1. Split.** Each winning decoy goes to *training* or *estimating* with probability ½ (seeded).

**2. Train.** Pseudo-targets are the targets plus the estimating decoys, with their labels hidden. There are 5 iterations, starting from the search-engine score (or, with `init_score=None`, the best single feature, Percolator-style). Each iteration:
- computes FDR̂_train = c_t·(D_train+1)/max(PT,1);
- takes as positives the pseudo-targets with q ≤ 1%, and as negatives all training decoys;
- fits the learner and rescores.

c_t = 3 for single-decoy and 2 for two-decoy. If there are fewer than 10 positives, training stops and the current score is kept, which is still valid.

**3. Report.** Drop the training decoys. Then FDR̂ = c_f·(D_est+1)/max(T,1), with c_f = 2 (single-decoy) or 1 (two-decoy).
- The **+1 is added to the estimating-decoy count before multiplying by c_f** (SeqStep+ with c/(1−c) = c_f).
- So at α = 0.1%, single-decoy RESET needs T ≥ 2000 before it can report anything.

**Two-decoy mode.**
- Each target competes against two decoys from independent shuffles. They must be exchangeable, so reverse + shuffle is not allowed.
- Winning decoys are split at random, not by decoy database.

**Learners.** `LDALearner` (the default; closed-form with ridge shrinkage) and `LinearSVMLearner` (sklearn `LinearSVC`).
- The SVM picks C+, C− from {0.1, 1, 10} by 3-fold CV on the training data only, as in Percolator.
- Any object with `fit(X_pos, X_neg)` and `decision_function(X)` can be plugged in.
- There are no CV folds over estimation data.

## 4. Entrapment (`entrapment.py`)

Wen, Freestone, Riffle, MacCoss, Noble, Keich, Nat Methods 22:1454 (2025). The formulas were checked against FDRBench `FDPCalc.java`, and the test suite includes a literal port of it.
- Lower bound: LB = N_E/(N_T+N_E). It is only a failure detector.
- Combined: N_E(1+1/r)/(N_T+N_E). This is an estimated upper bound.
- Paired (r = 1): (N_E + N_{E≥s>T} + 2·N_{E>T≥s})/(N_T+N_E). An exact E/T score tie counts as E > T.
- Protein level:
  - A group with ≥ 1 original target counts as a target.
  - For the paired estimator, an entrapment group is represented by its first member, as in FDRBench.

**Decision rule** (v0.3). The statistic is the mean FDP̂ over replicates (entrapment seed × decoy seed). Its 95% CI is mean ± 1.96·√(Var_boot(mean) + s²_between/J).
- Var_boot comes from a Poisson(1) bootstrap of target–entrapment pairs. One weight is drawn per pair key and shared across replicates.

| Verdict | Condition (checked in this order) |
|---|---|
| not-evaluable | α·N̄_T < 20, where N̄_T is the mean over replicates |
| violated | lower 95% CI of the mean LB > α |
| supported | mean ≤ α and upper CI ≤ 1.5α |
| supported-noisy | mean ≤ α and upper CI > 1.5α |
| inconclusive | otherwise |

- Gated levels are 1% and 5%; the others are reported only.
- `combine_arms` requires P+C for static low-capacity scorers and P+C+M for sequence-capacity models. Fixed scores need P or C. Any arm that is violated makes the overall verdict violated.

## 5. Protein FDR (`protein.py`)

Savitski et al., MCP 14:2394 (2015) (picked protein); The, Samaras, Kuster, Wilhelm, MCP 21:100437 (2022) (Picked Protein Group FDR; reference code `kusterlab/picked_group_fdr`).

**Procedure.**
1. **Grouping.** Subset grouping: proteins with identical peptide sets form one group, and subset proteins are merged into the largest superset protein. Razor and Occam assignment are never used.
2. **Group score.** Shared peptides are discarded. The group score is the best unique-peptide score.
3. **Picked competition.** Shuffle the groups (seeded), then sort by score. A group is dropped if any member's target/decoy pair was already seen; otherwise the pair keys of its leading proteins are marked as seen.
4. **q-values.** (D+1)/T. The reference implementation uses (D+1)/(T+1).
5. **Rescued subsets.** Filter peptides at the score of the 1%-protein-FDR cutoff and regroup. First-pass groups that do not appear in the new grouping are kept as "obsolete": they compete and count in the FDR, but are not reported. Then rescore and recompete.

**Not reproduced.** The connected-component regrouping of groups without unique peptides.

**Status.** Picked TDC on proteins is TDC with fair coins for absent proteins whose decoy mirrors them. Rescued subsets have only empirical support.

`classic_protein_fdr` is for comparison only.

## 6. Decoy-amino-acid FLR (`flr.py`)

Ramsbottom et al., JPR 21:1603 (2022); aggregation per Camacho et al., JPR 2023.

**Estimator.** FLR̂ = n_pA·(N_STY/N_A)/n_pSTY.
- N_STY and N_A count candidate residues over the distinct localized peptides in the reported set.
- Ala is the primary decoy residue. Leu is secondary, with I counted as L.

**Assumption.** A false localization picks an acceptor roughly uniformly. This requires residue-blind localization scoring (§4 L3).

`flr_curve` gives the estimate at each localization-score threshold, monotone by default. `flr_strata` labels each site mono/multi-phospho and "close", meaning another candidate acceptor lies within 2 residues.

## 7. Cross-link FDR (`xl.py`)

Walzthoeni et al. 2012; Fischer & Rappsilber 2017 (xiFDR); Lenz et al. 2021.
- At each threshold: FDR̂ = max(TD−DD, 0)/max(TT,1), always reported alongside the conservative TD/TT.
- Entities are merged to the claim level before FDR: CSM, peptide pair, residue pair, or PPI. The entity score is the best member score.
- FDR is computed separately for self and heteromeric links. A link to protein A's own decoy counts as self (`base_of`).
- Options:
  - `plus_one` adds 1 to the numerator;
  - `min_tt` is the minimum-count rule, which DESIGN leaves unspecified (default 0).
- Neither form has a finite-sample guarantee.

## 8. G2 null test (`g2.py`)

**Differential test** (v0.3). On the same null competitions, the learned pipeline's target-win rate may exceed the fixed pipeline's by ≤ 0.01 overall and ≤ 0.02 in the fixed-score top decile.
- The test is paired and one-sided: d̂ = mean(B−A), and the upper bound is UCB = d̂ + z·sd(B−A)/√n.
- The default confidence is one-sided 95%, which is the level implied by the ≈110k budget.
- It is run separately per instrument class.

**Planning helpers.** `g2_sample_size` and `g2_power` use a normal approximation for the joint power, taken as the product of the two clauses (a lower bound). Variance is worst case (independence) unless `rho` is given. `simulate_g2` is an exact multinomial Monte-Carlo.

**Diagnostics** (not gates):
- `g2_absolute`, the v0.2 Wilson-CI band [0.49, 0.51];
- `ks_winners`, comparing target-winner and decoy-winner scores;
- `leak_auc`, the G5 alarm at AUC > 0.55.

## 9. Family-max entrapment (`familymax.py`)

The search keeps the top-2 candidates per family per hypothesis. Families are T, E_s, dec_d(T), dec_d(E_s), and foreign together with its decoys.
- Any replicate database is a union of families. Its spectrum-level winner, and its runner-up score (needed for the label-free Δscore), are always among those rows.
- `reconstruct_competition` therefore reproduces a separate search exactly when ties are broken by a storage-order-free key. The test checks this against brute force on tie-heavy data, and also checks that top-1 is *not* enough for the runner-up.

## 10. Repository-scale (L6) global FDR (`global_fdr.py`)

DESIGN §4 L6 and §7 E9. Every run is scored with one frozen, label-blind score, and a single pooling group is fixed before estimation.
1. Spectrum-level competition within each run (`psm_winners`).
2. Best evidence per unit across all runs, kept separately for the target and the decoy side.
3. Global paired competition: target unit vs its own decoy unit.
4. (D+1)/T q-values (`global_competition`).

`run_competitions` gives every run's own double competition at once, so each run's list carries a run q. `naive_union` is the E9
comparator: the union of the per-run lists at q_run ≤ α. `global_outputs` builds the master list (unit, best run, score, global q), the
per-run lists (q_run plus the unit's q_global) and the expected false counts (q×N, and D+1 at the threshold).

**Validity.** Take a null pair whose target and decoy scores are exchangeable in each run, independently across runs. The maxima over
runs are then exchangeable too, so the coin stays fair for any N. Calibration across runs changes power, not validity.
`score_comparability` is a G4-style diagnostic on decoys only (dev/training data, §6.4). It gives each run's KS distance to the other
runs' decoys, together with a null KS from two random halves of the same run. `runner_up_stats`/`calibrate` is a label-blind per-run
monotone map built from runner-up candidate scores, so it keeps validity. On E9 it lowered power (see results/2026-10-27_E9_scale).

**Caveat (E9).** If a learned score tilts null competitions even slightly towards targets, best-of-N pooling amplifies the tilt as N
grows. The fixed hyperscore stays flat with N; check a learned score with entrapment at the N where it will be used.
`tests/stats/test_global_fdr.py` contains the Monte-Carlo check. Global control holds at N = 1 and 20, while the naive union inflates
(≈1% → 15% at N = 50 in the simulation). It also checks equivalence to the per-run procedure and permutation invariance.

## 11. Verification

Run the fast suite with `pytest tests/stats` (about 30 s). Add `--runslow` for the ≥10⁴-replicate versions.

`scripts/run_stats_sims.py` writes `results/2026-10-26_stats_sims/metrics.json`. A Monte-Carlo test passes when mean FDP − 2.326·SE ≤ α. The raw means are reported in metrics.json.
