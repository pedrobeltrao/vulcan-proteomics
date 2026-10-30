# Design issues found while implementing DESIGN §6 (v0.3) in `vulcan/stats`

Statistics engineer, 2026-10-26. DESIGN.md is not modified. Numbers come from `metrics.json` in this directory.

## Errors or impracticalities

### 1. G2 (v0.3) does not state its confidence level, and the 110k budget sits exactly on the 90%-power boundary
- The ≈110k budget per class only reproduces with a **one-sided 95%** upper bound, worst-case variance (fixed and learned outcomes independent), and a top-decile margin of 0.02. The top-decile clause dominates.
- I implemented exactly that. Monte-Carlo (10⁴ replicates, exact multinomial) at p_learned = p_fixed = 0.5, ρ = 0 (independent outcomes):

  | n per class | pass rate |
  |---|---|
  | 110,000 | 0.900 |
  | 107,570 (helper's n for 90% power) | 0.894 |
  | 116,620 (helper's n for 92% power) | 0.914 |

- So 110k gives ≈ 90% power with no margin; a replicate run can land below 90%.
- With realistic coupling between the two pipelines (ρ = 0.5), 110k gives 0.995, and 54k suffices.
- At p_learned = 0.52 and 0.55 the gate fires in 100% of replicates at every n ≥ 16.6k.
- **Suggestion:** state "one-sided 95%" in §6.4, and budget ≈ 117k per class (or the `g2_sample_size` output) for worst-case independence.

### 2. The v0.2 absolute G2 at n = 16.6k fails fair scorers (confirms round-2 B1)
- A fair scorer passes with probability 0.005 (Monte-Carlo, Wilson 99% CI within [0.49, 0.51]).
- It is kept only as `g2_absolute` (a diagnostic).

### 3. The entrapment verdict rule cannot rate a valid but tight method "supported"
- The paired estimator is an upper bound with positive bias.
- Simulation (`sim_entrapment_verdict`: augmented search, 3 dev entrapment seeds, TDC that is correct but nearly tight):
  - true mean FDP 0.0104 vs α = 1%; paired mean 0.0124, giving **"inconclusive"**;
  - at 5%: true 0.047, paired 0.054, also inconclusive.
  - 10⁴ independent simulations confirm the bias: true 0.0080 vs paired 0.0095 vs combined 0.0100 at α = 1%.
- Only a tool whose q-values are conservative (simulated as reporting at 0.6α) gets "supported".
- This is a property of upper-bound estimators, not a code issue, but it should be expected. Our TDC lists will often be *inconclusive* when the score is well calibrated.
- **Suggestion:** state this expectation in §6.3. Optionally pre-register that "inconclusive with a lower CI of paired FDP̂ ≤ α" is reported as "consistent with control".
- With these defaults, the 1.5α upper-CI criterion was easily met when N_T ≈ 5–9k.

### 4. "Seeded random" tie-breaking (§4 L4) conflicts with "ties broken by hash(sequence, seed)" (§4 L1)
- A seeded RNG stream depends on row order. Permuting input rows can therefore change winners, which is exactly the positional leak that round-2 [M4] warns about.
- Implemented both: `compete(..., seed=)`, and `compete(..., tiebreak=hash_tiebreak(sequence, seed))`. Production should use the hash. A test shows the hash is permutation-invariant.
- **Suggestion:** make §4 L4 say "hash(sequence, seed)".

### 5. RESET's training phase often cannot start on small or weak data
- With α_train = 1% and c_train = 3, the first positive set needs ≈ 300 pseudo-targets above the first training decoy.
- In a 2,000-peptide simulation with a moderate initial score, training never started in **87% (single) / 74% (two-decoy)** of replicates.
- RESET then falls back to the initial score with factor 2. This is valid but wastes power.
- The paper is silent on this fallback. I implemented "stop and keep the current score" (valid), plus an optional Percolator-style best-feature start (`init_score=None`), which uses training-decoy labels only.
- The headline RESET Monte-Carlo uses a 5,000-peptide regime in which training always starts. The weak regime is reported separately (`reset_*_weak_init`).
- **Suggestion:** pre-register the fallback, and possibly a relaxed α_train (e.g. 5%) when fewer than N positives are found.

### 6. Single-decoy RESET at α = 0.1% needs T ≥ 2000 before it can report anything
- FDR̂ = 2(D_est+1)/T ≤ 0.001 requires T ≥ 2000 even with D_est = 0.
- This is consistent with 0.1% being "report, not gate", but worth stating.

### 7. §6.7 "mean FDP ≤ α" needs a Monte-Carlo tolerance
- An estimator with FDR exactly α fails a literal test half the time.
- Tests use mean − 2.326·SE ≤ α.
- In all 10⁴-replicate runs the raw means were ≤ α anyway (see summary.md).

### 8. XL minimum-count rule is unspecified
- "Minimum-count rule following xiFDR conventions" has no number. It is implemented as `min_tt` (default 0).
- (TD−DD)/TT also gives q = 0 for short lists (no "+1"), so the conservative TD/TT must always be shown, as DESIGN says.
- A `plus_one` option is provided.

### 9. FLR "sites ≤2 residues apart" is ambiguous
- Implemented as: another candidate acceptor (S/T/Y or the decoy residue) lies within 2 residues of the localized site.
- The alternative reading is two *localized* phosphosites ≤ 2 apart (multi-phospho only).
- Please confirm which reading is intended.

### 10. Family-max: the "foreign" family needs its own decoy families
- §4 L1 lists T, E_s, dec_d(T), dec_d(E_s), foreign. The combined arm also needs dec_d(foreign); otherwise the foreign entrapment competes without decoys and q-values are about 2× liberal (round-1 M6).
- `ReplicateSpec` supports arbitrary families. Top-2 per family is needed and sufficient for the winner and the runner-up score: a test shows top-1 is not sufficient.

### 11. Group-walk merging with two crossed factors (closed: specificity × chimeric rank) has no order rule
- Implemented per factor: `pool_chimeric_rank` (≥3 pooled), `merge_adjacent` (offset bins by mass), `merge_rare_to_other` (PTM).
- For a crossed cell below 200 after these rules, §4 L4 does not say which factor merges first. Since min size is only a power rule, cells can simply stay small, but this should be stated.
- Note: G4-driven splitting must not be decided on estimation decoys (round-2 M1). The API only accepts label-blind counts (`label_blind_pair_counts`).

### 12. Picked protein group: data-dependent regrouping (rescue)
- The rescue cutoff is derived from decoy-based protein q-values, so the final grouping depends on decoy labels. There is no proof.
- Simulation (isoform pairs with shared peptides, 10⁴ replicates):

  | α | picked group, rescue on | classic, no rescue |
  |---|---|---|
  | 1% | 0.0060 | 0.0052 |
  | 5% | 0.0165 | 0.0340 |

- This is valid here, but it is the only estimator without a theorem, so protein-level entrapment (§6.3) remains necessary.
- The implementation uses (D+1)/T, not the reference's (D+1)/(T+1), and omits the reference's connected-component regrouping.

## Confirmed (no issue)
- RESET constants and +1 placement (single: train 3, report 2·(D_est+1)/T; two-decoy: 2 and 1) match Freestone et al.
- The FDRBench paired formula matches `FDPCalc.java` (`(n_E + 2·n_p_t_s + n_p_s_t)/(n_E + n_T)`), verified with a literal port.
- Biased null p = 0.55 in a tight TDC regime: FDP 1.19%/6.09% at α = 1%/5%, ratio ≈ 1.22 = p/(1−p). The differential G2 fires in 100% of replicates.
