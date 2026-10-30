# 3g decoy realism: analysis plan (written and committed BEFORE any new computation)

Tree item 3g; decisions D12, D14. Dev data only: E9/D8 outputs already on disk (`/tmp/vulcan_e9/out`), plus D8 runs
re-fetched by manifest order for the step-3 search. Excluded: `data/sealed`, T4, DN2 candidates (A2 pending), and any
non-human free-living-microbe PRIDE project as a spectral null source. Microbial **FASTA** sequence space may be used
(step 3), because no PRIDE project is consumed that way. Production code, metrics, labels, estimators and decoy
recipes are not changed. All new code lives in `scripts/decoy_realism/`. Any index this plan builds calls the
production builders (`vulcan.index.foreign.build_foreign`, `vulcan.index.decoys`) unchanged.

## Question
On arm-C null competitions (Arabidopsis entrapment E vs its pseudo-reverse decoy dec(E)), E wins 54.4% under the frozen
`static_gbm_pred` and 51.6% under hyperscore (E9, global 1%, N = 1000; constant in N). Three candidate explanations:
- **H1, homology / partial true match.** Human spectra share fragment ladders with homologous Arabidopsis peptides.
- **H2, naturalness.** Real sequences beat reversed ones through the features or predictors, independent of homology.
- **H3, something else.** Precursor mass, missed cleavages, length/charge mix, decoy collisions.

**Logical note, fixed before computing.** Pseudo-reverse preserves composition, length, mass, missed-cleavage count,
n_var and the C-terminal residue. So *marginal composition* cannot by itself move the E-vs-dec(E) coin. Only
order-dependent properties can: local motifs, terminal residues, ladder sharing with human sequences, and predictor
behaviour on order. "H2-composition" is therefore testable only in cross-family comparisons (step 3), never in the pair
coin itself. H3 effects on the coin can only act as moderators, or through the ~0.1% of decoys that are collision
shuffles rather than reversals.

## Common definitions
- **Coin.** c = N_E / (N_E + N_dec(E)) among **global pair winners** at or above a global threshold. Thresholds are the
  frozen L6 procedure exactly as in `scripts/e9/analyze.py`: spectrum winners, then best evidence per peptidoform,
  then paired competition, hash tie-break seed 11.
- **Thresholds.**
  - The **global 1% and 5%** score cutoffs are computed once on the *unfiltered* competition and held fixed for every
    filtered subset. Filtering removes whole E/dec(E) pairs only, so it cannot shift the human-target threshold.
  - **Primary decision threshold: global 5%** (about 22.6k E/dec(E) pairs at N = 1000, against about 4.1k at 1%). 1% is
    reported alongside.
- **Data, steps 1–2.** E9 arm C (`/tmp/vulcan_e9/out/C`, 982 runs, N = 1000 prefix). Hyperscore uses its own spectrum
  winners and global competition, as in E9.
- **Control.** E9 arm P (shuffled human entrapment E_2001 vs its pseudo-reverse decoy), from the same runs. Its global
  coin was 0.497 at 1%; it is recomputed with identical code.
- **Uncertainty.**
  - Wilson 95% CI for a coin.
  - For a difference of coins (or an excess ratio), a pair bootstrap with 2,000 replicates (seed 11). Pairs are the
    resampling unit because global pair winners are one per pair.
  - **Run/project cluster bootstrap** (resample projects) is a sensitivity analysis for the headline coins only.
- **Evaluable subset:** n_pairs ≥ 400 (SE ≤ 0.025). Smaller subsets are reported as "not evaluable" and never decide.
- **Excess:** e = c − 0.5. The **retained fraction** of a manipulated condition is e_cond / e_ref, with a bootstrap CI.

## Step 1: homology dose-response (H1 test)
Filters are applied to the **pair base sequence** E (I/L collapsed). A pair is removed or kept as a unit, whoever won,
so the filters are label-symmetric. The reference set is the human+contaminant target FASTA of the arm-C index
(`/tmp/vulcan_e2/human_contam.fasta`).

| Level (cumulative) | Removes a pair if E… |
|---|---|
| L0 | nothing beyond the current arm-C build (I/L-identical, single substitution, isobaric N↔GG/Q↔AG) |
| L1 | is within Hamming ≤ 2 of any same-length human target peptide (I/L collapsed) |
| L2 | shares any 6-mer with the human proteome |
| L3 | shares any 5-mer with the human proteome |
| L4 | has a Smith–Waterman local alignment (BLOSUM62, gap open 11 / extend 1) to any human protein with Karlin–Altschul E-value ≤ 10 (λ = 0.267, K = 0.041, n = proteome residues, m = peptide length); a numba SW over the whole proteome, so a "BLAST-like" exhaustive search |

- Each criterion is also reported **alone** (non-cumulative), together with the fraction of all foreign pairs it removes
  (from a random sample of 20,000 index pairs, seed 11).
- **Expectation, stated beforehand.** About 92% of foreign peptides share a 6-mer with human (E9), so L2 leaves about
  8% and L3 almost nothing. L3/L4 will probably not be evaluable at 1% and may not be at 5%. Strict k-mer filters also
  select short peptides, so every level is additionally reported **length-adjusted**: a Mantel–Haenszel pooled odds of
  E winning over length strata {7–8, 9–10, 11–13, 14–17, ≥18}, converted to a coin.
- **Mechanism-specific secondary (partial-ladder sharing; A2 draft (b) argues that shared b/y ladders, not internal
  k-mers, create partial matches).**
  - For a sequence x, h_y(x) = the longest C-terminal suffix of x (I/L collapsed) that occurs in the human proteome
    ending at a tryptic site.
  - h_b(x) = the longest N-terminal prefix of x that occurs in human starting after a tryptic site (or at a protein
    N-term).
  - Both are computed for x = E **and** x = the actual decoy sequence dec(E) (`fwd_seq` of the decoy row).
  - **Within-pair model:** logistic regression of 1[E wins] on Δ_y = h_y(E) − h_y(dec E) and Δ_b = h_b(E) − h_b(dec E),
    among threshold pair winners.
  - The intercept, converted to a coin at Δ = 0, is the asymmetry **not explained by ladder sharing**.
- **Decision rule, step 1 (at 5%; the length-adjusted coin decides, the raw coin is reported):**
  - **H1 supported:** the strictest evaluable cumulative level retains ≤ 50% of the L0 excess, with the upper 95% CI of
    the retained fraction < 1; **and** the Δ = 0 intercept coin has a 95% CI that includes 0.5.
  - **H1 rejected as the main cause:** the strictest evaluable level retains ≥ 50% of the excess, with a CI on its coin
    excluding 0.5; **and** the Δ = 0 intercept coin excludes 0.5.
  - **Partial / mixed:** anything else. The retained fraction is then reported as the homology-attributable share.

## Step 2: feature decomposition (which features carry the asymmetry)
On E9 arm C (and arm P as control), for the threshold pair winners at 5% (and 1%):
- **2a. Per-feature AUC.**
  - The comparison: E winners (positive) vs dec(E) winners, for every one of the 31 frozen-scorer inputs (22 base + 9
    predicted) and `rt_pred` (not a model input).
  - Hyperscore components: matched_b/y, longest_b/y, summed_intensity, matched_intensity_pct, poisson, peptide_len,
    missed_cleavages, charge.
  - Features are taken from the winning PSM (the unit's best evidence row).
  - The same table is computed on the hyperscore-tail pairs.
  - 95% CI by bootstrap. Reported, not gated. With equal scores on both sides, an AUC away from 0.5 shows *which
    features the E winners use* to reach the tail.
- **2b. Bulk.** The same AUCs among spectrum-level null winners (rank-1 is E or dec(E)) with no threshold, by
  frozen-score decile. This shows whether the tilt exists in the bulk or only in the tail.
- **2c. Within-spectrum paired.** Spectra whose top-5 contains both E and its own dec(E): the sign rate of
  (feature(E) − feature(dec E)) per feature, and of the frozen score. This isolates sequence-order effects at a fixed
  spectrum.
- **2d. Causal ablation through the frozen scorer** (the M3e "without RT" convention).
  - Feature group g is set to its per-run median over rank-1 PSMs for **all** candidates (both labels).
  - The frozen GBM is re-applied, spectrum winners, global competition and thresholds are recomputed, and the coin is
    recomputed. It is an ablation of the frozen model, not a retrain.
  - Groups:
    - MS2-pred: all 7 `pred_*`;
    - RT: rt_abs_resid and rt_resid_z;
    - each `pred_*` singly;
    - hyperscore-family: hyperscore, hs_runz, poisson;
    - ladder: longest_b/y/y_pct, matched_b/y/peaks;
    - mass/charge: ppm, isotope_error, z*.
  - Arm P under each ablation is the control, since its coin must stay ≈ 0.5.
- **2e. Spectrum-free predictor check** (does the predictor treat natural and reversed sequences differently with no
  spectrum?).
  - Samples: 5,000 Arabidopsis E and their dec(E), and 5,000 arm-P shuffled E and their dec(E) (seed 11).
  - Quantities: AlphaPeptDeep predicted fragment-intensity entropy, top-6 share, b/y balance at NCE 29, charge 2;
    predicted RT; and the M3e G1 sequence-classifier logit (`/tmp/vulcan_m3e/seqclf.pkl`, if present).
- **Decision rule, step 2.**
  - A feature group **carries** the asymmetry if ablating it retains ≤ 50% of the arm-C excess (CI upper bound < 1)
    while the arm-P coin stays within 0.5 ± 0.02.
  - If no single group does, the effect is called **distributed**.
  - If ablating predicted features (MS2-pred + RT) still leaves a coin whose CI excludes 0.5, the remainder is attributed
    to base or hyperscore features, which have no learned sequence model and so point to H1 or H3.

## Step 3: synthetic controls (new search; dev D8)
- **Index** `human_contam_3g`. Targets: human+contaminants as arm C. Foreign families are added as one foreign list to
  `build_foreign`; family membership is kept in a side table by sequence. All are pseudo-reverse decoys from the
  production recipe, seed 1001.
  - **A.** The existing arm-C Arabidopsis list, unchanged.
  - **B.** A natural, distant, unrelated prokaryotic proteome union: *E. coli* K-12 (UP000000625), *P. furiosus*
    (UP000001013) and the three bacterial DN FASTAs already in `data/fasta`. These are FASTA files only; no PRIDE
    project is used as a spectral source. Peptides go through the same arm-C neighbour filters vs human (the
    `arm_foreign` routine), and are additionally removed if I/L-identical to, or within single-substitution of, any A
    peptide, so that A and B are disjoint.
  - **S.** Shuffled A: each A peptide shuffled with the C-term fixed (`decoys.make_entrapment`, seed 2001), unique
    against T, A, B and all their decoys. These are non-natural with Arabidopsis composition.
- **Spectra.** D8 runs in manifest order 0–99, the same human dev runs as the E9 N = 100 prefix. Pipeline: the E9
  stream (`scripts/e9/stream.py --no-c --index-p <3g index> --work /tmp/vulcan_3g`), then `score_worker.py` with the
  frozen scorer. The prediction cache is copied, not shared. Nothing in `/tmp/vulcan_e9` is written.
- **Competitions** (global, N = 100; frozen score and hyperscore; thresholds 5% primary, 1% reported):
  - c_A (replicates E9 arm C in this search), c_B (natural vs reversed, distant), c_S (non-natural vs reversed:
    **control 1**);
  - **Natural vs natural (control 2):** the decoy-normalised share ratio
    R_AB = [N_A/(N_A+N_B)] / [N_decA/(N_decA+N_decB)], equivalently the odds ratio of c_A vs c_B;
  - **Natural vs shuffled-natural:** R_AS likewise (A vs S). This asks whether a *shuffled* decoy would also lose to
    natural sequences.
  - Each family also gets the step-1 cumulative homology filters (L1, L2) as a sensitivity analysis.
- **Pre-registered expectations:**

| | c_S (shuffled vs reversed) | c_B (distant natural vs reversed) | c_A − c_B | R_AS |
|---|---|---|---|---|
| H1 homology | ≈ 0.5 | smaller excess than A; falls with filters | > 0 | > 1 (A homologous) |
| H2 order-naturalness | ≈ 0.5 | ≈ c_A > 0.5, survives filters | ≈ 0 | > 1 |
| H2 composition-only | ≈ c_A (S has A's composition) | any | any | ≈ 1 |
| H3 | coin structure follows length/charge/mass/missed strata, not family | | | |

- **Decision rule, step 3 (5%, frozen score):**
  - Control 1 passes if the c_S 95% CI includes 0.5, or lies within 0.5 ± 0.015.
  - **H2-order supported** if the c_B CI excludes 0.5 **and** c_B − 0.5 ≥ 0.5·(c_A − 0.5), with the c_A − c_B CI
    including 0.
  - **H1 supported** if the c_B CI includes 0.5 (or c_B − 0.5 < 0.5·(c_A − 0.5)) **and** the c_A − c_B CI excludes 0.
  - Otherwise mixed.
  - Composition is supported only if c_S ≈ c_A.
- **H3 checks** (steps 1 and 3, descriptive): the coin stratified by length, charge, missed cleavages, precursor-mass
  decile, n_var, and whether dec(E) is a collision shuffle rather than a reversal. A Breslow–Day-type heterogeneity
  test is run across strata. H3 is flagged if one stratum carries more than 50% of the total excess while holding less
  than 25% of the pairs.

## Overall verdict rule
- **H1 / H2 / mixed** comes from the step-1 and step-3 rules together.
  - If they disagree, the verdict is "mixed", and both are reported with their numbers.
  - Step 2 names the carrier features.
- Step 3 decides H2 vs H1 only if control 1 passes. If it fails, H3 is flagged and investigated within strata.

## Step 4: implications and remedies (no adoption)
- Literature: decoy realism; natural vs reversed decoys; Wen/Noble and coworkers on entrapment and decoy choice.
- Remedies, each with a validity argument and the evidence from steps 1–3:
  - shuffled decoys;
  - Markov or natural-like decoy generation;
  - excluding or symmetrising naturalness-sensitive features;
  - entrapment-based calibration.
- **Nothing is adopted.** A decoy-recipe change needs a PI decision and a dated amendment before any claimed re-run
  (DESIGN §8).

## Outputs
`results/2026-10-27_3g_decoy_realism/{PLAN.md, metrics.json, NOTES.md, git_sha}`. Scripts go in
`scripts/decoy_realism/`. Timed sections are serialized with `flock /tmp/vulcan_bench.lock`. Long streams are not timed
claims.

## Deviations
Any deviation from this plan is recorded in NOTES.md with its reason and a statement of whether it was made before or
after seeing the affected result.
