# 3g decoy realism: why do natural foreign nulls beat their pseudo-reversed decoys? (dev data only)

Plan: `PLAN.md`, committed as 41288c0 before any new computation. Numbers: `metrics.json`. Code: `scripts/decoy_realism/`.
- No sealed, T4 or DN2 data was used. No microbial PRIDE project served as a spectral null; microbial FASTAs were used as
  sequence space only (step 3).
- Production code, metrics, labels, estimators and decoy recipes are unchanged. The 3g index was built with the
  production builders (`build_foreign`, decoy seed 1001).
- Scratch outputs are in `/tmp/vulcan_3g` (node gpu-node) and are not committed.

## Bottom line
- **The asymmetry comes from real fragment evidence on the foreign side (H1 in a broad sense), not from a naturalness
  bonus in the predictors (H2).** Under the pre-registered rules the formal verdict is **"mixed"**:
  - step 1 is partial;
  - step 3 control 2 is **not evaluable**, because the "distant natural" proteome turned out to be truly present in
    the samples.
- **Real evidence has two parts:**
  1. Homologous conserved proteins.
  2. Plant material that is actually present (post hoc).
- **H2 is not supported:**
  - neutralising all predicted features leaves the coin unchanged;
  - the predictor treats natural and reversed sequences almost identically when no spectrum is involved;
  - the shuffled control is a fair coin.
- **No remedy is adopted.** Any decoy-recipe or arm-C change needs a PI decision and a dated amendment.

## Reproduction of the E9 numbers
Re-derived with the identical L6 procedure on N = 1000 (981 arm-C runs). Frozen score:
- **1%: 0.5438 [0.529, 0.559]**, n = 4,097 pairs; project-cluster CI [0.525, 0.566].
- **5%: 0.5208 [0.514, 0.527]**, n = 22,649.

Controls and comparisons:
- Hyperscore: 0.516 (1%) and 0.505 (5%).
- Arm P (shuffled vs reversed): 0.497 (1%) and 0.501 (5%).
- **Bulk (no threshold), arm C:** 0.5015 over 1.3M null spectrum winners. By score decile the coin is 0.500–0.502 up to the
  8th decile, then 0.506 and 0.508. **The tilt lives in the high-score tail**, as expected for real partial matches and
  not for a per-PSM bias.

## Step 1: homology dose-response (frozen score)

| Filter (cumulative) | 5% pairs kept | coin 5% [95% CI] | coin 1% (n) |
|---|---|---|---|
| L0 (current arm-C build) | 22,649 | 0.521 [0.514, 0.527] | 0.544 (4,097) |
| L1 + Hamming ≤ 2 | 16,540 (73%) | 0.523 [0.516, 0.531] | 0.562 (2,242) |
| L2 + shared 6-mer | 979 (4.3%) | 0.513 [0.481, 0.544]; length-standardised 0.485 [0.44, 0.53] | 0.512 (123, n.e.) |
| L3 + shared 5-mer | 0 | not evaluable (100% of random foreign peptides share a human 5-mer) | – |
| L4 + SW E ≤ 10 | 0 | not evaluable | – |
| **SW E ≤ 10 alone** (removes 5.7% of tail pairs) | 21,364 | **0.514 [0.507, 0.521]**, excess retained 0.67 [0.34, 1.17] | **0.520 [0.504, 0.536]**, retained **0.46 [0.11, 0.96]** |

- **Rule outcome: partial/mixed.**
  - L2 is the strictest evaluable level. It retains 0.61 of the excess, but the CI [−0.88, 2.2] is uninformative.
  - The within-pair ladder model leaves a coin at Δ = 0 of 0.515 [0.509, 0.522], which excludes 0.5. That blocks
    "H1 supported".
- Dose-response is **not monotone at Hamming ≤ 2**: the excess grows (0.562 at 1%). Single-residue-level homology is not
  the driver.
- **Detectable homology does carry part of the effect.** Removing the 6% of pairs with a significant local alignment to
  human removes **about one third (5%) to one half (1%) of the excess**. Hyperscore behaves the same way: SW alone gives
  0.500 / 0.496.
- **Ladder sharing.** E shares longer human N-terminal (b) prefixes than its own decoy does: h_b 5.02 vs 4.85 at 1%. A
  b-ladder advantage predicts the winner (logit coefficient 0.19 [0.14, 0.25]).
  - For hyperscore, ladder sharing explains its tilt completely (coin at Δ = 0 is 0.499 [0.489, 0.506]).
  - For the frozen score, a residual remains (0.515).
- **H3 strata.** The pre-registered flag did not fire.
  - Length: p = 10⁻¹⁰ at 1%. Coin 0.507 for 7–8-mers vs **0.706 for ≥ 18-mers**. Long peptides are 9% of pairs and carry
    41% of the excess. The top mass decile behaves the same way (0.678).
  - Charge, missed cleavages and n_var are homogeneous.
  - Collision-shuffle decoys: 12 pairs, irrelevant.
  - Long natural peptides are where partial ladders accumulate the most matched ions.

## Step 2: which features carry it
- **Tail AUC, E winners vs dec(E) winners at 1%** (`metrics.json:features`).
  - **Higher for E: delta_next 0.564, peptide_len 0.556, matched_y 0.554, matched_peaks 0.546, hs_runz 0.541,
    pred_pearson 0.539, pred_obs_tic 0.537, pred_sa 0.535, longest_y 0.535.**
  - Lower for E: poisson 0.447, rt_resid_z 0.466 (E has *larger* RT residuals), pred_top6 0.474.
  - At 5% only delta_next (0.526) and peptide_len (0.516) depart from 0.5.
  - In the same table arm P is flat (all within 0.48–0.52).
  - Reading: the E winners look like *better matches on every evidence axis*, the profile of genuine partial matches.
    The RT channel points the other way.
- **Within-spectrum pairs** (E and its own decoy both in the top 5 of a spectrum where one of them wins; 10.7k pairs).
  - E has a longer y-ladder: longest_y 0.532 and matched_y 0.518.
  - E has fewer b ions: matched_b 0.474.
- **Causal ablation through the frozen model** (group set to its run median, model re-applied, global competition
  redone). Coin at 1% / 5%:
  - none: 0.544 / 0.521;
  - **all predicted (MS2 + RT): 0.555 / 0.523**;
  - MS2-pred: 0.555 / 0.520;
  - RT: 0.549 / 0.520;
  - each single `pred_*`: 0.541–0.558;
  - hyperscore-family, ladder, mass/charge, other: 0.540–0.550 / 0.520–0.523.
  - **No group carries it**: every retained fraction is about 1, with wide CIs of about [0.6, 1.6].
  - The verdict is **distributed**. The arm-P control stays within 0.48–0.51 throughout.
  - Only neutralising all base evidence features lowers the 5% coin (0.508 [0.501, 0.516]). That scorer is
    degenerate (14k targets at 1% instead of 361k), so it is only suggestive.
- **Spectrum-free predictor check** (AlphaPeptDeep; 5k pairs each). Arabidopsis E vs its decoy:
  - entropy sign rate 0.476, top-6 0.521, y-share 0.516;
  - predicted RT 0.512; expected RT residual against the human distribution 0.497, i.e. no RT advantage.
  - Human targets give the same pattern (0.49–0.51).
  - The G1 sequence classifier separates natural from reversed (logit sign 0.574; human 0.635), but **the scorer
    does not use that signal**, as the ablation shows.

## Step 3: synthetic controls (new search; D8 runs 0–99, 96 completed; frozen score; global; 5% / 1%)
- **Control 1 passes (shuffled A vs its reversed decoy): c_S = 0.504 [0.491, 0.516] / 0.494 [0.465, 0.523].**
  - Non-natural vs non-natural with Arabidopsis composition is a fair coin.
  - So **neither composition nor reversal per se** causes the effect.
- **Replication:** c_A = 0.517 [0.505, 0.529] / **0.549 [0.522, 0.576]**.
  - c_A − c_S = +0.014 [−0.003, 0.031] at 5% and **+0.055 [0.015, 0.093]** at 1%.
  - Hyperscore gives c_A 0.543 at 1% in this search.
- **Control 2 (B = E. coli, Pfu, B. cereus, S. aureus, Synechococcus): NOT EVALUABLE.**
  - c_B = 0.745 / **0.925**, and it survives the homolog filters: 0.787 after the 6-mer filter.
  - Post hoc (`diag_b`): the B target winners are **truly present proteins**: E. coli GroEL, DnaK, KatE, ArnA, CRP, ribosomal
    proteins, **BirA** (BioID), and *S. aureus* **protein A** (SAUSA300_0113, from IP beads). 50% of the excess comes from 2
    projects (PXD066795, PXD060271).
  - The rule's mechanical output ("mixed") is void, because the premise that B is a null is false.
  - Pfu alone gives 0.511 (n = 440; post hoc, uninformative).

## Post hoc diagnostics (after seeing step 3; no verdict changes)
- **Arm-C excess is protein-coherent.**
  - Permuting E_won across tail pairs (1,000×): E winners sitting in Arabidopsis proteins with ≥ 2 E winners exceed
    their decoy counterpart by **241 vs a permutation mean of 98** at 1% (p < 0.001). At 5% it is 884 vs 710 (p = 0.019).
  - So about 40% (1%) and 18% (5%) of the excess is *concentrated in specific proteins* beyond what a uniform tilt would
    produce. A per-PSM naturalness bias cannot do that.
- **Top E-winner proteins:**
  - **plant-specific: rbcL (RuBisCO, 13), atpB (chloroplast), LHCB, RBCS-1A, TUFA (chloroplast EF-Tu)**, which point to
    plant material actually present in some human samples;
  - **conserved: actin ×3, β-tubulin, HSP70, ClpB**, which are homologs.
- **Projects.** The 1% excess is heterogeneous across projects (χ² p = 9·10⁻⁵); half of it sits in 8 of 138 projects.
  The leave-one-project-out coin stays within 0.540–0.548, so the effect is **not one project**.

## Answer to 3g
| Hypothesis | Evidence | Status |
|---|---|---|
| H1, partial or true match | Tail-only; better on every evidence axis; SW-homologs carry ⅓–½ of the excess; b-ladder sharing predicts wins; protein-coherent (RuBisCO, actin, tubulin); worst for long peptides; hyperscore shows it too | **favoured (partial by the formal rule)** |
| H2, predictor naturalness | Ablating all predicted features: coin unchanged (0.555 / 0.523); spectrum-free predictor gap ≈ 0; RT residual favours the decoy | **not supported** |
| H2, composition | Shuffled-A control 0.504; reversal preserves composition anyway | **rejected** |
| H3, mass / missed cleavages / collisions | Homogeneous except length and mass (both follow from H1); 12 collision decoys | **not a separate cause** |

A residual is left after the measured homology/ladder covariates: Δ = 0 coin 0.515, 70% of the 5% excess. It is most
plausibly undetected homology plus truly present material (RuBisCO etc.), which the human-proteome ladder covariate
cannot see. That attribution is post hoc.

## Implications for TDC validity with pseudo-reverse decoys
- **The pseudo-reverse decoy is not shown to be unrealistic.**
  - Shuffled-vs-reversed is fair (arm P 0.497; S 0.504).
  - The scorer gives reversed sequences no "unnaturalness penalty" that a natural false target escapes.
- **The broken assumption is arm C's.** It assumes a natural foreign peptide is a pure null. On repository-scale human
  data it is not: homologous peptides earn genuine partial fragment evidence, and some foreign proteins are physically
  present (plant food or sample material, reagents).
  - E9's arm-C 1.21% at N = 1000 therefore overstates the FDP of a list whose false targets are true nulls.
  - The P/C gap (P 1.16%) is consistent with that.
- **What carries over to human absent targets** is the *neighbour* problem: absent human paralogs or isoforms of
  present peptides get partial matches that no decoy imitates.
  - This is the "neighbour" issue in Freestone, Noble and Keich's work on TDC correctness, and the "biased scoring"
    case in Danilova et al. 2019.
  - It is a property of what counts as "correct", not of reversal.
  - Its size for human targets is *not* measured here. Arm C overstates it, because foreign proteins add true presence.
- **Literature consulted:**
  - Wen et al., *Nat Methods* 2025 (entrapment framework; the combined estimator is valid only if entrapment ~ false
    targets);
  - Freestone, Noble and Keich, *J Proteome Res* 2024 (re-investigating TDC correctness; neighbours, chimeras);
  - Danilova et al. 2019 (bias in TDC FDR: target-favouring scores vs correlated spectra);
  - comparisons of decoy construction methods (reverse, shuffle, de Bruijn, DecoyPyrat).

## Candidate remedies, ranked (NONE adopted; each needs a PI decision and a dated amendment before a claimed re-run)
1. **Make arm C's null assumption true, not the decoys** (best validity-per-cost).
   - **What:** add an outcome-independent arm-C admission screen, fixed before any re-run, to the A1.3 neighbour
     filters:
     - an A2(b)-style ladder-identity screen: drop E if ≥ 20% of its b/y ions coincide with a human T peptide's in the
       window, or if its terminal ladders are shared;
     - an a priori taxon screen for plausibly present proteins (plastid/photosynthesis, food allergens);
     - a protein-level de-duplication rule.
   - **Validity:** it restores "E ~ false target" for the combined estimator (Wen et al.). Decoys and TDC are
     untouched.
   - **Caveat:** whatever the screen removes was chosen after this analysis, so it applies only to fresh data.
2. **Report arm P as the primary repository-scale validity check, and arm C as an upper bound until (1).**
   - **Validity:** arm P is exchangeable by construction (coin 0.497–0.501). The known share-effect overstatement only
     makes it conservative.
   - No recipe change is needed.
3. **Run the §6.2 shuffled-decoy variability arm (5 seeds) alongside pseudo-reverse.**
   - **Validity:** it detects any reversal-specific artefact.
   - It is **not expected to change the coin**: S vs reversed is already fair, and shuffled decoys lack homologs just
     as reversed ones do.
4. **Calibration via entrapment** (inflate D by κ estimated from an entrapment coin).
   - **Validity:** valid only if κ transfers from foreign to absent human targets, and these data show it does not
     (foreign presence and homology differ from human neighbour structure).
   - Not recommended.
5. **Markov or "natural-like" decoys.**
   - **Validity:** they target a naturalness gap that was not found. They risk decoys with partial matches to real
     peptides (conservative, lost power), and their exchangeability with absent targets is unproven.
   - Not recommended.
6. **Excluding or symmetrising predicted features.**
   - **Validity:** the ablation shows the predicted features do not carry the effect, so this would cost IDs (E9: the
     frozen scorer gives 1.55–1.7× hyperscore) for no validity gain.
   - Not recommended.

## Deviations from PLAN (all decided before the affected result unless marked)
- **Ablation neutral value:** the run median over **all** PSMs, matching the actual M3e code (`neutral_rt`), not over
  rank-1 PSMs as PLAN said. Decided before running.
- **Length adjustment:** direct standardisation to the L0 length mix, instead of a Mantel–Haenszel odds. MH reduces to
  the raw coin for one-sample proportions. Decided before running.
- **Step-3 run count:** 96 of 100 runs completed. Ingest or search failed for E9_0025, E9_0042 (PXD074511 reader), E9_0052 and E9_0092 (PXD072218, no PSMs); 3 of these 4 failed in E9 too. `cohort --skip-failed` was used.
- **Arm C has 981 runs.** One run (E9_0042) has no arm-C tables; the bulk/paired analyses skip it. E9's coins are
  reproduced exactly.
- **Post hoc (after step 3):** `diag_b`, `diag_a` and `diag_cluster_perm` are diagnostics only and are labelled as
  such in `metrics.json:post_hoc`.
- **SW:** an exhaustive numba Gotoh alignment over 11.58M residues (53k queries, 51 min on the shared node). It is
  "BLAST-like" in scoring and statistics, not BLAST seeding.

## Not done / next
- Measure the human-neighbour analogue directly: absent-target competitions for human paralogs vs their decoys, using
  a human-only sample with known absent isoforms.
- Build a principled control 2 that can be tested before seeing data: a proteome that cannot be physically present, for
  example a reconstructed ancestral or extinct-lineage set, or the A2 certification approach.
