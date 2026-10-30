# Statistical review of DESIGN.md v0.1 (FDR, entrapment, learned scoring)

Statistician review, 2026-10-26. Sources: DESIGN.md and lit/*; other reviews were not read.

## 0. Verdict

The skeleton is right: peptide/precursor double competition with stored pairs, paired plus combined entrapment, the lower bound used only to detect failure, a static model kept separate from RESET, Group-walk, picked protein groups, decoy-AA FLR, and XL merge-then-FDR.

**Formulas checked against primary sources:**
- Correct as written: paired (Wen 2025, r = 1), combined, lower bound, (D+1)/T and (TD−DD)/TT.
- **Missing:** the RESET constant and the decoy-AA FLR estimator.

**Where it falls short:** a feature leaks the label; the validation rule cannot see a learned sequence prior; underpowered gates pass by default; the same entrapment data gate merges and support claims.

---

## 1. Blockers

**B1. A label-dependent feature.**
- *Problem.* Stage-A "Δscore … to the best decoy" reads the decoy flag. Any score that uses it breaks equal chance by construction (the Winkelhardt 2026 near-perfect-separation failure).
- *Fix.*
  - Delete the feature. Δscore is computed only to the next candidate, whatever its label.
  - Enforce this in the structure: feature and scoring code never receives `is_decoy` or the pair side, which live in L4 only.
  - Add a unit test: swapping the labels of every pair leaves all features bit-identical.

**B2. The entrapment rule cannot see learned-prior leaks.**
- *Problem.* §5 accepts "paired (or combined)".
  - A learned peptide prior treats shuffled entrapment like decoys, because both look unnatural. False targets beat both, so the paired FDP̂ stays near α while the true FDP rises.
  - The *Arabidopsis* arm tests naturalness, not **memorization** of human peptides seen in PRIDE training.
- *Fix.* For any scorer with learned or predicted components (mode a, mode b with Stage-B, DIA tier 2), "supported" needs **all three** of:
  1. paired shuffled entrapment;
  2. combined natural foreign-proteome entrapment;
  3. **memorization entrapment**: a high-resolution sample from an organism absent from training, with *human* peptides as entrapment, subsampled (seeded) to r ≈ 1 and with I/L-identical peptides removed.
- S1 (Pfu vs Pfu+human) has the right structure, but it is one 2012 Velos run with r ≈ 15–30. Add an r ≈ 1 arm and a modern sample.

**B3. The §5 decision rule is statistically ill-posed.**
- *Problems.*
  - FDR is an expectation and FDP̂ is one noisy draw (decoy choice alone moves accepted PSMs by 4.7% at 1%; Keich 2019). "≤ α everywhere" fails correct methods by chance, and a single-list lower bound > α flags violations by chance.
  - Replicate aggregation is undefined. §7 adds a binomial tolerance that §5 lacks, and the binomial variance is wrong for the paired estimator.
  - Pseudo-reverse is deterministic, so the "3 decoy seeds" only change the replacement shuffles (<1% of peptides). That replication is illusory.
  - At α = 0.1% on S1 the expected false count is about 1–5, too few to evaluate.
- *Fix: pre-register this rule now, and use the same rule for the §7 gate.*
  - Statistic: the mean FDP̂ over replicates, with a 95% CI from bootstrapping target–entrapment pairs plus between-seed variance.

    | Verdict | Condition |
    |---|---|
    | Supported | mean paired/combined FDP̂ ≤ α |
    | Supported-noisy | mean ≤ α but the upper CI > α |
    | Violated | the lower CI of the mean lower bound > α |
    | Not evaluable | α·N_T < 20 |

  - Either add a real decoy-variability arm (C-term-fixed shuffles, 5 seeds) or state that only entrapment seeds are replicated.

**B4. Entrapment is used both for development and for claims.**
- *Problem.*
  - Tier S, including S1/S0 entrapment, gates every merge. After hundreds of merges, the surviving configuration has been *selected* to pass those sets.
  - The "≤1% peptide loss" gate sits inside decoy noise. It favours variants that carry more undetected false IDs.
- *Fix.*
  - Split the sets: *dev* (S tier, may gate merges) and *sealed claim sets*. List the sealed sets in DATA.md now, take them from other labs and instruments, and give them fresh entrapment seeds (2004+).
  - Evaluate sealed sets only at milestones. Log every look in a git ledger and report the count.
  - Gate IDs with a paired comparison over decoy shuffles and a margin larger than the seed SD.

**B5. The learned-score gates cannot work as written.**
- **G1 tests the decoys, not the scorer.**
  - AUC ≤ 0.6 is a property of the decoy recipe. Reversed decoys reach ≈ 0.64 (Reznikov 2026), so the default recipe fails by construction.
  - Any rich peptide embedding separates reversed sequences, so "no separation in the peptide branch" also fails.
  - Fix: keep the AUC as a diagnostic only. Test the scorer through G2, and require that null-spectrum scores do not correlate with the sequence-classifier logit in the tail.
- **G2 is an absence-of-evidence test.**
  - "0.5 inside the 99% CI" passes whenever n is small. It is also marginal, while priors act in the tail.
  - If null winners are targets with probability p, TDC underestimates FDR by p/(1−p): **p = 0.52 gives +8%, p = 0.55 gives +22%**.
  - Fix: use an equivalence test. The 99% CI of p must lie in [0.49, 0.51], which needs about 16.6k null competitions per stratum. Apply it overall **and** in the top score decile.
  - Also compare the score distributions of null target winners and null decoy winners.
  - Run G2 through the full production pipeline (retrieval, pruning, chimeric search, clustering) at the peptide unit, with a wrong-species-spectra-vs-human-DB arm.
- **Missing gates:**
  - near-perfect separation as a leak signal (lit R3c);
  - feature label-blindness (B1);
  - a consequence for failing the calibration check (R5).

---

## 2. Majors

**M1. Group-walk applicability.**
- *Problem.* Group-walk needs group membership that does not depend on labels, with fair-coin nulls inside each group.

  | Stratum | Problem |
  |---|---|
  | Chimeric rank | Computed among candidates that include decoys. |
  | Mass-offset group | The best target and best decoy of a peptide come from different spectra with different offsets. |
  | Charge | Undefined for a peptide (the DDA unit). |
  | XL class | (TD−DD)/TT is not a fair-coin competition. |

- *Fix.* Every stratum must be invariant within a pair:
  - charge at precursor (pep, z);
  - (stem, offset-bin) for mass offsets (CONGA's dynamic level);
  - length, PTM composition and specificity, which pseudo-reverse preserves;
  - chimeric stratum from the label-free rank of the precursor feature (MS1 intensity in the window).
- Pooled TDC remains the valid fallback. These invariants also keep each target and its entrapment in one group, as the paired estimator requires.

**M2. Chimeric search.**
- *Problem.* The competition unit is undefined. Secondary IDs are dominated by *neighbours*, targets that share fragments with the true peptide. Decoys rarely do this, so the competition tilts toward targets (Lin 2022 caveat).
- *Fix.*
  - Compete per (spectrum, precursor feature).
  - Apply CONGA's label-blind neighbour filter: drop any PSM, target or decoy, that shares >5% of b/y ions with a higher-scoring PSM of the same spectrum.
  - Run E7 as rank-stratified entrapment under B3.

**M3. Open search.**
- *Problem.* There is no narrow-plus-open structure, so unmodified peptides lose to random wide-window winners.
- *Fix.*
  - Adopt CONGA in full: narrow and open groups, neighbour filtering, stem-level competition and Group-walk.
  - Derive the offsets from a disjoint pilot set or a label-blind histogram.

**M4. Clustering and medoid propagation.**
- *Problem.* The medoid's top target was selected *because* it is a real peptide like the cluster; its decoy was not. Mis-clustered members therefore see targets winning with probability above ½.
- *Fix.*
  - Give each member its own candidate generation, or a pair-closed set (see M5) plus a full search below a similarity threshold.
  - Validate on vs off on identical data, using entrapment and G2.

**M5. Pruning, ANN, lazy prediction, routing and ties.**
- **Pair closure.** In every shortlist (top-K, XL α-first, ANN, lazy-prediction survival, isoform top-K), a scored target's paired decoy is scored with the same features, and vice versa.
- **Label-blind routing.** The bi-encoder side channel ("no good stage-1 hit") must route on the best score, not on q or target status. Otherwise treat it as cascade search with per-stage FDR. This channel is mostly nulls, so it has high leverage on FDP.
- **G3 on nulls.** Recall on correct targets says nothing about nulls. Measure ANN-vs-exhaustive recall of the best target *and* the best decoy on null spectra, binned by score.
- **Ties.** Use seeded random tie-breaking in competition, and count decoys first at equal scores in q-values.
- **Re-gating.** Re-run G2, G3 and dev entrapment whenever these parameters change.

**M6. Gaps in entrapment construction.**
- **Entrapment decoys.** Entrapment peptides need their own paired decoys: DB = T ∪ E ∪ dec(T) ∪ dec(E). L1 does not specify this. Without it, every q-value is about 2× liberal.
- **Protein level is missing.** Wen's protein entrapment uses paired shuffled proteins with shared peptides swapped consistently, and counts a group with ≥1 target as a target. Protein level is where DIA tools fail, and E9 needs it.
- **Peptidoform level** is undefined, so phospho is not covered.
- **r.** Compute r on the effective search space (unique, I/L-collapsed peptides after the length and mass filters) and report it per mass decile.
- **Near-identical peptides.** Remove single-substitution and isobaric (N↔GG, Q↔AG) peptides from foreign proteomes. Conserved families otherwise inflate N_E, which is conservative but costs power.
- **Seeds.** Each entrapment must be unique against the decoys of all decoy seeds.
- **Reporting.** FDP comes from the augmented search, sensitivity from the production search.

**M7. RESET is under-specified, with 2× bias at stake.**
- **Write into §5:**
  - FDR̂ = **2**(D_est+1)/T for single-decoy (factor 1 for two-decoy);
  - training constants 3 and 2 respectively;
  - the split applies to winning decoys after peptide competition;
  - estimating decoys enter training as unlabeled pseudo-targets.
- **No CV folds.** RESET has none; "peptide-grouped folds" comes from mokapot. If folds are used, group them by stem across the cohort.
- **Two-decoy mode** needs exchangeable decoys: two shuffles, not reverse plus shuffle.
- **Unproven combinations.** RESET + Group-walk and RESET + picked protein are unproven, because remaining nulls are targets with probability ⅔. Either adapt them (c = ⅔; decoy counts ×2) and verify by simulation and entrapment, or restrict mode (b) to pooled peptide FDR.

**M8. Training of the frozen scorer, and leakage.**
- **Negatives are unspecified.** Decoy negatives teach the model the decoy recipe.
  - Use incorrect *targets* as negatives: hard, wrong-species and shifted-precursor matches.
  - Never use decoy-recipe sequences as labeled negatives.
  - Consider a gradient-reversal head on the peptide encoder.
- **Accession-disjoint ≠ data-disjoint.** PRIDE contains re-deposits and MassIVE-KB/ProteomeTools overlaps. Deduplicate by raw sha256 and spectrum hash.
- **Third-party predictors.** Record the training accessions of the predictors in DATA.md.
  - Forbid per-run fine-tuning on confident targets in mode (a); in mode (b) it falls under RESET.
  - Allow only monotone RT alignment.
- **L6 fold assignment.** L6 must say which fold model scores each reprocessed project.
- **Peptide exclusion** is circular and creates an anti-prior. Prefer organisms held out at training.

**M9. PTM FLR.**
- *Problem.* A neural localization head learns residue priors (pS ≫ pT ≫ pY, never pA) and S/T-specific neutral losses. It then rarely places pA, which makes the decoy-AA FLR̂ anticonservative. This is why Ramsbottom found the estimate depends on which decoy residue is used.
- *Fix.*
  - Make localization scoring blind to residue identity given the masses: the same ions and losses for every acceptor, Ala included.
  - Make E8 (FLR̂ vs true FLR) a gate before any FLR claim.
  - Pre-register the estimator, e.g. FLR̂ = n_pA·(N_STY/N_A)/n_pSTY over candidate acceptors. Ala is primary, Leu secondary, with I counted as L.
  - Include Ala isoforms in the top-K isoforms.
  - Label "analogous choices for other PTMs" exploratory.
- **Dataset roles:**
  - PXD000138: dev/calibration;
  - PXD007058: sealed test;
  - PXD009449: exploratory.
- **Other checks.** No pool may contain positional isomers of one sequence. Identification FDR uses the localization-agnostic peptidoform.

**M10. Picked protein group.**
- The group score is unspecified. The 2022 used the best peptide PEP and showed PEP products are anticonservative. Use the best peptide score, or a PEP fitted on held-out decoys.
- Log the replacement shuffles and removed decoys that break the peptide→protein mirror.
- Rescued subsets have only empirical support. Validate with protein-level entrapment at 1, 100 and 1000 runs.

**M11. "No loss" has no statistical definition.**
- Pre-register non-inferiority on each sealed set. For example, IDs ≥ 0.99× the comparator, with a bootstrap CI over runs and decoy shuffles.
- Compare only methods whose control is "supported", and also at entrapment-matched FDP̂ = α.
- Apply the same to FLR (recall at matched true FLR), quant (ε and CV with CIs) and percolator-free (a ≥ b).

**M12. Repository scale.**
- **Calibration.** "Globally calibrated scoring" must not fit transforms on the *estimation* decoys; use a separate calibration-decoy seed. Define what happens when the decoy-quantile check fails.
- **Claims.** Say which claim "run-level ∩ global q" supports: the master list carries a global FDR and per-run lists carry run FDR.
- **Scale test.**
  - Define "FDP flat in N" as the B3 rule at every N plus a slope test.
  - Use heterogeneous runs.
  - Pre-register reduced replication at N = 1000, because 9 replicates are infeasible on 4 cores.

**M13. XL-MS.**
- (TD−DD)/TT has no finite-sample guarantee and is unstable at small counts. Report it with the conservative TD/TT, clip it at 0, and set minimum counts.
- Define self vs hetero for pairs containing decoys (A with decoy-A counts as self).
- Make α-first pruning pair-closed.
- Tune xiFDR-style "boost" only on dev data.
- Remove XL class from Group-walk (M1).

**M14. Missing hard regimes.**
- Wen's failures concentrate in single-cell DIA and at protein level. Tiers S and M lack both a single-cell DIA entrapment set (PXD023325/PXD049412) and an immunopeptidomics set.
- CF-MS has no interaction-FDR model. Use a gold standard held out at the complex level, or a shuffled-profile null.
- Validate MBR with PIP-ECHO and a two-proteome absence design (PXD014415).

---

## 3. Minors
- **§1 bar table.** Wen tested DIA-NN **1.8.1**, not 2.x, and showed its failure with the *lower bound*, not paired entrapment. The Casanovo-DB and DDA+ gains were not entrapment-verified. Label them as unverified.
- **PSM counts.** Counts at q ≤ 1% are not FDR-controlled under multiplicity. Label them as such.
- **Rigorous add-ons.** Add TDC-UB or FDP-SD bounds, and multi-decoy competition (Emery 2020, d = 2–4) for α ≤ 0.1%.
- **Monte-Carlo tests.** Simulate every estimator, including Group-walk, RESET and picked group: require mean FDP ≤ α over ≥10⁴ replicates. Add a biased (p = 0.55) simulation to show the gates fire.
- **Enforcement.** Freeze §5 in a git tag.
- **"Next seed".** Derive replacement shuffles as hash(peptide, seed) so they do not depend on DB order.
- **G2 offsets.** +7…50 Da includes +14, +16, +28, +42 and +44, where near-neighbour targets match real fragments. Exclude them or analyse them separately.
- **Quant.** Train the variance model on data disjoint from the ProteoBench sets.
- **PEP.** Define "PEP (where valid)" or drop it.

---

## 4. Keep
- **Statistical core:**
  - paired double competition with pair ids carried through every table;
  - correct formulas, with the lower bound kept as a failure detector only;
  - static model vs RESET, with "percolator-free" claimed only if (a) ≥ (b);
  - Group-walk rather than per-stratum TDC;
  - picked groups without razor peptides;
  - the global analyte constraint and q×N reporting.
- **Assays:**
  - XL merge-before-FDR with self/hetero separation;
  - decoy-AA FLR plus synthetic ground truth;
  - MBR transfers with their own q-value, never counted as IDs.
- **Reproducibility and symmetry:**
  - a content-addressed index whose hash includes the decoy recipe and seed;
  - fixed seeds and a frozen held-out list;
  - the G3 symmetry idea;
  - the Pfu+human structure, the right basis for the B2 memorization test.
