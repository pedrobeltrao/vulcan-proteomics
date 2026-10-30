# Amendment A3: implementation notes (A3.1.8-A3.1.10, A3.7)

Text: `reviews/amendment_A3_draft.md` rev 1 (DESIGN §6.0c, adopted in D22). Config: `config/a3.yaml`. Written
2026-10-28 by the A3 implementation agent (not the step-4b or step-5 agent). Nothing here read a DN2 per-label
quantity, anything under `/tmp/vulcan_a2`, any sealed file or any E7b output. Rule throughout (A2.9.5, A3.7): **where
the text is open, the stricter reading wins.** Each choice is also fixed in the docstring of the module that
implements it (`vulcan/gates/a3_stats.py`, `vulcan/gates/a3_asyn.py`, `scripts/a3/*.py`).

## What exists

| A3 item | Code | Status |
|---|---|---|
| A3.7 config | `config/a3.yaml`, `vulcan/gates/a3_config.py` | Model, offsets, AlphaPeptDeep and Sage hashes recorded and tested. Freeze fields (A3.1.3, A3.1.6, A3.1.7) filled by the A3 freeze agent in the commit after the freeze commit (notes 33-40; `results/2026-10-28_A3_freeze/`) |
| A3.7 offsets table | `config/a3_secondary_offsets.tsv` (+ `.json` sidecar), `scripts/a3/build_offsets.py` | 100,123 dev secondary groups (D0 2 runs, D1 1 run), sha256 `c53f6f0d...`; MS1 tables regenerated with the committed `vulcan.ms1` and bit-identical to E7's |
| A3.2.2 A-syn | `vulcan/gates/a3_asyn.py` (specification), `scripts/a3/step4b_chimeric.py` (`--asyn`) | Pure parts tested incl. the P1 swap test on synthetic T/D lists; GPU family screen not run |
| A3.2.2 step 4b | `scripts/a3/step4b_chimeric.py` | Phases `--ms1 --asyn --search --score --v12 --swap-test --record`. GPU/data paths not run on DN2 (they cannot be before 4b). The V12 exhaustive-scoring path (`exhaustive_v12_scores`) ran on dev (D7 yeast vs human) in `dev_truncation.py` |
| A3.2.3-A3.2.7 step 5 | `scripts/a3/gates_a3.py`, `vulcan/gates/a3_stats.py` | Pure layer tested on synthetic pools and end to end on synthetic labels (`tests/a3/test_a3_gates_e2e.py`). Data layer (`finalize`, `collect`, `collect_ce_s`, `g3_chimeric`) not run |
| A3.0 Q1 guard | `scripts/a3/v1_closure_guard.py`, `scripts/a3/v1_closure_allow.yaml` | Strengthened at the freeze (note 36). `gates_a2.py` closure: PASS (4 changed GPU files allowed by the 3c OT-HCD bit-identity evidence). `bench/t1_driver.py` closure: PASS at the freeze commit, every other changed file allowed by the A3.1.6 dry-run bit-identity records (note 37) |
| A3.1.9 power | `reviews/checks/a3_power.py` -> `a3_power.json` | see `results/2026-10-28_A3_preconditions/NOTES.md` |
| A3.1.8 dev preconditions | `scripts/a3/dev_g4.py`, `dev_armM.py`, `dev_mc_v2sep.py`, `dev_truncation.py` | Run; `results/2026-10-28_A3_preconditions/` |
| Tests | `tests/a3/` | synthetic pools for every V row, A-syn incl. P1 swap, no-per-label-count checks extended to `scripts/a3/`, closure guard, config pins, end-to-end synthetic-label run |

**Before step 5 (A3.6 step 3):** `gates_a3.py --finalize` and the label-free counts must run end to end on the real
DN2 inputs, plus a synthetic-label dry run on the real pool structure. The GPU wiring follows `scripts/a2/gates_a2.py`,
`scripts/e7/gates_v2.py` and `bench/t1/ours_chimeric.py`, but it has not been exercised on DN2.

## Ambiguity resolutions (stricter reading) and fixed degrees of freedom

### Freeze and guards (A3.1.3, Q1)
1. **Clean tree** = no modified tracked file anywhere, and no untracked file under `vulcan/ bench/ scripts/ config/
   tests/` (an untracked module could be imported). Untracked result files elsewhere are allowed.
2. **Closure guard covers the whole `bench/t1_driver.py` closure**, not only "its v1 path": any module the driver
   imports can change v1's T1 numbers. Each changed file is allowed only by the sha256 of its **exact** diff against
   `prereg-v0.4.2`, so any later edit voids the entry. Non-Python inputs are included (kernel source, build script,
   A2 config and tables).
3. **3c evidence (8e887f7)** is the committed statement in `results/2026-10-27_3c_itcid/NOTES.md` ("OT-HCD unchanged:
   D0 (2 runs) and D1 top-5 PSM tables re-searched bit-identical to M3"). There is no standalone artifact. This
   satisfies A3.0's "must cite its OT-HCD ppm bit-identity evidence". A reviewer may ask for an artifact.
4. **Not allowed yet:** `bench/t1_driver.py`, `bench/t1/{comparators,evaluate,ours_chimeric}.py`,
   `vulcan/chimeric/{__init__,compete,search}.py`, `vulcan/rescore/features_v2.py`. They are in the driver closure
   only, not in the `gates_a2.py` closure. Their evidence must come from the A3.1.6 dry run (v1 columns bitwise
   identical with and without the A3 additions). So the guard currently FAILS for the T1 path. That is intended.

### A-syn (A3.2.2)
5. **Offsets.** Δm/z = secondary mono m/z − the store's instrument precursor m/z (what "instrument precursor m/z"
   is on dev). Rows are sorted by (Δm/z, charge, dataset, run, spectrum, group), and that order is the empirical CDF.
   Δ is written with 6 decimals and read back from the file.
6. **Base m/z** for an A-syn hypothesis = the A spectrum's precursor m/z in the mixture, i.e. the certified mono m/z
   (A2.2(e)). The dev Δ is relative to an instrument precursor that is not always monoisotopic. The draw therefore
   reproduces offsets relative to the reported precursor. This is stated as part of the A-syn limitation.
7. **Draw.** u = int(sha256("{key}|{rank}|{seed}|{try}")[:16 hex]) / 2^64, row = floor(u·N), key =
   "<source run>|<store spectrum index>" (independent of the mixture), seed 47, tries 0..31.
8. **Envelope exclusion: both readings applied.** A draw is excluded if (i) any searched neutral mass (isotope errors
   0, 1 of the draw) is within 20 ppm of M + kN (k = −1..3), **or** (ii) its m/z is within 20 ppm of an envelope peak
   (M + kN)/zz + H⁺ at any zz = 1..6.
9. **W′ screens.** Isotope set {0, 1, 2} and ±10 ppm, the A2 screen convention, which is wider than the searched
   {0, 1}. (b) uses `a2_identity.max_f_in_window` unchanged. The Hamming clause does not depend on the window and
   every A spectrum already passed it, so it is not recomputed. (d) is the max exhaustive hyperscore over (T∪D) ∩ W′
   of the A spectrum's peaks ≥ hyperscore(p); ties are excluded.
10. **Certified p** (sequence, mask, hyperscore) is not stored by A2. Step 4b recomputes it with the unchanged A2
    own-species search (as `admit.py --swap-test` does) and stops unless every A spectrum's certified precursor equals
    the committed admission npz exactly.
11. **A-syn search.** Separate virtual spectra, isotope errors 0..1 (E7 `sec_isotope_errors`), report 5, group = rank
    − 1. A-syn units never enter the mixture's competition. Their features come from `run_v2` on the concatenation
    (real + A-syn). The real rows must be bit-identical to the real-only run (asserted), so A-syn never enters
    calibration.

### Step 4b search (A3.2.2)
12. **Primary group** = the committed step-4 PSM table (A2 hypotheses), reused as `primary`. A fresh search of the
    first mixture must reproduce it exactly.
13. **A1 offset copies** get their host spectrum's secondary groups shifted by the copy's mass offset. The offset is
    recovered as (mz_mix − mz_host)·z.
14. **No MS1 hills:** that source's spectra keep the closed search (no secondary groups), as E7. They are counted.

### Gate rows (A3.2.3, A3.2.4)
15. **τ_s** = min score of an accepted **target** pair winner (q_sep ≤ 1%) in merged stratum s, per pipeline and
    mixture (A1's `tail_threshold` convention). A stratum with no accepted target has τ = +∞, so no unit of that
    stratum is in the tail.
16. **V6 unit** = stage-1 winners of A∪B null spectra with **≥ 2 stage-1 group winners**. A group without any PSM
    cannot be filtered, so counting it would only dilute the drop rate. The top-decile edge is computed within the
    class on both labels pooled (label-free). A class is not evaluable unless **both** bins have n_T, n_D ≥ n_min.
17. **V6 n_min = 22,700** per label per bin. This is the smallest n (multiple of 100) at which a symmetric filter
    at the dev drop rate 0.118 passes one TOST bin with probability ≥ 0.90 (20,000 MC replicates per n, step-up
    after bisection so the target holds at n_min). The text's "about 23,000" agrees.
18. **V12 kept** = the exhaustive best candidate's peptidoform row is among the production top-5 rows of the spectrum,
    so the whole K1 top-K + top-5 chain is tested. Candidates are those with ≥ min_matched_peaks, as production.
    The arithmetic is `a2_stats.g3_gate_a2` with best_*_k = best_* if kept else −∞, i.e. identical A2.8 bins,
    merging and TOST.
19. **V12 n_min** (value in `config/a3.yaml`, = `reviews/checks/a3_power.json`) is taken at the dev miss rate of the
    **max-variance decile** (0.0981). The pooled rate (0.0889) would give a smaller n_min, so the decile rate is the
    stricter choice. Dev source (`truncation_dev.json`): a DN2-like mixture of the D0 host run b1906 and 12,000 D7
    yeast LM2480 spectra (single hypothesis, RT quantile-mapped), human prod index.
19b. **V12 calibration (review M1).** The exhaustive candidates are scored with the v2 calibration of the MIXTURE's
    production primary PSMs and the production NCE. They are appended as non-calibration rows (group 1), and the
    production rows' features are asserted bit-identical. The first implementation calibrated on the A candidates
    themselves; it was replaced before any use, and the dev rate was re-measured.
19c. **V12 ties (review minor 6):** kept = some production top-5 candidate of that label has a score equal to the
    exhaustive best (A2.8 score-equality convention). `dev_truncation.py` uses row identity; ties among
    exhaustive candidates are negligible on dev.
20. **V8:** n = null query units of the pool (G2-original's own unit). Not evaluable below `g2_sample_size()` = 107,570.
21. **V3:** recomputed G3(A) with **v1's** (A2) membership (V3 is "v1's G3(A)"). Byte-identity is tested as the
    sorted-key JSON of the gate dict versus `--v1-results`' `g3.A`. A mismatch gives "not-evaluable
    (implementation)".
22. **V4/V5** use v2_sep's (c-final) membership. Secondary units are the production virtual spectra (source-run MS1
    groups, isotope 0..1). The exhaustive search uses exactly those hypotheses.
23. **V10 (FS):** fail if any evaluable FS gate fails (counted as a tier-1 **fail**); not evaluable if all three are
    not evaluable (never blocks); pass otherwise.
24. **Tier outcome:** fail takes precedence over not evaluable. Tier 2 is only assessed if tier 1 passes.
25. **(c-final) for v2_sep:** the union over learned, reference and noRT of every unit's stage-1 winner **and** its
    `l1.partner`, real and A-syn units, against the project's A2 A1.3 reference. The training-overlap guard drops
    overlapping projects from both the v2_sep and the ce_s pools (label-free; the dropped list is recorded).
26. **ce_s:** A2's closed analysis (per-mixture paired double competition, threshold `tail_threshold`) with the ce_s
    logit of the step-4 top-5 tables. G2-original and G5 are run on A∪B **and** on A (A3.2.7). V12 uses the
    ce_s-best candidate. G3 = v1's G3(A) and G3(A∪B). The verdict uses the same tier rule; it licenses nothing.
27. **G1 logits** of A-syn and chimeric winners use the candidate's AlphaPeptDeep `rt_pred`, which step 4b stores, as
    the frozen classifier needs.

### Dev preconditions (A3.1.8) and power (A3.1.9)
28. **G4 partition** (frozen and committed in 8bfb431 before estimation): the 40 v2 training runs, stage-1 **unit**
    decoys, OT-HCD × HCD × chimeric stratum {1, 2, ≥3} × charge {2, 3, 4+}, with charge 1 kept as its own class.
    The statistic is `gates.stats.g4_gate` unchanged.
29. **Arm M:** A1's D7 design (`scripts/m3e/armM.py`, same indexes), with the v2_sep pipeline and the combined
    estimator, r = 1. Lists: overall and raw strata 1, 2, ≥3 at 1% and 5%. The decision uses CI convention (a)
    (A1.6/3f); (b) is reported.
30. **MC (§6.7):** the committed `chimeric_competition` and `picked_group_fdr` are called unchanged. There are 4
    scenarios × 10⁴ replicates, and the gated lists are each merged stratum, the union and pooled picked protein
    groups at 1% and 5%. Raw strata under a merged q are **reported only**, because the q_sep list that is controlled
    is the merged block. A3.3.6 also gates merged strata.
31. **Power projections** come from dev structure only (constants in the script docstring). B/A and the cap-scale
    |A| are assumptions, shown on a grid. Joint tier probabilities are products of separately simulated marginals
    (independence assumed and stated).
32. **Code state of the dev runs.** G4 (8bfb431), arm M (f258a9d) and MC (f293526) ran from clean worktrees at
    commits **before** E7b landed its speed paths (bd2cfb2: `vulcan/ms1/*`, `vulcan/rescore/fast.py`,
    `vulcan/chimeric/fast_compete.py`). The truncation run (31af0be) includes bd2cfb2, but none of the modules it
    uses was changed there. If the freeze adopts the E7b paths, A3.1.3 requires E7b's bitwise equivalence evidence.
    If that evidence holds, these dev results carry over unchanged.

### Freeze (A3.1.3, A3.1.4, A3.1.6, A3.1.7; A3 freeze agent, 2026-10-28)
33. **Freeze commit vs config commit.** `config/a3.yaml` cannot name the commit that contains it, so it is committed
    after the freeze commit. `a3_config.freeze_status` therefore accepts HEAD = the freeze commit **or a descendant** whose
    frozen paths have the recorded `git ls-tree` hashes (recorded hashes must also equal the freeze commit's). Any edit of
    a frozen path after the freeze changes HEAD's tree hash and the scripts refuse. The clean-tree rule (note 1) is
    unchanged. Frozen paths: `vulcan/`, `bench/`, `scripts/a3/`, `scripts/e7/`, and `scripts/e7b/` (added because
    `bench/t1/a3.py` imports `scripts/e7b/sage_rank_strata.py`).
34. **T1 driver, A3 arm** (`bench/t1/a3.py`, `bench/t1_driver.py --a3`; forced for T1, A3.3.8). Methods = A1 set +
    `vulcan_v2_sep` + `sage_chim2` + `sage_chim2+pred+mokapot` (C2). **R = 5 is not included**: E7b's pre-registered rule
    ran R = 5 beyond D0/D1 only if it beat R = 2 at matched FDP 1% on both, and it lost on D1 (6,238.3 vs 6,238.9; E7b
    NOTES). No wide-window variant was pre-registered (E7b config `wide_window` false; D25), so the claim wording is
    "Sage chimera at narrow precursor tolerance". Numerators `vulcan,vulcan_v2_sep`: `noninferiority` stays v1's;
    `noninferiority_by_numerator.vulcan_v2_sep` is family S's input. The sealed run refuses unless
    `require_freeze_commit()` passes and the model hashes match. Family S's fixed-sequence decision is not coded in the
    driver; `a3_secondary` supplies its inputs.
35. **Prediction-cache isolation.** The A3 arm (v2_sep chimeric scoring and C2's predicted features) predicts through its
    own cache `work/predcache_a3`. The first dry run shared `work/predcache`: the chimeric stage then predicted some
    Sage+pred peptides first, in different batches, and AlphaPeptDeep is not batch-composition invariant in the last
    bits (426 of 296,873 RTs and 8 MS2 vectors differed on D1), which moved Sage+pred+mokapot (an A1 comparator) by up to
    63 IDs at 0.1%. With the isolated cache the A1 comparators' prediction features and every mokapot output are
    byte-identical to the A1-only (prereg-v0.4.2) driver on the same Sage/Comet searches. v2_sep and C2 still share one
    cache (E7b fairness rule).
36. **Closure guard strengthened** (stricter only). Modules run as `python -m <module>` from a list literal are part of the
    closure (`vulcan.cohort` from `bench/t1/ours.py`, which ingests and searches for v1; `vulcan.ms1`), and `scripts/e7`,
    `scripts/e7b` are import roots. This adds `vulcan/cohort.py` (changed by E7b since the tag) and what it imports to the
    T1 closure; the previous guard did not see them.
37. **Allow-list entries at the freeze** are `kind: bit_identity` with evidence `results/2026-10-28_A3_freeze/
    dryrun_{D0,D1}.json`: the whole driver at prereg-v0.4.2 (A1 set) and at the freeze code (A3 set) give byte-identical
    v1 artefacts (every `ours/<rep>/<run>.{psms,scored}.parquet`), identical v1 metrics on every common key, and
    byte-identical A1 comparator artefacts from the same Sage/Comet searches. This is end-to-end evidence covering all
    changed files jointly; it is not per-hunk.
38. **Merged strata in `a3_secondary`.** For chimeric rank k, the merged stratum is the block of that replicate's
    label-blind `merge_adjacent` map that contains k (`chimeric_competition` `merge`/`merged`), so "merged stratum k" may
    differ between replicates only if the maps differ (they are reported per production replicate). Counts restricted to
    the block, paired target from the full list (`scripts/e7/strata.py`, unchanged). Raw strata are reported only.
39. **LOO backstop (A3.3.4).** For each run r: the point ratio with run r removed from both methods' run x decoy-seed
    curves, the same estimator as `evaluate.noninferiority` (own q, or the first crossing of the arm-P FDP of the
    remaining runs). The all-runs value is asserted equal to `noninferiority`'s point ratio. With one run (D1) it is None.
40. **Prediction cache for T1.** No pre-existing cache is read by the T1 look: `work/predcache` and `work/predcache_a3`
    are built inside the invocation from the pinned AlphaPeptDeep weights (sha256 in config). `prediction_cache_manifest`
    pins the caches behind the freeze evidence (dry runs, G0).

## Independent review (A3.7; fresh-context agent, 2026-10-28)

The reviewer read A3 rev 1, its review, DESIGN §6, `config/a2.yaml`, `gates_a2.py` and `a2_stats.py`, and reviewed
`gates_a3.py`, `a3_stats.py`, `a3_power.py` (+ json), `config/a3.yaml` and these notes. It ran `tests/a3` (53/53) and
modified nothing.

**Verdict: no BLOCKING findings, 4 MAJOR, 8 MINOR.** Rows V1-V14, τ_s, tail membership, merge maps, V3 byte identity,
V6/V12 statistics, tiers and ce_s were checked against the text and found correct.

| Finding | Outcome |
|---|---|
| M1 V12 scored with a calibration fitted on the A candidates, not production's | **Fixed** (note 19b; dev rate re-measured on a D0-host mixture: 0.0981 max-variance decile, was 0.0954; V12 n_min recomputed) |
| M2 ce_s V12 used v2_sep's membership | **Fixed**: ce_s V12 uses ce_s's own (c-final) membership |
| M3 power flags took a loose reading; 1.2M at 16 projects unreachable under A2.9.3's stopping rule; "can never pass" tested as evaluability | **Fixed**: the |A| grid is bounded by the stopping rule (J and |A| per per-run density); flags use pass probability < 0.05 at every reachable point, per B/|A| and overall; NOTES and json agree |
| M4 joint tier probability was a product of marginals, not simulated | **Partly fixed, disclosed**: the json reports bounds valid under any dependence (Fréchet lower, min upper) next to the independence product. A single synthetic DN2 per replicate was not built. This is a departure from "simulated joint", recorded here and in the NOTES. The flags use the dependence-free upper bound |
| minor 1 V13 tail fraction not from dev | Sensitivity at half the tail fraction is reported; the assumption is stated |
| minor 2 power gaps (FS size, C_CAND, MC reps, B/A) | |FS| = |A| is by construction (stated); V12 candidate share from dev stated; tail/G1 MC reps raised to 600; B/A is a stated grid; NOTES wording on margins corrected |
| minor 3 `--gates` without `--g3` | **Fixed**: refused |
| minor 4 finalize did not verify A2's membership | **Fixed**: `--a2-record` required, committed, and its sha256 must match the npz |
| minor 5 reported-only items missing; noRT reused v2's V6/V12 | **Fixed**: V13 per project, V6 class 2 by rank, V8 by stratum, B G1 secondary pooled, A1 offsets; V6 computed per learned pipeline; noRT V12 marked "not computed" |
| minor 6 V12 tie rule | **Fixed** (note 19c) |
| minor 7 training-overlap guard could be a silent no-op | **Fixed**: hard failure if the accession column is missing |
| minor 8 doc names; tol/conf in config unused; V11 inconsistency in the text | Names fixed; a test pins config tol/conf to the statistics constants; V11 stays "unused" (A3.2.4) and is reported as such |
