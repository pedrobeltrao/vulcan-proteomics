# Independent statistical review of Amendment A5 (draft revision 0; to become DESIGN §6.0d)

- **Reviewer:** independent statistician (agent), fresh context. I had no role in A1–A5, their reviews, the D28 review, the T1 driver, E7/E7b, the C9 pilot or any gate run.
- **Date:** 2026-10-28.
- **Draft reviewed:** `reviews/amendment_A5_draft.md` (rev 0, untracked at review time).
- **Read:**
  - mandate: RESEARCH_TREE D9–D28; `reviews/pi_proposal_D28_review.md` in full;
  - DESIGN §1, §6 (§6.0 A1, §6.0b, §6.0c, §6.1–§6.8);
  - `reviews/amendment_A2_draft.md` (P1–P5, A2.9.5–A2.11);
  - `reviews/amendment_A3_draft.md` (A3.0, A3.1, A3.2.5–A3.3.8, A3.5–A3.7);
  - `reviews/amendment_A4_draft.md` (A4.4–A4.6, A4.7.4, A4.8, A4.9);
  - code: `bench/sealed_run.py`, `bench/t1_driver.py`, `bench/t1/{a3,evaluate,ours}.py` (relevant parts), `vulcan/gates/a3_config.py` (freeze check), `vulcan/index/foreign.py` (header), `scripts/m3e/armM_prepare.py`, `scripts/m3e/armM.py`, `scripts/a3/dev_armM.py` (header and lists), `scripts/e9/stage_pred_env.sh`; env-var reads under `bench/`, `vulcan/`, `scripts/`;
  - `config/a3.yaml`; `.gitignore`;
  - results NOTES: `2026-10-27_T1driver_dryrun`, `2026-10-28_E7b_fair_speed`, `2026-10-28_A3_preconditions`, `2026-10-27_M3e_hardening` (gate table, to check A5.4's numbers), `2026-10-28_E7_chimeric_v2` (G0/G2/G5 lines);
  - metadata rows of `data/manifests/dev_sealed.tsv` and `data/manifests/train_m7.tsv` (file names, roles, species columns only).
- **Not read:** anything under `data/sealed/` or `data/sealed_derived/`; any `results/*_T1_sealed*` or `*_T5_sealed*` directory (none exists); T5's SDRF, `Sample_info.xlsx` or any truth file; any C9 pilot output (`results/2026-10-28_C9_pilot/`, `scripts/c9/`); any DN2 per-label quantity (none exists). I did not run `bench/sealed_run.py` or any driver.
- **T1 exposure: none. T5 exposure: none. C9 pilot exposure: none.**

---

## Verdict: **adopt with changes.** Do not adopt as written.

A5 does the hard parts well:
- it is written as an amendment, not as a timing note;
- the required first paragraph is verbatim;
- most clauses have quoted old and new text, and the quotes I checked are accurate;
- "entrapment-validated" and "contaminated" are gone, and every short form contains "FAILED";
- family S, the designations and arm M are closed;
- the S6 floor is mostly right;
- A5.4's gate numbers match their sources.

It is **not yet safe to execute**, for two reasons.

**First, the one-look machinery has two holes that cannot be repaired after the look (B1, B2):**
- a look can "succeed" with missing outputs;
- an interrupted look has no defined status and leaves plaintext readable.

**Second, several rules are still open or ill-defined at look time:**
- the claim-level arm-C void, and arm C at protein level, which does not exist (M1);
- the §1 conjunction (M2);
- the environment (M3);
- the GO branch (M6).

There are also regressions from A4.8's mechanics (M4), an incomplete T5 design (M5), a weakened S6 (M7) and missed clauses (M9).

Counts: **2 blocking, 10 major, 12 minor.**

---

## Blocking

### B1. A "successful" look can silently lack comparators, replicates or numerators, and family S's IUT then runs over fewer comparators
- **Problem:**
  - `evaluate.evaluate` catches missing outputs. The loader returns `None` on `FileNotFoundError`/`KeyError` (`bench/t1_driver.py` l. 389–391).
  - The method is then recorded only as `{"missing_replicates": [...]}` and skipped (`bench/t1/evaluate.py` l. 461–465).
  - `noninferiority` skips any method not in `curves`, and the driver exits 0.
  - So a mokapot, Comet or Sage-chimera output that is missing for one replicate removes that comparator from every ratio. "Against **every** comparator" (A3.3.3, A5.2 F1, P1) then becomes "against the ones that ran", which makes the test easier to pass.
  - A5.5.5 defines failure only as "a non-zero exit".
  - After the look nothing can be re-run (A5.5.3), so this cannot be repaired later.
- **Required fix:**
  1. Add a completeness rule, normative and in code. A look is **failed** unless:
     - every method of `A3_METHODS` has all 25 replicates and no `missing_replicates` key;
     - all 6 T1 runs and both T5 runs are present in every method's run list;
     - `meta.a3` lists every A3.3.3 comparator;
     - `fdp_reached` is reported;
     - `n_boot_nan` is recorded, with a fixed threshold for a not-evaluable ratio.
  2. Enforce it twice: the driver exits non-zero in sealed mode on any `missing output` log event, and `claims_a5.py` independently refuses incomplete inputs.
  3. Add a test with one deleted replicate output.

### B2. Interrupted looks are undefined, and their plaintext stays readable
- **Problem:**
  - On SIGTERM, `sealed_run.py` relocks and raises `SystemExit(143)` from the handler. The code after the `try/finally` never runs, so **no result line is written** (l. 125–135). SIGKILL, node loss and a scheduler time limit leave nothing at all.
  - The look is long: T1 alone was projected at ~4–6 h before A3's chimeric arm and Sage-chimera comparators (dry-run NOTES), and T5 comes first.
  - The default work dir is node-local `/tmp/vulcan_t1_<uid>/T1`, and the driver deletes only `stage`, `mzml` and `ms1` (l. 433–436). So PSM tables, mokapot outputs and partial `metrics.json` stay readable on that node, outside any move.
  - A5.5.5 covers only "non-zero exit". A5.5.5's retry bar ("refuse if any earlier look line has two exit-0 result lines") does not say what a look line with **no** result line is.
- **Required fix:**
  1. A look line without both result lines is a **failed look**. It is consumed and counts towards the retry limit.
  2. `sealed_run.py --look a5` must first check for an unterminated look id. It refuses to proceed until a `--recover <look-id>` step has written a failure line and moved or destroyed the work and out dirs.
  3. The work and out roots must live in one root under a known path, on shared storage or on a node that the recover step can reach, with a test.
  4. A deliberate interruption by the operator is a protocol violation and is reported.
  5. The look must run in a job whose time limit is at least 2× the measured dev `--look` dry run (M4), recorded in `config/a5.yaml`.
  6. Make the one-look check scan the **git history** of `results/SEALED_LEDGER.md` (`git log -p`), not only the current file, so that a rewritten ledger cannot hide an earlier T1 or T5 line.

---

## Major

### M1. The claim-level arm-C void (A5.7.1, A5.7.2) and the arm-C labels are ill-defined; arm C does not exist at protein level
- **Problem:**
  - **Arm C at protein level.** `evaluate.py` builds no protein table for arm-C replicates (`if kind != "C": tabs["protein"] = ...`, l. 443–445). Protein verdicts are arm P only (l. 485–486).
  - So A5.5.3 ("arm P (15) and C (5) … every unit") and A5.11.2.4 ("arm C … every unit") are factually wrong.
  - L-M ("arm C for the same list at the same threshold") cannot be filled for P2-protein, S3 or any protein-group number.
  - A5.7.1(1) silently reduces to arm P at protein level.
  - **The void itself** ("arm C at the same matched threshold, for the same list, has lower CI of LB > f") does not say:
    - (a) that the threshold is the matched **q** value `ours_q_at_match` (the interpolated point estimate from the 15-replicate mean arm-P curve), applied to each of the 5 arm-C replicates' own q;
    - (b) the r used (`arm_C_r` at that unit);
    - (c) the CI convention ((a) decides?);
    - (d) whether it applies to the numerator only or to both sides of the ratio;
    - (e) whether it applies to v2_sep per merged stratum or overall only;
    - (f) whether it applies at own-q thresholds.
  - Different readings give different voids. That is an option open at look time (A3 Q6).
- **Required fix:**
  - Define the void exactly as (a)–(f). Recommended: numerator only; overall list; convention (a); r = `arm_C_r[unit]`; peptidoform and sequence units; matched thresholds only.
  - State that at protein level arm C is "not computed (§6.3 arm C has no protein construction in the frozen driver)". Print this in L-E and L-M.
  - Correct A5.5.3 and A5.11.2.4.
  - Test the void in `claims_a5.py` against synthetic inputs at the boundary.

### M2. P1 drops half of §1's "no loss" test, and is weaker than v2_sep's S1
- **Problem:**
  - §1 computes no-loss "at each method's own q ≤ 1%, **and** at matched entrapment FDP". Own q is dropped only for comparators that are not "supported".
  - A3.3.5's S1 for v2_sep requires matched **and** own q.
  - A5.9.2 splits v1 into P1 (matched, the headline) and P1q (own q, a separate row). The headline (A5.9.4) cites P1 only.
  - So v1 can be headlined as "no loss" with P1q failing. That is an undeclared relaxation of §1, and the review said §1 is not to be touched.
- **Required fix:**
  - Define the §1 verdict "v1 no loss" as P1 ∧ P1q, where P1q is over supported comparators, as S1.
  - The headline cites that conjunction and prints both rows.
  - If A5 intends to relax §1, it must say so in A5.0 and in a quoted replacement of §1. I advise against that.

### M3. "Nothing to choose at look time" is false: the environment and external binaries are not pinned, and the proposed deny-list breaks the driver
- **Problem 1: the deny-list breaks the driver.** A5.11.2.2 denies `VULCAN_T1_*`, but the driver sets `VULCAN_T1_STAGED=1` itself when it re-execs (`staged_python`, l. 92–104). A prefix deny-list therefore refuses every sealed run, or has to be checked before the re-exec.
- **Problem 2: other variables that change results are not covered:**
  - `VULCAN_SAGE` (the Sage binary path, `vulcan/calib/searcher.py:28`);
  - `VULCAN_MOKAPOT_ENV` (`comparators.py:33`, default `/tmp/vulcan_m1e/envs/mokapot`, outside the repo);
  - `VULCAN_CORE_DIR` (the native core);
  - `VULCAN_TRFP_BIN`, `VULCAN_RAW_DLL`, `VULCAN_DOTNET_ROOT`;
  - `VULCAN_LOCAL_ENV` (`stage_pred_env.sh` reuses any stale staged env if `pred.stamp` exists);
  - `VULCAN_COHORT_TIMEOUT_SCALE`, `VULCAN_THREADS`, `VULCAN_T1_NO_LABEL_CACHE`.
- **Problem 3: files the freeze does not cover:**
  - `scripts/e9/` (the env stager) and `scripts/baselines/` (the mokapot env builder);
  - the Comet binary (`conda_envs/comet`), the TRFP version, the staged pred env and the mokapot package env.
  - A5.11.6's tree list omits all of these.
- **Required fix:**
  1. `sealed_run.py --look a5` builds the child environment from an **allow-list** in `config/a5.yaml`, with explicit values (paths and thread counts), and records it in the ledger line. Everything else is dropped.
  2. The driver asserts the allow-listed values after its re-exec.
  3. Add `scripts/e9/` and `scripts/baselines/` to the freeze.
  4. Record sha256 manifests of the Comet, TRFP, mokapot and staged pred envs in `config/a5.yaml`, verified at look time. Rebuild the staged env from scratch for the look.
  5. Test: a set `VULCAN_SAGE` is refused or overridden.

### M4. The look mechanics regress from A4.8 rev 3, which two reviews had already hardened
A4.8.2 put `--out` and `--work` under one root outside the repo. It also had a leak test (no file newer than the run start outside the root, under `/tmp/vulcan_t1_*`, `results/`, `data/cache/`, `work/`), a whitelist test of console lines, an end-to-end D0/D1 dry run of the mode, and a rule that no other agent of the project is active on the node or at that time. A5 keeps none of these.
- **Problems in A5 as written:**
  - The out dir defaults to `results/<date>_T1_sealed/T1` inside the repo, where other agents working in the repo can see it as untracked files while the run is in progress.
  - The T5 output location is unspecified.
  - There is no leak test of files written outside the declared dirs.
  - There is no test that the console shows only templates.
  - A5.11.1.10 tests use dummy files only. Nothing exercises `--look a5` end to end with the real driver and `t5_armM.py`.
  - "Nobody reads any output until both have result lines" is the only protection for the hours between T5 finishing and T1 finishing.
- **Required fix:**
  - Adopt A4.8.2.1 (one root, `--out` and `--work` under it, logs under it), A4.8.2.3 (leak test) and A4.8.1.4's console whitelist, minus encryption.
  - Add a dev-mode `--look` dry run (D7 in place of T5, then D0 or D1 in place of T1) that runs through `sealed_run.py` in a temp clone with dev manifests. Commit it before the freeze.
  - Require that no other agent of this project is active during the look, and log the node's process list at start and end.
  - The T1 out dir and the T5 out dir must not exist before the look.

### M5. Arm M on T5 is not yet a fixed, executable design
- **Problem 1: the dev code is not importable as-is.** `armM.py` hard-codes `RUNS = LM2486/2487/2490`, `SEEDS = 2001–2003` and the work paths. Its `evaluate()` **writes to `results/2026-10-27_M3e_hardening/armM.json`**, a committed dev result. `dev_armM.py` is the same with `$VULCAN_A3_WORK`. A5.11.3.1's "importing the dev arm-M code paths" therefore means a refactor.
- **Problem 2: v1's scoring path is ambiguous.** Dev arm M scores v1 through `scripts/m3/features.py`, while the T1 look scores v1 through `bench/t1/ours.py` `Scorer`. A5 does not say which one is v1 on T5, or require that they are equivalent.
- **Problem 3: the strata differ from dev.** Dev v2_sep arm M gated **raw** strata {1, 2, ≥3} (`dev_armM.py` docstring). A5.6 and the new A3.3.6 bullet use **merged** strata, yet A5.6 claims "exactly the dev design", and A5.11.3.2 demands that the dev tables are reproduced exactly.
- **Problem 4: details are left open:**
  - r is "on the effective search space" (dev hard-codes the peptide-level r = 1.0 for peptidoform units);
  - the dev decoy seed is unstated (1001 in the dev replicates);
  - the CI convention for the T5 void (A5.7.1(2)) is unstated.
- **Problem 5: nothing excludes a post-hoc excuse for a T5 violation.**
  - A T5 "violated" verdict would void v1. If hyperscore is also violated, the obvious post-hoc reading is "the sample contains human material".
  - T5's sample sheet is sealed, and the target list is fixed without it ("no spiked protein").
  - A5 has no clause excluding that excuse, which is exactly the D12 pattern.
- **Problem 6: T5 and T1 are not isolated.** AlphaPeptDeep is not batch-composition invariant (driver comment l. 197–200, 3i/E7b). If `t5_armM.py` shares any prediction cache or other state with the T1 run, running T5 first can change T1's v1 bits relative to the A3.1.6 dry run.
- **Required fix:**
  1. Write `bench/t5_armM.py` as a parametrised module that takes runs, seeds and output dir from its arguments, writes nothing outside its root, and is shared by the dev and sealed modes. The dev scripts call it.
  2. Fix v1's T5 scoring to the `bench/t1/ours.py` path, with a committed D7 equivalence against 3e's table, or explain the difference.
  3. State whether the strata are raw or merged; recommended: gate merged, report raw, and change the dev-equivalence check accordingly.
  4. Fix r = peptide-level effective r (as dev), convention (a) for voids, and dev decoy seed 1001.
  5. Add: "A T5 arm-M verdict is final. No explanation found after the look (sample composition, contamination, reference-pipeline behaviour) un-voids a claim. The hyperscore reference is reported, never used to excuse."
  6. Give T5 and T1 separate fresh caches. The dev `--look` dry run (M4) must show v1 columns bit-identical with and without the preceding T5-analogue step.

### M6. The GO branch is infeasible as ordered, and its trigger and record are not bound
- **Problem 1: A6's code has nowhere to land.** A5.11 and A5.12 put the A3 freeze (`prereg-v0.4.3`) and `prereg-v0.4.4` **before** step 6 (C9 GO, then A6). After the freeze nothing may land in `vulcan/`, `bench/`, `scripts/a3/`, `scripts/e7/` or the A5.11.6 trees. A6 is a removed-peptide null construction run through the frozen pipelines, and it will almost certainly need code under `vulcan/gates/` or `vulcan/index/` (or must import frozen code from outside). A5 gives no rule, so the GO branch either breaks the freeze check or has no legal implementation.
- **Problem 2: tagging before drafting delays A6.** A5.1.4 requires A5 to be "adopted **and tagged**" before A6 is drafted. Tagging needs all of A5.11, the freeze and the claim-script review, so A6 drafting waits days for no statistical reason. Only A5.10 needs to bind A6.
- **Problem 3: no time box after GO.** If A6 is not adopted, fails review, or never finishes, T1 waits indefinitely, which leaves open the option "never look".
- **Problem 4: GO's quantity is undefined.** "Simulated P(v1 minimum)" is not defined against A5.10.3's floor (A4 n_min, ≥ 10 clusters, outer envelope). A pilot that used A2 minima could declare GO on a looser basis.
- **Problem 5: the record is unbound.** `config/a5.yaml c9.outcome` is a mutable field that `sealed_run.py` trusts (A5.11.1.4), with no link to the committed pilot verdict.
- **Required fix:**
  1. Bind A5.10 at **adoption** (DESIGN §6.0d plus the RESEARCH_TREE note), not at `prereg-v0.4.4`.
  2. State that A6 code lives outside the frozen trees and imports frozen code unchanged. If it cannot, the GO branch takes the freeze after A6 code lands, and A5.11's order is re-stated for that branch.
  3. Time-box GO: if A6 is not tagged within N days of the GO record, or has not completed or stopped label-free within M days, record "A6 not completed (time box)" and proceed. N and M are PI points.
  4. Define "P(v1 minimum)" as the probability of reaching every v1 n_min of A5.10.3 (A4 floor values) under the pilot's measured yields.
  5. `c9.outcome` must cite the pilot verdict file by commit and sha256. `sealed_run.py` recomputes GO or NO-GO from that file (≥ 0.8) or from the time-box rule, and refuses on a mismatch.

### M7. S6 (A5.10) is weakened in two places and has three loopholes
- **Weakening 1: the cluster cap is dropped.** The review set "A4's margins, n_min floors, envelope **and IUT rules**" as the floor. A5.10.3 drops A4.3.6 (cluster cap and the capped/uncapped intersection–union) by classing it as foreign-organism selection, but clusters are labs, not organisms.
- **Weakening 2: "cluster" is left undefined.** A5.10.3 excludes A4.3.4 (the definition of a cluster) while requiring cluster-unit intervals and "fewer than 10 clusters ⇒ not evaluable".
- **Weakening 3: no stopping-rule requirement.** A label-free stopping rule committed before listing (A4.3.7-type) is not in the floor.
- **Loophole 1: draft shopping.** A T1-exposed PI can commission several T1-blind drafts, or re-commission after a draft fails review, and adopt the preferred one. The PI can also decide after T1 whether to commission any test. A5.10.4 forbids only choosing "among options" within one text.
- **Loophole 2: commit subjects reveal T1.** T1-blind agents in this repo see recent commit subjects in their session context (the harness shows `git log` subjects automatically), and they see the RESEARCH_TREE **status summary**, which is not a dated entry. A5.10.4's definition of T1-blind covers neither.
- **Loophole 3: no recorded fate of drafts.** Nothing requires that every commissioned null-set draft, and what happened to it, is reported.
- **Required fix:**
  1. Add A4.3.6's cap and capped/uncapped IUT, a definition of cluster (submitting lab or PI, as A4.3.4) and a label-free stopping rule to the floor.
  2. Every commissioned draft (prompt, text, reviews, fate) is committed and reported with T1. At most one null-set amendment per scorer may be open at a time. Abandoning one is reported with its reason.
  3. After the look, commit subjects and the RESEARCH_TREE status summary contain no T1 or T5 number or verdict; results are referenced by path only. T1-blind agents run in a clone or worktree at a pre-look commit, plus the specific files they need.
  4. Extend the definition of T1-blind to commit messages and the status summary.

### M8. D18's pre-registered secondary arm C is not computed, and A5.5.3 then makes it impossible
- **Problem:**
  - D18 pre-registered "a **secondary, reported-only** arm C … with an outcome-independent admission screen (≥ 4-ion b/y ladder vs human targets; plastid, photosynthesis or food-taxon proteins)".
  - `bench/t1/` has no such screen (no D18, ladder or screen code).
  - A5 makes arm C "the honest residual defence" (A5.13.3) and adds claim-level arm-C voids, yet omits the pre-registered diagnostic that separates real foreign matches from scorer excess.
  - After the look, A5.5.3 forbids new analysis.
- **Required fix:** either implement D18's screened arm C (reported only, never a verdict or a void) in `bench/t1/a5.py` before the freeze, with dry-run evidence, or record in A5 a dated deviation that it is not computed, with the reason. Do not leave it silent.

### M9. Clauses that are still normative and not amended (C1 requires explicit replacement text)

| Clause | Text | Conflict with A5 |
|---|---|---|
| A2 P3 | "Not evaluable blocks FDR claims, exactly like a fail." | A5 lets FDR-filtered numbers be reported, labelled; needs the same "FDR-*control*" edit as A1.2 |
| A2.11.6 | "P4 allows exactly one re-test. Any A3 must count A2 as a prior attempt" | Contradicts A5.10's test-2/test-3 sequence and the naming. Replace with a pointer to A5.10.1–A5.10.2 |
| A3.2.5, last bullet | "Family S (A3.3.5) is computed on exactly one list, determined by the tier." | Contradicts F1 |
| A3.0 Q5 | "The single production rule (A3.2.6) is a disclosed, pre-specified fallback …" | A3.2.6 is replaced; Q5 dangles |
| A3.3.2 and A3.6 step 6 | "the same `sealed_run.py` invocation and ledger line"; "whose ledger line names the v1 primary and the v2_sep secondary arm" | Compatible only if `look.purpose` for T1 names both arms. Say so and pin it in `config/a5.yaml` |
| A3.1.5, A3.1.10 | DN2 MS1 availability; `gates_a3.py` tested end to end | A5.12 step 4 lists only "freeze hashes, G0, E7b commit, driver commit". State explicitly which A3.1 preconditions are moot (DN2 not evaluable, D24) and which still gate `prereg-v0.4.3` |
| §1, "the bar" | "Compared at **entrapment-verified** FDP"; goal "at *empirically verified* FDR" | The same overclaim as "entrapment-validated". Add a §1 note that claims use the A5.9 labels and not this wording |
| §6.8 | PSM counts "labelled 'not FDR-controlled'" | Clashes with A5.11.4.4's string test, which bans "FDR-controlled". Allow the negated form |
| §6.3 required arms and A1.5 | v1 and v2_sep need P, C and M for "supported" | `evaluate.py` combines P and C under `static_low_capacity`. State that any "overall" verdict for v1 or v2_sep includes T5 arm M, or that only per-arm verdicts are printed |

**Required fix:** add a quoted old/new entry to A5.3 for each row.

### M10. Banned-wording enforcement covers only `claims.md`, and the list is incomplete
- **Problem:**
  - A5.11.4.4 tests generated claim text only. REPORT.md prose, slides and preprints (A5.9.4) are written by people or agents, and nothing checks them.
  - The banned list misses: "verified" (as in §1's "entrapment-verified"), "empirically verified", "calibrated", "FDR control", "controlled at", "at 1% FDR" when not followed by an L-Q or L-M tag, and "gate-passed".
- **Required fix:**
  - Add a committed lint (`tests/report/test_report_labels.py`) over REPORT.md and any release text. It checks that:
    - the banned phrases are absent outside quoted labels;
    - L-S is in the abstract and the headline;
    - every table caption with FDR-filtered numbers contains its L-labels;
    - the headline sentence does not start with a number.
  - Run it before any release.

---

## Minor

1. **F3 vs A5.10.6.**
   - F3 says "No T1 or T5 outcome changes either designation", but A5.10.6 requires "no void under A5.7.1/A5.7.2". Voids are T1/T5 outcomes.
   - Fix: reword F3 to "No T1 or T5 outcome can upgrade a designation; a T1/T5 void blocks any later upgrade."
2. **"Strongest comparator" in the headline (A5.9.4) is undefined.**
   - Fix: "the A1 comparator with the lowest P1 point ratio", fixed now.
3. **The GO-pass state (A5.9.3) keeps L-R implicitly.**
   - If A6 passes tail-G2/G1/G3, A1.8's precondition holds.
   - Fix: say whether L-R is dropped, or kept with "precondition met at test 2 (A6), before T1".
4. **L-M wording.**
   - "does not use either method's own q-values" is inaccurate: the matched threshold moves along each method's own q grid.
   - Fix: "uses each method's own q only to order its list, not as an error estimate".
5. **L-G "Passed:" list.**
   - Problem 1: "(A1 null set)" reads as applying to G0 and G4, which ran on D0/D1 and on training decoys.
   - Problem 2: G2-original was "nearly vacuous (winner changed on 0.08% of nulls)" (M3e NOTES).
   - Problem 3: dev arm M is missing.
   - Fix: give the data set per gate, add the vacuity note, and add "arm M (dev D7): inconclusive, not violated".
6. **Retry rules (A5.5.5, A5.11.5) contradict each other.**
   - (a) After a failure commit, the base commit has no `milestone-*` tag, so a retry needs a new tag. State it.
   - (b) `config/a5.yaml` says C9 is "the only field that may change after the tag", yet retries add deviations and may change hashes after a re-freeze. List all mutable fields, and require a `prereg-v0.4.4.x` tag after any re-freeze.
   - (c) A failure that cannot be reproduced on D0/D1/D7 has no rule. State "fix from the committed frames only, reviewed T1-blind", or "the look is permanently failed".
7. **`data/sealed_derived/` is not in `.gitignore`.**
   - The move makes the tree dirty, so `sealed_run.py` (and the freeze check, if any file lands under a checked path) refuses later sealed runs. It also invites a `git add` of T1 PSM tables.
   - The work dir also contains 25 replicate indexes (~108 M peptidoforms each on human; dry-run NOTES). These are not sealed-derived and are huge.
   - Fix: add the path to `.gitignore`, and tar only the outputs derived from sealed data.
8. **Derived files are described inconsistently.**
   - A5.5.2 tars "stores, mzML, MS1 tables". The driver deletes those (l. 433–436), and A5.11.2.1 forbids `--keep-sealed-derivatives`. A5.14.5 says "lock instead of deletion".
   - Fix: pick one and make the text match the code.
9. **C2's "the driver asserts both tags" is dropped.** Only `sealed_run.py` checks. Add the assertion to the driver (defence in depth).
10. **D28's cost figure.**
    - A5.0.4 cites 0.85–2.0 TB (A4.7.4, correct), but D28.1 says 0.75–1.9 TB.
    - Fix: say "A4.7.4's 0.85–2.0 TB supersedes the 0.75–1.9 TB quoted in D28.1 (A4.11.6)".
11. **The A1.8 exception covers T1 and T5 only, but A5.1.5 applies the labels to every sealed set.**
    - T2–T4 therefore still cannot use v1's re-ranking pipeline.
    - Fix: state that this consequence is intended (T4 stays blocked until a gate passes), or extend the exception by a later amendment.
12. **The C9 time box can be gamed by delay.**
    - A NO-GO by default favours the party who wants T1 now.
    - Fix: disclose in the adoption note who controls the pilot's commit, and that the pilot agent commits its verdict as soon as it exists, whatever its value.

---

## Checks that passed (no change needed)
- **Old-text quotes.** These are accurate against DESIGN and the A2/A3 texts, allowing for ellipses and flattened bullets: D12, §6.4, A1.2, A1.4, A1.5, A1.8, A2 P4, A2.9.6, A3.2.6, A3.3.1, A3.3.5, A3.3.6, D26.3, §6.6 item 2.
- **A5.4 numbers** match M3e NOTES, E7 §6, D23, D24 and A3_preconditions:
  - n_tail 751, UCB 0.742, diff UCB 0.124;
  - ρ 0.048 [0.019, 0.076];
  - G3 1.03 pp;
  - G5 0.506 / 0.504;
  - 894,579 nulls, 208,150 competitions;
  - τ 0.298 (v2_sep), arm M 1.05%.
- **A5.0.4** matches the A4 tables: 0.43 / 0.08, grid 0.01–0.98.
- **D25 / E7b caveats** match E7b NOTES: 4–13%, 0.999 [0.951, 1.039].
- **Wording:** "entrapment-validated" is removed; "contaminated" is removed as a qualifier; "FAILED" is in every short form; A4 is recorded as "drafted, reviewed twice, not adopted, not run"; the escrow-key argument is withdrawn.
- **Family S** is fixed on the full list, titled, and decided by the claim script. **Designations** are fixed. **v1 void-on-violated** is symmetric with A3.3.6 and stricter: it adds the sequence unit, T5 and A6.
- **A5.7.4's power statement is correct:** at r ≈ 1, LB ≈ FDP̂/2, so the void rules catch only gross violations.
- **T5 is in scope for training.** Yeast is absent from the training manifests (`train.tsv`, `train_m7.tsv`, no yeast rows). T5 is two HF-X runs, one lab.
- **`VULCAN_T1_FMATCH_0001`** is correctly excluded from the look.
