# Amendment A5 (DRAFT, revision 1): the T1 look without passed null-spectrum gates; closed options, literal labels, one look in code, arm M on T5, pre-committed rules for later null tests

- **Status:** draft for adoption as DESIGN.md §6.0d. It is not binding. DESIGN.md, the tags and every A1–A4 file are unchanged.
- **Author:** statistician (agent), fresh context. I did not write A1–A4, their reviews, the D28 review, the T1 driver, E7/E7b, the C9 pilot or any gate run. I will not run the look or any null-set step.
- **Date:** 2026-10-28.
- **Revisions:** revision 0 was reviewed by a fresh independent T1-blind reviewer (`reviews/amendment_A5_review.md`; verdict "adopt with changes": 2 blocking, 10 major, 12 minor). Revision 1 resolves every finding; the resolutions are listed in A5.15.
- **Mandate:** RESEARCH_TREE D28.3, implementing `reviews/pi_proposal_D28_review.md` changes C1–C8, C10, C11 and S6 (C9 is D28.2 and is not re-specified here).
- **Basis (read):**
  - RESEARCH_TREE D9–D28; `reviews/pi_proposal_D28_review.md`;
  - DESIGN §1 and §6 (§6.0 A1, §6.0b A2, §6.0c A3);
  - `reviews/amendment_A2_draft.md` (P1–P5, A2.9.5–A2.11); `reviews/amendment_A3_draft.md` (rev 1); `reviews/amendment_A4_draft.md` (rev 3);
  - code: `bench/sealed_run.py`, `bench/t1_driver.py`, `bench/t1/*.py`, `scripts/m3e/armM*.py` and `scripts/a3/dev_armM.py` (headers), `config/a3.yaml`, `.gitignore`;
  - results NOTES: `2026-10-27_T1driver_dryrun`, `2026-10-28_E7b_fair_speed`, `2026-10-28_A3_preconditions` (verdict tables), `2026-10-28_E7_chimeric_v2` §6 (gates), `2026-10-27_M3e_hardening` (gate table).
- **Not read:**
  - any sealed file (T1–T5, including T5's SDRF);
  - `results/SEALED_LEDGER.md` beyond its header (it has no entries);
  - any DN2 per-label quantity (none exists) and any DN3 listing (none exists);
  - any C9 pilot output.
- **T1 exposure: none. T5 exposure: none. C9 pilot exposure: none.**

---

## A5.0 What this amendment does (first paragraph, required verbatim by the D28 review C1)

**This amendment relaxes the conditions for FDR-based claims and permits the T1 look although v1 failed its only completed null-spectrum test. It is adopted knowing that failure, before any T1 or DN3 information exists.**

It is not a timing amendment of D12. It changes what may be claimed and in which words, and it closes every option that the relaxation opens. Specifically:

1. **The adopter's information set** (D28 review §A). D28 and A5 are decided knowing that:
   - v1 (`static_gbm_pred`) FAILED tail-G2, G1 and G3 under A1 (D12);
   - A2 ended not evaluable before any test (D24);
   - A4 (DN3) will not be run (D28.1).

   No T1, T5, DN3 or per-label DN2 quantity exists.
2. **What A5 does not do.** It makes no FDR-*control* claim available for v1 or v2_sep. The following wordings are not used for either scorer unless a null-spectrum test passes under A5.10: "validated FDR control", "controls FDR", "gate-validated", "entrapment-validated", "entrapment-verified", "empirically verified FDR". The full list is A5.11.7.
3. **A4 status.** A4 rev 3 (`reviews/amendment_A4_draft.md`, a158cf8) is recorded as **drafted, reviewed twice, not adopted, not run**. Nothing in A4 binds, except the parts that A5 adopts by reference (A5.5, A5.10).
4. **Provenance of "A4 not run" (C11), stated in full.** D28.1 was decided before any DN3 listing and before any T1 look, on label-free feasibility grounds only.
   - **Evaluability.** P(v1 reaches its minimum) is 0.43 (e = 0.15) and 0.08 (e = 0.35) in the reference cell, which is anchored on DN2's own listing counts. Across the sensitivity grid (6 to 24 usable files per project) it is 0.01–0.98. The number of usable RAW files per project is unknown before listing; the reference cell is the best estimate (`reviews/checks/a4_r3_walk.out`).
   - **Power.** P(clean v1 passes) ≈ 0.06–0.4. P(a margin leak is detected) ≈ P(evaluable) × P(fail | evaluable, margin leak) ≈ (0.08–0.43) × ≈0.9 ≈ 0.07–0.4.
   - **Cost.** 3–4 days and 0.85–2.0 TB (A4.7.4). This supersedes the 0.75–1.9 TB quoted in D28.1 (A4.11.6). A not-evaluable outcome would consume most of the remaining foreign-organism pool.
   - **The off-cluster escrow-key argument is withdrawn** as a reason. It argued only against running T1 in parallel with A4.

## A5.1 Scope and branch logic

D28.2 runs the C9 removed-peptide pilot first. A5 is written so that **its text is the same in both branches**:

| Part | Content | Applies |
|---|---|---|
| **I. Fixed choices** (A5.2–A5.4, A5.6–A5.11) | family S list and title, void rules, designations, arm M on T5 in the same look, literal labels, one look in code, every branch in one run, hashed claim script, full reporting, A1.8 statement, matched-FDP/own-q split, attempt ledger, S6 rules | **Either branch** |
| **II. The look without passed gates** (A5.5) | the T1 look may proceed although no null-spectrum gate has passed for v1 | Else-branch (C9 NO-GO). In the GO branch it applies once the C9 null-test amendment ("A6") has completed, stopped label-free, or timed out (A5.1.3), whatever its verdict |

1. **C9 outcome record (bound to the pilot's own file).**
   - **The quantity.** GO iff the pilot's committed P(v1 minimum) is ≥ 0.8. "P(v1 minimum)" is defined as the simulated probability that **every v1 n_min and cluster minimum of the A5.10.3 floor** is reached, under the pilot's measured label-free yields. The A4 floor values apply, not A2's.
   - **The record.** `config/a5.yaml c9` cites the pilot verdict file by commit and sha256. `sealed_run.py` recomputes GO or NO-GO from that file, or from the time-box rule below, and refuses on any mismatch.
   - **The time box.** If no committed pilot verdict exists 48 h after the D28 commit (3a71356), the outcome is NO-GO (D28.2: "one-day pilot").
   - **Disclosure.** The adoption note discloses who controls the pilot's commit. The pilot agent commits its verdict as soon as it exists, whatever its value. The pilot is label-free, and its value may be known to the adopter; this is disclosed, not hidden.
2. **NO-GO (else-branch).** A5 governs the look. At the look, every scorer's null-spectrum status is "no passed test" (A5.4).
3. **GO.**
   - A6 is drafted, reviewed, tagged (`prereg-v0.4.5`) and run under the A5.10 floor **before** the look.
   - **Time box (PI point A5.14.2; defaults written here):** A6 is tagged within 7 days of the GO record, and completes step 5 or stops label-free within 14 days of that record. Otherwise "A6 not completed (time box)" is recorded, the look proceeds, and A6 continues afterwards as a post-T1 test under A5.10.
   - The look proceeds **whatever A6 shows** (pass, fail, not evaluable, label-free stop, or time box). This closes the option "decide after A6 whether T1 happens". A6's verdict enters the labels and the void rules through A5.4, A5.7 and A5.9.3.
4. **When A5 binds.**
   - A5.10 binds **from adoption** (the DESIGN §6.0d text and the dated RESEARCH_TREE note). It therefore binds A6 even before `prereg-v0.4.4`. A6 may be drafted as soon as A5 is adopted.
   - **A6 code** lives outside the frozen trees (e.g. `scripts/a6/`, `config/a6*`) and imports frozen code unchanged. If A6 cannot be implemented that way, the GO branch moves the A3 freeze commit to after A6's code has landed. A5.11's order is then followed with "A6 code" inserted before the freeze, and the re-done A3.1.6 dry run covers it.
5. **Datasets.**
   - A5 authorises one look at **T1** and one at **T5** (arm M), in the same `sealed_run.py` invocation (A5.5, A5.6).
   - It authorises no look at T2, T3 or T4; those keep their own preconditions.
   - **The A1.8 exception (A5.3) covers T1 and T5 only.** So T2–T4 FDR-based uses of v1's re-ranking pipeline stay blocked until a gate passes or a later amendment extends the exception. This consequence is intended.
   - The label rules of A5.9 apply to every FDR-filtered v1 or v2_sep number anywhere in REPORT.md or in release text. Gate status is a property of the scorer, not of the data set.

---

## A5.2 Fixed choices (Part I; closed, so no option remains at look time)

| # | Choice | Source |
|---|---|---|
| F1 | **Family S runs on the full v2_sep list** (all merged strata, common `q_sep` cutoff), as A3.3.5's sequence S1–S5 against every A3.3.3 comparator, by LB and leave-one-run-out minimum. Its table is titled **"Secondary arm; not gate-validated"**. The stratum-1-only list (A5.5.5) is reported descriptively, never in its place | D28.3; review C3 |
| F2 | **v1 void-on-violated** (A5.7.1), symmetric with A3.3.6, plus arm M on T5 and any completed null-test fail | D28.3; review C4 |
| F3 | **Designations, fixed now:** v1 is "engine default, FDR not validated by null-spectrum tests"; v2_sep is "candidate". **No T1 or T5 outcome can upgrade a designation, and a T1/T5 void blocks any later upgrade.** Only a passed null-spectrum test can upgrade one, under A5.10.6, disclosed as before or after T1 | D28.3; review C5 |
| F4 | **Arm M on T5 is opened in the same look as T1**, for v1 and v2_sep (A5.6) | D28.3; review C6; DESIGN §6.6, A1.5 |
| F5 | **Literal per-arm verdict labels** (§6.3 words only), and the word **"FAILED"** (upper case) for A1. The contamination explanation is called **"unconfirmed"** and is not used | D28.3; review C.1–C.3, C10 |
| F6 | **One look, enforced in code** (`sealed_run.py`), for T1 and T5 | D28.3; review C7, B11 |
| F7 | **Every branch is computed in that one run**, with a completeness rule. There is no re-search and no new estimator afterwards (A5.5.4) | review C7 |
| F8 | **A hashed, committed claim script** maps outputs to claims and labels. It is reviewed T1-blind before the look (A5.11.4) | review C7 |
| F9 | **Full reporting:** the claim script's output goes into REPORT.md verbatim, including every void, every "inconclusive" and every comparator row where Vulcan loses. A committed lint checks REPORT.md and release text (A5.11.7) | review C8 |
| F10 | **A1.8 top-K re-ranking statement** on every number from a re-ranking pipeline (A5.3, A1.8) | review B4, C.5 |
| F11 | **Matched-FDP vs own-q split:** two label classes, L-M and L-Q (A5.9.1) | review C.7 |
| F12 | **S6 rules for any later null-set test:** v1's next completed test is test 2 at α/2; A4's margins and statistics are a floor; design and execution are by T1-blind agents (A5.10) | review S6 |
| F13 | **Finality of sealed verdicts:** no explanation found after the look (sample composition, contamination, reference-pipeline behaviour, comparator failures) un-voids a claim or changes a verdict (A5.7.5) | review D12 pattern; A5 review M5 |

---

## A5.3 Replacement text for every affected clause

Each entry quotes the clause now in force and gives its replacement. The "Old" quotes come from:
- DESIGN.md (§1, §6.0–§6.8);
- the normative A2/A3 texts, adopted verbatim by §6.0b/§6.0c;
- the A4 draft;
- RESEARCH_TREE.

Where a clause is only extended, the old text stays and the addition is marked.

### D12 (RESEARCH_TREE, 2026-10-27)
- **Old:** "The A1 gate failures for frozen static_gbm_pred **stand as recorded**: tail-G2, G1 and G3 FAIL. FDR claims and the T1 look stay blocked. The post-hoc explanation (human-contaminated S. aureus null project; offset nulls being shifted human spectra) is plausible but was found after seeing the failure, so it must not be used to re-score the *existing* nulls. Path: …"
- **New** (a dated note under the D-number of adoption): "The A1 gate failures for frozen static_gbm_pred stand as recorded: tail-G2, G1 and G3 FAILED, final for those nulls. **No FDR-control claim is made for static_gbm_pred unless a later null-spectrum test passes under DESIGN §6.0d A5.10.** The T1 look, together with arm M on T5, may happen **once**, after `prereg-v0.4.4`, under §6.0d (A5.5). Every FDR-dependent number carries the A5.9 labels. A post-hoc contamination explanation exists, is unconfirmed, and is not used, neither to re-score the existing nulls nor as a qualifier in any label."

### §1, "the bar" and the goal sentence
- **Old:** "Compared at **entrapment-verified** FDP" and "at *empirically verified* FDR".
- **Added note:** "Under §6.0d (A5), claims about v1 and v2_sep use the A5.9 labels and never the words 'entrapment-verified' or 'empirically verified'. Entrapment verdicts are quoted literally per arm. The §1 non-inferiority test itself (own q **and** matched FDP; own q only against 'supported' comparators) is unchanged (A5.9.2, P1)."

### §6.3, required arms (clarified)
- **Old:** "Models with sequence capacity (…any peptide-token model): **P, C and M**"; "Supported | mean FDP̂ ≤ α **and** upper 95% CI ≤ 1.5α, for each required arm".
- **Added:**
  - Under §6.0d, verdicts for v1 and v2_sep are printed **per arm only** (P, C, M).
  - The driver's combined "overall" field (`static_low_capacity`, P and C) is not quoted as a verdict for these scorers.
  - No "supported" statement is made for a list unless every required arm, including T5 arm M, is "supported" for it.
  - At protein level, arm C is **not computed**: the frozen driver has no arm-C protein construction. L-E prints this.

### §6.4, "Gates for learned or predicted components (all must pass before any FDR claim)"
- **Old:** "**Gates for learned or predicted components** (all must pass before any FDR claim):"
- **New:** "**Gates for learned or predicted components** (all must pass before any FDR-*control* claim, i.e. any statement that a list is FDR-controlled or that the method's FDR is validated). Under §6.0d (A5), FDR-*filtered* numbers of a scorer that has not passed every gate may be reported only with the A5.9 labels, which name each failed, not-evaluable or not-run gate."

### §6.8, PSM label
- **Old:** "PSM counts are labelled 'not FDR-controlled'."
- **Kept.** The A5.11.7 lint allows this negated form.

### A1.2 (§6.0 item 2), last paragraph
- **Old:** "With n_tail < 200 the gate is "not evaluable", and that blocks FDR claims for that scorer. The original PSM-level G2 is kept unchanged."
- **New:** "With n_tail < 200 the gate is "not evaluable". A fail or not evaluable blocks FDR-*control* claims for that scorer; under §6.0d (A5), FDR-filtered numbers may be reported only with the A5.9 labels. The original PSM-level G2 is kept unchanged."

### A1.4 (§6.0 item 4)
- **Old:** "**Gates G1, G3 and G4 are required** (§6.4) before any FDR claim. G3 is run on null spectra binned by score."
- **New:** "**Gates G1, G3 and G4 are required** (§6.4) before any FDR-*control* claim. Under §6.0d (A5), their status (pass, FAILED, not evaluable or not run, with the null set named) is stated in the A5.4 gate table, which is quoted with every FDR-filtered number. G3 is run on null spectra binned by score."

### A1.5 (§6.0 item 5)
- **Old:** "**Arm M is required for predicted-feature scorers.** … On dev, this means D7 yeast (PXD031616) searched against yeast + human entrapment at r ≈ 1. T5 is used only at the milestone."
- **New:** "… (unchanged up to "r ≈ 1"). **T5 is opened once, in the same `sealed_run.py` look as T1** (§6.0d A5.6), for v1 and v2_sep, with the hyperscore pipeline as reference. A "violated" arm-M verdict on T5 voids that scorer's FDR-based claims (A5.7). After this look, T5 is spent for every scorer, present or future."

### A1.8 (§6.0 item 8)
- **Old:** "**Re-ranking of the top-K by a static scorer** is allowed only if tail-G2, G1 and G3 pass."
- **New:** "**Re-ranking of the top-K by a static scorer** is allowed only if tail-G2, G1 and G3 pass.
  - **Exception (§6.0d A5), for the T1 look and the T5 arm-M run only.** These use the frozen pipelines as they are, although this precondition is not met:
    - v1: engine `report_psms 5`, rescored by `static_gbm_pred`;
    - v2_sep: chimeric search with top-5 per group and the frozen v2 model.
  - **Statement.** Every number from these pipelines carries: *"top-K re-ranking by the static scorer is used although A1.8's precondition (tail-G2/G1/G3 pass) is not met."*
  - **Reported.** The look reports, label-free, the fraction of spectra (all spectra, and among units accepted at own q 1%) whose rescored best candidate is not the engine's rank-1 candidate (A5.11.2).
  - **If the precondition is met at a test that completes before the look** (A6 pass in the GO branch), the statement becomes: *"top-K re-ranking by the static scorer; A1.8's precondition met at test 2 (A6), before T1."*"

### A2 P3 (§6.0b, adopted verbatim)
- **Old:** "minimum-n rules that make a test **not evaluable**. Not evaluable blocks FDR claims, exactly like a fail."
- **New:** "… Not evaluable blocks FDR-*control* claims, exactly like a fail. FDR-filtered numbers may be reported only with the §6.0d A5.9 labels."

### A2 P4 (§6.0b)
- **Old:** "**Not evaluable** keeps FDR claims blocked. Only a new amendment (A3) with yet another fresh set can follow, and it must report A2's result alongside its own."
- **New:** "**Not evaluable** keeps FDR-*control* claims blocked. Only a new amendment with yet another fresh set can follow. It must meet the floor of §6.0d A5.10, and it must report A1, A2, A4 (not run) and every later attempt alongside its own. (The name "A3" here is historical; A3 took that name for another purpose, A3.0.)"
- The rest of P4 ("A fail on DN2 is final…", "Nobody may re-draw…") is unchanged.

### A2.9.6 (§6.0b): pass wording and T1 precondition
- **Old:** "A **pass** satisfies tail-G2, G1 and G3 for `static_gbm_pred` and nothing else (G4, G5, arm M, B1/T1 preconditions still apply). It is quoted as "passed on fresh nulls under A2, after failing under A1 on contaminated nulls (D12)". The A1 failure is never omitted."
- **New:** "… It is quoted as "passed on fresh nulls under <amendment>, after A1 (2026-10-27): tail-G2, G1, G3 FAILED; final for those nulls; a post-hoc contamination explanation exists, is unconfirmed, and is not used". The A1 failure is never omitted. A pass is no longer a precondition of the T1 look (§6.0d A5.5)."
- The phrase "on contaminated nulls" is removed everywhere it qualifies the A1 result. That includes A4.9's draft wording, which is not adopted.

### A2.11.6 (§6.0b)
- **Old:** "**Repeat attempts.** P4 allows exactly one re-test. Any A3 must count A2 as a prior attempt and report both."
- **New:** "**Repeat attempts** follow the attempt ledger and level rule of §6.0d A5.10.1–A5.10.2. Every later attempt reports every earlier one, including label-free stops (A2) and not-run drafts (A4)."

### A3.0 Q5 (§6.0c)
- **Old:** "… The single production rule (A3.2.6) is a disclosed, pre-specified fallback, and both verdicts are quoted whenever it is invoked."
- **New:** "… Designations are fixed by §6.0d A5.2 F3 and change only under A5.10.6, which quotes every scorer's gate history whenever it is invoked."

### A3.1 preconditions (§6.0c): which still gate `prereg-v0.4.3`
- **Moot (D24: DN2 not evaluable before step 4; A3.2 is moot):**
  - A3.1.5 (DN2 MS1 availability);
  - the DN2 parts of A3.1.9 (already recorded, D23);
  - the DN2 end-to-end part of A3.1.10.
- **Still required before `prereg-v0.4.3`:**
  - A3.1.2 (review committed; met);
  - A3.1.3 (freeze, extended by A5.11.6);
  - A3.1.4 (G0 at the freeze);
  - A3.1.6 (driver plus dry run, re-done with A5 code);
  - A3.1.7 (E7b frozen);
  - A3.1.8 (met, D23);
  - A3.1.10's synthetic tests of `scripts/a3/` (already committed; re-run at the freeze).

### A3.2.5 (§6.0c), last bullet
- **Old:** "Family S (A3.3.5) is computed on exactly one list, determined by the tier."
- **New:** "Family S (A3.3.5) is computed on exactly one list, the full v2_sep list, fixed by §6.0d A5.2 F1."

### A3.2.6 (§6.0c), production rule
- **Old:** "**Production default for OT-HCD DDA:**
  - **v2_sep** iff all of: tier 2 passes; A3.1.8 (i)–(iii) hold with no "violated" arm M; and v2_sep has no "violated" verdict on T1 (A3.3.6);
  - otherwise **v1**, iff v1 passes A2;
  - otherwise no OT-HCD scorer has gate-validated FDR claims.

  … Until T1 has been looked at, v2_sep can be at most the **provisional** default."
- **New:** "**Designations, fixed by §6.0d (A5.2 F3):** v1 is "engine default, FDR not validated by null-spectrum tests"; v2_sep is "candidate". No T1 or T5 outcome can upgrade a designation. A designation changes only by the rule of A5.10.6 (a passed null-spectrum test), disclosed as before or after T1. Until then the status is "no OT-HCD scorer has gate-validated FDR claims", and it is stated wherever a designation is named."
- The tier definitions of A3.2.6 remain only as the meaning of "tier" in A5.10.6.

### A3.3.1 (§6.0c), first bullet
- **Old:** "The T1 look itself. It stays governed by D12 and A2.9.6: it requires v1 to pass A2.
  - If v1 fails or is not evaluable, there is no T1 look under A3, whatever v2_sep's result.
  - A3 does not pre-specify a v2-only route. Any such route needs a new amendment adopted before any T1 look, and that amendment must report v1's failure (closed; formerly O1)."
- **New:** "The T1 look itself is governed by §6.0d (A5.5), the new amendment that A3.3.1 required. It requires `prereg-v0.4.3`, `prereg-v0.4.4` and the C9 outcome record (A5.1). It does not require v1 to pass any null-spectrum gate. v1's failure under A1 is reported as A5.4 states."

### A3.3.1, second bullet ("v2_sep's claims on T1 need all of")
- **Old:** "tier 1 or tier 2 passed (A3.2.6); A3.1.8 (i)–(iii) hold, with no "violated" verdict; no "violated" verdict on T1 (A3.3.6)."
- **New:** "v2_sep's family-S results on T1 are reported as "secondary arm; not gate-validated" (A5.2 F1). They require all of:
  - A3.1.8 (i)–(iii) hold with no "violated" verdict (met, D23);
  - no "violated" verdict on T1 (A3.3.6) or on T5 (A5.6);
  - no completed null-spectrum test FAILED by v2_sep before the look (A5.7.2)."

### A3.3.2 and A3.6 step 6 (§6.0c)
- **Old:** "in **the same `sealed_run.py` invocation and ledger line**. There is no second look"; "one `sealed_run.py` invocation, whose ledger line names the v1 primary and the v2_sep secondary arm".
- **New:** "… in the single `sealed_run.py --look a5` invocation (§6.0d A5.5). That invocation writes one T1 line and one T5 line in one commit. The T1 line's purpose string, pinned in `config/a5.yaml look.purpose.T1`, names the v1 primary and the v2_sep secondary arm. There is no second look."

### A3.3.5 (§6.0c), family S list rule
- **Old:** "The list is the full v2_sep list if tier 2 passed, or the stratum-1 list if only tier 1 passed (A3.2.6). The list is fixed by gate outcomes under rules set here."
- **New:** "**The list is the full v2_sep list** (all merged strata, common `q_sep` cutoff). It is fixed by §6.0d (A5.2 F1), not by any outcome. The table is titled "Secondary arm; not gate-validated". The stratum-1-only list (A5.5.5) is reported descriptively in a separate table, never in place of the full list."
- The hypotheses, order, bars, comparators, LB and leave-one-run-out rules of A3.3.5 are unchanged.
- **The family-S decision is made by the claim script (A5.11.4), not by an analysis agent.** This replaces the `bench/t1/a3.py` docstring sentence "decided … by the analysis agent".

### A3.3.6 (§6.0c), last bullet
- **Old:** "**Arm M** for v2_sep on T5 at the milestone (A1.5)."
- **New:** "**Arm M** for v2_sep on T5, in the same look as T1 (§6.0d A5.6). Any "violated" verdict there (overall or any merged stratum, 1% or 5%, convention (a)) voids all v2_sep T1 claims, as the T1 void rule above does."
- The phrase "This applies whichever tier passed" is moot: the void rule applies to the full list.

### A3.3.7 (§6.0c), claim table
- **Old:** the five-row table beginning "fail / not evaluable | any | — | No T1 look under A3", and the required wording beginning "pre-registered secondary arm; gated on fresh nulls under A2/A3 (tier N); …".
- **New:** **replaced in full by the A5.9 claim and label tables.** v1's result is still reported first and in full in the same document (unchanged).

### A3.3.8 (§6.0c)
- **Unchanged.** It is extended by the completeness rule (A5.5.4): an incomplete look is a failed look, not a look with fewer arms.

### A4.6 (draft, not adopted): multiplicity across attempts
- **Old status:** draft text of an unadopted amendment; it binds nothing.
- **New:** A5.10.2 **adopts A4.6's definitions and geometric rule by reference as normative**, with the attempt ledger updated in A5.10.1. A4.6's DN3-specific sentences (simulated false-pass rates at the DN3 margins) are not adopted.

### D26 (RESEARCH_TREE, 2026-10-28)
- **D26.3, old:** "**T1 runs into escrow** after `prereg-v0.4.3` and the A4 tag. Outputs are hashed and sealed, and nobody (PI or agents) reads them until A4 step 5 or until the project nears its time limit, whichever comes first. … sealed_run.py must enforce one look. D12's timing is amended only in this escrow sense."
- **D26.3, new:** "**Superseded by D28/§6.0d (A5), as a reversal of a PI decision.** There is no escrow. T1 (with T5) is computed and read once, at the A5.5 look. The one-look rule of D26.3 is kept and implemented (A5.11.1), and A4.8.2's root, leak-test and console mechanics are adopted without encryption (A5.5.3)."
- **D26.1, D26.2 and D26.4** lapse with A4 (not adopted, not run). A later amendment that re-proposes foreign-organism nulls may cite D26.1's reasoning, but must decide it afresh under A5.10.

### D18 (RESEARCH_TREE): the secondary, reported-only screened arm C
- **Unchanged and implemented.** The screened arm C is computed in the look by `bench/t1/a5.py` (A5.11.2), reported only. It never changes a verdict, a void or a claim.

### D14, D23 and D25 caveats: kept, and added to, never replaced
- **D14** (arm C at 1.1–1.2% at α = 1%). On T1 the caveat is **made concrete**: every matched-arm-P-FDP number at peptidoform or sequence level is quoted with arm C's mean FDP̂ and LB lower CI for the same list at the same threshold (L-M). The dev statement is kept in the limitations section: "on dev, the natural-proteome (arm C) estimate sat at 1.1–1.2% at α = 1% (D14); 3g found no decoy-realism defect (D18); human paralog neighbours are unmeasured (D18)".
- **D23** (merged block). Any v2_sep per-stratum statement names the merged block (e.g. "merged stratum {2, ≥3}") and reports the raw stratum ≥ 3 figures inside it (dev: 1.18% / 8.0%).
- **D25** (narrow tolerance; peptide-level depth).
  - The wording "compared against Sage chimera **at narrow precursor tolerance**" is mandatory on every v2_sep comparison with C2 or `sage_chim2`.
  - v2_sep's headline is peptide-level depth. No protein-depth statement is made for v2_sep unless S3 is rejected (A3.3.5).
  - D8h's 0.999 [0.951, 1.039] is quoted as the dev expectation next to S3, whatever T1 shows.
- **D25 / E7b** (Sage rank ≥ 2 miscalibration). Quoted wherever C2 or `sage_chim2` ratios appear: "the comparator's rank ≥ 2 IDs were miscalibrated on dev (4–13% FDP at nominal 1%), so these ratios favour the comparator".

### §6.6 sealing mechanics, item 2 (extended)
- **Old:** "They are unlocked only by `bench/sealed_run.py`. The script requires a clean tree, a milestone tag and a ledger line written to `results/SEALED_LEDGER.md` before it runs."
- **Added:** "For T1 and T5 it also enforces one look (§6.0d A5.5, A5.11.1):
  - one joint invocation, with commands and environment taken from `config/a5.yaml`;
  - both prereg tags as ancestors;
  - refusal of any further T1 or T5 line, checked over the ledger's git history, except the A5.5.6 retry."

---

## A5.4 Gate-status record (quoted with every FDR-filtered number; review C.4)

`config/a5.yaml` holds these tables verbatim, and the claim script prints them. Every row names its data.

**v1 `static_gbm_pred` (sha256 31bfe685…)**

| Gate | Data | Result | Source |
|---|---|---|---|
| G0 label blindness | D0/D1 engine + scorer | PASS (bit-identical) | 3a |
| G2-original | A1 null set (6 projects, 894k species nulls, RT-mapped; offset nulls) | PASS, with and without RT; nearly vacuous for a re-scorer (winner changed on 0.08% of nulls) | 3e |
| tail-G2 | A1 null set, clean subset | **FAILED** (n_tail 751; UCB of target fraction 0.742 > 0.60; UCB of diff vs hyperscore 0.124 > 0.10) | 3e, D12 |
| G1 | A1 null set | **FAILED** (ρ 0.048 [0.019, 0.076]) | 3e, D12 |
| G3 | A1 null set, clean subset | **FAILED** (one bin 1.03 pp) | 3e, D12 |
| G4 | dev/training decoys | PASS | 3e |
| G5 | A1 null set | PASS | 3e |
| Arm M (dev) | D7 yeast + human entrapment, r = 1 | inconclusive, not violated (FDP̂_C 1.05%, equal to hyperscore) | 3e |
| §6.7 MC (estimators) | synthetic | PASS | 1c |
| A2 (DN2) | — | **not evaluable before any test** (6/10 projects) | D24 |
| A4 (DN3) | — | **not run** (provenance A5.0.4) | D28 |
| C9 / A6 | removed-peptide human nulls | NO-GO pilot, not run; or A6's verdict; or "A6 not completed (time box)" | A5.1 |

**v2_sep (model sha256 3f2c69cc…, frozen pipeline per A3.1.3)**

| Gate | Data | Result | Source |
|---|---|---|---|
| G0 | D1 + E2001, at the freeze commit | recorded at the freeze (A3.1.4; precondition) | A3.1.4 |
| G2-original | M3 mixture design through the chimeric search (208k null competitions) | PASS | E7 §6 |
| G5 | same | PASS (AUC 0.504) | E7 §6 |
| G4 | stratum × charge partition, training decoys | PASS (weak, τ = 0.298) | D23 |
| Arm M (dev) | D7, chimeric, raw strata | no "violated" in any list | D23 |
| §6.7 MC (full v2_sep competition, pooled picked groups) | synthetic | PASS | D23 |
| tail-G2, G1, G3, V-rows | — | **never tested** (DN2 not evaluable before step 4; A4 not run) | D24, D28 |
| C9 / A6 | removed-peptide human nulls | as for v1 | A5.1 |

---

## A5.5 The look (Part II)

### A5.5.1 Preconditions (all enforced in code by `sealed_run.py --look a5`, A5.11.1)
1. **Tags.** `prereg-v0.4.3` (A3) and `prereg-v0.4.4` (A5) exist and are ancestors of the base commit. In the GO branch, also `prereg-v0.4.5` (A6) and A6's step-5 commit or committed label-free stop, unless the GO time box has expired (A5.1.3).
2. **Milestone.** A `milestone-*` tag on the base commit (the unchanged §6.6 rule).
3. **Freeze.** Clean tree; the A3 freeze commit and tree hashes match (A3.1.3), with the tree list extended by A5.11.6.
4. **Hashes.** `config/a5.yaml` hashes match for the code files and environment manifests (A5.11.5). The C9 record matches the pilot file (A5.1.1).
5. **Ledger.** The ledger, including its git history, has no T1 or T5 line, or the A5.5.6 retry conditions hold.
6. **Verification.** `--verify` (the sha256 of every sealed file) is mandatory.
7. **Exclusive use.** No other agent of this project is active during the look. The node's process list is logged at start and end, and the look runs on an otherwise idle node so that speed is claims-grade (A5.13.9).
8. **Time.** The look runs in a job whose time limit is at least 2× the measured dev `--look` dry run (A5.11.8), recorded in `config/a5.yaml`.

### A5.5.2 One look
- **One invocation.** `python3 bench/sealed_run.py --look a5`. The commands **and the child environment** are read from `config/a5.yaml`, not from the command line or the caller's shell. Nothing is left to choose at look time.
- **Ledger.** One commit writes two lines (T5, then T1) with the same base commit, tag and shared look id, **before** any file is unlocked.
- **Order.** T5 (arm M) runs first, then T1. Both always run. Once the ledger lines exist, not running T1 after T5 is a protocol violation, reported with the results.
- **Isolation.** T5 and T1 have separate fresh work directories and separate fresh prediction caches; they share no state. The dev `--look` dry run (A5.11.8) shows that v1 columns are bit-identical with and without the preceding arm-M step.
- **Nobody reads** any output of either dataset until both have result lines.

### A5.5.3 Look root, console and leak test (A4.8.2.1, A4.8.2.3 and A4.8.1.4 adopted, without encryption)
- **One root, outside the repository.** `<LOOK_ROOT>/<look-id>/`, on shared storage reachable from any node, with path in `config/a5.yaml`, mode 0700. Each dataset has `out/`, `work/` and `logs/` under it.
  - Node-local scratch, if used for speed, is `/tmp/vulcan_look_<look-id>/`. It is deleted by the job epilogue, and by the recover step if reachable.
  - Neither the T1 nor the T5 out dir may exist before the look.
- **Console.** Stdout and stderr go to `<root>/<ds>/logs/`. The terminal shows only lines that match a committed whitelist of line templates (stage names, job counters, exit codes).
- **Leak test.** After the look (and in the dev dry run), a check asserts that no file newer than the look's start exists outside the look root and node-local scratch under `results/`, `data/cache/`, `/tmp/vulcan_t1_*`, `/tmp/vulcan_m3e*` or the repo's work directories. The ledger is the only exception.
- **Commit of outputs.** After both result lines, the operator copies `out/` of both datasets into `results/<date>_T1_sealed/` and `results/<date>_T5_sealed/`, each with its `manifest.sha256`. These are committed in one commit; the sha256 of each manifest is already in that dataset's ledger result line.
- **Derived files.** Sealed-derived intermediates (stores, mzML, MS1 tables, PSM tables, mokapot inputs and outputs, prediction caches of sealed peptides) are **deleted** at the end of each dataset's run: the driver's existing deletion, extended to the whole `work/` except the replicate indexes. The deletion is verified by the leak test. Nothing is kept for later re-analysis. Replicate indexes (built from FASTAs, not sealed data) may be kept.

### A5.5.4 Every branch in that one run, and the completeness rule (review C7; A5 review B1)
**The single run computes everything the claim script can need, whatever the outcome:**
- **v1 primary:**
  - every A1 comparator (Sage (RT on), Sage+mokapot, Sage+pred+mokapot, Comet+mokapot), plus v1 vs C2 and vs `sage_chim2` (A3 Q2);
  - at own q {0.1, 1, 5}% and at matched arm-P FDP {1, 5}%, for every unit (peptidoform, sequence, protein groups);
  - the RT-off Sage control (`--rt-off-control`, as in the 3f sealed command).
- **v2_sep:**
  - family-S inputs on the full list against every A3.3.3 comparator;
  - leave-one-run-out backstops;
  - per-merged-stratum verdicts and merge maps;
  - void flags, the Sage-rank split, and the runs without MS1 hills (all as `bench/t1/a3.py`).
- **A5 additions (A5.11.2):**
  - the stratum-1-only list (A5.5.5);
  - arm C at the matched arm-P threshold (A5.7.1);
  - the A1.8 rank-change fractions;
  - D18's screened arm C (reported only).
- **Entrapment verdicts:** for every method, arm P (15 replicates) at every unit; arm C (5 replicates) at peptidoform and sequence (**arm C is not computed at protein level**); both A1.6 conventions, with (a) deciding; α ∈ {0.1%, 1%, 5%}.
- **Speed and quant**, as the driver records them.
- **T5:** A5.6.
- **Not computed** (fixed now, so it is not an option later): matched arm-P FDP at 0.1%. `VULCAN_T1_FMATCH_0001` is excluded by the environment allow-list, because it was not part of the A3.1.6 dry run. 0.1% is reported at own q only (§6.3: reported, not gated).

**Completeness rule (normative; checked twice).** A look is **failed** unless all of the following hold:
- every method of `A3_METHODS` has all 25 replicates and no `missing_replicates` key;
- all 6 T1 runs are present in every method's run list, and both T5 runs in every arm-M list;
- `meta.a3` lists every A3.3.3 comparator;
- every ratio has `fdp_reached` recorded;
- every bootstrap's `n_boot_nan` is recorded.

A ratio with more than 5% NaN bootstrap draws is "not evaluable" for that comparator. For family S and P1, a not-evaluable comparator counts as **not beaten**: the IUT then fails, so it can never be passed by omission.
- **Check 1.** The driver exits non-zero in sealed mode on any "missing output" event.
- **Check 2.** `claims_a5.py` refuses incomplete inputs and prints "look failed: incomplete".

**After the look:**
- There is no re-search, no re-run of any stage, no new estimator, no new comparator and no re-weighting. Any further analysis of T1 or T5 data is a protocol violation.
- *Descriptive* tables computed from the committed `metrics.json` by code that did not exist at the tag are allowed only in a section titled "post-hoc, not pre-registered". They never appear in the abstract and never change a claim.

### A5.5.5 Stratum-1-only list (descriptive; definition fixed now, from A4.8.1.4)
- The list is the v2_sep accepted units of merged stratum 1 at a common `q_sep` cutoff.
- Matched arm-P FDP uses that list's own paired counts (target, entrapment and their decoys restricted to stratum-1 pair winners), with `evaluate.py`'s first-crossing rule.
- Ratios, bootstrap and the leave-one-run-out minimum are as in A3.3.4, against every A3.3.3 comparator.
- It is reported in a table titled "v2_sep stratum-1 list (descriptive; not family S)". It is never tested and never quoted as a claim.

### A5.5.6 Failure, interruption and retries (A5 review B2)
- **What counts as a failed look.** A look is **failed** if any of these holds:
  - either command exits non-zero;
  - the completeness rule fails;
  - the leak test fails;
  - **either ledger line lacks a result line** (SIGTERM, SIGKILL, node loss or time limit).

  A failed look is **consumed** and counts towards the retry limit. A deliberate interruption by the operator is a protocol violation and is reported.
- **Recover.** `sealed_run.py --look a5` first checks for an unterminated look id, and refuses until `sealed_run.py --recover <look-id>` has run. That step:
  - writes a "failed (unterminated)" result line;
  - deletes the look root and any reachable node-local scratch;
  - commits.

  The recover step reads nothing in the root except file names, which it does not print.
- **Failure summary.** It contains only the dataset, the stage name, the exception type and the traceback's file:line:function frames, with no messages and no locals. It is produced by the wrapper from the logs and committed. All outputs of a failed look, of both datasets, are deleted unread.
- **Retry.** At most two, each `sealed_run.py --look a5 --retry <look-id>`. All of the following are required:
  - an engineering fix that touches no analysis parameter, developed only from the committed frames and from reproductions on D0/D1/D7. If the failure cannot be reproduced on dev data, the fix may be made from the committed frames alone, with a T1-blind review stating that it does so;
  - a T1-blind review of the fix;
  - a dated deviation in `config/a5.yaml`;
  - the A3.1.6 bit-identity dry run and the dev `--look` dry run re-passed;
  - a new `milestone-*` tag on the fix commit;
  - if the fix touches a frozen tree, the A3 freeze is re-taken first: a new freeze commit, `prereg-v0.4.3.x`, and then `prereg-v0.4.4.x`.
- **What a retry runs.** A retry re-runs **both** datasets from fresh roots. **The last successful look is the result of record.** Every earlier ledger line and failure summary is reported.
- `sealed_run.py` refuses a retry if any earlier look has two exit-0 result lines, and refuses a third retry.

### A5.5.7 Who
- **Operator:** any agent that is not the author of A5, of the claim script or of the claim script's review. The operator commits outputs and runs the claim script once; its output is committed unedited.
- **The adopter and the claim-script reviewer** state their T1 exposure (expected: none) at the time of their act.

---

## A5.6 Arm M on T5 (same look; fixed now)

- **Data.** T5 = PXD052597, both sealed runs (Q Exactive HF-X), as one experiment. Yeast is absent from scorer and AlphaPeptDeep training (DESIGN §6.6; D4 manifest check).
- **Code.** `bench/t5_armM.py` (A5.11.3) is a parametrised module, shared by dev and sealed modes. It carries the dev arm-M design of `scripts/m3e/armM_prepare.py` and `armM.py` (v1) and of `scripts/a3/dev_armM.py` (v2_sep).
- **Databases:**
  - targets: `data/fasta/yeast_UP000002311.fasta` + `contaminants_frankenfield2022.fasta`. No spiked protein is added, and T5's sample sheet is not read;
  - entrapment: human Swiss-Prot peptides, neighbour-filtered as arm C (`armM_prepare.py`), restricted to the effective search space and subsampled to r = 1, one per **sealed** entrapment seed 2004–2006;
  - decoys: the 5 DECOY_SEEDS (C-term-fixed shuffles), giving **15 replicates** (5 decoy × 3 entrapment), the §6.3 milestone design. Dev used decoy seed 1001 only.
- **Estimator.** Combined, FDP̂_M = N_E(1+1/r)/(N_T+N_E), with **r = the peptide-level effective r** (as dev: 1.0 by construction) for peptidoform and sequence units. §6.3 verdict table. **CI convention (a) decides verdicts and voids**; (b) is reported.
- **Pipelines:**
  - **v1:** closed search, `report_psms 5`, scored through **`bench/t1/ours.py` `Scorer`** (the T1 path, not `scripts/m3/features.py`), paired double competition;
  - **v2_sep:** chimeric search, frozen v2 model through `bench/t1/ours_chimeric.ScorerV2`, `q_sep`;
  - **hyperscore** through each pipeline, as a reference (reported only).
- **Lists:**
  - v1: overall.
  - v2_sep: overall, plus each **merged** stratum (the frozen label-blind ≥ 200-pair rule, as on T1). **Merged strata decide voids**; raw strata {1, 2, ≥3} are reported (as dev).
  - Units: peptidoform decides voids; sequence is reported. There is no protein-level arm M.
  - α ∈ {1%, 5%} gated; 0.1% reported.
- **Not evaluable** (α·N_T < 20) is allowed and labelled.
- **Void:** A5.7. **A T5 arm-M verdict is final** (F13). The hyperscore reference is reported, never used to excuse a violation.
- **Dev equivalence** (A5.11.3): run in dev mode on D7 with seeds 2001–2003 and decoy seed 1001, `bench/t5_armM.py` reproduces:
  - the committed v2_sep raw-stratum table (`results/2026-10-28_A3_preconditions/armM.json`) exactly;
  - v1's 3e arm-M table. Exactly if the `ours.py` path is bitwise equal to `scripts/m3/features.py` on D7; otherwise the difference is committed and explained before the freeze.
- **Consequence.** T5 is spent for every scorer after this look. A future scorer (e.g. v3) needs a new sealed arm-M set, added by a dated note before its first look.

---

## A5.7 Void rules (outcome-independent; each check can only remove claims)

### A5.7.1 v1
All FDR-based v1 claims are **void, recorded as final**, if any of the following holds. FDR-based claims are: non-inferiority at matched FDP or own q, IDs at any q, ratios, and quant on q-filtered lists.
1. Any "violated" verdict on T1 under convention (a), at 1% or 5%:
   - arm P in the peptidoform, sequence or protein-group list;
   - arm C in the peptidoform or sequence list.
2. Any "violated" arm-M verdict on T5 (peptidoform, 1% or 5%, convention (a)).
3. v1 FAILED a completed null-spectrum test after A1 and before the look (in the GO branch, an A6 fail). A1's own failure does not void: it is the premise of A5.

**Claim-level void (in addition; defined exactly):**
- **(a) Which ratios.** Every matched-arm-P-FDP ratio with v1 as **numerator** at peptidoform or sequence level and f ∈ {1%, 5%}.
- **(b) The threshold.** The point-estimate matched threshold `ours_q_at_match`: the first crossing of the 15-replicate mean arm-P FDP̂ above f, interpolated on the code's grid. It is applied to each of the 5 arm-C replicates' own q.
- **(c) The quantity.** Arm-C LB = N_E/(N_T+N_E) per replicate at that threshold, with mean and 95% CI under convention (a). r = `arm_C_r[unit]` is used only for the reported FDP̂_C.
- **(d) The rule.** The ratio is **void** if the lower CI of the mean LB is > f.
- **(e) Lists.** The overall list only, and only at matched thresholds. Own-q thresholds are covered by rule 1.
- **(f) What is printed.** The same quantities are printed for the comparator side of every ratio, but they never void (comparators carry their own verdicts).

### A5.7.2 v2_sep
- **Unchanged from A3.3.6:** any "violated" in any arm, in any peptidoform list (overall or merged stratum), at 1% or 5%; a protein-group "violated" (arm P) stops family S at S3.
- **Plus:**
  - any "violated" arm-M verdict on T5 in any v2_sep merged-stratum or overall list;
  - v2_sep FAILED a completed null-spectrum test before the look;
  - the claim-level matched arm-C void of A5.7.1, applied to v2_sep as numerator on the overall list.

### A5.7.3 What a void leaves
- The voided numbers are still printed in full, in a table titled "VOID under §6.0d A5.7 (reason)".
- Speed and the entrapment verdicts themselves are not void.
- Comparator rows are never void. Comparators have their own verdicts, which are reported.

### A5.7.4 Power of the void rules (stated so they are not over-read)
"Violated" needs the lower CI of LB above α. At r ≈ 1, LB ≈ FDP̂_C / 2. So the void rules fire only against a **gross** violation: true FDP of roughly ≥ 2α in arm C, and similarly for arm M. They are a backstop, not a validation. This sentence is part of the limitations text.

### A5.7.5 Finality (F13)
- No explanation found after the look un-voids a claim or changes a verdict. That includes sample composition, contamination of T1 or T5, behaviour of the hyperscore reference, and comparator failures.
- Such explanations may be written only in a "post-hoc, not pre-registered" section, with the word "unconfirmed".

---

## A5.8 Family S (summary of the fixed rule)
A3.3.5 is unchanged, except:
1. the list is the full v2_sep list (F1);
2. the title is "Secondary arm; not gate-validated";
3. it is decided by the claim script;
4. a not-evaluable comparator counts as not beaten (A5.5.4);
5. the level is 0.025 one-sided within S, as A3.3.5. The study-wise error across the primary and S families is up to the sum of the levels (disclosed, unchanged);
6. every S row carries L-M or L-Q as appropriate (S1 carries both; S5 carries L-Q).

---

## A5.9 Claim and label tables

### A5.9.1 Label texts (exact; held in `config/a5.yaml`; inserted by the claim script)

**L-G (gate status, long form, per scorer).** For v1:
> "Null-spectrum gates: A1 (2026-10-27): tail-G2, G1, G3 FAILED; final for those nulls. A post-hoc contamination explanation exists, is unconfirmed, and is not used. A2 (DN2): not evaluable before any test (D24). A4 (DN3): drafted and reviewed twice, then not run by PI decision D28 (2026-10-28), made before any DN3 listing and before any T1 look, on label-free feasibility grounds: P(v1 evaluable) 0.08–0.43 in the reference cell (grid 0.01–0.98). {C9 clause} Passed: G0 (D0/D1), G4 (dev decoys), G2-original and G5 (A1 null set; G2-original nearly vacuous for a re-scorer), §6.7 estimator simulations; arm M (dev D7): inconclusive, not violated."

- **{C9 clause}** is one of:
  - "C9 removed-peptide pilot: NO-GO (P(v1 minimum) = x < 0.8, commit c); not run."
  - "C9 pilot: no verdict within the time box; not run."
  - "A6 (removed-peptide human nulls, test 2 at α/2): PASS | FAILED | not evaluable | not completed (time box) (commit c), before T1."
- **For v2_sep:**
  > "v2_sep: no null-spectrum test completed (DN2 not evaluable before any test; A4 not run) {C9 clause for v2_sep}. Passed: G0 at freeze (D1), G2-original and G5 (M3 mixtures, chimeric), G4 (weak, training decoys), §6.7 MC; arm M (dev D7): no violated list."

**L-E (entrapment, per list).**
> "arm P: <verdict>; arm C: <verdict | not computed at protein level>; arm M (T5): <verdict | not applicable to this unit> (<unit>, α = <α>)"

Verdict words are only those of §6.3: supported, supported-noisy, inconclusive (+ "consistent with control" when the paired lower CI ≤ α), violated, not evaluable. "Consistent with control" never ranks scorers.

**L-M (matched-entrapment-FDP comparison).**
> "compared at matched arm-P (shuffled paired entrapment) FDP̂ = f; each method's own q is used only to order its list, not as an error estimate; the comparison rests on the entrapment assumptions; arm P cannot detect a preference for natural sequences; arm C for the same list at the same threshold: mean FDP̂ = x%, LB lower CI = y% [protein level: arm C not computed]."

**L-Q (own-q number).**
> "at the method's own q; estimator validity for this scorer is not established by null-spectrum tests (A1 FAILED)."

**L-R (A1.8).**
> "top-K re-ranking by the static scorer is used although A1.8's precondition (tail-G2/G1/G3 pass) is not met (rank-1 changed in p% of spectra, q% of accepted units)."

In the GO-pass state, A1.8's met-precondition text is used instead (A5.3).

**L-S (short form for prose; the minimum).**
> "FDR-filtered at the method's own threshold and checked by entrapment (verdicts as stated); not validated by null-spectrum tests: the only completed test (A1) FAILED."

In the GO branch with an A6 verdict, the last clause becomes "A1 FAILED; A6 <verdict>". **Every short form contains the word "FAILED".**

**L-2 (v2_sep).**
> "pre-registered secondary arm; not gate-validated; compared against every comparator in the look, including Sage chimera at narrow precursor tolerance; per-stratum statements name the merged block."

**L-D (designations).**
> "v1: engine default, FDR not validated by null-spectrum tests. v2_sep: candidate. No OT-HCD scorer has gate-validated FDR claims."

**L-P (post-T1 pass, A5.10.6).**
> "null-spectrum gates passed after T1 was seen (<amendment>, test k at α/2^(k−1))"

### A5.9.2 Claim table (one row per claim; every row is printed whatever its outcome)

| # | Claim | Statistic and bar | Labels | Void if |
|---|---|---|---|---|
| P1 | **v1 no loss (§1)**: the conjunction of (i) and (ii) | (i) peptidoform, matched arm-P FDP 1%: LB ≥ 0.97 **and** leave-one-run-out minimum ≥ 0.97 vs **every** A1 comparator; (ii) peptidoform, own q 1%: LB ≥ 0.97 vs every A1 comparator whose own-q verdict is "supported" (§1). If none is supported, (ii) prints "no supported comparator" and P1 rests on (i), which is stated. Both rows are always printed | L-S (headline), L-M on (i), L-Q on (ii), L-E, L-R, L-G | A5.7.1 |
| P2 | v1 vs each comparator at 5%, at sequence and protein-group level, and at 0.1% (own q) | reported ratios with CI | L-M or L-Q, L-E, L-R | A5.7.1 |
| P3 | v1 IDs | counts at own q {0.1, 1, 5}% per unit | L-Q, L-E, L-R | A5.7.1 |
| P4 | v1 vs C2 and `sage_chim2` | reported in the primary section (A3 Q2), not a claim | L-M, the narrow-tolerance phrase, the E7b rank caveat | A5.7.1 |
| P5 | v1 entrapment verdicts | per arm, unit, α, both conventions | literal | — |
| P6 | v1 arm M (T5) | per §6.3 | literal | — |
| P7 | v1 quant (ε, CV, DE pAUC) | on q-filtered lists | L-Q | A5.7.1 |
| P8 | speed (core-s, runs/h) | §2 | none needed; if expressed per ID, L-Q | — |
| S1–S5 | v2_sep family S | A3.3.5 on the full list | title "Secondary arm; not gate-validated"; L-2; L-M and/or L-Q; L-E; L-R; L-G (v2_sep) | A5.7.2 |
| R1 | v2_sep stratum-1-only list | descriptive (A5.5.5) | "descriptive; not family S", L-2 | — |
| R2 | v2_sep per-merged-stratum verdicts | A3.3.6 | literal; merged block named; raw stratum ≥ 3 shown (D23) | — |
| R3 | v2_sep/v1 ratio, per-stratum IDs, sequence level | reported | L-M or L-Q, L-2 | A5.7.2 |
| R4 | D18 screened arm C | reported only | "secondary arm C (D18 screen); reported only" | — |
| D | designations | fixed (F3) | L-D | — |
| G | gate-status tables | A5.4 | verbatim | — |

The leave-one-run-out minimum in P1(i) mirrors A3.3.4's backstop for S. It is a stricter addition for v1, fixed now.

### A5.9.3 State-dependent wording (closed list)

| State at the look | v1 headline clause | v2_sep title |
|---|---|---|
| NO-GO, or no pilot verdict in the time box | L-S as written | "Secondary arm; not gate-validated" |
| GO, A6: v1 PASS | L-S's last clause becomes "passed null-spectrum gates at its second completed test (A6, α/2) before T1, after A1 FAILED (…)"; L-Q's last clause becomes "validated by A6 (test 2), after A1 FAILED"; L-R uses the met-precondition text; designation per A5.10.6 | unchanged, unless v2_sep also passed A6 on the full list (A5.10.6). Then "Secondary arm; null-spectrum gates passed under A6 (first test)" (PI point A5.14.3) |
| GO, A6: v1 FAILED | all v1 FDR-based claims void (A5.7.1.3); L-S becomes "…: A1 FAILED; A6 FAILED" | per v2_sep's own A6 status |
| GO, A6 not evaluable, label-free stop, or time box | as NO-GO; the C9 clause names A6's status | as NO-GO |

### A5.9.4 Placement (review C.8)
- **Where labels go.**
  - The L-S label (or its A5.9.3 variant) is in the REPORT.md abstract **and** in the headline sentence.
  - The full L-G, L-E, and L-M or L-Q labels are in **every table caption** with FDR-filtered numbers, including protein groups and quant metrics computed on q-filtered lists.
  - Labels never appear only in a footnote or a limitations section.
- **The headline sentence does not start with a gain figure or any number.**
  - The first T1 sentence states the gate status.
  - The second states P1's verdict with its CI against **the A1 comparator with the lowest P1(i) point ratio** (fixed now), including a loss if there is one (e.g. the D17 dev pattern, 0.990× vs Sage+pred+mokapot).
- **Order.** v1 is reported first and in full; v2_sep follows under its title.
- **Release.** External release of T1 numbers (slides, abstracts, preprints) carries L-S at minimum and passes the A5.11.7 lint.

---

## A5.10 Future null-set tests (S6; binding from adoption, while the PI is T1-blind)

### A5.10.1 Attempt ledger (replaces A4.6's ledger)

| Scorer | Completed tests | Not tests (label-free stops / not run) | Next completed test |
|---|---|---|---|
| v1 | 1: A1 (FAILED) | A2 (not evaluable, label-free), A4 (not run) | **test 2, level α/2** |
| v2_sep | none | A2/A3 tiers (not evaluable), A4 (not run) | test 1, level α |
| ce_s | none | as v2_sep | test 1, level α |

A completed test is one where step 5 ran for that scorer, so that some per-label quantity was computed, whatever the verdict (the A4.6 definition).

### A5.10.2 Level rule (A4.6, adopted by reference)
- **Levels.** Completed test k uses level α/2^(k−1) on every re-tested row (tail-G2 one-sided; G1 two-sided; TOST rows two-sided). The rule is the same for all scorers.
- **Sequencing.** After a completed test without a pass, the next test uses the next level. After a label-free stop, the level is unchanged.
- **Gates not re-tested.** G2-original and G5 keep their levels.
- **No discounting of A1.** Ignoring A1 because of its post-hoc explanation is not permitted.

### A5.10.3 Statistical floor (a later amendment may only tighten it)
Every later null-set amendment (A6, or any after T1) uses at least:
1. every A1 margin and n-minimum; A2 P1–P5; A2's validity controls (hyperscore tail fraction ≈ ½); A3's V-rows for v2_sep (with tiers as in A3.2.6); and intersection–union within a scorer;
2. **A4.5.2:** the outer envelope of the analytic, project-unit and cluster-unit intervals for every bound-based row (B = 10,000, `default_rng(83)`);
3. **A4.5.3:** n_min under the A4 rule. Final values are re-derived with ≥ 4,000 replicates, p̂ ≥ 0.92 and base-rate heterogeneity of logit-sd 0.5, and are never below the A4.5.3 provisional table. Fewer than 10 clusters means every clustered row is not evaluable;
4. **a cluster definition at least as coarse as A4.3.4**: the submitting lab or PI, with order-insensitive person keys, extended to any shared-source key the design introduces (e.g. the originating project of human donor spectra);
5. **A4.3.6:** the cluster cap and the capped/uncapped intersection–union rule;
6. **a label-free stopping rule** of the A4.3.7 type, committed before any listing;
7. **A4.5.5:** the subset check. For any instrument or source subset that reaches n_tail ≥ 200, the point estimates must not exceed the tail-G2 margins;
8. **A4.5.6 and A4.9:** reporting (per cluster, per tier, both interval components, leave-one-cluster-out);
9. exclusions at least as broad as A4.3.3's list, plus: T1–T5 projects and runs; train and dev sets; predictor-training accessions; every earlier null set; and the C9 pilot's runs;
10. a label-free evaluability simulation committed before any listing, quoting the reference cell and the grid;
11. a closed text (A3 Q6): no option left open after adoption.

Selection rules that are specific to foreign organisms (A4.3.2 lineage classes, marker screens, tier 0) are not part of the floor. A new design must justify its own null construction and simulate its validity risks.

### A5.10.4 Who designs and runs it
- **T1-blind** means the agent has not read, and has not been shown:
  - any file under `results/*_T1_sealed*` or `results/*_T5_sealed*`, or any look root;
  - the T1/T5 sections of REPORT.md, or the claim script's output;
  - any RESEARCH_TREE entry or status-summary text written after the look's ledger line;
  - any commit message after that line;
  - any message containing T1/T5 numbers or verdicts.
- **Running T1-blind agents.** T1-blind agents run in a clone or worktree at a **pre-look commit**, plus only the specific later files they need, listed in their commissioning prompt.
- **Repository hygiene after the look.** Commit subjects and the RESEARCH_TREE status summary contain no T1/T5 number or verdict; results are referenced by path only.
- **Which roles must be T1-blind.** The **author, reviewer(s), selection/walk agent, step-4/4b/5 agents and deviation reviewers** of any null-set amendment adopted after the look. Each artefact carries "T1 exposure: none".
- **Commissioning prompts.** The prompt from a T1-exposed PI to these agents is committed verbatim with the artefact. It contains no T1/T5 information.
- **Adoption.** A T1-exposed PI may adopt such an amendment **only verbatim, or with changes reviewed by a T1-blind reviewer**. The PI may not choose among options after T1; the text must arrive closed.
- **No draft shopping.**
  - At most **one** null-set amendment per scorer may be open at a time.
  - Every commissioned draft (its prompt, text, reviews and fate) is committed and reported with T1.
  - Abandoning a draft is reported, with its reason.
  - A rejected draft may be followed by a new one only after the rejection, and the reason for it, are committed.

### A5.10.5 What a later result does
- A **fail** is final for that scorer under P4, and is reported with T1.
- **Not evaluable** keeps FDR-control claims blocked. The next attempt follows A5.10.2.
- A **pass** relabels T1 numbers only by adding L-P ("null-spectrum gates passed after T1 was seen"). It never deletes L-G's history.

### A5.10.6 Designation rule after a pass (the only route that changes F3)
- **v1** becomes "engine default; null-spectrum gates passed at test k (<amendment>)" iff:
  - it passes tail-G2, G1 and G3 at the A5.10.2 level; and
  - it has no void under A5.7.1.
- **v2_sep** becomes the OT-HCD production default iff all of:
  - it passes the tier-2-equivalent rows on the full list (in the A3.2.4/A3.2.6 sense, under the A5.10.3 floor);
  - A3.1.8 holds with no "violated";
  - it has no void under A5.7.2.

  A tier-1-only pass licenses statements on the stratum-1 list only, and does not change "candidate".
- Whenever a designation is named, the gate history (L-G) is quoted, and L-P too if the pass came after the look.

---

## A5.11 Implementation (everything in frozen trees lands **before** the A3 freeze commit)

**Ordering** follows A4.4's rule. The A3 freeze hashes `vulcan/`, `bench/`, `scripts/a3/` and `scripts/e7/`, and A5.11.6 adds further trees. So:
1. every A5 change below that lives in those trees is committed and tested first;
2. then the A3.1.6 dry run and the dev `--look` dry run are re-done;
3. then the freeze commit and `prereg-v0.4.3`;
4. then `config/a5.yaml` and `prereg-v0.4.4`.

The D27 GPU MS1 merge waits until the look is committed.

### A5.11.1 `bench/sealed_run.py`: one look in code
1. **`ONE_LOOK = {"T1", "T5"}`.** The one-look check reads the current ledger table **and** `git log -p -- results/SEALED_LEDGER.md`. Any row ever written for a ONE_LOOK dataset makes a new run **refused**, except an A5.5.6 retry.
2. **`--look a5`:**
   - reads `config/a5.yaml` (`look.commands.{T5,T1}`, `look.purpose.{T5,T1}`, `look.env`, `look.root`, `look.time_limit`);
   - refuses any trailing command;
   - builds the child environment **only** from `look.env` (A5.11.5), dropping every other variable;
   - writes both ledger lines in one commit with a shared look id;
   - runs T5, then T1, each under its own root;
   - filters the console through the whitelist;
   - writes result lines carrying the exit code, wall time and the sha256 of the output `manifest.sha256`;
   - runs the leak test.
3. **Plain runs refused.** `--dataset T1` or `--dataset T5` without `--look a5` is refused.
4. **Tag checks.** The base commit (after skipping ledger-only commits) carries a `milestone-*` tag **and** has `prereg-v0.4.3` and `prereg-v0.4.4` as ancestors (`git merge-base --is-ancestor`). In the GO branch it also needs `prereg-v0.4.5`, unless the time box has expired. Both are decided from the C9 record, which is recomputed from the pilot file (A5.1.1).
5. **Hash checks** of every file and environment manifest listed in `config/a5.yaml`.
6. **`--verify`** is forced for ONE_LOOK datasets.
7. **Look id.** Exports `SEALED_LOOK_ID` to the commands.
8. **Failure handling** (A5.5.6):
   - the summary writer (frames only);
   - deletion of the roots on failure;
   - `--recover <look-id>` for unterminated looks;
   - `--retry <look-id>` with its conditions, a count ≤ 2, and a new milestone tag.
9. **The SIGTERM handler** writes a "failed (terminated)" result line before exiting.
10. **Tests** (`tests/bench/test_sealed_run_one_look.py`, in a temp git repo via `--root`, with dummy sealed files and dummy commands):
    - a second look is refused, including after the ledger file has been rewritten (history check);
    - plain `--dataset T1` is refused;
    - ledger-only commit skipping still finds the tag;
    - a missing prereg tag is refused;
    - a hash mismatch is refused;
    - a C9 record that mismatches the pilot file is refused;
    - a caller-set `VULCAN_SAGE` is dropped;
    - both ledger lines are written before unlock;
    - an unterminated look blocks a new look until recover;
    - SIGTERM writes a failure line;
    - a retry is allowed only after a failed look, is refused after a successful one, and is refused at the third attempt;
    - the leak test detects a planted file.

### A5.11.2 `bench/t1_driver.py`: sealed hardening, and `bench/t1/a5.py`
1. **Refused arguments in sealed mode:**
   - non-default `--stages`, `--methods`, `--ours`, `--mzml-cache` and `--keep-sealed-derivatives`;
   - any `--out` or `--work` that is not under `$SEALED_LOOK_ROOT`, or that exists and is non-empty (**no resume in sealed mode**).

   `SEALED_LOOK_ID` must match the open ledger line, and both prereg tags must be ancestors of HEAD (asserted in the driver too, as defence in depth).
2. **The environment** is asserted **after** the re-exec in `staged_python()`. Every `VULCAN_*` variable must equal `config/a5.yaml look.env`, except `VULCAN_T1_STAGED=1`, which the driver sets itself. Analysis switches such as `VULCAN_T1_FMATCH_0001` and `VULCAN_T1_NO_LABEL_CACHE` are absent from the allow-list and therefore refused.
3. **`--rt-off-control`** is forced in sealed mode (it is in the pre-registered 3f command).
4. **The completeness rule** (A5.5.4): a sealed run exits non-zero on any missing output.
5. **`bench/t1/a5.py`** (additive, called after `a3.secondary`; nothing is re-searched) computes:
   - the stratum-1-only list (A5.5.5);
   - the claim-level arm-C quantities of A5.7.1 (a)–(f) for every numerator and comparator;
   - the A1.8 rank-change fractions, label-free on the prod replicates:
     - v1: the best rescored PSM's engine rank > 1;
     - v2_sep: accepted units whose PSM had engine rank > 1 within its group;
   - **D18's screened arm C** (reported only): entrapment peptides that share a ≥ 4-ion b/y ladder with any human target, or that come from plastid, photosynthesis or food-taxon proteins (3g remedy 1 lists), are excluded. The arm-C estimate is then recomputed on the 5 arm-C replicates.
6. **Deletion.** At the end of a sealed run, the driver deletes the whole `work/` except the replicate indexes (A5.5.3), and writes `<out>/manifest.sha256`.
7. **Dry run.** The A3.1.6 dry run on D0 and D1 is re-done with A5's code. v1 columns must be bitwise identical with and without the A3+A5 additions, and Sage-based comparators must be within the documented non-reproducibility. Committed.
8. **Test:** a dev run with one replicate output deleted must exit non-zero, and the claim script must refuse it.

### A5.11.3 `bench/t5_armM.py` (new, parametrised) and a T5 entry
1. **Interface.** It takes runs, seeds, decoy seeds, root and mode as arguments, writes nothing outside its root, and uses its own fresh prediction cache. The dev scripts `scripts/m3e/armM.py` and `scripts/a3/dev_armM.py` are refactored to call it (their committed outputs are not overwritten). Sealed-mode refusals are as for the T1 driver.
2. **Contents.** A5.6 exactly: v1 through `bench/t1/ours.py Scorer`, v2_sep through `ours_chimeric.ScorerV2`, merged strata deciding, raw reported, r peptide-level, convention (a) deciding, 15 replicates in sealed mode.
3. **Dev equivalence** (A5.6), committed.
4. **Outputs:** `armM.json` and `manifest.sha256`.

### A5.11.4 `bench/t1/claims_a5.py`: the claim script
1. **Inputs.** It is a pure function of:
   - T1 `metrics.json` (+ manifest);
   - T5 `armM.json` (+ manifest);
   - `config/a5.yaml` (the A5.4 gate tables, the C9/A6 record, the label texts);
   - D23's committed A3.1.8 records.

   It verifies every input's sha256 against the ledger result lines, and it applies the completeness rule.
2. **Outputs.** `claims.json` and `claims.md`, containing:
   - the A5.9.2 table, with every void and its reason;
   - family S decisions (sequence, comparators, LB and leave-one-run-out);
   - own-q eligibility, the S3 stop, and the P1 conjunction;
   - the gate tables, with the labels inserted verbatim and the A5.9.3 state;
   - the headline comparator choice (A5.9.4);
   - every reported-only table.
3. **No new decisions.** It decides nothing that is not written in A5 or A3. Where code and text disagree, the stricter reading wins, and this is recorded as a dated deviation **before** the look.
4. **Tests.**
   - Synthetic `metrics.json` and `armM.json` inputs cover every branch:
     - each void trigger, and the claim-level void at its boundary;
     - S stopping at each S_k, and the S3 skip;
     - a not-evaluable comparator;
     - no supported comparator;
     - not-evaluable strata;
     - each A5.9.3 state;
     - incomplete inputs.
   - An end-to-end run on the D0/D1 dry-run outputs and the D7 arm-M output.
   - The lint of A5.11.7, applied to `claims.md`.
5. **Review.** A T1-blind agent reviews the script and its tests before `prereg-v0.4.4`. The review is committed, and its sha256 is recorded in `config/a5.yaml`.

### A5.11.5 `config/a5.yaml` (outside frozen trees; committed before `prereg-v0.4.4`)
- **`look`:**
  - commands for T5 and T1 (full argument lists);
  - purpose strings (T1's names the v1 primary and the v2_sep secondary arm);
  - `env`: the complete allow-list with explicit values, covering `PATH`, the thread counts, `VULCAN_SAGE`, `VULCAN_MOKAPOT_ENV`, `VULCAN_CORE_DIR`, `VULCAN_TRFP_BIN`, `VULCAN_RAW_DLL`, `VULCAN_DOTNET_ROOT`, `VULCAN_LOCAL_ENV`, `VULCAN_THREADS` and `VULCAN_COHORT_TIMEOUT_SCALE`;
  - root, time limit, and the console whitelist.
- **`hashes`:** `sealed_run.py`, `t1_driver.py`, `bench/t1/*.py`, `t5_armM.py`, `claims_a5.py`, the claim-script review, and **environment manifests**. The manifests are sha256 lists of the Sage binary, the Comet binary (`conda_envs/comet`), TRFP, the mokapot package env, and the staged pred env. The staged env is rebuilt from scratch for the look, and the look refuses on any mismatch.
- **`gates`:** the A5.4 tables.
- **`labels`:** the A5.9.1 texts.
- **`c9`:** the pilot verdict file commit and sha256.
- **`deviations`:** the retry deviations (empty at the tag).
- **Mutable after `prereg-v0.4.4`** (each by a dated commit that touches nothing else): `c9`, `deviations`, and the hashes that a re-freeze changes. A re-freeze also requires `prereg-v0.4.4.x`.
- **Test:** a test pins the label texts to this document.

### A5.11.6 Freeze extension (stricter)
`config/a3.yaml freeze.tree_hashes` adds:
- `scripts/m3/`, `scripts/m3e/`, `scripts/entrapment/`, `scripts/e7b/` (T5 and the comparators import from them);
- `scripts/e9/` (the env stager);
- `scripts/baselines/` (the mokapot env builder).

`vulcan.gates.a3_config.freeze_status` checks them.

### A5.11.7 Wording lint (`tests/report/test_report_labels.py`)
It runs over `claims.md`, REPORT.md and any file under `release/`. It is run before any release, and a failing lint blocks the release. It checks:
- **Banned phrases**, outside quoted label text: "entrapment-validated", "entrapment-verified", "empirically verified", "validated FDR", "FDR control" (except in "FDR-control claim" inside A5.9 labels), "controls FDR", "controlled at", "FDR-controlled" (except "not FDR-controlled"), "gate-passed", "calibrated" (when applied to v1/v2_sep FDR), "contaminated nulls".
- **"validated"** appears only in the negated forms "not validated" and "not gate-validated", or in the exact A5.9.3 pass phrases.
- **Labels present:**
  - L-S is in the abstract and in the headline sentence;
  - every table caption with FDR-filtered numbers carries its L-labels;
  - "at 1% FDR"-style phrases are followed by an L-M or L-Q tag.
- **The headline sentence** does not start with a number.

### A5.11.8 Dev `--look` dry run (before the freeze)
- **Setup.** `sealed_run.py --look a5 --dev` runs in a temp clone with dev manifests. D7 stands in for T5, and D1 (and, timed, D0) for T1. It goes through the real `t5_armM.py`, the real driver and the real claim script.
- **It shows:**
  - the console whitelist holds;
  - the leak test passes;
  - outputs are complete;
  - v1 columns are bit-identical with and without the preceding arm-M step;
  - the recover and retry paths work.
- **Record.** The run is committed, and its wall time sets `look.time_limit`.

---

## A5.12 Order of operations
1. This draft, its independent review and the revision log are committed.
2. **Adoption.** The PI adopts A5 as DESIGN §6.0d, with a dated RESEARCH_TREE note. The note contains:
   - "A5 is decided knowing that v1 failed its only completed null-spectrum test (A1)";
   - the adopter's exposure to T1, T5, DN3 and the C9 pilot;
   - who controls the pilot's commit;
   - D26.3's reversal;
   - the PI's answers to A5.14.

   **A5.10 binds from here.** In the GO branch, A6 drafting may start now.
3. **Implementation and dry runs.** A5.11.1–A5.11.4 and A5.11.6–A5.11.7 are committed and tested. The dry runs (A5.11.2.7, A5.11.3.3, A5.11.8) and the claim-script review are committed. In the GO branch, A6 code that cannot live outside the frozen trees lands here (A5.1.4).
4. **The A3 freeze.** The remaining A3.1 preconditions (A5.3): freeze hashes, G0 at the freeze, E7b commit, driver commit. Then the A3 freeze commit and `prereg-v0.4.3`.
5. **The A5 tag.** `config/a5.yaml`, then `prereg-v0.4.4`.
6. **The C9 record** (A5.1.1). In the GO branch, A6 runs under A5.10 with tag `prereg-v0.4.5`, to completion, label-free stop, or the time box.
7. **The look.** A `milestone-*` tag, then the look (A5.5). The outputs are committed, the claim script runs once, and its output is committed unedited and copied verbatim into REPORT.md. The lint passes.
8. **Decision.** A RESEARCH_TREE decision records every claim row and void, referenced by path. The commit subject carries no numbers (A5.10.4).

---

## A5.13 What is lost, and residual risks (stated so A5 is not read as stronger than it is)
1. **No null-spectrum evidence for v1's tail at the look** (else-branch). The only completed test FAILED. A1's margins already left a +0.05 tail blind spot. Tail FDR now rests on entrapment arms P and C, whose void rules fire only on gross violations (A5.7.4).
2. **v2_sep has never had a null-spectrum test.** Its chimeric secondary PSMs are checked only by per-stratum entrapment, and stratum ≥ 3 was marginal on dev (D23). "Candidate" is permanent until a pass.
3. **Arm P cannot see a natural-sequence preference.** Arm C can, but only at r ≈ 1 with list-level power, and not at protein level. Human paralog neighbours are unmeasured (D18).
4. **T1 and T5 are spent** as confirmatory sets for every scorer, present or future, whatever any later test shows (as A4.8.0 had accepted).
5. **Every later null-set design is T1-exposed at the PI level.** A5.10 limits this through pre-commitment, T1-blind agents at pre-look commits, and the no-draft-shopping rules. It cannot remove it: the PI still decides whether to commission a test at all.
6. **Controls are procedural, not cryptographic.** The following are rules, auditable through the ledger and commit history, not technical guarantees:
   - "nobody reads outputs until both result lines exist";
   - the operator's console whitelist (the files are readable by the same uid);
   - T1-blindness of agents.
7. **Small run counts.** T1 has 6 runs and T5 has 2, so CIs over runs are coarse. The leave-one-run-out check is conservative but coarse, and T5's arm M rests on 2 runs of one lab.
8. **Study-wise error** across the primary and S families is up to the sum of their levels (unchanged, disclosed).
9. **Speed** from the look is claims-grade only if the node was otherwise idle (A5.5.1.7). Otherwise it is labelled "not claims-grade".
10. **The C9 time box favours NO-GO by delay.** This is mitigated only by the disclosure and immediate-commit rule (A5.1.1).

---

## A5.14 Points for the PI (decide at adoption; defaults are written into the text)
1. **GO-branch continuation (A5.1.3).**
   - A5 lets T1 proceed after A6 whatever A6 shows.
   - An A6 fail voids v1's FDR-based T1 claims.
   - D28 states neither explicitly. Both are needed to close options. Confirm or reject.
2. **Time boxes:**
   - C9: 48 h after 3a71356, otherwise NO-GO;
   - GO: A6 tagged within 7 days, and completed or stopped within 14 days, of the GO record.

   Confirm or change these before adoption. They cannot change afterwards.
3. **v2_sep title if A6 passes v2_sep before T1** (A5.9.3). D28 fixes "secondary; not gate-validated". A5 keeps that title unless v2_sep has passed a null-spectrum test on the full list before the look, because otherwise the title would be false. Confirm.
4. **Retries re-run both datasets** (A5.5.6). This costs a T5 re-run, but avoids a partially read look.
5. **Sealed-derived intermediates are deleted, not kept** (A5.5.3). No audit re-analysis is possible. That matches "no re-analysis" and avoids storing T1 PSM tables.
6. **P1 adds a leave-one-run-out backstop for v1** (A5.9.2), which is stricter than §1 as written. Confirm, or drop it to keep §1 exactly.
7. **Exclusive node** (A5.5.1.7). The look waits until no other agent of this project is active.

---

## A5.15 Revision log

**Revision 0:** first draft.

**Revision 1** answers `reviews/amendment_A5_review.md` (2 blocking, 10 major, 12 minor). Every finding is accepted.

| Finding | Resolution |
|---|---|
| B1 incomplete look can "succeed" | Completeness rule, normative and checked twice (driver exit + claim script); a not-evaluable comparator counts as not beaten (A5.5.4, A5.8, A5.11.2.4, A5.11.2.8, A5.11.4) |
| B2 interrupted look undefined; plaintext left | Unterminated look = failed and consumed; `--recover`; SIGTERM writes a failure line; one root on shared storage; time limit ≥ 2× the dev dry run; one-look check over the ledger's git history (A5.5.1.8, A5.5.3, A5.5.6, A5.11.1) |
| M1 claim-level arm-C void ill-defined; no arm C at protein level | Void defined by (a)–(f) (A5.7.1); arm C "not computed at protein level" in §6.3 note, L-E, L-M, A5.5.4; A5.7.1(1) restricted to arm P at protein level |
| M2 P1 drops §1's own-q half | P1 = (i) matched ∧ (ii) own q over supported comparators; §1 note; LOO added as a stricter backstop (PI point 6) |
| M3 environment not pinned; deny-list breaks driver | Allow-list environment built by `sealed_run.py`, asserted after re-exec (exempting `VULCAN_T1_STAGED`); environment manifests (Sage, Comet, TRFP, mokapot env, staged pred env) hashed; staged env rebuilt; `scripts/e9/`, `scripts/baselines/` frozen (A5.11.2.2, A5.11.5, A5.11.6) |
| M4 regressions from A4.8 mechanics | A4.8.2.1 root outside repo, A4.8.2.3 leak test, A4.8.1.4 console whitelist adopted; dev `--look` dry run; exclusive node; out dirs must not pre-exist (A5.5.1.7, A5.5.3, A5.11.8) |
| M5 T5 design not executable/ambiguous | Parametrised `bench/t5_armM.py` shared by dev and sealed; v1 via `ours.py Scorer`; merged strata decide, raw reported; r peptide-level; convention (a); dev decoy seed 1001; finality clause; separate fresh caches and bit-identity check (A5.6, A5.7.5, A5.11.3, A5.11.8) |
| M6 GO branch infeasible/unbound | A5.10 binds at adoption; A6 code outside frozen trees, or the freeze moves after A6 code; GO time box (PI point); P(v1 minimum) defined on the A5.10.3 floor; C9 record bound to the pilot file and recomputed by `sealed_run.py` (A5.1, A5.12) |
| M7 S6 weakened; loopholes | Cluster definition, cap + capped/uncapped IUT, and label-free stopping rule added to the floor; one open amendment per scorer, every draft reported; T1-blind definition covers commit messages and the status summary; T1-blind agents at pre-look commits; commit hygiene (A5.10.3–A5.10.4) |
| M8 D18 screened arm C missing | Computed in the look, reported only (A5.3 D18 entry, A5.11.2.5, R4) |
| M9 missed clauses | Added: A2 P3, A2.11.6, A3.2.5, A3.0 Q5, A3.3.2/A3.6 step 6, A3.1 moot/still-required list, §1, §6.8, §6.3 required arms, A3.3.8 |
| M10 lint only on claims.md; list incomplete | `tests/report/test_report_labels.py` over claims.md, REPORT.md and release/, with the extended banned list and the placement checks (A5.11.7) |
| m1 F3 vs A5.10.6 | F3 reworded: no outcome can upgrade; a void blocks any later upgrade |
| m2 "strongest comparator" undefined | The A1 comparator with the lowest P1(i) point ratio (A5.9.4) |
| m3 L-R in the GO-pass state | Met-precondition text defined (A1.8 entry, A5.9.3) |
| m4 L-M wording inaccurate | "own q is used only to order its list, not as an error estimate" |
| m5 L-G "Passed" list | Data named per gate; G2-original vacuity noted; dev arm M added (A5.4, L-G) |
| m6 retry contradictions | New milestone tag on retry; mutable config fields listed; `prereg-v0.4.4.x` after a re-freeze; rule for irreproducible failures (A5.5.6, A5.11.5) |
| m7 `data/sealed_derived/` not ignored; indexes huge | Lock directory dropped: sealed-derived intermediates deleted, look root outside the repo; indexes may be kept (A5.5.3; PI point 5) |
| m8 derived files described inconsistently | One rule: delete (A5.5.3, A5.11.2.6, PI point 5) |
| m9 driver should assert both tags | Added (A5.11.2.1) |
| m10 D28 cost figure | A4.7.4's 0.85–2.0 TB stated to supersede D28.1's 0.75–1.9 TB (A5.0.4) |
| m11 A1.8 exception scope | Stated as intended: T2–T4 stay blocked for v1's re-ranking pipeline (A5.1.5) |
| m12 C9 time box gameable | Disclosure of who controls the pilot commit, and an immediate-commit rule (A5.1.1, A5.13.10) |
