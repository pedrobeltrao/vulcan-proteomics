# T1-blind review of the A5 claim script (A5.11.4.5)

- **Reviewer:** fresh agent, T1-blind. I did not write the claim script, its helpers, its tests or `config/a5.yaml`.
- **Date:** 2026-10-28.
- **Reviewed at:** HEAD 3b15bbb.
  - `bench/t1/claims_a5.py`, `bench/t1/wording_lint.py`, `bench/t1/completeness.py`.
  - `tests/bench/test_claims_a5.py`, `tests/bench/test_a5_config.py`, `tests/report/test_report_labels.py`.
  - `config/a5.yaml`, sections `gates`, `labels`, `a6`, `c9`, `claim_script_stricter_readings` and `comparator_nonreproducibility`. These are committed. The uncommitted working-tree diff touches only `look`, which is out of scope.
- **Normative basis:** `reviews/amendment_A5_draft.md` (A5.4, A5.5.4, A5.5.5, A5.7, A5.8, A5.9, A5.11.4, A5.11.7); A3.3.3–A3.3.6; DESIGN §6.0d/§6.0e; A6.2 "Status file"; RESEARCH_TREE D29–D33.
- **Supporting code read for key names only:** `bench/t1/evaluate.py` (`noninferiority`), `bench/t1/a3.py` (`loo_ratios`, `void_flags`), `bench/t1/a5.py` (output schema).

## Verdict: **adopt with changes**

The core decision logic is faithful to the text:
- the v1 void rules 1–3;
- the v2_sep A3.3.6/A5.7.2 voids, with T5 merged strata deciding and raw strata not deciding;
- the claim-level arm-C rule, which is strict `>` f, applies to the numerator only, and never applies to comparators;
- the family-S fixed sequence, an IUT over every A3.3.3 comparator by LB and LOO minimum;
- not evaluable = not beaten;
- the S3 stop on a protein-group "violated";
- the S5 skip;
- the P1 conjunction, and "no supported comparator";
- the A6 status read from the look base commit, with absent, unparseable or out-of-enum → "stopped: node time limit";
- the completeness refusal (exit 2, "look failed: incomplete");
- deterministic output.

The banned-phrase list matches A5.11.7. The label texts and gate tables are pinned verbatim to the draft.

There is one blocking defect: in the GO/A6-PASS state, v2_sep numbers get v1's pass wording. There are also five major items. All are small code or test changes, and they should land before `prereg-v0.4.4`.

## Tests run

- `pytest tests/bench/test_claims_a5.py tests/bench/test_a5_config.py tests/report/test_report_labels.py`: **103 passed, 1 skipped**. The skip is the REPORT.md/release parametrisation, because no REPORT.md or `release/` file exists.
- **Synthetic runs.** These used the test fixtures, wrote to a scratch dir only, and were all built on dev data:
  - the GO/A6-PASS state (temporary git repo);
  - a v2_sep stratum void;
  - a v1 claim-level void at peptidoform/1%;
  - no supported comparator.
- **End-to-end run on the committed D0 dry-run metrics** (`results/2026-10-28_A3_freeze/D0/with_a3_metrics.json`, dev, 2 runs), with a synthetic `a5` block and armM.
  - Result: exit 0 and lint clean.
  - P1 is met, with headline comparator sage+pred+mokapot 0.990 (LB 0.989, the D17 pattern). S1 and S2 are rejected; S3 is not rejected, so S4 and S5 are not tested.
  - The script works on real driver structure beyond D1.

## BLOCKING

**B1. GO/A6-PASS state: v2_sep tables carry v1's A6 "validated" L-Q and the A1.8 "precondition met" L-R.**
- **Where:** `bench/t1/claims_a5.py:216-231`, `:237-240`, `:916`, `:718-720`, `:776-781`.
- **Mechanism.**
  - `render_labels` makes **one** state-dependent `L-Q`. In the PASS state it becomes "…validated by A6 (test 2), after A1 FAILED", and `L-R_template` becomes None.
  - That single `lab['L-Q']` is then used in the captions of Table S, Table R3 and Table R3-strata.
  - `res['lr_v2'] = render_lr(lab, …)` therefore renders `L-R_pass` for v2_sep.
- **Verified on a synthetic PASS run:**
  - the Table S caption contains "validated by A6 (test 2)" and "precondition met at test 2";
  - the Table R3 caption contains "validated by A6";
  - `lr_v2` = "top-K re-ranking by the static scorer; A1.8's precondition met at test 2 (A6), before T1."
- **Why it is wrong.**
  - A6 gates v1 only (D32.1; A6.6.2: v2_sep "not run").
  - A5.9.3 changes only the v1 clauses; the v2_sep column reads "per v2_sep's own A6 status" or "unchanged".
  - The script's own recorded reading says the same (`STRICTER`, `:87`; `config/a5.yaml:256`): "v2_sep wording is always the no-A6 wording".
  - The output would therefore print a null-spectrum validation label on v2_sep numbers, which is exactly what F3/A5.0.2 forbid. The lint cannot catch it, because the pass phrase is whitelisted.
- **Fix.**
  - Render separate `L-Q_v1`/`L-Q_v2_sep` and `L-R_v1`/`L-R_v2_sep`. v2_sep always gets the no-A6 text.
  - In the PASS-state test (`tests/bench/test_claims_a5.py:393-396`), assert that the Table S, R3 and R3-strata captions and `lr_v2` do **not** contain "validated by A6" or "precondition met".

## MAJOR

**M1. The S1 "no supported comparator" reading is recorded as "stricter", but it is the permissive reading.**
- **Where:** `claims_a5.py:505-510`, `:83`; `config/a5.yaml:254`.
- **The text.** A3.3.5 gives an explicit rule only for S5 ("If none is supported, S5 is not tested"). The P1 "rests on (i)" rule (A5.9.2) is written for P1 only.
- **The code.** It lets S1 be rejected on its matched-FDP part alone, so S1–S4 stay claimable.
- **Why this needs a change.** A5.11.4.3 allows only the stricter reading where the text is open. The stricter reading is "S1's own-q part cannot be met, so S1 is not rejected and the sequence stops". The vacuous-IUT argument is defensible, but then it is a reading of the text and not a "stricter reading".
- **Fix (either):**
  - switch to the stricter reading; or
  - have the PI record it before the tag as a dated deviation explicitly labelled **non-stricter**, mirroring A5.9.2 P1(ii).
- Either way, `claims.md` must state it on the S1 row; it does today, via the notes.

**M2. A claim-level-void P1(i) ratio is quoted in the headline without the word VOID, and P1 is reported as "not met".**
- **Where:** `claims_a5.py:459-467`, `:581-586`, `:946-951`.
- **Verified on a synthetic run** (v1 arm-C LB lower CI 2% > f = 1%). The headline reads: "…P1 (v1 no loss) is not met; against … comet+mokapot, v1's peptidoform ratio at matched arm-P FDP̂ 1% is 1.2 (95% CI [1.1, 1.1] …)".
- **Why it is wrong.**
  - The quoted ratio is void under A5.7.1(d).
  - A5.7.3 requires voided numbers to be shown under a VOID title.
  - "not met" reads as a loss.
- **Fix.**
  - Give P1 a distinct outcome, e.g. "VOID (claim-level, A5.7.1: peptidoform/1%)".
  - Make the headline say VOID with the reason whenever the quoted ratio is void.
  - Update `test_claim_level_void_boundary` (`tests/bench/test_claims_a5.py:204`).

**M3. The stratum-1 list is put through the family-S decision arithmetic, and the outcomes are written into `claims.json`.**
- **Where:** `claims_a5.py:910-914`.
- **What happens.** `R1_descriptive_S_arithmetic` holds "rejected (claim holds)" / "not rejected" outcomes for the stratum-1 list. A D0-based synthetic run gives S1/S2 "rejected (claim holds)".
- **Why it is wrong.** A5.5.5: "It is never tested and never quoted as a claim". Test outcomes for that list in a committed output invite exactly that misuse.
- **Fix.** Drop the key. Table R1 already prints the descriptive ratios, CIs and LOO minima, which is all A5.5.5 asks for.

**M4. Voided numbers are not shown "in a table titled 'VOID under §6.0d A5.7 (reason)'" (A5.7.3).**
- **v2_sep.** When v2_sep is void, `family_s` returns S rows with `rows: []` (`:487-489`), so Table S shows only "VOID —" cells, with no ratios. The numbers appear only in Table R3, whose caption does not say VOID.
- **v1.** When v1 is void, the title is a level-3 heading (`:630-632`), followed by unrelated `### P1` / `### P2` headings. No table caption carries the VOID title.
- **Fix.**
  - Still compute and print the S rows when void, with outcome "VOID".
  - Prefix every affected table caption (P1(i), P1(ii), P2, P3, S, R3, R3-strata) with the VOID title when the corresponding void holds.

**M5. Tests do not meet A5.11.4.4.**
- **Gap.** The end-to-end test uses D1 only, with a synthetic armM (`tests/bench/test_claims_a5.py:495-505`). A5.11.4.4 requires "the D0/D1 dry-run outputs and the D7 arm-M output".
  - D0 (`results/2026-10-28_A3_freeze/D0/with_a3_metrics.json`) is committed and works; I ran it.
  - No committed `bench/t5_armM.py` dev output exists yet. `bench/t5_armM.py` and `bench/t1/a5.py` are still untracked.
- **Fix.**
  - Add D0.
  - Add the D7 arm-M end-to-end once the A5.6 dev-equivalence run is committed.
  - Add the B1 assertions.
  - Add a P1(ii)-fails-while-(i)-passes case.
  - Add a v2_sep T5 **overall** violated case; only merged is tested now.
  - Add a test that the code `STRICTER` list equals `config/a5.yaml claim_script_stricter_readings` (minus the date prefixes).

## MINOR

1. **Claim-level void when the arm-C LB CI is not finite.** It is treated as void (`claims_a5.py:442`). This is stricter, but it is not in `STRICTER` or the config list. Record it.
2. **Headline comparator ties and None ratios.** A None ratio counts as "lowest", and ties are broken by name (`:462-463`). Neither is recorded.
   - A None ratio means the comparator has 0 IDs, which is the opposite of a low ratio.
   - Prefer "lowest finite ratio; None never headlines", and record it.
3. **A6 FAILED is kept even when unverifiable** (`:171-182`; `config/a5.yaml:255`).
   - This departs from A6.2's explicit text ("absent, unverifiable … → 'stopped: node time limit'"). It is not a reading of open text.
   - It is stricter (it can only void), but it should be confirmed by the PI as a deviation.
   - Non-PASS/FAILED statuses are not verified at all. This affects only the C9-clause wording.
4. **The A6 time box is not checked.** `config/a5.yaml a6.time_box` (18:00, D31.3) is recorded but not checked: a PASS committed after 18:00 but before the look base commit would be accepted. Either enforce it for PASS (stricter), or state in `claims.md` that A6's own clock governs.
5. **Sealed-input handling.**
   - Sealed mode is detected only from `meta.role` (`:898`).
   - `run()` does not itself require ledger hashes for sealed inputs; only `main()` does (`:1060-1061`).
   - In `ledger_look`, a look without a matching line leaves `short=None`, so `rev-parse` crashes (`:1034-1038`).
   - Manifest entries are matched by basename only (`:871-873`).
6. **Captions.**
   - Table R1 prints own-q rows under an L-M-only caption, with no L-Q and no L-R (`:790-791`).
   - The Table S caption's L-E and L-M are for peptidoform 1% only, although S3 is protein (L-M protein variant; arm M "not applicable") and S4 is at 5% (`:718-720`).
   - v2_sep's L-R renders "—" for the spectra fraction, because `a5.py` gives v2_sep only the accepted-units fraction (`:916`; `bench/t1/a5.py:314`).
7. **Arm-M "consistent with control".** It uses arm M's own CI lower bound (`:283-284`). A5.9.1 defines it by the *paired* lower CI, and arm M has no paired control. This reading is unrecorded; record it or drop the suffix for arm M.
8. **Cosmetic label rendering.**
   - The v2_sep C9 clause renders as "A4 not run) ; C9 pilot: …" (stray space before ";", `:205-209`).
   - The A5.9.3 "(…)" in L-S pass is filled with invented text, "(A1, 2026-10-27: tail-G2, G1, G3)" (`:224`). Record it or keep the clause minimal.
9. **Lint scope.** It accepts "FDR-control claim" anywhere (`wording_lint.py:38`, `:73-76`); A5.11.7 limits that exception to the A5.9 labels, and quoted labels are already stripped. Remove the exception outside labels.
10. **Completeness rule gaps** (`completeness.py:10`, `:46-51`, `:75-79`).
    - Own-q ratios carry no `n_boot_nan`, and `evaluable()` skips the NaN check for them (`claims_a5.py:309-314`). The argument is that the 1e-9 floor in `evaluate.py:533` makes NaN impossible. It is sound, but A5.5.4 says "every bootstrap's n_boot_nan". Record it as a dated reading.
    - "All 6 runs in every method's run list" and "both T5 runs in every arm-M list" are checked only at top level. The LOO list length partly covers the first.
11. **Missing input for one v2_sep void.** "v2_sep FAILED a completed null-spectrum test before the look" (A5.7.2) has no input; it relies on D32 (v1-only A6). State this in `claims.md`. Likewise, the A3.1.8 record is checked by hash only; say "A3.1.8 (i)–(iii) hold (D23)" in the output.

## Checks that passed (no finding)

- **Rows.** Every A5.9.2 row is always printed: P1–P8, S1–S5, R1–R4, D and G. Incomplete inputs print only "look failed: incomplete" with the list of violations, as A5.5.4 requires.
- **v1 voids.**
  - Rule 1: arm P at peptidoform, sequence and protein; arm C at peptidoform and sequence; 1% and 5%; convention (a) only. 0.1% and convention (b) do not void (tested).
  - Rule 2: T5 arm M v1 peptidoform. The hyperscore reference never voids.
  - Rule 3: A6 FAILED.
  - Claim-level (a)–(f): the boundary at exactly f does not void (tested), and the comparator side is printed and never voids.
- **v2_sep voids.** Overall plus merged strata, arms P and C. T5 overall plus merged, with raw strata not deciding (tested). The union with `a3_secondary.void_flags` is recorded. The protein-group "violated" stop at S3 is recorded.
- **Family S.** Order, bars (≥ 0.97; > 1.00 strict), the "every comparator" IUT and LB + LOO for every S part. An undefined LOO minimum does not clear. Not evaluable (> 5% NaN of 2000, exactly 5% evaluable; `fdp_reached` not both True; CI not finite) = not beaten.
- **P1.** (i) LB and LOO ≥ 0.97 vs every A1 comparator. (ii) LB only, against "supported" (deciding convention) A1 comparators. "no supported comparator" is printed, and P1 then rests on (i).
- **A5.9.3 mapping.** NO-GO, GO+PASS, GO+FAILED and GO+other. PASS needs a verifiable gate output and a step-5 commit that is an ancestor of the base commit. A file that is uncommitted at the base commit maps to "stopped: node time limit" (tested).
- **Labels.** They are verbatim from `config/a5.yaml`, which is pinned to the draft by `test_a5_config.py`. Every L-S variant contains "FAILED". The L-G clause is filled from the C9/A6 state. The D33 non-reproducibility values are reported only.
- **Lint.** The banned list equals A5.11.7 (entrapment-validated, entrapment-verified, empirically verified, validated FDR, FDR control, controls FDR, controlled at, FDR-controlled except "not FDR-controlled", gate-passed, calibrated in an FDR/scorer context, contaminated nulls). "validated" is allowed only negated or in the A5.9.3 pass phrase. Quoted-label exemption applies only to exact rendered labels. The headline must start with L-S, and so cannot start with a number.
- **Determinism.** Sorted keys, no timestamps and no randomness; byte-identical across repeated runs (tested).

## T1 exposure

**None. T5 exposure: none. A6/C9 output exposure: none.**
- I did not read `data/sealed/`, any T1 or T5 output, or `results/SEALED_LEDGER.md`; I did not open the ledger at all.
- Data I read or ran: the committed dev dry-run metrics (D0 and D1, `results/2026-10-28_A3_freeze/`) and synthetic fixtures, all in a scratch directory.
- I modified no code.

---

## Addendum: re-check of 27fee10 (2026-10-28, 13:10 CET)

- **Same reviewer, same rules:** T1-blind, no code edits, and dev or synthetic inputs only, in a scratch directory.
- **What I re-checked:**
  - Commit 27fee10, as of HEAD b3b7e11. The only later change in scope is `config/a5.yaml look.root`, which is outside the sections reviewed here.
  - The three test files: **112 passed, 2 skipped**.
    - The REPORT.md/release parametrisation is skipped as before.
    - `test_end_to_end_d7_armm` is skipped because no committed dev `armM.json` from `bench/t5_armM.py` exists yet.
  - I re-ran my synthetic scenarios against the fixed script.

### Status of the original findings

- **B1: resolved.** `render_labels` now returns separate `L-Q_v1`/`L-Q_v2_sep` and v1/v2_sep L-R templates (`claims_a5.py:250-268`).
  - Synthetic GO/A6-PASS run: no v2_sep caption contains "validated by A6" or "precondition met". These are Tables S, R1, R2*, R3, R3-strata and the v2_sep gate table.
  - `lr_v2` keeps the "is used although … not met" text. `lr_v1` gets the pass text. The lint passes.
  - The PASS-state test now asserts all of this (`tests/bench/test_claims_a5.py:417-425`).
- **M1: resolved.** S1 with no supported own-q comparator is now "not rejected" and the sequence stops; verified synthetically: S2–S5 are "not tested (S1 not rejected)". The S5 "not tested" rule is unchanged, and the stricter reading is recorded in the code and the config.
- **M2: resolved.** A claim-level void at peptidoform/1% gives P1 the outcome "VOID (claim-level, A5.7.1: peptidoform/1%)".
  - The headline reads "P1 … is VOID under §6.0d A5.7 (claim-level …: arm-C mean LB lower CI 2% > f = 1%)".
  - The quoted ratio is marked "(VOID under §6.0d A5.7; printed in full, not a claim)".
- **M3: resolved.** `R1_descriptive_S_arithmetic` has been removed, and no stratum-1 decision outcome appears in `claims.json` (tested). Table R1 is descriptive only.
- **M4: resolved.**
  - A void v2_sep still computes every S row and prints it in full, with outcome VOID (synthetic: 10/6/6/6/4 rows).
  - The Table S, R3 and R3-strata captions carry "VOID under §6.0d A5.7 (reason)".
  - A v1 scorer-level void titles the P1(i), P1(ii), P2 and P3 captions. A v1 claim-level void titles P1(i) and P2, and the affected rows are marked VOID.
- **M5: resolved in code; one condition remains (R1 below).** New tests cover:
  - D0 end to end;
  - P1(ii) failing while (i) passes;
  - a v2_sep T5 overall violation;
  - STRICTER = config;
  - no stratum-1 arithmetic;
  - a None ratio never headlining;
  - sealed inputs requiring ledger hashes;
  - an A6 PASS after the deadline.
- **Minors 1–11: resolved or recorded.**
  - Recorded as dated readings:
    - claim-level void when the CI is not finite;
    - the headline uses the lowest finite ratio;
    - an unverifiable A6 FAILED is kept, explicitly flagged for PI confirmation as a deviation from A6.2;
    - own-q `n_boot_nan`;
    - the "(…)" fill;
    - arm M without "consistent with control";
    - "v2_sep FAILED a completed test" cannot fire, which `claims.md` states.
  - Fixed in code:
    - an A6 PASS must be dated before the A6 step-5 deadline (17:55, `config a6.pass_deadline`);
    - `run()` itself refuses sealed inputs without ledger hashes;
    - `ledger_look` checks for a missing base line;
    - the lint no longer accepts "FDR-control claim" outside quoted labels (now tested as banned);
    - caption labels are completed: R1 has L-Q/L-R, and v2_sep's L-R says "not defined for this pipeline" instead of "—".

### Remaining

- **R1 (MAJOR, a pre-tag condition, not a code defect): the D7 arm-M end-to-end test (A5.11.4.4) has not been run.** `test_end_to_end_d7_armm` skips while no committed dev `armM.json` exists, so a green suite at the tag would not show that this test ran.
  - **Condition:** before `prereg-v0.4.4`, commit the A5.6 dev-equivalence `armM.json`.
  - **Check:** `test_end_to_end_d7_armm` must then run and pass rather than skip.
  - **Suggestion:** turn the skip into a failure once the freeze is reached.
- **R2 (MINOR): the A6 PASS deadline uses the step-5 commit's committer date (`git log %cI`), which the committer can set.**
  - It only ever downgrades PASS, so it cannot loosen anything. But it is not tamper-proof evidence of timing.
  - **Suggestion:** have the A6 runner also record its step-5 time in its ledger or status file, and say so in the A6 NOTES.
- **R3 (MINOR, for the PI): confirm the recorded deviation that keeps an unverifiable A6 FAILED** (A6.2 literally maps "unverifiable" to "stopped: node time limit"). The script flags it for confirmation; the PI decision should be logged before the tag.

### Final verdict: **adopt**

Adoption is conditional on R1: the D7 arm-M end-to-end test must run and pass before `prereg-v0.4.4`. No BLOCKING finding remains.

**T1 exposure for this re-check: none. T5: none. A6/C9 outputs: none.** I read no sealed data and no ledger, and edited no code.
