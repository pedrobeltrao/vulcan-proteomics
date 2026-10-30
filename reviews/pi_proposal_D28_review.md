# Independent statistical review of the PI proposal D28 (do not run A4/DN3; T1 without null-spectrum gates)

- **Reviewer:** independent statistician (agent), fresh context. I had no role in A1–A4, their reviews, the T1-sequencing opinion, DN2/DNIT, E7/E7b, the T1 driver or any gate run.
- **Date:** 2026-10-28.
- **Proposal reviewed (D28, not adopted):**
  1. Do not run A4/DN3, citing P(v1 evaluable) 0.08–0.43, P(clean v1 passes) ≈ 0.06–0.4, a cost of 3–4 days and 0.75–1.9 TB, the consumption of the remaining pool by a not-evaluable outcome, and the lack of a human holder for the off-cluster escrow key.
  2. Amend D12 so that the single T1 look goes ahead after `prereg-v0.4.3`, with A3.3's v2_sep arm and no escrow. FDR-based claims would carry the label "entrapment-validated (arms P/C/M per §6.3); null-spectrum gates not passed: A1 fail on contaminated nulls, A2 not evaluable, A4 not run". The proposal also bans "validated FDR control" wording, reports everything in full, and requires any future null-set amendment to disclose that a T1-unblinded PI designed it.
  3. v2_sep stays "candidate".
- **Read (committed files only):** DESIGN.md §1 and §6 (§6.0, §6.0b, §6.0c). Because §6.0c and §6.0b adopt those texts verbatim, I also read the normative A3 text (`reviews/amendment_A3_draft.md`: A3.0, A3.2.5–A3.3.8, A3.5–A3.6) and A2 (`reviews/amendment_A2_draft.md`: P4, A2.4, A2.9.5–A2.11). Also read: RESEARCH_TREE D9–D27; `reviews/amendment_A4_draft.md` (rev 3); `reviews/amendment_A4_review.md`; `reviews/amendment_A4_review2.md`; `reviews/a4_t1_sequencing_review.md`; `reviews/checks/a4_r3_walk.out`; `results/2026-10-27_T1driver_dryrun/NOTES.md`; `results/2026-10-28_E7b_fair_speed/NOTES.md`; PROMPT.md.
- **Facts checked:**
  - `results/SEALED_LEDGER.md` has no entries.
  - Tags exist up to `prereg-v0.4.2`; `prereg-v0.4.3` does not exist yet.
  - `bench/sealed_run.py` contains no one-look rule (grep).
- **Not read:** anything sealed, escrowed or per-label from any null set, and any T1 output (none exists).

---

## Verdict: **adopt with changes.** Do not adopt as written.

**Part 1 (do not run A4) is defensible.** The evidence is label-free and predates any listing. A test whose chance of being evaluable is 0.08–0.43 in the reference cell, and whose failure to evaluate burns most of the scarce foreign-organism pool, is a poor use of 3–4 days. Two corrections are needed:
- **The escrow-key argument is not a reason to skip A4.** It only argues against running T1 *in parallel* with A4. "A4 first, then T1" needs no escrow at all. The decisive reasons are low evaluability, cost and pool consumption, and D28 should say so.
- **The cited figures must include the sensitivity grid**, not only the reference cell (see M5).

**Part 2 is a real relaxation of the pre-registered protocol, but D28 presents it as a D12 timing amendment.** It breaks at least eight other normative clauses (section B). Its proposed label is itself an overclaim: "entrapment-validated" and "arm M" are wrong, and "contaminated nulls" uses the post-hoc explanation that D12 forbids. D28 can be adopted honestly only if it:
- is written as a dated, tagged amendment that says plainly it relaxes the conditions for FDR claims, not just the timing of the look;
- closes every option it opens (v2_sep family S, the production default, the conditions that void a claim);
- pre-commits the statistics of any future null-set test now, while the PI is still T1-blind.

**Before T1 is spent, a cheaper and probably far better-powered null-spectrum test exists:** in-proteome "removed-peptide" human nulls (section E). It avoids the foreign-organism bottleneck entirely. I recommend making it a bounded, label-free go/no-go step before the look (change C9). This is not a precondition for the honesty of D28.

---

## A. Is the proposal sound, honest and non-outcome-driven?

### What is in its favour
- **No T1 information exists.** The ledger is empty and no `prereg-v0.4.3` look has happened.
- **No DN3 information exists.** No DN3 candidate has been listed.
- **DN2 has no per-label quantity.** Every number D28 cites comes from label-free walk simulations (`a4_r3_walk`) committed before this proposal.
- So D28 meets the letter of §6's rule ("a dated note … written before the affected results are seen") for T1.

### Where it is outcome-driven in the wider sense, which must be disclosed
- **The PI's information set includes a completed gate test that v1 failed** (A1: tail-G2, G1 and G3 FAIL, D12). D28 removes the requirement that this failure blocked.
- "The only completed test failed, the next one is probably not evaluable, so we stop requiring the gate" is the pattern pre-registration exists to prevent. It does not matter that each step has a label-free justification.
- This is acceptable only if the claim language really does not upgrade, that is, no FDR-*control* claim is made for v1 or v2_sep. And the D28 note must state the PI's information set: "D28 is decided knowing that v1 failed its only completed null-spectrum test."
- **What can falsify the honesty of D28 later** is how it is read: whether labelled claims are presented in REPORT.md as the headline. Section C fixes this mechanically.

### The escrow argument is misattributed
- Escrow (A4.8, D26.3) was needed only because T1 was to be computed *before* DN3 step 5.
- Running A4 and then T1 needs no key and no escrow. It costs the same 3–4 days, plus the delay of T1.
- So the key problem argues "if A4 runs, T1 waits". It does not argue "do not run A4".
- **Required:** D28 drops the key as a reason for item 1, or demotes it to "it also removes the option to overlap T1 with A4".

### The prospective numbers are honest but must be quoted completely
- The reference cell, which is anchored on DN2's own listing counts, gives P(v1 minimum) = 0.43 at e = 0.15 and 0.08 at e = 0.35.
- The sensitivity grid runs from 0.01 (6 files per project) to 0.98 (24 files per project).
- **The number of usable RAW files per project is unknown before listing.** The reference cell is the best estimate, but D28 should quote the grid too. Otherwise a reader can fairly suspect the pessimistic end was chosen.
- "P(clean v1 passes) ≈ 0.06–0.4" is **conditional on v1 being clean**. The relevant quantity for what is lost is the probability of *detecting* a leak at the margin. That is about P(evaluable) × P(fail | evaluable, margin leak) ≈ 0.08–0.43 × ≈ 0.9 ≈ **0.07–0.4**. Quote both figures.
- A4.7.4 gives a volume of 0.85–2.0 TB, while A4.11.6 gives 0.75–1.9 TB. Cite one source.

---

## B. What running T1 without gates breaks in §6, beyond D12
Each clause below is normative, and each needs explicit replacement text in the amendment. A D-note saying "D12 amended" is not enough.

| # | Clause | What it says | Effect of D28 |
|---|---|---|---|
| B1 | §6.4, "Gates for learned or predicted components" | all must pass "before any FDR claim" | D28 makes FDR-based claims (labelled) without G1–G3. This is a relaxation of the claim condition, not the timing |
| B2 | A1.2 (§6.0 item 2) | fail or not evaluable "blocks FDR claims for that scorer" | same as B1 |
| B3 | A1.4 | G1, G3 and G4 required "before any FDR claim" | same as B1 |
| B4 | **A1.8** | "Re-ranking of the top-K by a static scorer is allowed only if tail-G2, G1 and G3 pass" | The T1 driver searches with `report_psms 5` and rescores with the frozen static scorer. Wherever that rescoring can change the rank-1 candidate, every v1 number on T1 comes from a pipeline mode §6 does not permit. v2_sep is the same (top-5 truncation, which is the A3 V12 row). This is not just a labelling issue: the amendment must permit the mode explicitly and name it in the label |
| B5 | **A3.3.1 bullet 1; A3.3.7 row 1** | "If v1 fails or is not evaluable, there is no T1 look under A3 … Any such route needs a new amendment adopted before any T1 look, and that amendment must report v1's failure" | A3 itself requires a **new amendment before the look**, not a D-note. D28 must become that amendment |
| B6 | **A3.3.5** (family S) | the list is "the full v2_sep list if tier 2 passed, or the stratum-1 list if only tier 1 passed" | Without any tier, family S has no defined list. As written, D28 leaves open whether S1–S5 are tested and on which list. That is an analysis option open at look time (A3 Q6 forbids this) |
| B7 | A3.2.6 production rule | default is v2_sep only with tier 2; otherwise v1 iff v1 passes A2; otherwise "no OT-HCD scorer has gate-validated FDR claims" | D28 item 3 covers v2_sep only. It must also say what v1 is called (the de facto default since D7), and that T1 cannot change the default |
| B8 | A1.5 and §6.3 arm M | arm M required for predicted-feature scorers; on dev D7; **T5 at the milestone** | T1 does not compute arm M (T1 is not the arm-M design). D28's label "arms P/C/M" claims something T1 does not provide. Either T5 is opened in the same milestone with its own one-look rule, or the label says "arm M: dev D7 only (D23)" |
| B9 | A2 P4, A2.9.6; A4.6 (not adopted) | the attempt ledger and α spending | If a null-set test ever happens after T1, its level and statistics would be chosen by a T1-exposed designer. A4.6 is a draft, so nothing binds this today |
| B10 | D26.3 (escrow) | T1 only into escrow | Superseded. D28 must say so explicitly, as a reversal of a PI decision |
| B11 | §6.6 sealing; the sequencing opinion's safeguard 2 | one look | `sealed_run.py` allows several sealed runs per milestone (checked). The one-look rule exists only in the unadopted A4.8.2 |
| B12 | D14, D23, D25 caveats | arm-C 1.1–1.2% caveat at α = 1%; "name the merged block" (stratum ≥ 3: 1.18%/8.0%); "Sage chimera at narrow precursor tolerance" | Already binding wording. D28's label must not replace it, only add to it |

What is not broken:
- the §1 non-inferiority test;
- §6.3's verdict table;
- seeds;
- the estimators;
- the A3 freeze (A3.1). That freeze becomes simpler, because no A4 code has to land before it, and D27's merge can follow the freeze.

---

## C. Is the label sufficient? No. Required wording changes
1. **Delete "entrapment-validated".**
   - §6.3 has no such category. Its expected outcome for a calibrated method at T1 size is "inconclusive", reported as "consistent with control".
   - "Validated" contradicts D28's own ban on "validated" wording.
   - **Replace it with the literal §6.3 verdicts per arm, level and unit**, for example: "arm P: inconclusive, consistent with control; arm C: supported (peptidoform, α = 1%)".
2. **Arm M:** write "arm M: dev D7 only, no violation (D23); not run on sealed data", unless T5 is opened (B8).
3. **Remove "on contaminated nulls" as a qualifier of the A1 result.**
   - D12: the contamination explanation "was found after seeing the failure" and is never used.
   - A2.9.6 allows the phrase only inside a *pass* statement on fresh nulls.
   - **Required text:** "A1 (2026-10-27): tail-G2, G1, G3 FAILED; final for those nulls. A post-hoc contamination explanation exists, is unconfirmed, and is not used."
4. **Add the gates that passed, and say on what.** G0 strict; G4 pass (weak, τ = 0.298); G2-original and G5 on the A1 nulls (state the null set); §6.7 MC pass. The reader should see the whole gate table, not only the failures.
5. **Add the pipeline-mode statement (B4):** "top-K re-ranking by the static scorer is used although A1.8's precondition (tail-G2/G1/G3 pass) is not met."
6. **Add the decision provenance:** "A4 was drafted and reviewed twice, then not run by PI decision D28 (2026-10-28), made before any DN3 listing and before any T1 look, on label-free feasibility grounds: P(v1 evaluable) 0.08–0.43 in the reference cell (grid 0.01–0.98)." Do not write "A4 not run" bare, which reads as "pending".
7. **Separate the two kinds of FDR-dependent number.**
   - (a) **Matched-entrapment-FDP comparisons.** These do not trust either method's own q, but rest on the entrapment assumptions. Arm P's shuffled entrapment cannot see a natural-sequence preference; arm C can.
   - (b) **Own-q numbers.** These trust Vulcan's TDC, which is exactly what the gates were meant to protect.
   - Label (b) "own-q; estimator validity on this scorer not established by null-spectrum tests". Where (a) is quoted at matched arm-P FDP, quote arm C's mean FDP̂ and verdict for the same list next to it.
8. **Placement.**
   - The label goes in the REPORT.md abstract or headline sentence and in every table caption with FDR-filtered numbers. That includes protein groups and quant metrics computed on q-filtered lists.
   - It never goes only in a footnote or a limitations section.
   - The headline sentence must not start with a gain figure.
9. **Short form for prose:** "FDR-filtered at the method's own threshold and checked by entrapment (verdicts as stated); not validated by null-spectrum tests: the only completed test (A1) failed." This is the minimum, and every short form must contain the word "failed".

---

## D. What is lost relative to running A4
1. **Detection of a gross sequence-level leak, at PSM-tail level and independent of the entrapment assumptions.**
   - If DN3 were evaluable, tail-G2 at the A4 levels fails a margin leak about 0.86–0.94 of the time (`a4_r3_walk` tail_r3). Unconditionally, about 0.07–0.4.
   - The A1 margins already leave a +0.05 tail blind spot (A2.11.2: "tail FDR control rests on the entrapment arms"). So A4's incremental protection was always limited to large leaks.
   - On T1, arm C has list-level power against the same failure mode. That is the honest residual defence, and it is why change C7 and safeguard S4 matter.
2. **Any v2_sep null-spectrum evidence.**
   - v2_sep has never had a completed null test.
   - The V-rows (secondary hypotheses, A-syn, V12 top-5 truncation) were its only check on chimeric secondary PSMs other than per-stratum entrapment, and the stratum ≥ 3 entrapment is already marginal (D23).
   - v2_sep can therefore never be more than "candidate" under D28. Item 3 is correct and must be stated as permanent until a future gate pass.
3. **ce_s:** a scientific answer (does a sequence-capacity model leak on clean nulls?). It licenses nothing, so the loss is minor.
4. **T1 as an untouched confirmatory set** for any future validated scorer. A4.8.0 had already accepted this. It is unchanged, but D28 should restate it.
5. **A T1-blind design for any future null set.** This is the largest *procedural* loss, and it is mostly avoidable (S6, C9).
6. **What D28 gains:** 3–4 days; 0.85–2.0 TB; the remaining foreign-organism pool, kept unspent; and no dependence on a human key holder.

---

## E. A cheaper null-spectrum test with real power: in-proteome removed-peptide human nulls
Every attempt so far (A1, A2, A4, DNIT) has been bottlenecked by **foreign-organism supply**: about 82–135 eligible projects, few labs, and marker contamination. The test that would remove the bottleneck uses **human** spectra.

### Construction (a sketch, for a new amendment; it is not specified here)
1. **Source.** Seeded draw from the T4 frame (14,190 projects) of human OT-HCD DDA projects. Exclude:
   - the T4 selection, T1–T5, train and dev, E9/D8;
   - predictor-training accessions;
   - DN2/DNIT/A1.

   Take 2–3 runs per project from ≥ 30 distinct labs. The pool is thousands of projects, so a label-free stopping rule can simply keep drawing; a not-evaluable outcome burns nothing scarce. Volume is about 0.1 TB.
2. **Certification by a fixed score.** Certify spectra with hyperscore or Comet at q ≤ 0.1%, never with v1 or v2. Require a single precursor hypothesis, as A2.2(e), using the MS1 feature finder.
3. **Symmetric per-run removal (P1).** For each run, remove from **both** target and decoy sides the following, together with their paired decoys:
   - every peptide identified in that run or project at 1% by any fixed-score engine, including chimeric secondary IDs;
   - their I/L-collapsed Hamming ≤ 2 neighbours;
   - all modforms and nested or extended peptides that share a ≥ 4-ion b/y ladder (the D18 screen);
   - peptides sharing a 6-mer (A1.3).

   A certified spectrum is then a null competition by construction: its true peptide, and everything a fixed score could recognise in that run, is absent from T ∪ D.
4. **Run the production pipeline unchanged.** RT is the real human RT, with no quantile mapping. Compute tail-G2, G1 and G3, plus v2_sep's primary-unit rows, with A2/A4 statistics.

### Why it tests the right thing
The feared leak is that predicted features, trained on human data, give remaining *human* targets an advantage over decoys regardless of the spectrum. Here every candidate is a real human peptide that is absent from the spectrum, scored against human-trained predictors. That is closer to deployment than foreign spectra.

### Why A2.4's objection to human donors does not apply as stated, and what remains
- A2.4 rejected human donors because "keeping [the precursor] would leave the true peptide as a candidate", and because nested targets and correlated features favour targets. Removing the true peptide and its ladder-sharing relatives answers both of those.
- **What remains:**
  - unidentified co-isolated human peptides, and sub-threshold paralogs;
  - population shift towards good spectra, as in A2's certification.
- **These residuals favour targets for every scorer**, and a stronger scorer exploits real partial evidence more. So they push towards *fail* or *invalid pool*, not towards a false pass. The A2 validity control (hyperscore tail fraction ≈ ½) then turns them into "not evaluable".
- Their main risk is to evaluability, not validity. A new amendment must still argue and simulate this; I have not.

### Power and cost
- A rough estimate: about 5–10k certified single-precursor spectra per run, and a DN2-like tail rate of 0.09%. Then n_tail ≥ 200 needs about 25–45 runs, and G1 n_top and G3 are easily met at about 80 runs.
- Clusters are not binding with thousands of labs, so the A4 envelope rules work as intended.
- **These are guesses.** A **label-free pilot** on already-local dev human runs (D0, D8 stream) can measure certified yield, the tail proxy and the removal fraction in hours, with no per-label quantity. That pilot is the right go/no-go instrument.
- Expected total: design plus independent review about 1 day, implementation (masked replicate DBs already exist in the L1 index) under 1 day, and run hours.

### What does not work
- **Decoy-vs-decoy (two shuffles)** nulls have exact exchangeability but no power against a natural-sequence preference, since both sides are shuffled. That is G0's territory.
- **Simulated spectra from the same predictors** are circular.
- **A1-recipe offset nulls** were correctly rejected (A2.4).
- **FS nulls** are expected to be not evaluable.

### Status if adopted
- It is v1's **completed test 2**, at the A4.6 level α/2. It is v2_sep's test 1.
- It must report A1, A2 and A4-not-run.
- It is **not** a fallback that rescues a failure: a fail is final under P4.

---

## F. Required changes (the conditions for "adopt")

### Form
- **C1. Make D28 a dated amendment (§6.0d, "A5"; do not reuse the A4 name).**
  - It replaces, with quoted old and new text: D12; §6.4's gate sentence; A1.2; A1.4; A1.8; A2.9.6's T1 precondition; A3.3.1 bullet 1; A3.3.5's list rule; A3.3.7 row 1; A3.2.6's production rule; D26.3.
  - Its first paragraph says: "This amendment relaxes the conditions for FDR-based claims and permits the T1 look although v1 failed its only completed null-spectrum test. It is adopted knowing that failure, before any T1 or DN3 information exists."
  - It records A4 rev 3 as "drafted, reviewed twice, not adopted, not run".
- **C2. Have the amendment text independently reviewed T1-blind**, then tag it (`prereg-v0.4.4`) after `prereg-v0.4.3` and before any `sealed_run.py` T1 line. The driver asserts both tags.

### Close every option before the look
- **C3. Family S for v2_sep is fixed now.** My recommendation: run S1–S5 on the **full v2_sep list** as pre-specified, under the C-label, titled "secondary arm; not gate-validated". Report the stratum-1 list descriptively, never in its place.
  - The alternative is to run family S not at all, with a table titled "exploratory". Either choice is acceptable if made now. What is not acceptable is leaving it open.
- **C4. Void rule for v1, symmetric with A3.3.6.** Any "violated" verdict, in arm P or C, at 1% or 5%, in peptidoform or protein-group lists, voids all FDR-based v1 claims as final.
  - This is stricter, outcome-independent, and needed, because entrapment is now the only line of defence.
  - The A3.3.6 void rule for v2_sep stays.
- **C5. Production wording, fixed now:**
  - "v1: engine default, FDR not validated by null-spectrum tests";
  - "v2_sep: candidate";
  - **no T1 outcome changes either designation.** Only a future passed null-spectrum test can (per A3.2.6), and that would be disclosed as post-T1.
- **C6. Arm M:** either open T5 in the same milestone (one look, pre-specified, A1.5) or use the dev-only label (C.2). Decide now.

### Enforce one look and report everything
- **C7. One look in code.**
  - `sealed_run.py` refuses a second T1 ledger line, with a test.
  - The single run computes every branch: both v2_sep lists, every comparator, the LOO backstops, all verdicts, and matched arm-C figures if they can be added before the freeze without touching the v1 path.
  - A hashed, committed claim script maps outputs to the claim table and labels. It is reviewed before the look. No re-search and no new estimator afterwards.
- **C8. The report-regardless commitment is extended to the claim script's output verbatim**, including every void, every "inconclusive" and the comparator rows where Vulcan loses (e.g. the D17 0.990× against Sage+pred+mokapot, if it recurs).

### Pre-commit future null tests while the PI is still T1-blind
- **C9. Recommended: run a removed-peptide pilot before T1.**
  - Before the T1 look, commission a T1-blind agent to draft the removed-peptide null amendment (section E) and run its **label-free** pilot, with a time box of about 1 day.
  - **Pre-specified rule, fixed in D28:**
    - if the pilot's simulated P(v1 minimum) is ≥ 0.8, adopt, tag and run the test before T1 (no escrow needed);
    - otherwise proceed to T1 under C1–C8.
  - Either way the outcome is recorded.
  - This keeps the one cheap test that could still be T1-blind from being designed after T1.
- **S6 (required even if C9 is declined). Fix the rules for any later null-set test now, in D28:**
  - v1 = completed test 2 at α/2, with A4.6's geometric rule;
  - A4's margins, n_min floors, envelope and IUT rules are a **floor** that any later amendment may only tighten;
  - its author, reviewer, walk agent and step-4/5 agents must be **T1-blind agents**, which are always available in this project, unlike a human key holder;
  - a T1-exposed PI may adopt it only verbatim or with changes reviewed T1-blind;
  - a later pass relabels T1 numbers only with the text "null-spectrum gates passed after T1 was seen".

  "Disclosed as designed by a T1-unblinded PI" alone is too weak. Disclosure does not remove the degrees of freedom; pre-commitment does.

### Language and release
- **C10. Labels:** every change in section C, in particular dropping "entrapment-validated", the arm-M correction, removing the "contaminated" qualifier, the A1.8 statement, and quoting the D14/D23/D25 caveats.
- **C11. Quote the full evaluability evidence in D28:** the reference cell, the grid, and the detection probability in A. Drop or demote the escrow-key argument for item 1.

---

## G. Bottom line for the PI
- **Not running A4 is a reasonable, label-free call.**
- **Running T1 without gates is permissible only as an openly declared relaxation**, done by amendment, with closed options and a label that says "failed", not "validated".
- **As written, D28's label overclaims in three places** ("entrapment-validated", "arm M", "contaminated nulls"). It also leaves v2_sep's family S undefined, and it protects future null-set design only by disclosure.
- **Make the changes above, or at minimum C1–C8, C10, C11 and S6, and D28 is honest.**
- **Consider C9 seriously.** A human removed-peptide null is the one test whose supply does not run out. It is worth a one-day label-free pilot before T1 makes every later null-set design T1-exposed.
