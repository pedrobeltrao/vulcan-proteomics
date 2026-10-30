# T1 sequencing relative to A4: independent statistical opinion

- **Reviewer:** independent statistician (agent), the same reviewer as `reviews/amendment_A4_review.md`. I had no role in A1–A4, the T1 driver or E7/E7b.
- **Date:** 2026-10-28.
- **Question (verbatim, from the project lead):** "D12 blocks the T1 look until null gates pass. Scorers are frozen, so T1 numbers don't depend on gate outcomes; only claim labels do. Is it statistically sound to run the T1 look now (after the A3 freeze), with all FDR claims labelled 'not gate-validated, A4 pending', and then A4 later? What risks or leakage paths (T1 informing A4 design, etc.) exist and what safeguards are required?"
- **Facts checked:**
  - `results/SEALED_LEDGER.md` has no entries.
  - `prereg-v0.4.3` (the A3 freeze) is **not** tagged yet; tags exist up to v0.4.2.
  - `bench/sealed_run.py` allows several sealed runs per milestone (its docstring, rule 2). Nothing technically enforces a single T1 look.
  - `bench/t1/a3.py` computes family-S inputs for the **full** v2_sep list (`vulcan_v2_sep`), plus per-stratum entrapment verdicts and IDs. I found no output for the **stratum-1-only** list's matched-FDP ratios against each comparator, which is the list family S uses if only tier 1 passes (A3.3.5). I found this by grep and have not verified it exhaustively.

## Verdict: **sound only with safeguards, and not sound "now".**
A T1 look before DN3 step 5 is acceptable once A4 is revised, adopted and tagged, and DN3's label-free walk (steps 2–3) is committed with the v1 minimum met. It is not acceptable before then.

**The premise is mostly right.** The scorers, comparators, estimators, seeds and the driver are frozen (A3.1.3, A3.1.6, A3.1.7). So the T1 *numbers* are a fixed function of frozen code and sealed data, and the gate outcome changes no number.

**But the premise is not complete in either direction.**

### 1. The gate outcome selects what is computed or reported on T1
- **Family-S list.** v2_sep's tier selects which list family S runs on: stratum-1 only, or the full list (A3.2.6, A3.3.5, A3.3.7).
  - If T1 is looked at before the tier is known, both branches must be computed in the single look.
  - Otherwise a second, post-hoc analysis of sealed outputs is needed, after the tier and after the analysts have seen T1.
- **Already pre-specified, provided nothing is recomputed:**
  - the production default (A3.2.6), which needs tier 2 and no T1 "violated";
  - the claim table and wording (A3.3.7).

### 2. T1 can inform A4, and that is the main risk
- **A4 design.** A4 is still a draft, and my review asks for material changes:
  - the organism tiers;
  - the keyword handling;
  - the cluster cap;
  - the stopping rule.

  T1 contains no null labels. It does contain arm-C entrapment, which is *Arabidopsis* natural sequences against decoys, and per-stratum verdicts. That is direct information on whether the scorers prefer natural sequences, which is exactly what DN3's tier II (land plants) and the V-rows test. Anyone who sees T1 could, even unconsciously, choose:
  - the organisms;
  - the cap;
  - the n_min values;
  - whether to keep tier II.

  These choices shift the probability that DN3 passes.
- **A4 execution.** These steps involve judgement:
  - the DN3 walk itself (DN2 needed 21 dated deviations);
  - step-4b engineering fixes;
  - the ce_s 24-h rule;
  - "investigate" branches such as a host-threshold mismatch;
  - V3 "not evaluable (implementation)".

  A T1-exposed approver can steer any of them.
- **Outcome-dependent persistence.**
  - A disappointing T1 creates an incentive to deprioritise or abandon A4. That would leave "A4 pending" numbers permanently.
  - A good T1 creates pressure to "fix" a failing A4.
  - My A4 review finds a substantial chance that DN3 ends **not evaluable** on label-free counts (B2). The *next* null-set amendment would then be designed after T1 had been seen. That is the strongest version of T1 informing gate design.
- **Spent sealed set.**
  - The rationale of A3.3.1 ("no T1 look if v1 fails or is not evaluable") includes not spending the one-shot T1 on a scorer without validated FDR.
  - If T1 is looked at now and A4 then fails, T1 is spent. A future validated scorer (v3) cannot use it as an untouched confirmatory set.
  - This is a cost, not a bias, but it must be accepted knowingly.
- **Published-in-effect numbers.** If A4 later fails, T1 numbers will already exist. The matched-FDP ID claims use arm P, whose shuffled entrapment does not detect a "natural-sequence" leak (arm C does). So they are FDR-dependent claims in D12's sense, not label-independent facts. They must stay labelled for good, not just until A4 reports.

**Conclusion.** A T1 look does not change any error rate of the A4 gate tests, *provided* nothing in A4's design or execution can depend on T1. The safeguards below make that checkable. Before A4 is tagged, they cannot hold.

## REQUIRED safeguards (all checkable)
1. **Order.** No T1 look before all of the following are committed:
   1. `prereg-v0.4.3` (A3 freeze, A3.1 preconditions);
   2. A4 revision ≥ 1 with its review, adopted, implemented and tagged as `prereg-v0.4.4`;
   3. DN3 steps 2–3 (membership hashes, label-free counts, swap test) with the A4.3.8 v1 minimum **met**.

   If the v1 minimum is not met, T1 stays blocked until a further null-set amendment is tagged and its walk reaches the same point.
   - **Check:** the `sealed_run.py` ledger line for T1 is dated after those commits, and the driver refuses unless they exist. Add this assertion to the driver.
2. **One look, every branch computed.**
   - Exactly one `sealed_run.py` T1 invocation, enforced in code: refuse if the ledger already has a T1 line.
   - That invocation must output everything any gate outcome could need:
     - the v1 primary;
     - v2_sep family-S inputs for **both** the full list and the stratum-1-only list;
     - all comparators, LOO backstops and entrapment verdicts.
   - After the look, the analysis is a hashed, committed, deterministic script that maps (A4 verdicts, T1 outputs) to the claim table. No re-search and no new estimator.
   - **Check:** driver and script hashes in `config/a3.yaml`/`a4.yaml` before the look; a dry run on D0/D1 shows both lists present.
3. **A4 is blind to T1.** These people must not have seen any T1 output:
   - the DN3 step-4b and step-5 agents;
   - every reviewer of an A4 deviation or engineering fix;
   - the author and reviewer of any later null-set amendment.

   Any A4 deviation approved by a T1-exposed person needs a co-sign from a T1-blind reviewer. Each A4 artefact (NOTES, deviations, adoption notes) states its author's T1 exposure.
   - **Check:** a T1 access list (who opened `results/<T1>/`), and exposure statements in every A4 commit after the look.
4. **No T1-contingent change to A4.**
   - After the look, A4 may change only through engineering fixes that touch no statistic, pool, threshold or selection rule (A3.6 rule), reviewed T1-blind.
   - A4 is run to completion whatever T1 shows: steps 4–5 run whenever the v1 minimum holds. Stopping or abandoning A4 for any reason other than its own label-free rules is a recorded protocol violation, and it is reported with the T1 numbers.
5. **Labelling and reporting regardless of outcome.**
   - Before A4 step 5, every T1 number carries: "not gate-validated; A4 pending; v1 FAILED tail-G2/G1/G3 under A1 (D12); A2 not evaluable (D24)".
   - After A4, the A3.3.7 claim table applies unchanged, with "A4/DN3" in place of "A2/DN2".
   - If v1 fails or is not evaluable on DN3, the T1 numbers are still reported in full, under "FDR not gate-validated (v1 failed / not evaluable under A1 and A4)". No FDR, ID-at-FDR or superiority claim is made for either scorer, and v2_sep is not production.
   - T1 numbers are never withheld because of the A4 outcome.
6. **Production default and wording only after A4 step 5** (A3.2.6 unchanged). Until then, nothing is called "default" or "validated".

## Recommended (optional)
- **Escrow instead of reading.** Run T1 when it is convenient, but write its outputs to a sealed location (`chmod 000` with a hash committed, as for `data/sealed/`), and open them only after A4 step 5 is committed. This removes every T1→A4 path at no statistical cost. It is the cleanest option if the goal of running T1 now is scheduling rather than reading.
- **No external release before A4 step 5.** No preprint, talk or slide shows T1 numbers before A4's step 5, even with labels.
- **Pre-commit a deadline and resources for A4** (and for one further attempt at α/4 if DN3 is not evaluable), so that persistence cannot depend on T1.

## Text that must change (by dated amendment, before the look)
- **D12.** Replace "FDR claims and the T1 look stay blocked until gates pass" with: "FDR claims stay blocked until v1 passes the null gates. The T1 look may occur only under the safeguards of `reviews/a4_t1_sequencing_review.md` (A4 tagged; DN3 steps 2–3 committed with the v1 minimum met; one look computing every branch; A4 agents blind to T1; reporting regardless of outcome)." State plainly that this relaxes the timing of the look, not any claim condition.
- **A3.3.1 and the first row of A3.3.7.** Replace "If v1 fails or is not evaluable, there is no T1 look under A3" with:
  - "The T1 look follows DESIGN §6.0d (A4.8). If v1 fails or is not evaluable on DN3, the T1 numbers are reported as 'FDR not gate-validated' and no v1 or v2_sep claim is made."
  - Add the stratum-1-list outputs to A3.3.2 and A3.3.8 as part of the single look.
- **A4.8.** Make the six required safeguards normative A4 text at the tag. Do not include them "by reference to revision 1". Add to A4.4:
  - an adoption-note line stating the adopter's T1 exposure (none, at adoption);
  - the T1-blindness requirement for the 4b and step-5 agents and for deviation reviewers.

  Add to A4.9: "A4 is reported whatever T1 showed; T1 is reported whatever A4 shows."
