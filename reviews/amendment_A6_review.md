# Review of Amendment A6 draft (revision 0): HN1 human removed-peptide null gate

- **Reviewer:** fresh T1-blind agent (statistical reviewer). I wrote no amendment and ran no gate step.
- **T1 exposure:** none. T5 exposure: none. I read no file under `data/sealed/`, no `*T1_sealed*`/`*T5_sealed*` path, and nothing from `results/SEALED_LEDGER.md`. I computed no per-label quantity.
- **Date:** 2026-10-28, about 12:15–12:40 CET. Time-boxed to about 20 minutes, so I report only blocking and major items.
- **Read:**
  - `reviews/amendment_A6_draft.md` (all);
  - `reviews/amendment_A5_draft.md` (A5.0–A5.2, A5.10);
  - `reviews/amendment_A3_draft.md` (A3.2.2–A3.2.6);
  - `reviews/amendment_A4_draft.md` (A4.5.2, A4.5.5);
  - `results/2026-10-28_C9_pilot/NOTES.md` (section 5, the cost section, limitations);
  - `reviews/checks/a6_power.py` (all), `a6_power.out`, `a6_cert_yield.json` (summary keys);
  - `scripts/a3/gates_a3.py` (V6 path, around l.137–231);
  - `scripts/c9/c9_power.py tail_iut` (contamination model);
  - function lists of `vulcan/gates/a2_stats.py` and `scripts/c9/pilot_yield.py`;
  - RESEARCH_TREE D29–D31; E9 NOTES (throughput line).
- **Not read:** A2 in full, A4.3.x, A4.6 text (I relied on A5.10.1–A5.10.2, which adopt A4.6), `a6_cert_yield.py`, `pilot_yield.py` bodies.

**Verdict: ADOPT WITH CHANGES.** B1–B3 must be fixed in the text before the tag. If the B3 feasibility criterion cannot be met, record "not tagged: infeasible" (a label-free stop) rather than tagging anyway.

**Checks that passed.**
- **Quoted numbers match.** Every number quoted in A6.10 matches `a6_power.out`/`.json`: walk quantiles, P(v1 minimum) by budget, every TAIL cell at κ = 1 and 0.5 for both ×1 and ×0.5 yields, and G1/G3. The certification retention figures (98.7%, 98.5%, 118/137 = 86%) match `a6_cert_yield.json`.
- **Levels follow A5.10.2.** α/2 rows are 0.975 one-sided, 0.995 two-sided, and TOST at 0.95. G2-original and G5 keep their levels. An invalid pool counts as a completed test.
- **v2_sep tier-1 conclusion is correct.** `gates_a3.py` l.141 selects V6 units with `nwin >= 2`. A3.2.4 puts V6 class 1 in tier 1, and under A3 "not evaluable blocks the tier". "Not run, therefore not a test" matches A5.10.1's definition: a completed test means step 5 ran for that scorer and some per-label quantity exists.
- **No floor margin, n_min, envelope, cap or IUT rule is relaxed** that I could find. The exceptions are M2 and M4 below, which are gaps rather than relaxations.

---

## Blocking

**B1. The swap test cannot detect the asymmetry it exists to catch** (A6.4.6.2, A6.4.3; pilot condition 2).
- **Problem.** The swap test recomputes I "from the frozen pair winners" with the certification outputs frozen. The only things that differ under a label swap are the labels, and I as defined does not use them. So the test compares a frozen file with itself.
- **What it misses.** It cannot detect the bug A6.4.3 was written to prevent: an implementation that writes only target-won pairs into the frozen I file, as the pilot did. Such a file passes the swap test bit for bit.
- **Coverage.** The test also runs on only 2 projects.
- **Why it matters.** False-pass direction (R6, pro-decoy). Pilot condition 2 asks for a swap test that can fail. As written, this one cannot.
- **Fix (text).**
  - The swap test recomputes I and R_P from the **raw** full-index rank-1 PSM table and the pair map, not from the frozen winners. It runs on **every** admitted project; the cost is set operations only.
  - Add a synthetic-label unit test in `tests/a6/`: a decoy-won pair above the threshold must be in I, and a target-won pair just below it must not be.
  - Add a label-free assertion per project: |I| equals the number of distinct pairs whose rank-1 pair winner is ≥ the run threshold. The count is computed by the certification module without splitting by label.

**B2. Step 4 computes per-label quantities on null candidates** (A6.5.2 vs A6.4.6.3, A6.1.1 P2, A6.8.1 role S).
- **Problem.** A6.5.2 sits under "Searches and mixtures (step 4; label-free)". It says each pipeline runs "exactly as on an A2 mixture … and each pipeline's own per-mixture 1% threshold from `double_competition` (as `gates_a2.collect`)".
- **Why that is a label read.** That threshold is a TDC quantity. It is computed over winners that include the null copies, so it uses decoy counts among null candidates. Agent S would compute it before the memberships are committed.
- **Why it matters.** This violates P2 and the label guard, and agent S's "must not" list. It is an invalid-procedure route that is not caught by any test, because the label guard exempts only the certification module and `mixtures_hn1.py` would import `double_competition`.
- **Fix (text).**
  - Every mixture threshold, tail membership and target/decoy split is computed only in `gates_hn1.py` at step 5.
  - Step 4 commits only label-free winners, scores and features.
  - The certification module may read labels only for **host** rows of the full-index search (rank 1 plus the pair partner). A test asserts that no label array is indexed at rows mapped from any reduced-index output.

**B3. The clock is not feasible as specified, and the pre-tag feasibility criterion ignores the clock** (A6.8.2, A6.8.3, A6.3.6, A6.10; pilot condition 6; D31.2).
- **Step 2 has 45 minutes (14:15–15:00) for:**
  - six new modules;
  - a reduced index with an exactness test;
  - the TRFP + Sage certification path;
  - self-host featurisation;
  - the dev dry run;
  - the `a6_power` re-run (about 300 s) and the cert-yield re-run (about 220 s).
- **The A4 statistics do not exist as tested code.** A4 was never implemented. `vulcan.gates.a2_stats` has no project or cluster jackknife-t for G1, no project or cluster bootstrap for TOST or G2-original, no cap subsample, no IUT precedence and no subset check. The only versions are simulation code in `reviews/checks/a4_power.py`. Untested gate statistics in a one-shot test is a validity risk as well as a feasibility one.
- **The stop rule needs full admission inside the walk.** It needs A_cand and the tail proxy, so each project must go through download, TRFP, Sage, certification, I, R_P and the reduced-index search. All of that must happen within the 60-minute walk for about 100 runs (p90 130).
  - E9's 580 runs/h is download-only.
  - The pilot did not measure TRFP, Sage or decode cost (C9 NOTES, cost section).
- **The pre-tag criterion uses the wrong budget.** It is "P(min within 400 runs) ≥ 0.8", but the binding constraint is the 16:00 deadline. P(min within 100 runs) is only 0.47 (×1 yields) and 0 (×0.5 yields).
- **Node contention.** A6 steps 4–5 (15:00–17:45) overlap the A5 freeze + G0 (by 17:00), the speed track's bench-lock timings (to 17:30) and the "quiet node" T1 look (18:00).
- **Why it matters.** Feasibility. A likely label-free stop burns about 40 projects of pool. A rushed stats module is a validity risk.
- **Fix (text).**
  1. Commission the implementer now, on parts not affected by the review.
  2. Gate v1 only (see M5).
  3. **Tag criterion.** Measure end-to-end throughput r (runs/h, including download, TRFP, Sage and GPU) in the dry run. Tag iff P(v1 minimum within r × (16:00 − tag time) runs) ≥ 0.8 at ×1 yields and ≥ 0.5 at ×0.5; otherwise "not tagged".
  4. **`a6_stats.py` acceptance test.** On fixed synthetic inputs it must reproduce the reference functions in `a4_power.py`/`a4_r1_walk.tail_iut` (verdict and bounds to 1e-12) before the tag.
  5. **Resource clause.** A6 uses a named GPU or lock that is disjoint from G0, the bench and the look. A6 work stops at 17:55 whatever its state, and the completed-test rule applies.

## Major

**M1. Integration with the A5 look is undefined** (A6.2 tag preconditions, A6.12, A5.1.3, A5.9.3).
- **Tags.** Under D31.2, `prereg-v0.4.3` and `prereg-v0.4.4` cannot exist at 14:30. So `prereg-v0.4.5` would precede them, and the A5 freeze (17:00) would contain `scripts/a6/`.
- **A6's verdict reaches the claim script after the freeze.** The verdict exists only at 17:45–18:00, after the hashed claim script is frozen at 17:00, and the look starts at 18:00. Nothing says how the claim script learns A6's state. Setting the A5.9.3 state by hand at look time is a selection route.
- **Fix.**
  - A6 writes `results/2026-10-28_A6_hn1/a6_status.json` with a closed enum: {pass, FAILED, not evaluable, invalid pool, not completed, label-free stop, not tagged, stopped: node time limit}, plus the step-5 commit sha.
  - The A5 claim script reads this file by a fixed path, committed before the freeze. If the file is absent or unverifiable when the look starts, the state is "stopped: node time limit".
  - Name the tag `prereg-a6-hn1`, or let the PI confirm the out-of-order `v0.4.5` explicitly with this sequencing.

**M2. The tail-G2 difference clause can pass falsely under differential contamination** (A6.11 R1/R4 "never a false pass", A6.10 contamination cells).
- **Problem.** UCB_diff = learned − hyperscore.
  - Residual homologs and false certifications that carry shared fragment ions but wrong retention time inflate hyperscore's tail target fraction δ_h more than v1's δ_l. v1 uses retention time.
  - The validity control only catches δ_h above roughly 0.03 at n_tail ≈ 1,250 (95% Wilson half-width ≈ 0.028).
  - So the difference clause is loosened by up to δ_h − δ_l.
- **Not simulated.** `c9_power.tail_iut` "contam" sets learned = 0.5 + 1.5δ (δ_l > δ_h), which is only the benign direction. A6.11 asserts safety only for the absolute clause and stays silent on the difference clause.
- **Fix (stricter, writable now).** Compute the difference clause against min(p̂_hyperscore, 0.5), so it is never more lenient than A2. Also report the δ_h = 0.02, δ_l = 0 cell with a true learned excess of 0.10. Correct R1/R4's "direction" column.

**M3. The power summary is optimistic, and the simulation differs from the stated rules** (A6.10 vs `a6_power.py`).
- **What matches.** The quoted values match the outputs. Pool exclusions (items 1–6), the person-key components and cluster-full skips match A6.3. Stopping at the v1 minimum is equivalent to the simulation's snapshot.
- **Differences.**
  - **(a)** "P(clean v1 passes) = P(min) × tail × G1 × G3" omits G2-original, G5, the FS rows and the subset check, and multiplies them as if independent.
  - **(b)** The deadline and throughput are not modelled (see B3).
  - **(c)** The simulation's hard cap is 2,000, and there is no 400-run budget. "Within 400 runs = 1.00" is inferred from the 300-run value.
- **Fix.** Label 0.86–0.93 as an "upper bound (tail, G1, G3 only)". Quote the 300-run figure. Add the deadline-based P(min) from B3.

**M4. The subset check is narrower than floor item 7** (A6.6.1 subset row vs A5.10.3 item 7).
- **Problem.** Item 7 covers "any instrument **or source** subset". In HN1 the source is the project. A6 checks clusters and instrument classes only, so a single project with n_tail ≥ 200 inside a larger cluster escapes the check.
- **Why it matters.** Floor compliance, and the false-pass direction for a project-specific leak diluted inside its cluster.
- **Fix.** Add "each project with n_tail ≥ 200" to the subset row, evaluated on both pools; an exceedance is a fail. Update A6.14 item 7.

**M5. Close v2_sep and ce_s outright, and spend HN1 for every scorer** (A6.2, A6.6.2, A6.6.3, A6.16.1/6).
- **v2_sep and ce_s will be dropped anyway.** The conditional "gated iff tier-1 evaluable" is false by construction, and D31.2 guarantees that the 14:30 tag precondition fails. Keeping the machinery (ce_s readiness at 16:30, V-row mapping) costs implementation time and leaves an option open, contrary to A3 Q6 / A5.10.3 item 11.
- **HN1 is not declared spent for v2_sep and ce_s.** Only "all others" are spent. A later multi-group amendment could reuse HN1 after v1's labels on the same candidates have been seen.
- **Fix.**
  - State: "A6 gates v1 only. v2_sep and ce_s are not run; A6 is not a test for them (A5.10.1). HN1 is spent for every scorer after v1's step 5 starts."
  - Keep the V6 class 1 reasoning as the stated reason. A5.9.3 state for v2_sep: "A6: not run (V6 class 1 not evaluable by construction)".

---

## Minor (one line each)

1. A6.10 "pilot f < 0.05 cell 64,507 / 23,527 / 139" differs from `a6_cert_yield` `cert_hs` (64,507 / 23,406 / 137, contaminants dropped). The 86% retention is relative to the latter; say so.
2. A6.7/A6.8.2: generalise the rule "per-label quantity exists → completed test, verdict 'not completed'" to any stop during step 5 (guard failure, crash), not only the time limit.
3. R12: give the closed rule now. Any v1 feature that counts host PSMs or proteins of a candidate is reported, with direction pro-target (fail/invalid), never a reason to drop spectra.
4. A6.8.3: name the exact reuse targets: the E9 download module (`scripts/e9/stream.py`), the A2 step-4 featurisation function, and the Sage comparator config path. Otherwise the implementer spends the 45 minutes searching.
5. A6.3.7: state that "cluster full" projects skipped before listing are *not* burnt (metadata only). This is consistent with condition 6 and saves pool.
