# Independent review of Amendment A3 draft

- **Date:** 2026-10-28
- **Reviewer:** independent statistician (agent). I had no role in E7, E7b, M7, A2's implementation, the DN2 walk or the A3 draft.
- **Draft under review:** `reviews/amendment_A3_draft.md` (untracked in the working tree, 268 lines).
- **Read:** DESIGN §1 and §6 (incl. §6.0 A1, §6.0b A2); `reviews/amendment_A2_draft.md`; `config/a2.yaml`; `scripts/a2/IMPLEMENTATION_NOTES.md`; `scripts/a2/gates_a2.py` and `build_mixtures.py` (header); RESEARCH_TREE D9–D21; E7 `NOTES.md`, `config.yaml` (incl. E7-A1) and `gates_v2.json`; M7 `NOTES.md`; 3c follow-up `NOTES.md`; T1 dry-run `NOTES.md`; E7b `config.yaml`; `vulcan/chimeric/compete.py`; `bench/t1/evaluate.py`, `bench/t1/ours_chimeric.py`; `git diff` of the working tree and of `prereg-v0.4.2..HEAD`.
- **Not read:** anything under `data/sealed`, `/tmp/vulcan_a2` or `results/2026-10-27_A2_dn2/`. I did not open the contents of the E7b outputs (`results/2026-10-28_E7b_fair_speed/t1/D1/`). I only noted that the directory exists. I ran no gate, search or sealed script.

**Hash check.** All 11 sha256 values in A3.1 match the files on disk. The A2 text and `config/a2.yaml` are byte-identical to `prereg-v0.4.2`. The hashes are correct. What is wrong is the *set* of objects they cover (B1).

## Verdict

**Do not adopt as written.** The draft's architecture is sound, and §6 below lists what should not change: per-scorer roles fixed in advance, v1 committed first, B never used to pass, a single T1 invocation, and an intersection–union test (IUT) within each scorer. But:
- six problems are **blocking**: an incomplete freeze, a gate that cannot pass, a hidden relaxation of A1.2, an unenforced v1 isolation, a T1 inclusion that can be steered after DN2, and open options in normative text;
- fourteen are **major**.

Most fixes are a sentence or a test. Two need new pre-tag work: B2 (V6 redesign with power), and B3 (a clean secondary-unit tail test, or the stricter O3 reading).

---

## 1. Blocking

### B1. The frozen object set does not cover the v2_sep pipeline, and part of that pipeline is being changed right now (A3.1)
- v2_sep's secondary precursor groups come from `vulcan.ms1` (feature finder and hypothesis table; A3.1 definition, A3.2.2). **`vulcan/ms1/isotopes.py` is modified and uncommitted in the working tree.** The change is the E7b "fast exact path". Its own docstring says it is the **default** whenever hills are contiguous.
- A3.1 hashes none of `vulcan/ms1/*`. It also does not hash:
  - the GPU kernels and engine (`vulcan/gpu/*`), which K1/K2, per-group top-K and G3 depend on;
  - `vulcan/stats/tdc.py` and `vulcan/stats/groupwalk.py`, which hold double competition, `merge_adjacent` and blocked TDC (= `q_sep`);
  - `vulcan/ptm/conga.py`, the neighbour filter;
  - `vulcan/rescore/features.py`, the v1 features that `features_v2` wraps;
  - the AlphaPeptDeep weights and prediction cache;
  - the index builder.

  Any of these can change v2_sep's DN2 or T1 output while every A3.1 hash still matches.
- The E7b replacement clause ("may replace *these files* for T1 only if … bitwise-identical … on D0, D1 and D8h") has three problems:
  - (a) it does not cover `vulcan/ms1`, which is exactly what E7b is changing;
  - (b) it allows **the gated pipeline (DN2) and the claimed pipeline (T1) to differ**, so equivalence rests on three dev sets that cannot prove equivalence on T1;
  - (c) it is silent on step 4b, which will run while E7b is in flight.
- The E7 G0 was produced at `dac5ec3` (`git_sha`), not at the commit A3 would freeze.

**Fix.**
- Freeze the pipeline by a **commit plus a clean-tree requirement**, not by a list of five files. `config/a3.yaml` records:
  - the commit;
  - `git ls-tree` hashes of `vulcan/`, `bench/t1/`, `scripts/a3/` and `scripts/e7/`;
  - the predictor weights and prediction-cache manifest;
  - the Sage binary hash.
- Step 4b, step 5 (A3 part) and the v2_sep T1 numerator must all run from that commit (`git status` clean; scripts refuse otherwise).
- E7b speed changes either:
  - land **before** the A3 tag, with the E7b equivalence evidence committed (score columns and MS1/hypothesis tables bitwise on D0/D1/D8h), and G0 is re-run at the tag commit; or
  - are excluded from both DN2 and T1.

  The draft's mixed regime (exact path on DN2, fast path on T1) is not allowed.

### B2. V6 (neighbour-filter symmetry) rests on a false premise and, as specified, fails a clean scorer almost surely (A3.2.4 row V6)
- **Factual error.** "secondary units (only they can be dropped)" is wrong. `vulcan/chimeric/compete.py` runs the CONGA filter "in descending score order over the group winners" of a spectrum, whatever their chimeric rank. A **rank-1 group winner is dropped** whenever a higher-scoring secondary winner shares > 5% of its ions. The draft's stratum handling and its argument that A is "structurally not evaluable" both need restating. (The latter conclusion is still right, because A spectra have one group.)
- **Power.** V6 copies A2.8's n_min = 2,000 per label per top-half decile and the ±0.01 TOST margin. A2.8's n_min was sized for recalls of 0.998–0.999, where the binomial variance is tiny. The filter drops 11.8% of group winners on D0 (E7 NOTES §1). With p ≈ 0.12 and n = 2,000 per label:
  - the 90% interval half-width is 1.645·√(2·0.12·0.88/2000) ≈ **0.017 > 0.01**;
  - **the interval cannot fit inside [−0.01, +0.01] even at Δ̂ = 0.** The bin is "evaluable" and fails, and under Q4/P4 **a fail is final**;
  - fitting at Δ̂ = 0 needs ≈ 5,700 per label per bin. 90% power at a true Δ = 0 needs ≈ 22,900 per label per bin (≈ 10,300 at p = 0.05).
- **Fix** (decided now, with a synthetic check):
  - (i) Define the unit correctly: every group winner of a ≥ 2-group null spectrum, stratified by the winner's chimeric class (1 vs ≥ 2).
  - (ii) Derive n_min from power at the drop rate observed on **dev** (D0/D1 or the M3 mixtures, not DN2), so that a symmetric filter passes with ≥ 0.9 probability. Below n_min the row is **not evaluable** (blocks, not final).
  - (iii) Check the rule with a simulation committed with the A3 implementation, as A2.10 did.
  - If this cannot be done before the tag, replace V6 by a reported-only statistic and list "neighbour-filter symmetry untested" as a residual risk that blocks strata ≥ 2 claims (see B3).

### B3. Strata ≥ 2 have no tail test: this is a relaxation of A1.2, not a structural fact (A3.2.4 "Consequence", A3.5.1, O3)
- A1.2 defines tail-G2 on null spectra "whose best target or decoy passes **the method's own 1% threshold**". For v2_sep, the method's own thresholds are per stratum (`q_sep`). The draft gates tail-G2 and G1 only on stratum 1 and leaves strata 2 and ≥ 3 to G2-original, G5 and entrapment.
- These are **≈ 30% of v2_sep's IDs**: on D8h, 24,258 of 80,038 at matched FDP 1% (E7 NOTES §5). A2.5 itself says G2-original and G5 do not bound a +0.05 tail bias.
- The draft calls this "the same kind of residual risk as A2.11 item 1". It is not. A2.11.1 is a narrower *population* of the gated list. Here a whole *sub-list with its own threshold* is never tail-tested.
- Secondary units are lower-evidence, prior-dominated matches: exactly the regime where a sequence-capacity leak in predicted features would show (A2.1 C(3)). Each stratum also has its own threshold in a different score region.
- Reason (i) in the draft ("`chim_rank` is not a feature") does not address this. The leak acts through the prior on features, not through `chim_rank`.
- Under A2.9.5 ("stricter reading wins") and P3, the draft's choice of O3 option 1 is the non-strict reading.

**Fix.** Choose one before the tag:
- (a) **O3 option 2** (strict): v2_sep FDR claims cover the stratum-1 list only. Strata ≥ 2 are reported as "not tail-gated", and family S is computed on the stratum-1 list. The draft says this removes most of the chimeric gain. That is a cost, not a reason.
- (b) A **pre-specified clean secondary-unit tail test**: give each A spectrum k label-free *synthetic* secondary hypotheses.
  - Draw (m/z, z) deterministically (fixed seed) from the dev distribution of secondary-hypothesis offsets relative to the instrument precursor.
  - Exclude masses within the certified precursor's isotope envelope ± tolerance.
  - Apply the A2.2(b) identity screen to each synthetic window W′ as well.
  - These are not MS1 features, so A2.2(e)'s co-isolation route cannot arise.
  - Gate tail-G2 and G1 on them at the stratum-2/≥ 3 thresholds. They are added to the mixture *after* the real `q_sep` thresholds are fixed, so they cannot move them.
  - This must be specified and simulated before the tag. If it cannot be argued a priori, take (a).

### B4. "v1 untouched" (Q1) is asserted but not enforced (A3.0 Q1, A3.7 tests)
- A3.7's guard diffs only `scripts/a2`, `vulcan/gates/a2_*` and `config/a2.yaml`. `gates_a2.py` also imports:
  - `vulcan/gates/stats.py`, `vulcan/stats/{g2,tdc,_util}.py` and `vulcan/gates/homolog.py`;
  - `vulcan/rescore/static.py` and the v1 feature builder;
  - `vulcan/gpu/*` and `vulcan/index/build.py`;
  - `scripts/m3/*`, `scripts/m3e/{gates_a1,seqclf}.py` and `scripts/m2/e2_correctness.py`.
- **`vulcan/gpu/kernels.cu`, `engine.py` and `host.py` have already changed since `prereg-v0.4.2`** (8e887f7, the 3c Da-tolerance path, claimed bit-identical for ppm). The proposed test would pass anyway. This shows the guard is too narrow.
- A3 implementation, E7b and the T1-driver edits are being written now in the same tree. Any shared-module edit changes v1's step 5 without tripping the guard.

**Fix.**
- The guard is `git diff prereg-v0.4.2 <step-5 commit>` over the **transitive import closure** of `gates_a2.py`, computed by a committed script.
- Every hunk must be either a walk-agent dated deviation, or a change with a committed bit-identity test on the v1 path (for example the 3c kernel change, with its OT-HCD bit-identity evidence cited).
- Step 5 for v1 runs on a clean tree.

### B5. T1 inclusion of v2_sep can still be chosen after seeing DN2 (A3.3.8 vs A3.3.1 items 4–5, A3.6 step 4)
- A3.3.8 makes inclusion unconditional "once A3.3.1 items 4–5 hold". But items 4–5 (driver committed, dry run, E7b freeze) are scheduled **after** step 5 (A3.6 step 4, "in any order").
- If v2_sep fails DN2, nobody is obliged to complete items 4–5. The T1 look then runs without v2_sep, which is exactly the "v2 only when it looked good" table A3.3.8 says it prevents.

**Fix.** Items 4 and 5 are **adoption preconditions** (before the tag, hence before step 5). Alternatively, state that no T1 look may occur until items 4–5 hold, whatever the DN2 outcome for any scorer. The driver commit and the E7b freeze hash go into `config/a3.yaml` at the tag.

### B6. Normative text with open options; adoption not blinded to DN2 label-free information (A3.2.1 "Optional; recommended", A3.2.7, A3.8 O1–O4, preamble)
- A pre-registration cannot contain "optional", "recommended" or open questions. O1–O4 and ce_s inclusion must be closed in the adopted text.
- The walk is in progress. Its committed records (notes 17–19 and the running label-free counts) are visible to whoever adopts A3. Those counts predict evaluability of V5, V6, V8 and V9 and of tail-G2 under the stratum-1 threshold. A PI who knows |A|, the tail proxy or the B yield can choose O2, O3 or O4 (or V6's form) by expected power.
- The draft's own statement that it "depends on no DN2 quantity" binds the author, not the adopter.

**Fix.**
- Close O1–O4 and ce_s in the text.
- The RESEARCH_TREE adoption note records which DN2 label-free quantities (if any) the adopter had seen, with their values, and that no per-label quantity existed.
- If that cannot be stated, the choices default to the stricter option in each case: O3 strict, O4 two-sided, O1 closed, ce_s excluded (ce_s licenses nothing, so excluding it costs no claim).

---

## 2. Major

### M1. The production-default rule is a "best-of-two" selection, and A3.2.6 misdescribes it (A3.2.6, A3.2.5, A3.3.7)
- "v2_sep iff v2_sep passes …; otherwise v1, if v1 passes" means **the production scorer has gate-validated FDR claims iff at least one of two scorers passes the same nulls**. That is exactly the claim Q5 forbids.
- The probability that the production default falsely carries validated claims is up to α_v1 + α_v2 (roughly, given correlated tests). The sentence "This rule depends only on v2_sep's own verdict" is false: the fallback depends on v1's.
- The rule is also inconsistent with A3.3.7:
  - A3.2.6 sets production before T1, from A3.2 plus A3.3.1 items 1–4. Item 4 is the driver dry run, which is irrelevant to production, and item 5 is missing.
  - A3.3.7 says a T1 "violated" means "Not production".

**Fix.**
- State the production rule once, in one place, with the full condition (A3.2 pass, A3.3.1 items 1–3, and no T1 "violated").
- Any statement about "the Vulcan default" carries both scorers' DN2 verdicts and the note that the default was chosen by a pre-specified fallback over two gated scorers.
- Drop the "depends only on" sentence.

### M2. The comparator for the gain claim can be steered, and the "fair" comparator is structurally weak (A3.3.3)
- As configured in E7b (`chimera: true`, `wide_window: false`, ±10 ppm, isotope errors 0–2), Sage's chimera mode re-scores after peak removal **within the same precursor window**. It cannot report co-isolated peptides at a different precursor m/z, which are exactly the units v2_sep adds through MS1 features. "Fair comparator" overstates it.
  - A1.9 asks for comparators "at their strongest".
  - §1 lists DDA+ chimeric (MSFragger) gains as the literature reference.
- Choosing C2 as the R with more matched-FDP IDs **on T1** (point estimate) and then computing the ratio CI against it ignores the selection. It is also unspecified for S3–S5: is the same R used for protein level, 5% and own q?
- `sage_chim{R}` (Sage's own LDA, all ranks) is "reported only". If it has more matched-FDP IDs than C2 on T1, the gain is being quoted against a weaker comparator. At matched FDP its calibration is irrelevant.
- E7b is already running (`results/2026-10-28_E7b_fair_speed/t1/D1/` exists). A3.3.5 says "Nothing is known yet vs C2". That must be true *at adoption*, and verifiable.

**Fix.**
- Make every S hypothesis an IUT over **all** comparators run in the look, at matched FDP: C1, C2 at R = 2, C2 at R = 5 if it qualified, `sage_chim{R}`, and the other A1 comparators. No R selection is then needed, and steering is impossible. Own-q comparisons keep the §1 "supported" condition per comparator.
- Pre-register, in an E7b amendment dated before any E7b readout, a wide-window Sage chimera variant if it is admissible for DDA. Otherwise qualify the claim as "vs Sage chimera at narrow precursor tolerance".
- The A3 adoption commit must predate any E7b ratio or verdict output, or the adoption note must state which E7b outputs existed and whether anyone involved had seen them.

### M3. Adding C2 to the T1 look interacts with v1's primary analysis (A3.0 Q2, A3.3.3)
- A1.9: "Gains are quoted against the strongest comparator." If C2 out-performs C1 on T1, is C2 the "strongest comparator" for **v1's** primary? A3 does not say.
- Either reading is a problem:
  - if yes, A3 has changed v1's primary analysis;
  - if no, the same document shows v1 against a weaker comparator than one in the adjacent table.

**Fix.** State that v1's primary comparator set is the A1-tagged set, unchanged. v1 vs C2 (and vs `sage_chim{R}`) is pre-specified as reported, with CIs, **in the primary section**. Any v1 "no loss" sentence cites it.

### M4. S1 is weaker than §1's non-inferiority (A3.3.5)
§1 computes non-inferiority "at each method's own q ≤ 1%, **and** at matched entrapment FDP". S1 uses matched FDP only, yet cites §1.

**Fix.** S1 = LB ≥ 0.97 at matched arm-P FDP 1% **and** at own q 1% against every comparator whose own-q verdict is "supported". It is an IUT across all of these.

### M5. Superiority claims on a 6-run percentile bootstrap (A3.3.4, A3.5.5)
- The draft acknowledges under-coverage and does nothing about it. For a non-inferiority primary that is conventional. For **secondary superiority claims**, under-coverage directly inflates the family-S error rate. The T1 dry run notes that with 2 runs the run bootstrap has 3 distinct outcomes. With 6 runs, the run dimension of the bootstrap is very coarse.
- Expected effect sizes vs C1 are large (≈ 1.2×). A conservative rule therefore costs little power, except vs C2, where nothing is known.

**Fix.** Pre-specify that each S hypothesis also requires one of:
- the **minimum** of the six leave-one-run-out ratios to clear the bar; or
- a t-interval on run-level log-ratios with 5 df to clear it.

Keep the percentile LB as well (an IUT, so no splitting).

### M6. No power analysis for the new rows; the joint-power statement is unsupported (A3.2.5, A3.5.2)
- "Joint pass ≈ 0.5–0.85" multiplies A2 marginal pass rates as if independent. It ignores:
  - the 4–7% invalid-pool rate of the tail-G2 control (IMPLEMENTATION_NOTES MC);
  - V5, V6, V8 and V9-secondary, whose n are unknown;
  - V6, whose pass probability is ≈ 0 as specified (B2).
- V8 (secondary-unit G2-original) has no n minimum, although §6.4 budgets ≈ 120k null competitions per class for G2. An underpowered top-decile UCB fails, and fail is final.
- A2 shipped `reviews/checks/a2_power.*`; A3 ships nothing comparable.

**Fix.**
- Add `reviews/checks/a3_power.py`, synthetic only, with sizes projected from **dev** structure (e.g. E7's 35% secondary share on the M3 mixtures). It gives the pass probability of each V row under a clean scorer and its firing rate at the margins.
- For every new row (V5, V6, V8 and V9-secondary), set an n_min below which the row is not evaluable, so it blocks without a final fail. V8's n_min follows §6.4's `g2_sample_size`.
- Report the simulated joint pass probability honestly.

### M7. The tail membership and threshold for v2_sep are ambiguous (A3.2.3 "Tail membership")
- `q_sep` assigns a stratum by the chimeric rank **of the pair winner after double competition**. The draft's tail rule uses the stratum of the *unit*. They differ whenever an A unit's peptidoform pair is won by a PSM from another unit, e.g. a host secondary group.
- The merge map (`merge_adjacent`, ≥ 200 pairs) is computed per competition call. The draft does not say:
  - whether it is per mixture;
  - whether learned, reference and noRT pipelines each use their own map;
  - what happens when stratum 2 is merged into 1, in which case "the stratum-1 threshold" is a merged-stratum threshold.
- It is not stated whether tail and G2-original winners are taken before or after the neighbour filter.

**Fix.** Write the definitions exactly:
- τ_s = the minimum pair-winner score with `q_sep` ≤ 1% in merged stratum s of that mixture, per pipeline, with each pipeline's own merge map.
- A unit is in the tail iff its stage-1 winner score is ≥ τ of the merged stratum containing **the unit's own** chimeric rank (1 for A and FS).
- G2-original, G5 and G1 use stage-1 (pre-filter) unit winners, as `gates_v2.py` did.
- Report the merge maps.

### M8. G3 does not test the pruning stage that matters for a learned rescorer (A3.2.4 V3–V5; A1.8)
- A2's G3, and so V3 and V4, measures K1 top-K retrieval of the **hyperscore-best** target or decoy. v2_sep, like v1 and ce_s, then rescores only the top `report_psms = 5` per group, ranked by hyperscore.
- Whether the **v2-best** candidate of the exhaustive set survives that truncation equally for T and D is a pruning symmetry (§6.4 G3: "any pruning … stage") that no row tests. A1.8 conditions top-K re-ranking on G3.
- This is inherited from A2 for v1, and A3 must not change v1. But A3 can add the row for v2_sep, or must at least list it as a residual risk.

**Fix.** Either:
- add V5b: on A, recall_T − recall_D of the exhaustive v2-best candidate within the production top-5, with A2.8 bins, TOST and n_min power-checked; or
- add it explicitly to A3.5 and to the required claim wording.

### M9. Missing estimator verification and calibration preconditions (DESIGN §6.7, §6.4 G4; A3.3.1, A3.3.4)
- §6.7 requires a Monte-Carlo verification (≥ 10⁴ replicates) for each estimator. `q_sep` is a composite:
  - stage-1 winners;
  - the score-ordered neighbour filter;
  - paired double competition;
  - **strata assigned from the pair winner**, which is a function of the competition outcome;
  - a label-blind merge;
  - blocked TDC.

  Its validity is plausible by symmetry but unverified. `tests/stats` has only a unit test of blocked TDC.
- **Protein groups** (S3) use picked-group FDR on sequence scores **pooled across strata**. This is the cross-stratum pooling whose peptide-level analogue (Group-walk/pooled) E7 rejected because it moves error into strata ≥ 3.
- §6.4 G4's partition (instrument, activation, charge) does not include the chimeric stratum, which is exactly the new partition.

**Fix.** Add these A3.3.1 preconditions:
- (i) a §6.7 MC of the full v2_sep competition, including null pairs whose T and D best PSMs come from units of different rank, and of pooled picked protein under stratum-heterogeneous nulls;
- (ii) the G4 partition extended to {1, 2, ≥ 3} chimeric strata, frozen before estimation;
- (iii) a protein-group entrapment verdict (arms P and C, overall) added to the A3.3.6 void rule, with S3 not tested if it is "violated".

### M10. DN2's freshness after step 5 (A3.0 Q4)
- "One shot per scorer" implies that DN2 remains available for scorers not yet built, such as a hybrid CE + static or a v3. After step 5, per-label DN2 results for three scorers are known.
- A scorer designed after that, and then gated on DN2, is gated on a set that informed its design. That is the A1 → D12 problem again.

**Fix.** "DN2 is valid only for scorers whose frozen hash is in `config/a3.yaml` at the `prereg-v0.4.3` tag. After step 5, DN2 is spent for every other scorer." If ce_s is included, say that any future hybrid needs a new fresh set.

### M11. No rule for fixing `gates_a3.py` after v1's result is committed (A3.6 3b→3c)
- v1's per-label results on the same nulls are committed (3b) before `gates_a3.py --gates` runs (3c). They are strongly correlated with v2's results (V3 is identical by construction).
- A crash or bug in new code paths (V5/V6 are new; D16 shows bugs happen) would be fixed with that knowledge.

**Fix.**
- Before 3b, run `gates_a3.py` end-to-end on real DN2 label-free inputs (finalize, counts) plus a synthetic-label dry run.
- Any change after 3b is limited to engineering that touches no statistic, parameter, pool or unit. It is written as a dated deviation and reviewed by an agent who has not seen 3b.
- Alternatively, run both `--gates` in one sealed session and commit v1's file first. The order of commit matters less than the impossibility of editing after seeing.

### M12. Missing adoption preconditions
1. **DN2 MS1 availability.** Step 4b needs each DN2 source run's MS1 (for `vulcan.ms1`). The mixtures (`*.spectra.pkl`) hold MS2 only, and the walk notes mention keeping RAWs only "until walk end". Before the tag, require a label-free confirmation from the walk agent that the MS1 scans of every admitted run are retained.
2. The D21 requirement "independently reviewed before step 5": this review and the author's responses are committed before the tag.
3. G0 (label swap) of the full v2_sep pipeline at the frozen commit (B1), not at `dac5ec3`.
4. The B6 disclosure, and the E7b status record (M2).

### M13. A3 is not "only additions" (A3.0; A2 P4 naming)
- A1.1 says "**Exactly one scorer goes into T1**." A3 puts a second, *claimable* scorer into the T1 look. That changes A1.1, not just "every other variant is exploratory".
- It is the first amendment that is not purely stricter. A1's and A2's preambles claim that property, and readers will assume it carries over.
- A2 P4 reserves the name "A3" for "a new amendment … with yet another fresh set" after a not-evaluable v1. This A3 has no fresh set.

**Fix.** Say plainly that A3 relaxes A1.1 by admitting one pre-registered secondary arm into T1, and why that is acceptable (one look, a separate family, primary first). Rename, or state that P4's "A3" refers to a later amendment.

### M14. Provenance (A3.4) is incomplete
- "D8h … nothing was chosen on it" is inaccurate. Under E7-A1, D8h was the **go/no-go set** for `sep` ("sep becomes the v2 default iff … AND (D8h) no stratum … violated"), and D21's candidate status rests on it. It was used for a decision, even if not for tuning.
- The pre-registered D8h method `vulcan_v2` (Group-walk; E7 `config.yaml` `evaluation.d8h`) was not run. A first, all-method D8h run was lost to a memory limit. State whether any partial output of that run was viewed.
- D8h has now served v1's A1.1 confirmation, E6 (ce_s), E7-A1 (sep), D21, and will serve E7b. "Reused holdout" understates this: it is a four-times-used decision set.
- Missing: what DN2 label-free information and what E7b outputs existed at drafting and at adoption (B6, M2).

**Fix.** Correct the D8h sentence and add the items above.

---

## 3. Minor
1. **Matched FDP description (A3.3.4).** The draft says "largest common cutoff at which FDP ≤ f". `bench/t1/evaluate.py` uses the **first crossing** of FDP above f with linear interpolation, on `geomspace(1e-4, 0.25, 160)` ∪ the α grid. Cite the code rather than paraphrase it.
2. **"bench/t1/evaluate.py, unchanged" (A3.3.4).** It, `comparators.py` and `t1_driver.py` are currently modified and uncommitted (the E7b `all_ranks`, chimera and `report_psms` arguments). They appear default-preserving, but say "unchanged for default arguments, verified by item 4".
3. **Sage version (A3.3.3, "0.14.x").** Pin it: tag v0.14.7, commit 99407db6, binary prints 0.14.6 (`env/ENVIRONMENT.md`). E7b says "0.14.6 binary" and DESIGN says 0.14.7. Record the binary sha256.
4. **A3.3.1 item 4, "v1 primary outputs bit-identical".** Sage is not bit-reproducible on augmented searches (dry-run finding 4). Require bitwise identity for Vulcan v1 columns and the documented tolerance (≤ 1e-4 arm-P, ±0.0004 ratio) for Sage-based comparators.
5. **"Frozen at HEAD b2c3e01" (A3.1).** The files were last changed at `d64664a`/`dac5ec3`. "A2 + D16" for `config/a2.yaml` is misleading: the file is byte-identical to `prereg-v0.4.2`.
6. **"About 10 must-pass rows" (A3.5.2).** Counting V9's three alarms and V10's three FS gates, it is ≈ 15 checks. State the count.
7. **IUT level (A3.2.5).** The argument holds for gates with a confidence bound. G5 is a point rule (AUC ≤ 0.55) with no controlled level, which is inherited. Say so.
8. **V3 assertion failure.** If V3 does not reproduce v1's G3(A) byte for byte, the outcome must be **not evaluable (implementation)** for v2_sep, not "fail".
9. **Training-overlap guard (A3.2.3).** Dropping a project "by accession" is a pruning of DN2 for v2. That is acceptable only because it is pre-specified and label-free. Say so, and report the change in the number of clusters. Also name `train_v2_folds.tsv`, the actual v2 training list.
10. **A3.3.6.** The four lists and sixteen checks do not name their unit; state peptidoform. Add protein-group overall (M9). Say how a stratum that is "not evaluable" at 1% on T1 is labelled in the claims.
11. **V4/V5 hypotheses.** Exhaustive search for G3 must use the chimeric virtual-spectrum hypotheses for B secondary groups, which `gates_a2.g3` does not provide. Specify this in `gates_a3` docstrings before the tag.
12. **ce_s (if kept).** Its G3 is v1's (same top-5 tables), so M8 applies to it too.

---

## 4. Multiplicity: summary assessment
- **Within scorer (IUT):** correct. Adding rows cannot raise the false-pass rate. No Bonferroni is right. The cost is power, which the draft has not quantified (M6), and one row is mis-powered (B2).
- **Across scorers:** per-scorer licensing and fixed roles are right. The production fallback rule re-introduces "at least one passes" (M1).
- **Family S:** fixed-sequence testing at full α and the NI → superiority order on one endpoint are sound. "vs C1 and C2" as an IUT is sound. Remaining problems:
  - the R selection on T1 (M2);
  - the comparator set too narrow (M2);
  - S1 weaker than §1 (M4);
  - interval coverage (M5).
- **Primary vs secondary:** a study-wise error rate up to the sum of the two levels is the usual convention and is disclosed. That is acceptable once M3 is fixed.
- **Per-stratum entrapment, 16 uncorrected checks:** correct, because they can only void claims.

## 5. Post-hoc and peeking: summary
- **D0/D1 tuning:** acknowledged correctly.
- **D8h:** understated (M14).
- **DN2:** no per-label leak in the design, but adoption is not blinded to label-free yields (B6). DN2 freshness after step 5 is not protected (M10).
- **T1:** no peek route found. Inclusion can be steered after DN2 (B5).
- **E7b:** running concurrently with drafting and adoption, with no timestamped separation (M2).

## 6. What is fine and should not change
- Q1's structure: v1's (c-final) membership is computed first and A3 never enters it. v2 pools are label-free subsets of v1's. v1's result file is committed before any A3 statistic. (This needs B4 to be enforced.)
- Not gating tail-G2 or G1 on B, per stratum or otherwise (A3.2.4 "Why …"). This is A2's rule, correctly kept.
- Keeping A2.2(e) (single hypothesis) for A and FS spectra. Not inventing a "clean secondary" pool from real MS1 features.
- Using hyperscore through the *identical* chimeric pipeline as the reference, on the same null competitions.
- One shot per scorer, P1–P5 inherited, and not-evaluable blocks, including the secondary-unit G5 (stricter than FS).
- `itcid_linear_pred` not gated on DN2 (D20).
- Step 4b separated from step 5 by agent. P1 swap test of the chimeric search. The no-per-label-count test extended to `scripts/a3/`.
- A single `sealed_run.py` invocation and ledger line, with no second T1 look for v2. "v2_sep never becomes T1 primary."
- The per-stratum "any violated voids all v2_sep claims" rule, mirroring E7.
- ce_s licensing nothing, if included.
- No quant or speed claim for v2_sep.
- The required claim wording in A3.3.7 (extend it per B3 and M8).
- All A3.1 sha256 values (verified).

## 7. Checklist for the revised draft
- [ ] B1: freeze by commit + clean tree + closure hashes; settle E7b before the tag; one pipeline for DN2 and T1; G0 at the tag commit.
- [ ] B2: V6 unit corrected; n_min from dev-based power; simulation committed.
- [ ] B3: O3 strict, or a pre-specified and simulated synthetic-secondary tail test.
- [ ] B4: closure-diff guard for v1.
- [ ] B5: A3.3.1 items 4–5 before step 5.
- [ ] B6: O1–O4 and ce_s closed; DN2 label-free exposure disclosed.
- [ ] M1–M14 as above; minor items 1–12.
