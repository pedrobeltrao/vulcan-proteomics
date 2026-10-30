# Amendment A3 (DRAFT, revision 1): v2_sep and ce_s in the single DN2 step-5 gate run; v2_sep as a pre-registered secondary T1 arm

- **Status:** draft for adoption as DESIGN.md §6.0c. It is not binding. DESIGN.md, the tags and every A2 file are unchanged.
- **Author:** statistician (agent). I did not build E7, E7b, M7, the A2 implementation or the DN2 walk, and I will not run step 5.
- **Date:** 2026-10-28. Revision 1 answers the independent review `reviews/amendment_A3_review.md`. Every change is listed in A3.10.
- **Basis:** DESIGN §1 and §6, including §6.0 (A1, `prereg-v0.4.1`) and §6.0b (A2, `prereg-v0.4.2`); RESEARCH_TREE D9–D21; `reviews/amendment_A2_draft.md`; `config/a2.yaml`; `scripts/a2/` (incl. IMPLEMENTATION_NOTES.md); `results/2026-10-28_E7_chimeric_v2/` (NOTES, config.yaml incl. E7-A1, gates_v2.json); `results/2026-10-28_M7_encoders/NOTES.md`; `results/2026-10-27_3c_itcid/followup/NOTES.md`; `results/2026-10-27_T1driver_dryrun/NOTES.md`; `results/2026-10-28_E7b_fair_speed/config.yaml` (E7b pre-registration, 7a2d6ee).
- **Not read by the author:**
  - any sealed file (T1–T5);
  - anything under `/tmp/vulcan_a2`;
  - any DN2 record or count (`results/2026-10-27_A2_dn2/*`);
  - any E7b output (`results/2026-10-28_E7b_fair_speed/t1/` exists, untracked, and was not opened).

  What I know of DN2: the walk (steps 2–4) is in progress, step 5 has not run, and commit titles mention "notes 18–19". No A3 rule uses a DN2 or E7b quantity.

---

## A3.0 Scope and principles

**What A3 changes.** A3 is **not** purely stricter, unlike A1 and A2. It **relaxes A1.1** ("Exactly one scorer goes into T1") by admitting one pre-registered, claimable **secondary** arm (v2_sep) into the T1 look. That is acceptable for four reasons:
- there is still one T1 look;
- the secondary arm is a separate error family (A3.3.5);
- the primary is analysed and reported first and is unchanged;
- no outcome can promote the secondary to primary.

Everything else in A3 adds analyses or makes rules stricter.

**Naming.** A2 P4 used "A3" for a possible later amendment with a fresh null set after a not-evaluable v1. This document takes the name A3 because D21 assigned it. Any such fresh-set amendment will take the next free number, and P4's meaning is otherwise unchanged.

- **Q1, v1 untouched, and enforced.** The A2 step-5 analysis of v1 (`static_gbm_pred`) is run by `scripts/a2/gates_a2.py` exactly as at `prereg-v0.4.2`, plus the walk agent's dated deviation notes.
  - Its pipelines (`fixed`, `gbm`, `gbm_noRT`), its (c-final) membership, its pools, its verdict logic and its comparator sets are unchanged. No A3 pipeline enters the v1 (c-final) union.
  - **Enforcement** (A3.7): a committed script computes the transitive import closure of `gates_a2.py` and `bench/t1_driver.py`'s v1 path and diffs it against `prereg-v0.4.2`. Every changed hunk must be either a walk-agent dated deviation or a change with committed bit-identity evidence on the v1 path. One example is the 3c Da-tolerance change in `vulcan/gpu/*` (8e887f7), which must cite its OT-HCD ppm bit-identity evidence.
  - Step 5 runs on a clean tree.
- **Q2, T1 primary untouched.** v1 is the only T1 primary. Its comparator set for every primary claim is the A1-tagged set: Sage (RT on), Sage+mokapot, Sage+pred+mokapot, Comet+mokapot. v1 vs the chimeric comparators is reported, with CIs, **in the primary section** (A3.3.3). Any v1 "no loss" sentence cites those numbers next to the primary ones.
- **Q3, inherit A2.** P1–P5 of A2 apply to every A3 analysis. Every A1 and A2 margin, n-minimum, seed, bootstrap and validity control is reused unchanged.
- **Q4, one shot per scorer; DN2 spent afterwards.** DN2 is the first and only tail-G2/G1/G3 attempt for v2_sep and for ce_s. A fail is final under this protocol. Not evaluable keeps FDR claims blocked. **DN2 is valid only for scorers whose frozen hash is in `config/a3.yaml` at the `prereg-v0.4.3` tag.** After step 5, DN2 is spent for every other scorer, including any hybrid of ce_s and static or any v3.
- **Q5, roles fixed now.** Each verdict licenses claims for its own scorer only. No claim of the form "a Vulcan scorer passes" is made. No scorer takes another's role after seeing results. The single production rule (A3.2.6) is a disclosed, pre-specified fallback, and both verdicts are quoted whenever it is invoked.
- **Q6, no open options.** Every choice below is closed. Where the author's first draft offered options, revision 1 takes the stricter one, or specifies an a-priori test (A3.10).

---

## A3.1 Adoption preconditions (all before the tag `prereg-v0.4.3`, hence before step 5 and before any T1 look)

1. **Text.** The PI adopts A3 as DESIGN §6.0c, with a dated RESEARCH_TREE note. That note **discloses** two things:
   - which DN2 label-free quantities (walk notes, running counts) the adopter had seen, with their values, and that no per-label DN2 quantity existed;
   - which E7b outputs existed, and whether the adopter or the author had seen any E7b ratio or verdict.

   If the adopter cannot make the disclosure, A3 is adopted as written anyway. Every choice in it is already the strict one, or a pre-specified test.
2. **Review.** This draft, `reviews/amendment_A3_review.md` and A3.10 are committed (D21: "independently reviewed before step 5").
3. **Pipeline freeze (by commit, not file list).** `config/a3.yaml` records:
   - the freeze commit;
   - `git ls-tree` object hashes of `vulcan/`, `bench/`, `scripts/a3/`, `scripts/e7/`;
   - the AlphaPeptDeep weights hash and prediction-cache manifest hash;
   - the Sage binary sha256 (`tools/bin/sage`: tag v0.14.7, commit 99407db6, prints 0.14.6);
   - the mokapot wrapper hash;
   - the model hashes of A3.2.

   Step 4b, the A3 part of step 5, and the T1 v2_sep numerator and comparators all run **from the freeze commit, on a clean tree**. The scripts refuse otherwise. **DN2 and T1 use the same v2_sep implementation.** The E7b speed work (currently uncommitted changes to `vulcan/ms1/isotopes.py`, `vulcan/rescore/fast.py` and a new `vulcan/chimeric/fast_compete.py`) must do one of two things:
   - land before the freeze, with E7b's equivalence evidence committed: bitwise-identical MS1 feature tables, hypothesis tables and v2_sep score columns, and identical ID sets, on D0, D1 and D8h; or
   - stay off both the DN2 and T1 paths. That is enforced by the clean-tree rule and by `config/a3.yaml` naming the code path.
4. **G0 at the freeze commit** for the full v2_sep pipeline: label swap on D1 + E2001, as in `gates_v2.py`, with bit-identical chimeric PSM tables, MS1 hypotheses, features and v2 scores. E7's G0 was produced at dac5ec3 and does not count.
5. **MS1 availability.** The walk agent confirms, label-free, that the MS1 scans or stores of every admitted DN2 run are retained, so that step 4b can run `vulcan.ms1`. If they are not, V5, V6, V8, V9-sec, V13 and V14 (A3.2.4) are not evaluable.
6. **T1 driver.** It is committed with the v2_sep numerator and every comparator of A3.3.3. A dry run on D0 and D1 shows:
   - v1 columns **bitwise identical** with and without the A3 additions;
   - Sage-based comparators within the documented non-reproducibility (≤ 1e-4 in arm-P means, ±0.0004 in ratios; dry-run finding 4).
7. **E7b comparator frozen.** The configuration is committed and its commit recorded in `config/a3.yaml`. If E7b pre-registers a wide-window Sage chimera variant, that must happen in a dated E7b amendment **before any E7b readout**. Otherwise no such variant exists, and the claim wording carries "Sage chimera at narrow precursor tolerance".
8. **Dev preconditions for v2_sep claims** (results may be pass or fail; inclusion does not depend on them, A3.3.8):
   - (i) **G4** for v2_sep with the partition extended to chimeric stratum {1, 2, ≥3} × (instrument, activation, charge), frozen before estimation, on training/dev decoys;
   - (ii) **arm M** for v2_sep on D7 (yeast + human entrapment, r ≈ 1, seeds 2001–2003, chimeric search). Rule: no "violated" at 1% or 5%, overall or in any stratum;
   - (iii) **§6.7 Monte-Carlo verification** (≥ 10⁴ replicates, pass if mean − 2.326·SE ≤ α) of the full v2_sep competition. It covers stage-1 winners, the score-ordered neighbour filter, double competition, strata assigned from the pair winner (including null pairs whose target and decoy best PSMs come from units of different rank), label-blind merge and blocked TDC. It also covers pooled picked protein groups under stratum-heterogeneous nulls. Each stratum list and the union must control mean FDP.
9. **Power check** `reviews/checks/a3_power.py` (synthetic data only, seed 131, sizes projected from **dev** structure: E7's G2 mixtures, 35% secondary units, D0 neighbour-filter drop rate 11.8%). It must report:
   - the pass probability of every V row under a clean scorer and its firing rate at the margins;
   - the n_min values of V6, V11, V12 and V14, derived by the rules in A3.2.4;
   - the simulated **joint** pass probability for each claim tier.

   If the simulation shows a row can never pass at any n reachable within A2's DN2 cap (16 projects), that row is still kept, and the tier depending on it is recorded in advance as "expected not evaluable".
10. **A3 implementation** (A3.7) is committed and tested on synthetic data, including an end-to-end `gates_a3.py` run on synthetic labels.

If step 5 for v1 starts before the tag, A3.2 lapses: v2_sep and ce_s are not gated on DN2, and A3.3 lapses.

---

## A3.2 Step 5 (part a)

### A3.2.1 Scorers and roles (closed)

| Scorer | Frozen identifier | Role | What its verdict licenses |
|---|---|---|---|
| v1 `static_gbm_pred` | `results/2026-10-27_M3_scoring/models/static_gbm_pred.pkl`, sha256 `31bfe685830b91f4403d8a91c78426e84eca3cdd84157a1d1a9c5b0e6c9126c0` | Per A2, unchanged (Q1) | As A2.9.6 |
| **v2_sep** | model `results/2026-10-28_E7_chimeric_v2/models/static_v2.pkl`, sha256 `3f2c69cc6cad002e42c51475dc1ad52b9bfbbc84f3dd51997b22cb24a70f0186` (frozen 6607ecd); pipeline pre-registration `results/2026-10-28_E7_chimeric_v2/config.yaml`, sha256 `de018e957deceec541b2e2c853a4072ece04c4cd6a9ed587d001f897014ce4ed` (8be8f45 + E7-A1 c7bdaae); code at the freeze commit (A3.1.3) | **Gated** (A3.2.3–A3.2.4) | Tier-1 or tier-2 FDR claims for v2_sep (A3.2.6), subject to A3.1.8 |
| **ce_s** | `results/2026-10-28_M7_encoders/models/ce_s.pt`, sha256 `aa19c8f3fbec07d7532a005132efc689244d6e75bbd8ac69f9df24209f9a97c6` (0d9d795) | **Gated, reported-only role** (A3.2.7) | Nothing. Final for ce_s; ce_s stays a NO-GO research arm (D19) |
| `itcid_linear_pred` | — | Not gated (DN2 is OT-HCD; D20 stands) | — |
| v2 (Group-walk), chim_v1, chim_v1_sep, RESET arms | — | Not gated; DN2 is spent for them (Q4) | — |

For reference, the code files as of the draft (all last changed at d64664a, except `ours_chimeric.py` at dac5ec3 and `static.py` at e3650c4) are:
- `vulcan/chimeric/search.py` `d1a2af54…`
- `vulcan/chimeric/compete.py` `454e2f97…`
- `vulcan/rescore/features_v2.py` `1fbfff07…`
- `vulcan/rescore/static.py` `de4b1fd9…`
- `bench/t1/ours_chimeric.py` `10324a9f…`

These are informational only. The binding freeze is the commit plus tree hashes (A3.1.3). The sequence classifier is `4ec134eb…`, as in `config/a2.yaml`, which is byte-identical to `prereg-v0.4.2`.

**v2_sep** means:
1. E7 chimeric search;
2. stage-1 group winners;
3. the CONGA neighbour filter, over **all** group winners of a spectrum in descending score order, so a rank-1 winner can also be dropped;
4. paired double competition (peptidoform units);
5. strata = chimeric rank of the pair winner {1, 2, ≥3}, merged label-blind at < 200 pairs;
6. q = TDC within each merged stratum (`q_sep`).

### A3.2.2 Step 4b: label-free chimeric search of the DN2 mixtures

This is a new step between step 4 and step 5. It is run by an agent that does not run step 5, from the freeze commit.
- **Hypotheses.**
  - Host (D0) and B-candidate spectra get E7 chimeric hypotheses. The primary group is A2's production hypothesis set. Secondary groups come from each spectrum's **source** run's `vulcan.ms1` hypothesis table (E7 defaults), as in `scripts/e7/gates_v2.py`.
  - **A and FS spectra keep A2.2(e)**: one hypothesis, no real secondary groups.
  - A1-recipe offset copies shift every group (reported only).
- **Synthetic secondary hypotheses on A spectra (pool "A-syn"; used only by V13–V14; see A3.2.4).**
  - Each A spectrum gets two synthetic hypotheses, assigned chimeric rank 2 and 3. Their (Δm/z, z′) pairs are drawn deterministically from the **frozen dev table** `config/a3_secondary_offsets.tsv`. That table is built before the tag, label-free, from all secondary groups of the D0 and D1 MS1 hypothesis tables: Δm/z = secondary mono m/z − instrument precursor m/z, and the secondary charge.
  - The draw is `hash(spectrum key, rank, seed 47)` into the table's empirical CDF.
  - **Exclusions.** A draw is replaced by the next draw (same hash chain, up to 32 tries; if none qualifies, the spectrum has no synthetic hypothesis of that rank) if:
    - the drawn neutral mass lies within ±20 ppm of the certified precursor's isotope envelope (k = −1…+3 isotopes at any z ≤ 6); or
    - the drawn window W′ fails A2.2(b) (identity screen of p against (T∪D) ∩ W′, f ≥ 0.20, or Hamming ≤ 2) or A2.2(d) (family screen: max hyperscore over (T∪D) ∩ W′ ≥ hyperscore(p)).
  - These are label-free and symmetric (P1). They are not MS1 features, so A2.2(e)'s co-isolation route cannot arise. Every human candidate in W′ is a chance match to the certified non-human fragments of p.
  - A-syn units are searched as separate virtual spectra. **They do not enter the mixture's competition**, so they cannot move any real threshold.
- **Search and scoring.** `vulcan.chimeric.search` with the E2 production SearchParams and `report_psms` 5 per group, on the A2 prod gate index. Then v2 features (primary-only calibration, computed on real primary PSMs only, so A-syn units never enter calibration). Scores: v2, hyperscore, and v2 with `neutral_rt` (noRT, reported).
- **For V12 only:** on A spectra, the v2 and ce_s scores of **every** candidate in the single-hypothesis window W (exhaustive), not just the top 5.
- **Commit.** Per-mixture sha256 of the chimeric PSM tables, the MS1 hypothesis tables of every DN2 source run, the A-syn hypothesis table, label-free unit counts by rank and pool, and a P1 swap test of the chimeric search and of A-syn construction on the first mixture. The no-per-label-count test (`tests/a2/test_no_label_counts.py`) is extended to `scripts/a3/`.

ce_s needs no new search beyond V12. It scores the step-4 closed top-5 tables, as in M7.

### A3.2.3 Pipelines, reference, thresholds and membership for v2_sep (exact definitions)

- **Null query unit:** (null spectrum, precursor group). On A and FS every real unit is rank 1.
- **Pipelines, per mixture:**
  - learned: v2_sep;
  - reference: hyperscore through the **identical** pipeline (stage-1 winner by hyperscore, neighbour filter ordered by hyperscore, double competition, strata, merge, `q_sep`);
  - noRT: v2 with `neutral_rt` (reported).
- **Merge map:** per mixture and **per pipeline**, each pipeline's own `merge_adjacent` map. All maps are reported.
- **Thresholds:** τ_s(pipeline, mixture) = the minimum pair-winner score with `q_sep` ≤ 1% in merged stratum s, on that mixture.
- **Tail membership:** a null unit is in the tail iff its **stage-1 (pre-filter) unit winner's** score is ≥ τ of the merged stratum that contains **the unit's own chimeric rank**. That rank is 1 for A and FS, and 2 or 3 for A-syn.
- **Unit winners for G1, G2-original and G5:** stage-1 (pre-filter) unit winners, as in `gates_v2.py`.
- **(c-final) membership for v2_sep:**
  - start from the committed A2 (c-final) membership (v1's);
  - additionally drop a null spectrum if, for any of its units (real or A-syn), the stage-1 winner or its pair partner under the learned, reference or noRT pipeline is A1.3-near the project's reference (the same `HomologFilter` reference set as A2).

  This is label-free. v2 pools are subsets of v1's, and v1's are unaffected. The membership file and its sha256 are committed before any v2 statistic, behind a guard like `gates_a2.py`'s.
- **Training-overlap guard.** Assert that no DN2 accession is in `data/manifests/train.tsv`, `data/manifests/train_m7.tsv` or `results/2026-10-28_E7_chimeric_v2/train_v2_folds.tsv`. An overlapping project is dropped from the v2_sep and ce_s pools by accession. That is a pre-specified, label-free pruning. The change in the number of clusters is reported. None is expected.

### A3.2.4 Gate rows for v2_sep

Every row uses A1/A2 rules unchanged unless stated. Not evaluable blocks the tier the row belongs to, except FS rows (A2.4 rule). The unit is the null query unit.

| # | Gate | Pool / units | Rule | Tier |
|---|---|---|---|---|
| V1 | tail-G2 | A, rank-1 | A2.5: n_tail ≥ 200; UCB_abs ≤ 0.60; UCB_diff ≤ 0.10 vs reference; cluster bootstrap B = 10,000, seed 83; reference validity control | 1 |
| V2 | G1 | A | A2.6: n_top ≥ 7,200; 99% CI ⊂ (−0.05, 0.05); reference control CI ∋ 0; frozen classifier | 1 |
| V3 | G3 | A | A2.8. Scorer-independent: it is v1's G3(A), recomputed by `gates_a3` and asserted byte-identical. **A mismatch makes V3 "not evaluable (implementation)", never "fail".** | 1 |
| V4 | G3 | A ∪ B, all query units | A2.8 TOST and n_min 2,000. Exhaustive search uses the same chimeric virtual-spectrum hypotheses as production (specified in the `gates_a3` docstring before the tag) | 1 |
| V5 | G3 | A ∪ B, secondary units (rank ≥ 2) | as V4 | 2 |
| V6 | G3-NF, neighbour-filter symmetry | every stage-1 group winner of A ∪ B null spectra with ≥ 2 groups, split into class 1 (rank-1 winners) and class 2 (rank ≥ 2) | Δ = P(kept \| target winner) − P(kept \| decoy winner). Two bins per class: top decile of winner v2 score, and the rest. 90% Newcombe interval ⊂ [−0.01, +0.01] (two-sided, as A2.8). **n_min per label per bin = the smallest n at which a symmetric filter at the dev drop rate passes with probability ≥ 0.90** (A3.1.9; about 23,000 at p = 0.12). Below n_min: not evaluable | class 1 → 1; class 2 → 2 |
| V7 | G2-original | A ∪ B, all units | §6.4: one-sided 95% UCB ≤ 0.01 overall, ≤ 0.02 top decile, vs reference | 1 |
| V8 | G2-original | A ∪ B, secondary units | same margins. **n_min = `g2_sample_size`** (§6.4 budget); below it, not evaluable | 2 |
| V9 | G5 | (a) A ∪ B all units; (b) A; (c) A ∪ B secondary units | AUC ≤ 0.55 (a point rule with no controlled level, inherited from §6.4). (c) is evaluable iff ≥ 200 target and ≥ 200 decoy winners | (a), (b) → 1; (c) → 2 |
| V10 | FS gates | FS (rank-1) | A2.4 rule 3 unchanged: tail-G2, G1, G5 if evaluable; an evaluable fail blocks; not evaluable does not block | 1 |
| V11 | — | (number intentionally unused) | | |
| V12 | top-5 truncation symmetry | A, rank-1 | Per A spectrum, the exhaustive **v2-best** target and decoy in W, and whether each is inside the production top-5 (hyperscore-ranked). recall_T − recall_D, A2.8 bins (deciles of the exhaustive v2-best score), 90% Newcombe ⊂ [−0.01, +0.01]. n_min from A3.1.9 by the V6 rule. This tests the pruning stage that a top-K rescorer adds (A1.8, §6.4 G3) | 1 |
| V13 | tail-G2, secondary | **A-syn** (rank 2 and 3 pooled) | A2.5 rule and margins. Tail = unit winner ≥ τ of its rank's merged stratum. Reference validity control: hyperscore's A-syn tail fraction Wilson CI ∋ 0.5. n_tail ≥ 200 | 2 |
| V14 | G1, secondary | A-syn | A2.6 rule; n_top ≥ 7,200; reference control | 2 |

**Reported, never used to pass:**
- tail-G2 and G1 on B per stratum (1, 2, ≥3, secondary pooled);
- V13 and V14 split by rank 2 vs 3;
- per-project and leave-one-project-out for V1 and V13;
- every row with noRT;
- V7, V8 and V6 split by stratum;
- A1 offsets.

**Why B is not used to gate the tail.** A2.2 says B is never used to pass. Gating on B would either pass a scorer on a contaminated pool, or block it on a pool A2 calls invalid.

**Why A-syn is a valid clean secondary pool.** Each A-syn unit is a certified non-human spectrum (A2.2(a)) whose candidates at the synthetic window W′ passed the same identity and family screens as A (A2.2(b), (d)). So every human candidate is a chance match, exactly A's argument. It tests the question that matters for strata ≥ 2: whether the v2 score prefers targets among chance matches at the stratum-2/≥3 thresholds. **Limitation:** real secondary units come from MS1 features of co-isolated peptides, which shape which peaks remain after the primary's ions. A-syn reproduces the mass offsets and charges of real secondaries, but not their co-isolation. This is stated in every tier-2 claim.

### A3.2.5 Multiplicity

- **Within a scorer: no adjustment.** Every row with a confidence bound tests the null "a leak or asymmetry of at least the margin exists". A tier passes only if all its rows reject (intersection–union test). The false-pass probability for a scorer that violates one margin is at most that row's level. This holds for confidence-bound rows; G5 is a point rule with no controlled level, as already in §6.4.
- **Power** is the cost. It is quantified by A3.1.9 and reported. A3 makes no unsimulated claim about the joint pass probability. Tier 1 has about 13 checks (V1–V4, V6 class 1, V7, V9a/b, V10 × 3 if evaluable, V12). Tier 2 adds 7 (V5, V6 class 2, V8, V9c, V13, V14, and the stratum-2/≥3 split being reported only).
- **Across scorers.** Each verdict licenses only its own scorer. ce_s licenses nothing. The production rule (A3.2.6) is a pre-specified fallback and is disclosed as such: its probability of carrying falsely validated claims is up to α_v1 + α_v2.
- **Tiers.** Tier 2 ⊃ tier 1, so the claim set is nested and monotone. The tier is determined by gate outcomes under rules fixed here. Family S (A3.3.5) is computed on exactly one list, determined by the tier.

### A3.2.6 What a v2_sep result licenses; production rule (stated once)

- **Tier 2 passes** (all tier-1 and tier-2 rows): tail-G2/G1/G3 are satisfied for the **full v2_sep list (all strata)**.
- **Only tier 1 passes:** they are satisfied for the **v2_sep stratum-1 list only**. This is the list of rank-1 pair winners under the stratum-1 `q_sep`, which is controlled at α by its own TDC. Strata ≥ 2 are labelled "not gate-validated".
- **Tier 1 fails:** final for v2_sep under P4. No v2_sep FDR claim.
- **Tier 1 not evaluable:** no v2_sep FDR claim. Only a later amendment with a new fresh set can follow, and it reports this result.
- **Production default for OT-HCD DDA:**
  - **v2_sep** iff all of: tier 2 passes; A3.1.8 (i)–(iii) hold with no "violated" arm M; and v2_sep has no "violated" verdict on T1 (A3.3.6);
  - otherwise **v1**, iff v1 passes A2;
  - otherwise no OT-HCD scorer has gate-validated FDR claims.

  Whenever the production default is named, both scorers' DN2 verdicts are quoted, together with the sentence "default chosen by a pre-specified fallback over two scorers gated on the same nulls". Until T1 has been looked at, v2_sep can be at most the **provisional** default.

### A3.2.7 ce_s

- **Decision: included.** It licenses nothing, so its inclusion cannot be steered towards a claim. It also answers a question the D12 nulls could not: does a sequence-capacity model leak on clean nulls?
- **Analysis:** exactly A2's for v1, with the ce_s logit as the score: closed pipeline, paired double competition per mixture, pooled 1% threshold. Rows: tail-G2(A), G1(A), G2-original and G5 on A ∪ B and on A, FS if evaluable, and V12 (top-5 truncation symmetry for the ce_s-best candidate).
- **G3:** v1's G3(A) and G3(A ∪ B), which are scorer-independent with the same units.
- **(c-final) membership:** v1's, minus spectra whose ce_s winner or pair partner is A1.3-near.
- **Verdict:** final for ce_s (P4). No noRT ablation, since ce_s uses no RT.

### A3.2.8 Reporting
`results/<date>_A2_gates/NOTES.md` gets one section per scorer: v1 first, unchanged in form; then v2_sep, with every row, its tier, the merge maps and the A-syn construction counts; then ce_s. Everything is reported whatever the outcome. The RESEARCH_TREE decision records every verdict. **Residual risk for v1 and ce_s, unchanged by A3:** the top-5 truncation symmetry (V12) is not a gate for v1. For v1 it may be reported if computed, but it never changes v1's A2 verdict.

---

## A3.3 v2_sep as a pre-registered secondary T1 arm (part b)

### A3.3.1 When v2_sep claims are possible
- **The T1 look itself.** It stays governed by D12 and A2.9.6: it requires v1 to pass A2.
  - **If v1 fails or is not evaluable, there is no T1 look under A3, whatever v2_sep's result.**
  - A3 does not pre-specify a v2-only route. Any such route needs a new amendment adopted before any T1 look, and that amendment must report v1's failure (closed; formerly O1).
- **v2_sep's claims on T1 need all of:**
  - tier 1 or tier 2 passed (A3.2.6);
  - A3.1.8 (i)–(iii) hold, with no "violated" verdict;
  - no "violated" verdict on T1 (A3.3.6).

### A3.3.2 Arm and data
- **Numerator:** v2_sep, from the freeze commit, on T1 (ProteoBench DDA-QExactive, 6 runs).
- **Look:** the same 25 replicate DBs as the primary (5 prod; 15 arm-P with sealed seeds 2004–2006; 5 arm-C), in **the same `sealed_run.py` invocation and ledger line**. There is no second look.
- **Runs without MS1 hills** are kept; their chimeric search equals the closed search. Their number is reported.

### A3.3.3 Comparators
- **All comparators run in the look, and family S tests v2_sep against every one of them** (an IUT; no selection):
  - C1 = Sage+pred+mokapot (A1.9);
  - C2 = Sage chimera + pred + mokapot as frozen by E7b (A3.1.7). This uses R = 2, and R = 5 as well if R = 5 qualified under E7b's own dev rule (both are then run and both are comparators);
  - `sage_chim{R}` (Sage's own LDA, all ranks, for each R run);
  - Sage (RT on), Sage+mokapot, Comet+mokapot.
- Because every comparator must be beaten, there is no choice of "the strongest" and nothing to steer.
- **Own-q comparisons** use only comparators whose own-q verdict is "supported" (§1).
- **v1 vs C2 and v1 vs `sage_chim{R}`** are reported in the primary section (Q2). They are not v1 primary claims.
- If A3.1.7 is not met, the adoption precondition fails and there is no A3.

### A3.3.4 Units, estimators, CIs
- **Estimator.** `bench/t1/evaluate.py` at the freeze commit. Its default-argument behaviour must be unchanged, as verified by A3.1.6. It is cited, not paraphrased:
  - the matched-FDP threshold is the **first crossing** of mean arm-P FDP̂ above f, with linear interpolation, on the code's grid (`geomspace(1e-4, 0.25, 160)` ∪ the α grid);
  - for v2_sep the grid variable is a common `q_sep` cutoff applied to every stratum.
- **Units.** Peptidoform (I/L, Met-ox positions) is primary; sequence and picked protein groups are secondary.
  - For v2_sep, protein groups use picked-group FDR over sequence-level scores **pooled across strata**, as E7 did. Its validity is a precondition (A3.1.8 iii), and it has its own entrapment void rule (A3.3.6).
- **Ratio.** Run-level IDs summed over runs and averaged over 5 decoy shuffles. 95% percentile bootstrap over runs × shuffles (2000 draws, seed 11). "LB" means the lower end of this interval.
- **Few-cluster backstop** (6 runs; required for every S hypothesis in addition to the LB, as an IUT): the **minimum of the six leave-one-run-out point ratios** must also clear the bar.

### A3.3.5 Family S (FWER 0.025 one-sided within S; fixed sequence)

The list is the full v2_sep list if tier 2 passed, or the stratum-1 list if only tier 1 passed (A3.2.6). The list is fixed by gate outcomes under rules set here. Each hypothesis is tested only if all earlier ones were rejected. Each must clear the bar against **every** comparator of A3.3.3, by both the LB and the leave-one-run-out minimum.

| Order | Hypothesis | Unit, threshold | Bar |
|---|---|---|---|
| S1 | non-inferiority | peptidoform, matched arm-P FDP 1% **and** own q 1% (own q: "supported" comparators only) | ≥ 0.97 (§1) |
| S2 | more IDs | peptidoform, matched arm-P FDP 1% | > 1.00 |
| S3 | more protein groups | protein groups, matched arm-P FDP 1% (not tested if v2_sep's protein-group verdict is "violated", A3.3.6) | > 1.00 |
| S4 | more IDs at 5% | peptidoform, matched arm-P FDP 5% | > 1.00 |
| S5 | more IDs at own q | peptidoform, own q 1%, "supported" comparators only. If none is supported, S5 is not tested | > 1.00 |

- **Reported, not tested:** sequence level, 0.1%, v2_sep / v1, per-stratum IDs, speed (§2), quant (ε, CV). **There is no quant or speed claim for v2_sep.**
- **Family S is controlled within itself.** The study-wise error across the primary and S families is up to the sum of the two levels. This is disclosed, and every v2_sep claim is labelled secondary.
- **Dev expectation, not a criterion:** on D8h, v2_sep vs C1 was 1.222 [1.138, 1.279]. Nothing is known to the author vs C2.

### A3.3.6 Entrapment verdicts on T1 (peptidoform unit unless stated)
- **Computed for:** v2_sep at α ∈ {1%, 5%}; arms P (15 replicates) and C (5 replicates); CI convention (a) (A1.6); for the overall list and each merged stratum (1, 2, ≥3), plus the **protein-group overall list**.
  - Stratum FDP̂ uses counts restricted to the stratum, with the paired target's acceptance taken from the full list, as in `scripts/e7/strata.py`.
  - The merge map is the frozen label-blind ≥ 200-pair rule and is reported.
- **Void rule:**
  - **Any "violated"**, in any arm, in any peptidoform list (overall or any stratum), at 1% or 5%: all v2_sep T1 claims are void, recorded as final.
  - **"Violated" in the protein-group overall list:** S3 is not tested, and the sequence stops at S3.
  - This applies whichever tier passed: a violation in a stratum ≥ 2 voids tier-1 claims too.
- **"Not evaluable"** (α·N_T < 20) in a stratum is allowed. That stratum is labelled "not evaluable at α" in every claim, and no stratum-specific statement is made for it.
- **Otherwise** each list carries its §6.3 label: supported, supported-noisy, or inconclusive with "consistent with control" when the paired lower CI is ≤ α. "Consistent with control" never ranks scorers.
- **Multiplicity:** 2 arms × 5 lists × 2 levels = 20 "violated" checks, uncorrected. That is deliberate: each check can only void claims.
- **Comparators:** verdicts as in the primary. The arm-P split by Sage rank for `sage_chim{R}` and C2 is reported (E7b).
- **Arm M** for v2_sep on T5 at the milestone (A1.5).

### A3.3.7 Claim table and wording

| v1 on DN2 | v2_sep tier + A3.1.8 | v2_sep T1 "violated" | Allowed |
|---|---|---|---|
| fail / not evaluable | any | — | No T1 look under A3 |
| pass | tier 1 failed or not evaluable, or A3.1.8 failed | — | v1 primary as A1. v2_sep is run and reported in a separate table titled "exploratory, not gate-validated". No v2_sep claim. Not production |
| pass | tier 1 or 2 | yes | v1 primary as A1. v2_sep violation reported as final. No v2_sep claim. Not production |
| pass | tier 1 only | no | v1 primary. Family S on the **v2_sep stratum-1 list**. Strata ≥ 2 reported as "not gate-validated". Not production |
| pass | tier 2 | no | v1 primary. Family S on the full v2_sep list. Production default per A3.2.6 |

**Required wording for any v2_sep claim:** "pre-registered secondary arm; gated on fresh nulls under A2/A3 (tier N); tail tests on certified non-human primary units [and, tier 2, on synthetic secondary hypotheses that reproduce secondary mass offsets but not co-isolation]; top-5 truncation symmetry gated; compared against every comparator in the look, including Sage chimera at narrow precursor tolerance [unless a wide-window variant was pre-registered]". v1's result is always reported first and in full in the same document.

### A3.3.8 Inclusion is not outcome-dependent
A3.1.6–A3.1.7 are adoption preconditions, so the driver and comparators exist before step 5. If the T1 look happens, v2_sep and every A3.3.3 comparator run in it, whatever v2_sep's DN2 and dev-precondition outcomes. Only the claims depend on those outcomes.

---

## A3.4 Provenance of the v2_sep dev evidence (part c)

- **D0 and D1 informed v2_sep.**
  - The `sep` competition (E7-A1) was added **after** the D0/D1 per-stratum results of `chim_v1` were seen, because Group-walk had violated stratum ≥3 there.
  - D0/D1 also served as the first half of E7-A1's go/no-go rule.
  - All D0/D1 numbers for v2_sep and chim_v1_sep (e.g. 1.211× and 1.171× vs C1, per-stratum verdicts) are therefore **tuned selection-set numbers**. They are never quoted as evidence.
- **What was not chosen on D0/D1:**
  - the v2 model, chosen by leave-project-fold-out CV on 40 training runs and frozen before any dev evaluation (6607ecd);
  - the chimeric search parameters, neighbour filter, strata and merge rule, fixed before any chimeric run (8be8f45).
- **D8h is the only dev set on which v2_sep was not tuned, but it is not unused:**
  - It was the **go/no-go set** for `sep` under E7-A1 ("sep becomes the v2 default iff … AND (D8h) no stratum … violated"). D21's candidate status rests on it. So D8h was used for a decision, though not for tuning.
  - It is a reused decision set: v1's A1.1 confirmation, E6 (ce_s), E7-A1 (sep), D21, and now E7b (comparator evaluation).
  - The pre-registered D8h method `vulcan_v2` (Group-walk) was not run. A first all-method D8h run was lost to a job memory limit. E7 NOTES do not say whether any partial output of that lost run was viewed. **The E7 author is asked to state this in the adoption note.**
  - The decision to propose A3 came after D8h, so D8h's 1.222× **motivated** the secondary arm and carries a follow-up optimism bias. **Only T1 is confirmatory.**
- **DN2 hosts are D0 runs.** `sep` was chosen after seeing D0, and D0 spectra fill every DN2 mixture's host part. This does not leak null labels, but the real per-stratum thresholds on the mixtures are set partly by a run the pipeline was tuned on. The same holds for v1 under A1/A2.
- **Disclosure at adoption:** the DN2 label-free and E7b exposure (A3.1.1).

---

## A3.5 Residual risks
1. **Tier-2 tail validity rests on A-syn.** A-syn matches secondary mass offsets and charges, but not co-isolation (A3.2.4).
2. **Power.** The rows are many, and some (V6, V12, and V14 at n_top ≥ 7,200 on A-syn) may be not evaluable at DN2's size, which was sized for v1 (A2.9.3). A3.1.9 quantifies this before the tag. The outcome may well be tier 1 only, or nothing.
3. **B-based rows** (V4–V9) use the A1-type pool. Contamination there makes the rows stricter, and can block a clean scorer.
4. **V6 and V12 are two-sided** (as A2.8). A pro-decoy, FDR-conservative asymmetry fails them too (closed; formerly O4).
5. **T1 has 6 runs.** The leave-one-run-out backstop is conservative but coarse.
6. **Family S is controlled within itself only** (A3.3.5).
7. **Top-5 truncation (V12) is gated for v2_sep and ce_s, but not for v1** (Q1). v1's A2 verdict stands without it.
8. **D8h is a reused decision set** (A3.4).

## A3.6 Order of operations
1. Adoption preconditions A3.1.1–A3.1.10, then the tag `prereg-v0.4.3`.
2. Walk agent: steps 2–4 per A2 (unchanged), then step 4b (A3.2.2), each committed.
3. Before any per-label statistic: `gates_a3.py` runs end-to-end on the real DN2 label-free inputs (finalize and counts), plus a synthetic-label dry run on the real pool structure. Engineering bugs are fixed now.
4. Step-5 agent (not the step-4b agent; not E7's or A3's author), **in one session**:
   - a. `gates_a2.py --finalize` (v1, unchanged), commit; `gates_a3.py --finalize` (v2_sep, ce_s), commit.
   - b. `gates_a2.py --gates --g3` (v1), then **commit v1's results**.
   - c. `gates_a3.py --gates --g3` (v2_sep, ce_s), commit.
   - **After 4b, the only permitted change to `gates_a3.py` is an engineering fix that touches no statistic, parameter, pool, unit or threshold.** It is written as a dated deviation and reviewed, before it is run, by an agent who has not seen 4b's results.
5. NOTES and RESEARCH_TREE decision (A3.2.8).
6. T1 look: one `sealed_run.py` invocation, whose ledger line names the v1 primary and the v2_sep secondary arm.

## A3.7 Implementation (before the tag)
- **`config/a3.yaml`:** the freeze (A3.1.3), the model hashes, the E7b commit, every V row (pool, unit, margin, n_min from A3.1.9), seed 47 for A-syn, and the tiers.
- **Data table:** `config/a3_secondary_offsets.tsv` (dev, label-free) with its sha256.
- **Scripts:**
  - `scripts/a3/step4b_chimeric.py` and `scripts/a3/gates_a3.py`. They reuse `vulcan.gates.a2_stats`, `vulcan.gates.stats` and `vulcan.stats.g2` unchanged, and do not modify `scripts/a2/`;
  - `scripts/a3/v1_closure_guard.py` (Q1).
- **Tests:**
  - synthetic pools for every V row;
  - P1 swap tests for A-syn construction;
  - the no-per-label-count test extended to `scripts/a3/`;
  - the closure guard;
  - an end-to-end synthetic-label run.
- **New statistics** (V6, V12, V13/V14 thresholds, exhaustive chimeric G3): their docstrings are the specification. Where the text is open, the stricter reading wins (A2.9.5).

---

## A3.8 Points for the PI

The open options of the first draft are closed (Q6, A3.10). What remains for the PI is adopting or rejecting A3 as a whole, plus the disclosure in A3.1.1. Two consequences need to be accepted knowingly:
1. **The adoption preconditions are substantial:** G4, arm M, the §6.7 MC, the power simulation, A-syn, the frozen E7b, the frozen driver and G0 at the freeze commit, all before step 5. **Step 5 for v1 is delayed until they are done.** The alternative is to run v1's step 5 now and let A3 lapse, in which case v2_sep is never gated on DN2 and needs a new fresh null set.
2. **The likely outcome may be tier 1 only, or not evaluable,** because DN2 was sized for v1. The A3.1.9 simulation will say so before any data are seen.

---

## A3.9 What did not change from the first draft
- Per-scorer licensing; v1 committed first; B never used to pass.
- A2.2(e) kept for real hypotheses on A and FS.
- The reference is hyperscore through the identical pipeline.
- One shot per scorer; `itcid_linear_pred` not gated.
- Step 4b is run by a separate agent.
- One T1 invocation; v2_sep is never the primary.
- The "any violated voids" rule; no quant or speed claim.
- The A3.1 hashes, which the reviewer verified.

## A3.10 Changes after the independent review (`reviews/amendment_A3_review.md`)

| Review item | Change in revision 1 |
|---|---|
| **B1** incomplete freeze; E7b fast path; DN2 and T1 could differ | Freeze by commit + clean tree + tree hashes + predictor, Sage and mokapot hashes (A3.1.3). One implementation for DN2 and T1. E7b speed either lands before the freeze with bitwise evidence or stays off both paths. G0 re-run at the freeze commit (A3.1.4). The file hashes are now informational |
| **B2** V6 premise wrong; cannot pass at n = 2,000 | Unit corrected: every group winner, rank 1 included, classes 1 and ≥ 2 (A3.2.1 definition, V6). Two bins per class. n_min from dev-based power (≈ 23,000 per label per bin at p = 0.12), simulated in A3.1.9. Below n_min: not evaluable, which blocks the tier but is not final |
| **B3** strata ≥ 2 never tail-tested (relaxation of A1.2) | Two nested tiers: tier 1 (stratum-1 list only, the strict reading) and tier 2 (full list), which requires new clean secondary tail rows V13 and V14 on pre-specified synthetic secondary hypotheses on A spectra (A-syn: dev offset table, seed 47, A2.2(b)/(d) screens at W′, kept out of the competition) (A3.2.2, A3.2.4, A3.2.6) |
| **B4** v1 isolation not enforced | Transitive import-closure diff guard against `prereg-v0.4.2`. Each hunk must be a dated deviation or have bit-identity evidence (e.g. 8e887f7). Clean tree for step 5 (Q1, A3.7) |
| **B5** T1 inclusion steerable after DN2 | Driver and E7b freeze are adoption preconditions (A3.1.6–7), so they exist before step 5. Inclusion is unconditional (A3.3.8) |
| **B6** open options; adoption not blinded | O1 closed (no T1 look without v1 pass), O2 closed (ce_s included, licenses nothing), O3 replaced by the tiers, O4 closed (two-sided). Adoption note discloses DN2 label-free and E7b exposure (A3.1.1, Q6) |
| **M1** production rule is best-of-two; misdescribed | Stated once, with full conditions, disclosed as a fallback with error up to α_v1 + α_v2. Both verdicts quoted. "Provisional" until T1. The "depends only" sentence removed (A3.2.6, A3.2.5) |
| **M2** comparator steerable or weak | Family S is an IUT over **all** comparators run (C1, C2 at R = 2 and, if qualified, R = 5, `sage_chim{R}`, Sage, Sage+mokapot, Comet+mokapot); no R selection. Wide-window variant only via an E7b amendment before any E7b readout; otherwise the claim wording names the narrow tolerance. E7b exposure disclosed (A3.3.3, A3.1.7, A3.1.1) |
| **M3** C2 vs v1's primary | v1 comparator set = A1 set, unchanged. v1 vs chimeric comparators reported in the primary section (Q2) |
| **M4** S1 weaker than §1 | S1 = matched FDP **and** own q (supported comparators) |
| **M5** 6-run bootstrap for superiority | Leave-one-run-out minimum required in addition to the LB, for every S hypothesis (A3.3.4) |
| **M6** no power analysis | `reviews/checks/a3_power.py` is an adoption precondition. n_min for V6, V8 (`g2_sample_size`), V12 and V14. The unsupported 0.5–0.85 figure is removed (A3.1.9, A3.2.5) |
| **M7** tail threshold ambiguous | Exact τ_s per pipeline, per mixture, per merged stratum, each pipeline with its own merge map. Tail by the unit's own rank. Stage-1 winners for G1/G2/G5. Merge maps reported (A3.2.3) |
| **M8** top-5 truncation untested | New gated row V12 for v2_sep and ce_s. Listed as residual risk for v1 (A3.5.7) |
| **M9** no §6.7 MC; pooled protein; G4 partition | Adoption precondition A3.1.8: §6.7 MC of the full competition and pooled protein; G4 extended to chimeric strata; protein-group entrapment void rule for S3 (A3.3.6) |
| **M10** DN2 freshness after step 5 | DN2 valid only for scorers hashed at the tag; spent afterwards (Q4) |
| **M11** fixing gates_a3 after v1 seen | End-to-end run on real label-free inputs and synthetic labels before step 5. One session. After 4b only engineering fixes, as a dated deviation reviewed by an agent who has not seen 4b (A3.6) |
| **M12** missing preconditions | MS1 retention, review committed, G0 at the freeze, disclosures (A3.1.1–5) |
| **M13** not "only additions"; A3 name | A3.0 says plainly that it relaxes A1.1, and why that is acceptable. Naming clarified |
| **M14** D8h understated | D8h described as the E7-A1 go/no-go set and a reused decision set. Lost-run viewing to be stated (A3.4) |
| minor 1 matched-FDP description | Cites `evaluate.py`: first crossing, interpolation, grid |
| minor 2 evaluate.py modified | "default-argument behaviour unchanged, verified by A3.1.6" |
| minor 3 Sage version | Pinned: v0.14.7, 99407db6, prints 0.14.6, binary sha256 in config |
| minor 4 Sage not bit-reproducible | v1 bitwise; Sage comparators within the documented tolerance |
| minor 5 commit attributions | Corrected (d64664a / dac5ec3 / e3650c4). `config/a2.yaml` described as byte-identical to the tag |
| minor 6 row count | Counts given per tier (A3.2.5) |
| minor 7 G5 has no controlled level | Stated (V9, A3.2.5) |
| minor 8 V3 mismatch | "Not evaluable (implementation)", not fail |
| minor 9 overlap guard | Names `train_v2_folds.tsv`; described as pre-specified label-free pruning; cluster change reported |
| minor 10 A3.3.6 unit, protein, not-evaluable labelling | Peptidoform stated; protein-group overall added; not-evaluable strata labelled |
| minor 11 exhaustive chimeric G3 | Specified: same virtual-spectrum hypotheses; docstring before tag |
| minor 12 ce_s and M8 | V12 also for ce_s |
| **Not adopted** | None. Every objection was accepted. Where the reviewer offered alternatives (B3 a vs b), revision 1 takes both in nested form: strict tier 1, plus a pre-specified tier 2 |
