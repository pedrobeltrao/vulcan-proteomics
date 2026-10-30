# Amendment A2 (DRAFT): null admission, offset nulls, G1 and G3 power, fresh-null selection

- **Status:** draft for adoption as DESIGN.md §6.0b. It is not yet binding and DESIGN.md is unchanged.
- **Author:** independent statistician (agent). I have had no prior role in Vulcan, M3, M3e or A1.
- **Date:** 2026-10-27.
- **Basis:** DESIGN §6 including §6.0 (A1, `prereg-v0.4.1`); RESEARCH_TREE D9, D12 and items 3a, 3b and 3e; reviews/m3_independent_review.md; results/2026-10-27_M3e_hardening (NOTES, config, `vulcan/gates/stats.py`, `scripts/m3e/*`).
- **Checks:** `reviews/checks/a2_power.py`, with output in `a2_power.json` / `a2_power.out`. They use synthetic data only.

**Adoption preconditions (binding once adopted):**
1. The text below goes into DESIGN as §6.0b with a dated RESEARCH_TREE note.
2. The implementation (admission scripts and config) is committed.
3. The tag `prereg-v0.4.2` is created.

All three happen **before any file of the fresh null set (DN2) is listed, downloaded or searched**. The DN2 selector itself runs after the tag.

---

## A2.0 Scope and principles

D12 stands. The A1 failures of the frozen `static_gbm_pred` (tail-G2, G1, G3) are final for the existing nulls, which are never re-scored, re-filtered or re-pooled under A2. A2 governs one re-test on a fresh, never-seen null set (DN2). The scorer, the hyperscore reference, the sequence classifier and every A1 margin stay unchanged.

Every A2 rule satisfies all of the following.

- **P1, label-swap invariance.** An admission decision is a function of the spectrum, the source-species data and the *unordered* gate database T∪D (targets plus their paired decoys). It never depends on which label a candidate carries. Swapping the T and D files must leave the admitted set bit-identical. This is tested like G0.
- **P2, decided before outcome.** Every admission quantity comes from metadata, the own-species search, or label-collapsed quantities such as max over T∪D. The admitted spectrum lists are hashed and committed before any gate statistic, per-label count or tail is computed.
- **P3, no relaxation.** All A1 margins (n_tail ≥ 200, 0.60, 0.10, ±0.05 at 99%, 1 pp, 0.55) and all A1 clauses are kept. A2 only adds:
  - tighter null definitions;
  - stricter tests (G3 TOST, cluster-robust bounds);
  - minimum-n rules that make a test **not evaluable**. Not evaluable blocks FDR claims, exactly like a fail.
- **P4, one shot.** There is exactly one fresh set, drawn by the seeded rule in A2.9.
  - A **fail** on DN2 is final for this scorer under this protocol.
  - **Not evaluable** keeps FDR claims blocked. Only a new amendment (A3) with yet another fresh set can follow, and it must report A2's result alongside its own.
  - Nobody may re-draw, extend or prune DN2 after the gate run.
- **P5, report regardless of outcome** (A2.9.6).

**Why A2 is a validity correction and not a re-roll after failure.**
- It fixes two defects that are visible without any gate result:
  - (i) A1.3 filters partial true matches only against the *source* proteome plus contaminants. By construction it cannot remove real human material in a non-human sample.
  - (ii) A1 offset nulls are transforms of real human spectra (A2.4).
- Neither fix is a choice of which projects to keep, and both are applied to fresh data only.
- Nothing in A2 was tuned on the M3e nulls. The one place where M3e informed a threshold (G3 n_min, from the observed ~100 per bin) makes the test stricter, not looser.

---

## A2.1 Evaluation of the candidate admission ideas

| Idea | Verdict | Main flaw(s) |
|---|---|---|
| **M. Metadata criteria** (non-human source, culture-only, no human-tissue contact, SDRF/PRIDE fields) | **Adopt as a necessary first filter, never as sufficient** | <ul><li>PXD023811 already passed the A1 metadata filter (single organism, LFQ, keyword screen).</li><li>Free text is incomplete, and SDRF is absent for most frame projects.</li><li>Metadata cannot see handling contamination (skin and dust), LC carry-over from human runs on a shared instrument, or animal-derived media (FBS, BHI, peptone).</li><li>An organism-class rule (free-living, axenic, non-host-associated) is the only robust metadata lever. It is a priori and does not depend on any gate outcome.</li></ul> |
| **H0. Pre-screen against a human-free DB that *removes* spectra explained by the own proteome** | **Reject as stated; adopt its inverse (certification)** | <ul><li>Removing own-explained spectra *enriches* contamination: human spectra are exactly the ones the own proteome cannot explain.</li><li>It also throws away the most valid nulls, whose true identity is known to be non-human.</li><li>The inverse is sound: keep only spectra *positively certified* as own-species peptides (A2.2a).</li></ul> |
| **C. Certification** (keep only spectra confidently identified as own-species peptides) | **Adopt as the primary null definition (stratum A)** | <ul><li>(1) False certifications. At most a fraction q_cert of certified spectra are wrong, and a wrongly certified human spectrum is a strong target match. Mitigated by q_cert = 0.1% plus the symmetric screen (d).</li><li>(2) Chimeric co-isolation of a human peptide. Mitigated by (e), a single certified precursor hypothesis.</li><li>(3) The population shifts towards good spectra. That is acceptable: false positives in the 1% tail come from good spectra of absent peptides, not from noise, and the human candidates of a certified spectrum are all chance matches, which is the prior-dominated regime where a sequence leak acts.</li><li>(4) Own peptides homologous to human ones give partial true matches. Handled by the identity screen (b), which is outcome-independent.</li><li>(5) Yield: |A| is roughly 10–15% of MS2. This is an *assumption*; A2.9 sizes DN2 by a label-free count, not by this guess.</li></ul> |
| **X. "Remove spectra whose best hit in human+own is human"** | **Reject, because it is label-asymmetric** | <ul><li>The human decoys are not in the screen, so chance *target* hits are removed while chance *decoy* hits are kept.</li><li>In simulation (S4) the learned-tail target fraction falls from 0.50 to **0.04**, and a 0.25-SD target leak still reads **0.06**. That is a guaranteed pass for a leaking scorer.</li></ul> |
| **X′. Symmetric version: remove if the best hit belongs to the T∪D *family*** (a human target or its paired decoy) rather than the own proteome | **Adopt only as a second line on certified spectra** (A2.2d) | <ul><li>It is label-swap invariant (P1). In S4 it keeps a clean null at 0.49 ± 0.01 and keeps (even sharpens) a 0.25-SD leak (0.65 vs 0.58 unscreened).</li><li>Flaw 1: on *uncertified* spectra the human family nearly always beats the own proteome, so it deletes almost the whole tail (S4: n_tail 4,000 → ~110). It has no power unless it is applied after certification.</li><li>Flaw 2: it needs the T∪D search, so the script must emit only the label-collapsed boolean (P2).</li></ul> |
| **F. Protein-family variant: remove if the matched pair's protein (target or its decoy) is on a human-exposure panel** | **Adopt as an extension of the A1.3 contaminant list** | <ul><li>Pair-level, so symmetric.</li><li>It only catches listed proteins. PXD023811's tail hits (DSP, DSG1, JUP, TGM3, PKP1, DCD, DSC1, CASP14, PIP) are all panel members, but contamination from cultured human cells or carry-over would not be.</li></ul> |
| **K. Project-level marker detection** from the own-species search (keratin, skin, sweat, serum, bovine-serum panels with their own decoys) | **Adopt as a project exclusion rule** (A2.3) | <ul><li>Label-blind for the gate: the scorer under test, the gate DB outcomes and T-vs-D wins on the human proteome are never used.</li><li>Flaws: the thresholds are detection-level, not quantitative, and low-level carry-over escapes.</li><li>Hence it is defence in depth on top of C, not a substitute for it.</li></ul> |

**Not permissible as admission or exclusion criteria (they would select on the gate outcome):**
- any per-project or per-source statistic computed from the gate search, such as T/D win counts, tail fractions, the hyperscore target excess by project, G1 ρ by source, or the "source dominates the tail" rule;
- any threshold set after DN2 is listed;
- any re-weighting of projects after the gate.

The pool-validity controls in A2.5/A2.6 do use gate output, but they can only turn a verdict into *not evaluable*. They never admit, drop or re-weight data, and so they cannot help a scorer pass.

---

## A2.2 Spectrum admission: stratum A (gated) and stratum B (reported)

A DN2 spectrum enters **stratum A** iff all of (a)–(e) hold. Each step is computed once, by `scripts/a2/admit.py`, for all compared pipelines jointly. The output is a per-spectrum boolean list plus label-free counts; its sha256 is committed before the gate run.

- **(a) Certified own-species identity.**
  - Own search: production closed search (E2 SearchParams; precursor ±10 ppm; fragment ±20 ppm; max fragment charge 1; min matched peaks 4).
  - Database: own UniProt reference proteome (upid recorded at selection) ∪ cRAP ∪ Panel H ∪ Panel B (Appendix), with pseudo-reverse decoys of all entries.
  - Score: **hyperscore**, the fixed, G0-verified score. The scorer under test is never used for admission.
  - Estimator: paired double competition.
  - Certification requires the rank-1 PSM to be an *own-proteome* target with PSM-level q ≤ 0.001 and peptide-level q ≤ 0.01. Its peptide p must not occur (I/L collapsed) in cRAP ∪ Panel H ∪ Panel B and must not be A1.3-near them (6-mer or Hamming ≤ 2).
- **(b) Identity screen against the gate DB (symmetric, outcome-independent).**
  - Let W be the production precursor windows of the certified precursor (±10 ppm, with the production isotope-error set).
  - For every candidate c ∈ (T∪D) ∩ W, let f(c) be the fraction of p's theoretical b2…b(L−1) and y2…y(L−1) ions (z = 1) that coincide within ±20 ppm with c's theoretical b/y ions.
  - Remove the spectrum if max over c of f(c) ≥ 0.20, or if p is within Hamming ≤ 2 (I/L collapsed, N↔GG, Q↔AG) of any T∪D peptide.
  - Why this works:
    - It depends only on p and the DB, not on the observed spectrum or on which candidate wins.
    - It targets the mechanism of a partial true match, which is shared b/y ladders. Internal 6-mers do not create shared b/y ions unless the flanking masses agree.
    - A 6-mer test of p against all of T∪D is *not* used: with ~2·10⁷ distinct I/L-collapsed 6-mers in T∪D, it would remove nearly every peptide.
  - Implementation: score a unit-intensity synthetic spectrum of p's b/y ions with the K2 kernel in exhaustive mode, keeping only the max over T∪D.
- **(c) A1.3 winner filter, retained** (P3).
  - The contaminant list is extended by Panel H ∪ Panel B.
  - The filter applies to the matched pair of every compared pipeline (union), as in M3e.
- **(d) Symmetric family screen (X′).**
  - Remove the spectrum if, under hyperscore on the same spectrum and precursor, max over (T∪D) ∩ W of the hyperscore ≥ the hyperscore of p.
  - This removes falsely certified human spectra, whose true human peptide beats p. It is symmetric in T and D.
- **(e) Single precursor hypothesis.** In the mixture, an A spectrum carries only its certified precursor (mono m/z, z). Secondary MS1 hypotheses, which could be co-isolated human peptides, are dropped for A spectra. This applies to both labels and all pipelines.

**Stratum B** = DN2 species spectra from admitted projects that are not in A but pass A1.3 (c).
- B is the A1-type pool. B can be contaminated (pro-target), and its winner-conditioned filter leaves sub-threshold partial matches.
- **B is reported for tail-G2 and G1 and is never used to pass.**
- B enters G2-original, G5 and G3 "all admitted nulls" exactly as A1's clean pool did (A2.7, A2.8).

**Residual-contamination arithmetic (why (d) is required).**
- The worst case for (a) alone: every false certification is a human spectrum and lands in the learned tail as a target. With |A| = 120,000 and q_cert = 0.001, that is ≤ 120 spectra.
- Against n_tail ≈ 800, that is a 15% contaminated share. S2 shows this would push a non-leaking scorer towards failing (pass 0.73 at 10%; 0.03 at 20%).
- (d) removes exactly these cases, because their true human target beats p.
- Realistic false certifications are low-quality spectra, which do not reach the tail.

---

## A2.3 Project admission

A project enters DN2 only if it passes, in order, the three steps below. Steps 1–2 use metadata only. Step 3 uses the own-species search (A2.2a) and is run **before** the project is ever searched against the human DB.

1. **Frame filter (metadata).**
   - Source: `data/manifests/t4_frame.tsv`.
   - The project must have a single organism whose genus is on the ALLOW list (A2.9.2), and be DDA, Thermo, LFQ, enrichment none, on Q Exactive / Q Exactive Plus / HF / HF-X / Exploris 240/480 (the OT-HCD class; identical to A1).
   - The title, keywords and experiment types must not match BAD_KW. BAD_KW is A1's `select_nulls.py` list plus: `serum|FBS|FCS|fetal|blood|plasma|patient|clinical|skin|wound|swab|sputum|tissue|biopsy|macrophage|cell line|HeLa|HEK|K562|spike|UPS|standard|benchmark|mixed proteome|entrapment|human`.
2. **Protocol text (metadata).**
   - Fetch the PRIDE project record (`description`, `sampleProcessingProtocol`, `dataProcessingProtocol`) and the SDRF if present.
   - Exclude the project if the text matches BAD_KW or `brain heart|BHI|blood agar|peptone from (meat|casein)|Homo sapiens|human (database|proteome)|co-?culture|host`.
   - Exclude it if an SDRF lists more than one organism or an organism part other than whole organism or cell culture.
   - Missing text or SDRF is *not* a reason to exclude (most projects lack an SDRF). The rule only removes.
3. **Marker screen (own-species search; label-blind for the gate).**
   - Search all selected runs with A2.2a.
   - Count panel proteins detected at peptide q ≤ 1% with ≥ 2 distinct peptides, excluding peptides A1.3-near the own proteome.
   - **Exclude the project if ≥ 2 non-keratin Panel H or Panel B proteins are detected**, pooled over its runs.
   - Keratins alone (KRT1/2/5/9/10/14, ubiquitous in dust) do not exclude, because their spectra are already removed by (a) and (c).
   - Excluded projects are logged with their counts and replaced by the next project in the seeded walk. Replacement is allowed because the decision cannot see any gate outcome.

---

## A2.4 Offset nulls

**Are n + 0.5 Da shifted human spectra valid nulls when the fragment ions still match? No.** The A1 recipe takes a D0 MS2 spectrum and moves its precursor by n + 0.5 ± 0.1 Da (n ∈ {−2, −1, 1, 2}) while keeping its fragments. There are three structural target-favouring routes:

- **(i) Precursor mis-assignment.** If the original mono-isotopic mass, charge or isotope peak was wrong (D0 has 53–55% chimeric-candidate spectra), the shifted mass can land within ppm of the true peptide or a variant of it. M3e's diagnosis found exactly this.
- **(ii) Nested and extended targets.** Missed-cleavage, semi-specific and variant peptides that contain the true sequence share its b or y ladder. They exist only among *targets*; reversed decoys do not share it.
- **(iii) Correlated features.** RT and the predicted fragment pattern stay consistent with the true human peptide.

Each route lets a better scorer legitimately prefer targets. An offset null therefore tests "does the scorer recognise real human spectra", not "is it label-blind". They are also almost empty by design (595 of 177k had a candidate).

**A2 rule:**
1. A1-recipe offset nulls are still generated (P3: nothing is silently dropped) but are **reported only**. They enter no gate pool.
2. They are replaced by **fragment-shifted nulls (FS)**. Donors are stratum-A spectra (non-human, certified, identity-screened). Every MS2 peak m/z is shifted by **δ = +2/7 Th (0.2857)**, and the certified precursor and the RT mapping are kept.
   - Production matches z = 1 fragments at ±20 ppm. A shift of δ moves every fragment ≥ 0.26 Da from any physical z = 1 relation (isotope, NH₃, H₂O, CO, deamidation).
   - Even if fragment charges up to 4 are enabled later, z·δ ∈ {0.29, 0.57, 0.86, 1.14} Da stays ≥ 0.12 Da from the 13C and deamidation spacings (S7).
   - So FS spectra keep realistic peak statistics while their fragment evidence for *every* candidate is chance. This removes even sub-threshold homology.
   - Human spectra are never used as donors: shifting the precursor would leave (i)–(iii), and keeping it would leave the true peptide as a candidate.
3. FS gating:
   - Gate pools: G1, G5 and tail-G2, under the same margins. A fail on an *evaluable* FS pool blocks.
   - With min matched peaks = 4, most FS spectra are expected to yield no PSM. In that case FS is "not evaluable" and, **as an addition beyond A1**, does not block.
   - The prior-dominated regime that FS was meant to probe is already gated in stratum A, where every human candidate is a chance match.

---

## A2.5 Tail-G2 (A1.2 margins unchanged)

- **Gate pool:** stratum A. The same spectra, mixtures and host thresholds are used for the learned scorer and for hyperscore (each pipeline's own per-run 1% threshold, as M3e).
- **Pass requires all of:**
  1. n_tail ≥ 200; otherwise not evaluable (blocks).
  2. UCB_abs ≤ 0.60. UCB_abs = max(exact CP one-sided 95% bound, project-cluster bootstrap 95th percentile) of the tail target fraction. The bootstrap uses 10,000 resamples of projects with seed 83 (SEEDS[3]).
  3. UCB_diff ≤ 0.10. UCB_diff = max(the M3e delta-method bound, the cluster-bootstrap bound) on (learned − hyperscore).
  4. **Pool-validity control (new).** The *hyperscore* tail target fraction on A must have a two-sided 95% Wilson CI containing 0.5.
     - Hyperscore is fixed and label-blind (G0). On a valid exchangeable null it is ½ in expectation.
     - If the CI excludes 0.5 in either direction, tail-G2 is **not evaluable (invalid pool)** and blocks.
     - A pro-decoy pool would otherwise mask a leak in clause 2. A pro-target pool means residual contamination.
- **Why the cluster term:**
  - The claim is about null spectra in general, not just the 10–16 sampled projects. M3e showed source-level heterogeneity (one source held half of the tail).
  - S3: CP covers the realised-pool estimand (0.95), and the max rule is conservative (0.96–0.99). The term can only make passing harder.
- **Reported:** per project, leave-one-project-out, with and without RT (A1.3), stratum B, FS and the A1 offsets.

**What tail-G2 can and cannot detect with the unchanged A1 margins (S1, contamination-free null, p_fixed = 0.5):**

| n_tail | P(pass), true tail target fraction 0.50 | 0.55 | 0.60 | 0.65 |
|---|---|---|---|---|
| 200 (mean; half are < 200, so not evaluable) | 0.42 | 0.14 | 0.005 | 0 |
| 400 | 0.99 | 0.47 | 0.013 | 0 |
| 750 | 1.00 | 0.80 | 0.015 | 0 |
| 1500 | 1.00 | 0.98 | 0.011 | 0 |

- A +0.10 tail bias fails at ≥ 98.5% from n_tail ≈ 400.
- **A +0.05 tail bias passes 80–98% of the time: tail-G2 under A1 margins is not designed to detect it, and A2 does not change that** (P3 forbids loosening and this draft does not tighten).
- A +0.05 bias at PSM level is caught by G2-original (margin 0.01, ≥ 120k competitions). In the tail it is bounded only by the entrapment arms.
- A tail target fraction of 0.60 corresponds to decoys under-counting tail nulls by 1.5×. The team should not describe a tail-G2 pass as "no tail leak". The correct wording is "tail null target excess ≤ 0.10 (95%)".

---

## A2.6 G1

**Is ±0.05 meaningful under contamination? No.** The gate statistic is Spearman(score, sequence-classifier logit) over the top-decile null winners, which mixes target and decoy winners.

- Any process that makes target winners score higher creates ρ > 0 through the label path. The classifier separates targets from decoys at AUC ≈ 0.70 even for a scorer that ignores sequence entirely, so truly-human spectra produce this effect.
- S5: a non-leaking scorer reaches ρ ≈ 0.015 / 0.028 / **0.054** when 5 / 10 / 20% of top-decile winners are real human spectra. M3e's 0.048, and hyperscore's own 0.028 [−0.001, 0.056], fit this.
- Conversely, G1 is a weak detector of a pure label-preference leak: a 0.3-SD target bonus gives ρ ≈ 0.049. G1's real job is to detect *use of sequence information*, and that is only interpretable on a valid pool.

**Valid fresh-null design:**
1. **Pool:** stratum A. FS is also gated if evaluable (A2.4). B is reported.
2. **Minimum n.** n_top ≥ **7,200** top-decile winners, so |A winners| ≥ 72,000.
   - At ρ = 0 the 99%-CI-inside-±0.05 rule then passes with probability ≥ 0.90 (S5: n_min = 7,150). At ρ = 0.05 it still passes only 0.5% of the time.
   - Below 7,200 the gate is not evaluable (blocks).
   - The rule itself (two-sided 99% Fisher-z CI with Bonett–Wright variance, inside (−0.05, 0.05)) is unchanged.
3. **Negative control (new).** Hyperscore's ρ on the same pool must have a 99% CI containing 0. Otherwise G1 is **not evaluable (invalid pool)** and blocks.
4. **Frozen classifier.** The classifier, `scripts/m3e/seqclf.py` with its training data and features, is frozen:
   - its pickle sha256 is recorded in the A2 config before DN2 is listed;
   - if the /tmp copy is lost, it is regenerated deterministically by the committed script and hashed *before* DN2 is listed;
   - it is never retrained afterwards.
5. **Reported, not gated:** ρ within target winners and within decoy winners (Fisher-z pooled), which removes the label path; and ρ on B and on the A1 offsets.

---

## A2.7 G2-original, G4, G5

- **G2-original and G5:** unchanged rules and margins, on **A ∪ B** (the A1-equivalent "clean" pool from admitted projects). G5 is also run on A alone and on FS if evaluable, and an alarm on any of these blocks. The A1 offsets are reported only.
- **G4:** unchanged (training decoys; not a null-pool gate).

---

## A2.8 G3 power rule (TOST per bin, minimum n)

The A1 rule was a point estimate |recall_T − recall_D| ≤ 0.01 in every decile bin. At ~100 spectra per bin, one miss (1 pp) decides the outcome.
- A symmetric pruner with recall 0.999 passes all 10 bins only **14%** of the time (S6).
- A pruner that loses 1 pp more decoys in the top bin passes 7% of the time. In other words, the rule decides on noise.

**A2 rule (stricter for any given data, P3):**
- **Pools:** (i) stratum A; (ii) all admitted nulls, A ∪ B. Both must pass, as A1 required for both of its pools.
- **Data:** G3 runs exhaustive and production top-K search on **all** DN2 mixtures, not just the first two.
- **Bins:** deciles of the exhaustive best score within the pool, as in A1.
  - Each of the **top five** deciles must contain ≥ **2,000** spectra with an exhaustive best target and ≥ 2,000 with an exhaustive best decoy. Otherwise G3 is not evaluable (blocks).
  - A bottom-half decile below 2,000 is merged with the next higher bottom-half decile until it reaches 2,000. Merging never crosses into the top half.
- **Test per bin:** Δ = recall_T − recall_D, and a 90% two-sided Newcombe hybrid-score interval (the two one-sided 95% bounds of TOST).
  - **Pass iff the interval ⊂ [−0.01, +0.01] in every bin.**
  - Both directions are kept. A1's |Δ| ≤ 1 pp was two-sided, and dropping the conservative direction would be a relaxation.
  - The interval treats T and D as independent samples. They are paired per spectrum with positive correlation, so the interval is conservative.
- **Why this does not loosen the 1% margin in expectation:**
  - TOST passes ⇒ |Δ̂| ≤ 0.01 − half-width, so on identical data it never passes when the A1 rule fails.
  - At a true 1 pp decoy deficit, its pass probability is ≤ 5.4% at every n tested; A1's point rule passed ~45%.
  - At a 2 pp deficit it passes ~0% (A1's point rule passed 0.2–4%).
  - Pooling across bins was rejected. It would dilute a top-bin asymmetry, which is the one that matters for the 1% tail, by nine uninformative bins.
- **Power (S6, all 10 bins must pass):**

  | n per bin per label | symmetric, recall 0.999 | 0.995 | 0.99 |
  |---|---|---|---|
  | 1,000 | 0.97 | 0.07 | 0.00 |
  | 2,000 | 1.00 | 0.89 | 0.20 |
  | 5,000 | 1.00 | 1.00 | 0.99 |

  The 2,000 minimum gives ≥ 0.89 power at the recalls E2 measured (0.998–0.999). With |A ∪ B| ≈ 10⁵ spectra that have a candidate, bins hold ~10⁴, so the binding case is the lower recall.

---

## A2.9 Fresh-null selection (DN2), procedure, and reporting

### A2.9.1 Exclusions (mechanical)
- The union of the `accession` columns of **every** `data/manifests/*.tsv` except `t4_frame.tsv`, as of the A2 tag. Today that is 300 accessions, covering:
  - dev and sealed sets, tier S, D7b, D9–D12, D10cf;
  - the T4 selection, train.tsv and the predictor-training accessions;
  - the E9/D8 stream and DN (g2null).
- The ProteoBench (non-PXD) sets.
- All six previously used null species: *P. furiosus*, *S. cerevisiae*, *T. brucei*, *S. aureus*, *Synechococcus* sp., *B. cereus*.
- Training species (human, mouse, *E. coli*) and *Arabidopsis* (arm C).

DN2 accessions are added to the T4 selector's exclusion list and to `bench/sealed_run.py`'s T4 assertion.

### A2.9.2 ALLOW list (organism class fixed a priori)
The class is free-living, non-host-associated, axenically culturable unicellular organisms in animal-free media. Human commensals and pathogens, parasites, anything cultured on host cells or serum, plants and animals are all excluded. The list is matched on genus:

*Chlamydomonas, Synechocystis, Sordaria, Phaeodactylum, Streptomyces, Myxococcus, Cupriavidus, Caballeronia, Frankia, Abortiporus, Bacillus* (not *cereus*/*anthracis*/*thuringiensis*), *Komagataella, Fervidobacterium, Corynebacterium glutamicum, Emiliania, Propionibacterium freudenreichii, Deinococcus, Chlorella, Aeropyrum, Acaryochloris, Aspergillus niger/nidulans, Amycolatopsis, Sorangium, Neurospora, Pseudomonas putida, Rhodobacter, Haloferax, Halobacterium, Shewanella, Geobacter, Trichoderma, Schizosaccharomyces, Caulobacter, Sulfolobus/Saccharolobus, Thermococcus, Methanosarcina, Thermus, Zymomonas, Sinorhizobium, Thalassiosira*.

A metadata-only count on the frame, run while drafting (no PRIDE file listing), gives **82 eligible projects in ~35 genera**.

### A2.9.3 Seeded walk
- **Order:** sort the candidates by accession and shuffle with `random.Random(23)` (SEEDS[1]; DN used 11).
- **Per project:**
  - At most **one project per genus**.
  - It needs ≥ 3 usable RAW files (as A1: `.raw`, 0.3–1.5 GB, file name not matching BAD_FILE).
  - Pick 3 files with `random.Random(23 + accession digits)`.
  - The protocol check (A2.3 step 2) and the marker screen (A2.3 step 3) must both pass.
- **Mixtures:** as A1 (D0 hosts, 20,000 per mixture, each spectrum used once, RT quantile-mapped), permutation seed 23. FS is generated with δ = 2/7 Th.
- **Stopping rule.** Uses only label-free counts from `admit.py`. Stop when **all** of the following hold, or at a cap of **16 admitted projects**, whichever comes first:
  - ≥ 10 projects admitted;
  - |A| ≥ 120,000;
  - A spectra with an exhaustive best candidate ≥ 72,000 (G1 n_top);
  - the hyperscore tail proxy ≥ 400. The proxy counts A spectra whose max over T∪D hyperscore ≥ the host's hyperscore 1% threshold, so it is label-collapsed.

  If the cap is reached first, the gates run anyway and may come out not evaluable (P4).

### A2.9.4 Order of operations (each step is committed before the next)
1. Tag `prereg-v0.4.2`.
2. Run the selector and download, recording sha256 and the selection log.
3. Run the own-species searches, the marker screen and `admit.py`. Commit the per-spectrum A/B/FS membership hashes, the label-free counts and the swap-test result (P1).
4. Build the mixtures.
5. Run the gates once with the frozen scorer, hyperscore and the frozen sequence classifier.

No code under `scripts/a2/` may print or store per-label counts before step 5.

### A2.9.5 Who runs it
- Steps 2–4 should be run by an agent that does not run step 5, or at least in a separate commit with the gate scripts untouched.
- Any deviation from this text is written as a dated note before the step it affects. Where something is ambiguous, the stricter reading wins.

### A2.9.6 Reporting commitment
- `results/<date>_A2_gates/NOTES.md` reports every gate (pass, fail or not evaluable) on A, B, FS and the A1 offsets. It includes the per-project tables and the pool-validity controls, and is written whatever the outcome.
- A RESEARCH_TREE decision records the verdict.
- A **pass** satisfies tail-G2, G1 and G3 for `static_gbm_pred` and nothing else (G4, G5, arm M, B1/T1 preconditions still apply). It is quoted as "passed on fresh nulls under A2, after failing under A1 on contaminated nulls (D12)". The A1 failure is never omitted.
- A **fail** is final under P4.

---

## A2.10 Simulation and analytic checks (summary of `reviews/checks/a2_power.*`)

All checks are synthetic and use seed 131.

| Rule | Check | Result |
|---|---|---|
| Certification + symmetric screens keep leak power (A2.2) | S4: learned-tail target fraction under no screen / asymmetric X / symmetric X′ | No leak: 0.50 / **0.04** / 0.49. 0.25-SD target leak: 0.58 / **0.06** / 0.65. The asymmetric screen hides a leak; the symmetric one keeps it. The cost of X′ on uncertified spectra (n_tail 4,000 → 110) is why it is only used after certification |
| Contamination drives tail-G2 failure regardless of leakage (motivates A) | S2: non-leaking scorer, contaminated share c of the learned tail | P(pass) 0.93 / 0.92 / 0.87 / 0.73 / 0.03 / 0.00 at c = 0 / 1 / 5 / 10 / 20 / 40%. M3e's *S. aureus* share was ~50% |
| Tail-G2 power with A1 margins (A2.5) | S1 | Fires (≥ 98.5%) at ≥ +0.10 from n_tail 400. Does **not** fire at +0.05 (by design of the A1 margins) |
| Cluster term (A2.5) | S3: coverage of the one-sided 95% UCB, 10 projects | CP 0.947–0.958; max(CP, cluster) 0.958–0.994. Stricter, never looser |
| G1 n_min, contamination (A2.6) | S5 | n_min 7,150 for 90% pass at ρ = 0; boundary ρ = 0.05 passes 0.5%. Contamination induces ρ ≈ 0.054 at 20% |
| G3 TOST + n_min (A2.8) | S6 | See the A2.8 table. At a 1 pp decoy deficit: TOST ≤ 5.4% pass; A1 point rule ~45% |
| FS δ (A2.4) | S7 | δ = 2/7 Th: z·δ ≥ 0.12 Da from the 13C/deamidation spacings for z ≤ 4. Production uses z = 1 (gap ≥ 0.26 Da) |

---

## A2.11 Residual risks (stated so nobody reads A2 as stronger than it is)

1. **A narrower gated population.** Tail-G2 and G1 are now gated only on certified, identity-screened spectra.
   - A leak that acts *only* on low-quality, uncertifiable spectra would be tested in the tail only through B (reported) and FS (likely not evaluable).
   - It would still be tested at PSM level through G2-original and G5 on A ∪ B.
   - This is the price of a valid null, and the main point a skeptic should press.
2. **The +0.05 tail blind spot** remains under the A1 margins (A2.5). Tail FDR control rests on the entrapment arms.
3. **Yield uncertainty.** If |A| or the tail proxy falls short at 16 projects, the result is not evaluable and the scorer stays blocked. The stopping rule protects validity, not the schedule.
4. **Pool-validity controls can block a clean scorer** if residual sub-threshold homology makes hyperscore's tail fraction differ from ½. That is a deliberate trade: block rather than pass on an invalid pool.
5. **Degrees of freedom in the implementation** (panel resolution, f(c) details, merge order). These are fixed in code docstrings and config before DN2 is listed. The ambiguity rule (stricter wins) applies.
6. **Repeat attempts.** P4 allows exactly one re-test. Any A3 must count A2 as a prior attempt and report both.

---

## Appendix: panels

The panels are resolved to UniProt reviewed accessions at tag time. The sha256 of the panel FASTA is recorded in the A2 config.

- **Panel H (human exposure)**, non-keratin:
  - desmosome and cornified envelope: DSP, DSG1, DSC1, DSC3, JUP, PKP1, PKP3, EVPL, PPL, TGM1, TGM3, CASP14, FLG, FLG2, HRNR, LOR, IVL, SPRR1A, SPRR1B, SPRR2A, KRTDAP, SBSN, DMKN, CDSN, S100A7, S100A8, S100A9, ARG1, SERPINB3, SERPINB4, SERPINB12, CSTA, CST6, CALML5, ALOX12B, KLK5, KLK7, LYPD3;
  - sweat, saliva and tears: DCD, PIP, SCGB1D2, SCGB2A2, LYZ, LTF, AZGP1, ZG16B, MUCL1, LCN1, STATH, HTN1, AMY1A, CST1, CST4, SMR3B;
  - blood and serum: ALB, HBA1, HBB, HBD, APOA1, APOA2, SERPINA1, TF, HP, A2M, FGA, FGB, FGG, IGHG1, IGKC, CP, C3.
- **Panel H (human exposure)**, keratins (do not trigger A2.3 step 3): KRT1, KRT2, KRT5, KRT6A, KRT9, KRT10, KRT14, KRT16, KRT17, KRT77.
- **Panel B (animal-derived media):**
  - bovine: ALB, TF, APOA1, SERPINA1, AHSG, FETUB, HBA, HBB, LALBA, CSN1S1, CSN2, CSN3, LGB;
  - porcine: ALB, TF, HBA, HBB;
  - chicken: OVAL.
