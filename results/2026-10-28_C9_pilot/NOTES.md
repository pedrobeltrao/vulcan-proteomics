# C9 pilot: human "removed-peptide" nulls, label-free feasibility (D28.2)

- **Author:** fresh T1-blind agent, 2026-10-28. No role in A1–A4, DN2, DNIT, E7/E7b or the T1 driver.
- **Read:** RESEARCH_TREE D9–D28, `reviews/pi_proposal_D28_review.md` (C9, section E), DESIGN §6.0–6.0c, A2 draft (A2.0–A2.8, esp. A2.4), A4 draft rev 3 (A4.0, A4.3.6–8, A4.5–A4.7), `scripts/a2/*`, `vulcan/gates/*`, `reviews/checks/a4_power.py` and `a4_r1/r3_walk.py`, LESSONS.
- **Facts checked:**
  - `results/SEALED_LEDGER.md` has a header only, with no entries.
  - No sealed file, T-set, escrow or `bench/sealed_run.py` run was touched.
  - Only DEV data was used: D0 (PXD001468, 2 runs) and D8h, the first 10 D8/E9 runs (10 projects).
  - The frozen models and `scripts/a2`, `scripts/a3` are unmodified.
- **Label discipline:**
  - No target/decoy label of any null candidate was read. The exhaustive null table is never joined to `labels.arrow`.
  - Labels of the full index were read only for:
    - certification of the human spectra, which is an ordinary identification search;
    - each run's hyperscore 1% threshold (host threshold, as `scripts/a2/host_threshold.py`);
    - the identified set I.
  - Every null quantity below is label-collapsed: counts, max over T∪D, and the ladder overlap with the certified peptide.

**Bottom line: GO under the D28 rule.** Simulated P(v1 minimum) is ≥ 0.95 within 200 runs in every cell, including halved yields and the strictest screens. It reaches 1.0 at 300 runs. The cost is about 0.1–0.2 TB and hours of compute.

GO has conditions. The design must use a **stricter ladder screen (f < 0.05)** than A2's 0.20 (section 3c), and must **forbid removal driven by hits on the null spectrum itself** (section 1.3).

The main remaining uncertainty is not the minimum size but the **validity control**. A residual target excess of ≥ 0.02 in the hyperscore tail makes "not evaluable (invalid pool)" likely. The pilot cannot measure that excess without labels. Every route to it pushes towards invalid or fail, never towards a false pass (section 2).

---

## 1. Design (for the amendment; precise enough to implement, not normative)

### 1.1 Source and certification
- **Source runs.** A seeded walk over the eligible human pool. The frame is the E9/D8 one: *Homo sapiens* only, DDA, LFQ, no enrichment, Thermo QE / QE Plus / HF / HF-X / Exploris 240/480, with the BAD_KW screen.
  - **Excluded:**
    - every accession in `data/manifests/*.tsv` (dev, sealed T1–T5, the T4 selection, tier S, train, train_m7, train_itcid, predictor training, DN2, DNIT, D8/E9);
    - every project E9 tried (166);
    - every project in any T4 draw under `results/**/t4_selector/` (T4-reachable).
  - **Pool** (`pool_count.json`, metadata only): 1,666 frame projects, of which 201 are excluded, leaving **1,465 projects from 969 labs**. 753 labs have a single project; the largest has 25.
  - **Per project:** R = 3 runs (2 simulated as well), drawn from usable RAWs of 0.1–2 GB. Cluster = normalised lab, ≤ 5 projects per cluster.
  - **Role:** every walked project becomes `dev-null` and is excluded from future T4 redraws.
- **Certification score.** Hyperscore (fixed, G0-verified, not under test), from the production closed search of the run against the **full** human + contaminant index.
  - Estimator: paired double competition per run.
  - **Certified** means all of:
    - the rank-1 PSM is a human target;
    - PSM-level q ≤ 0.001 and peptidoform q ≤ 0.01;
    - the peptide is not A1.3-near the contaminant FASTA.
  - This is A2.2(a) with "own proteome" = human (`vulcan.gates.a2_admission.certify_run`).
  - Per-run q ≤ 1e-4 is **not resolvable** by TDC (the +1 floor): only 8,501 of 64,507 certified spectra survive, and 11 of 12 runs give zero (`certify_strict.json`). A stricter certification needs project-pooled TDC or a two-engine agreement (section 2, risk R4).
- **Single precursor (A2.2(e)).** A certified spectrum carries exactly one hypothesis: certified mono m/z, z, isotope error 0.
  - The MS1 co-isolation flag (≤ 1 feature in the isolation window) was measured as a variant. It is not recommended as a filter (section 3d).

### 1.2 Removal set (label-free, closed under pairing, the same function for both labels)
For a certified spectrum s with peptide p in project P, let I be the I/L-collapsed stripped target peptides accepted at peptidoform q ≤ 1% by the same fixed-score search in any run of P. I is built from **rank-1 PSMs only**.

A T∪D peptidoform c in the precursor window of s is removed iff **c or its pair partner** meets either condition:
- **(i) Near I.** It is I/L-identical to, within Hamming ≤ 2 of any window of, or shares a 6-mer with, any peptide in I or an isobaric variant of it (N↔GG, Q↔AG).
  - This uses `HomologFilter` over I, which covers all modforms, nested and extended peptides, and missed-cleavage and semi-specific relatives.
- **(ii) Shares ions with p.** It shares ≥ f_max of p's b2…b(L−1) and y2…y(L−1) ions with its own b/y ladder.
  - Ions are z = 1 at ±20 ppm, cross-type matches are allowed, and all mod states count (A2.2(b) f(c), exact numpy).
  - **Recommended f_max = 0.05** (section 3c). A2 used 0.20.

The paired decoys of removed targets go too, and so do the paired targets of removed decoys. The null competition of s is the exhaustive or production search against T∪D \ R. In production this is a masked index per project.

Recommended extra: drop contaminant-FASTA pairs from the gate DB. This is label-free and removes non-human co-eluters.

### 1.3 What the removal must not use (P1)
- **No secondary or chimeric hits of the null spectrum itself, and no rank ≥ 2 PSMs of it, may drive removal.**
  - Such hits are selected because a candidate scored high on s. Removing their neighbourhoods is a label-asymmetric, X-type removal (A2.1).
  - It removes target pairs that scored well by chance on s, which is a **pro-decoy, false-pass** direction.
  - Co-isolated peptides are handled only through I, which is built from the rank-1 IDs of other spectra, and optionally a label-free MS1 screen.
- **The removal decision** is a function of (p, I, the unordered T∪D with its pair map) only.
  - G0-style swap test: rebuild the index with `swap_labels=True`. Certification then uses the swapped full index, so it must be **frozen** (the certified list and I hashed) before the gate DB is touched.
  - With certification frozen, the removal and admission membership must be bit-identical.

### 1.4 Gate search and statistics (A2/A4 reused; nothing relaxed)
- **Mixture = the null's own run.**
  - Certified spectra are searched against their reduced DB with a single hypothesis.
  - All other spectra of the run keep the full DB and production hypotheses.
  - Each pipeline's own per-run 1% threshold is computed on that combined PSM list, as in the A2 mixtures.
  - RT is the real human RT, with no quantile mapping.
  - The pilot's tail proxy used each run's own full-search hyperscore threshold, with the D0-host max (A2 convention) as a sensitivity.
- **Rows** (v1 = completed test 2 at A4.6 levels α/2; v2_sep = test 1):
  - tail-G2 (n_tail ≥ 200, UCB_abs ≤ 0.60 and UCB_diff ≤ 0.10, one-sided 0.975), with the Wilson validity control on hyperscore;
  - G1 (n_top ≥ 11,100, 0.995, hyperscore control);
  - G3 (TOST ±0.01, ≥ 2,000 per label per bin);
  - G2-original and G5 as A2.7;
  - v2_sep rows V1–V14 per A3 and A4 (n_min V6 30,600, V12 26,800).
- **Inference** follows A4.5: the outer envelope of the analytic, project-bootstrap and lab-cluster-bootstrap intervals, the 10% cluster cap (minimum of the proxy, A_cand and |A| factors), and the intersection–union rule over the capped and uncapped pools.
- **Stopping rule** (A4.3.7, label-free, evaluated after each admitted project):
  - v1 minimum: J ≥ 20, J_eff ≥ 10, capped proxy ≥ 1,200, capped A_cand ≥ 120,000;
  - full stop: the same with A_cand ≥ 450,000.
  - The amendment may tighten J (C9 suggested ≥ 30 labs), but must re-simulate.
- **Reporting:** with and without RT, per lab, leave-one-cluster-out, and the A1/A2/A4 history (C9: "must report A1, A2 and A4-not-run").

---

## 2. Against A2.4's rejection of human donors

**What A2.4 said.** A2.4 argued that n + 0.5 Da offset nulls built from human spectra are invalid through three target-favouring routes:
- (i) **precursor mis-assignment** puts the shifted mass on the true peptide or a variant;
- (ii) **nested and extended targets** share the true ladder, and only targets do;
- (iii) **correlated features**: RT and the predicted pattern stay consistent with the true peptide.

So an offset null tests whether the scorer recognises real human spectra, not whether it is label-blind. For FS nulls, A2.4 added: "Human spectra are never used as donors: shifting the precursor would leave (i)–(iii), and keeping it would leave the true peptide as a candidate."

**Why removal answers each route:**
- **(i)** The precursor is not shifted. The window is centred on the certified mono mass, and the true peptide, all its modforms and their decoys are gone, so no mis-assignment can put it back.
- **(ii)** Every nested, extended or missed-cleavage relative shares a 6-mer with p (all peptides are ≥ 7 aa) and is removed. Isobaric relatives that share ≥ f_max of the ladder are removed, both labels.
- **(iii)** The competitors are now unrelated human peptides and their pseudo-reverse decoys. A pseudo-reverse decoy has the same composition as its target, so composition-driven RT and pattern agreement is symmetric within each pair. A *sequence-learned* preference for real human peptides is exactly the leak the gate is for. Human nulls test it in the deployment regime (human-trained predictors, human candidates), which the D28 reviewer called the right target.
- **"Keeping the precursor leaves the true peptide as a candidate"** is the objection removal answers by construction.

**Residual risks and their direction:**

| # | Residual | Who gains | Direction | Mitigation |
|---|---|---|---|---|
| R1 | Unremoved homologs (Hamming ≥ 3, no shared 6-mer, f < f_max) carrying partial true evidence | targets only (biological homology) | hyperscore tail > ½ → **invalid pool**; learned excess → **fail** | f_max = 0.05; I built per project |
| R2 | Unidentified co-isolated peptides (chimeras) | targets that share sub-ladders with the co-isolate (rarely the co-isolate itself, since the ±10 ppm single hypothesis excludes it) | invalid / fail | I removes identified co-isolates' neighbourhoods; MS1 screen optional (no label-free signal, section 3d) |
| R3 | Contaminants (keratins are human and partly identified; trypsin, BSA in the contaminant FASTA) | targets | invalid / fail | drop contaminant pairs from the gate DB; certification excludes contaminant-near p |
| R4 | False certification: the true peptide stays in the DB, unrelated to the false p | targets (strong) | invalid / fail | Bound: ≤ 0.1% of \|A\|, i.e. ≤ ~330 at the f < 0.05 v1 minimum (|A| ≈ 120k / 0.365) against n_tail ≈ 600–1,250 (worst case if all land in the tail). A2's family screen (d) is vacuous here (the reduced DB ⊂ the full DB). **Needed:** project-pooled TDC at q ≤ 1e-4, or two-engine agreement, or a label-free score-gap screen hs(p) − max_reduced ≥ Δ (symmetric, but it censors the tail). Its yield cost must be measured before the tag |
| R5 | Population shift towards good spectra | none | changes the estimand (as A2) | report |
| R6 | Removal driven by hits on s itself, or any label-dependent pruning | **decoys** | **false pass** | forbidden by design (section 1.3); swap test |

- Every residual except R6 favours targets for every scorer.
  - The hyperscore validity control turns this into "not evaluable (invalid pool)".
  - A learned scorer that exploits partial evidence more strongly turns it into "fail".
- Neither can produce a false pass of the absolute clause (UCB ≤ 0.60). Contamination only raises the learned tail target fraction.
- Contamination does inflate the hyperscore reference, so it loosens the difference clause by about δ_hs. But the absolute clause still binds, and δ_hs ≳ 0.02 is flagged invalid at these n (section 4).
- **The only false-pass route is design-side (R6) and is excluded by construction.**
- Cost of the fail-safe: an invalid pool counts as a completed test (step 5 ran, A4.6), so v1's next attempt would run at α/4. That is a reason to adopt the strict screens before the run, not after.

---

## 3. Label-free yield pilot (12 DEV runs, 11 projects; `per_run/*.json`, `metrics.json`)

### 3a. Per run
"Proxy" is the label-collapsed count of A_cand spectra whose max surviving hyperscore is ≥ the host 1% threshold.

| project | run | MS2 | certified | A_cand (f < 0.20) | proxy own-thr | A_cand (f < 0.05) | proxy (f < 0.05) |
|---|---|---|---|---|---|---|---|
| PXD011447 | E9_0000 | 134,146 | 11,126 | 4,878 | 40 | 2,200 | 3 |
| PXD018572 | E9_0001 | 20,275 | 5,594 | 3,208 | 122 | 1,469 | 21 |
| PXD014823 | E9_0002 | 20,519 | 1,440 | 760 | 16 | 271 | 2 |
| PXD016022 | E9_0003 | 50,257 | 5,362 | 2,793 | 98 | 1,324 | 21 |
| PXD027980 | E9_0004 | 59,625 | 2,349 | 2,161 | 1 | 1,690 | 0 |
| PXD057360 | E9_0005 | 66,278 | 3,525 | 1,970 | 40 | 846 | 6 |
| PXD069854 | E9_0006 | 25,220 | 3,899 | 3,553 | 56 | 2,609 | 17 |
| PXD011751 | E9_0007 | 10,411 | **0** | 0 | 0 | 0 | 0 |
| PXD026819 | E9_0008 | 65,101 | 3,848 | 3,029 | 8 | 1,675 | 3 |
| PXD007543 | E9_0009 | 44,644 | 6,874 | 4,458 | 62 | 2,156 | 15 |
| PXD001468 | D0 01A | 41,820 | 10,644 | 8,080 | 166 | 4,900 | 23 |
| PXD001468 | D0 01B | 46,816 | 9,846 | 7,396 | 199 | 4,387 | 28 |
| **total** | 12 | 585,112 | **64,507** | **42,286** | **808** | **23,527** | **139** |

### 3b. Rates

| Quantity | Value |
|---|---|
| Certified per MS2 | 11.0% |
| Certified per run | 5,376 mean (0–11,126; one zero run) |
| Certified spectra with ≥ 1 surviving candidate (≥ 4 matched peaks) | 65.6% (f < 0.20); 36.5% (f < 0.05) |
| Tail proxy per A_cand, own threshold | 1.91% (f < 0.20); 0.59% (f < 0.05) |
| Tail proxy per A_cand, D0-host max threshold | 1.45% (f < 0.20) |
| DN2 reference (foreign certified spectra, D0 host) | 0.23% |

- The tail proxy per run is ≈ 67 at f < 0.20 and ≈ 12 at f < 0.05, against DN2's ≈ 1.2.
- Yields vary by about 10× between projects: log-sd of certified per run 0.63 between projects; tail rate log-sd 1.25.

### 3c. Ladder diagnostic (label-free, the reason for f < 0.05)
The table bins every A_cand spectrum by the ladder overlap f between its **winning surviving candidate** and its own certified peptide p.

| f of winner vs p | A_cand winners | tail winners | tail rate |
|---|---|---|---|
| [0, 0.05) | 9,378 | 110 | 1.17% |
| [0.05, 0.10) | 9,691 | 146 | 1.51% |
| [0.10, 0.15) | 11,314 | 190 | 1.68% |
| [0.15, 0.20) | 11,903 | 362 | **3.04%** |

- The tail is enriched 2.6× in winners that share 15–20% of the true peptide's ions. Nearly half of the tail (45%) sits in that bin.
- This is partial true-match evidence that A2's f < 0.20 screen (designed for *foreign* spectra, where overlap with human peptides is rare) leaves in.
- Whether those winners are targets or decoys is unknown and must stay unknown.
- But biological homology exists only on the target side, so this is the R1 risk in concrete form.
- At **f < 0.05** (≤ 1 shared ion), this route is closed by construction. The cost is A_cand × 0.56 and tail proxy × 0.17.
- It is a design choice made on label-free counts before any label exists (P2). It applies to both labels identically, and a stricter symmetric removal cannot create a pro-decoy pool.
- Under the f < 0.05 design the tail rate is 0.59% per A_cand, still about 2.6× DN2's 0.23%. That is plausibly spectrum quality (high-resolution human runs with rich fragment ladders). It is not identifiable label-free, and it is the reason the validity-control risk is quoted below.

### 3d. Co-isolation diagnostic (label-free)
- Tail rate per A_cand is 1.914% in spectra with a single MS1 feature in the isolation window (14,942 A_cand) and 1.909% with several (27,344).
- So co-isolation leaves no label-free trace in the tail, and the single-feature filter (yield × 0.38) is not recommended as a requirement.
- E9_0001's MS1 table reports one feature for every spectrum, which looks like an instrument or method artefact. Excluding it, the rates are 1.4% vs 1.9%.

### 3e. Cost per run
Measured in one warm process on the shared node (EPYC 9354 + A100, load 7–10; `per_run/_meta.json`).

| Item | Cost |
|---|---|
| Full search, GPU | 0.06–0.6 s wall |
| Certification (Python) | ≈ 4 s |
| Exhaustive null search | ≈ 0.1 s |
| Homolog removal | 1–33 s (`HomologFilter` Hamming loops; a masked index per project would remove this) |
| Ladder removal | ≈ 0.4 s |
| **Mean per run** | **11.7 s wall, 11.5 process-CPU s**; index setup 28–47 s once |

Not measured here: decode (≈ 5 core-s per run, D27), the gate scoring of v1/v2_sep, and download (≈ 1.08 GB per run, E9 mean).

---

## 4. Simulation (`c9_power.py` → `c9_power.json`; seed 131; A4 machinery reused)
**Model:**
- The walk runs over the **real** 1,465-project pool in seeded order, with real lab keys and ≤ 5 projects per lab.
- Usable RAWs per project come from E9's label-free listing counts (166 projects: 24 zero, median 15).
- Run decode failure is 0.2%; dead projects are 1.4% (E9).
- Per-run yields are log-normal or logit-normal with the pilot's between-project spreads, plus the pilot's 1/12 zero-yield runs.
- A4 10% cap (min of three factors), J ≥ 20, J_eff ≥ 10, capped proxy ≥ 1,200, capped A_cand ≥ 120,000 (full stop 450,000).
- 400 walks per cell. Tail-G2: 200 states per scenario, A4 levels, IUT over capped and uncapped pools (`vulcan.gates.a2_stats.tail_g2_a2`).
- G1 uses `a4_power.g1_pass_prob` with the jackknife hull at the walk's J. G3 uses `a4_power.tost_pass_prob`.

### 4a. Evaluability: P(v1 minimum) by run budget, and size at the v1 minimum

| screen | R/proj | yields | P at 60 / 100 / 150 / 200 / 300 runs | runs p10/50/90 | projects p50 | J p50 | runs at full stop p50 |
|---|---|---|---|---|---|---|---|
| f < 0.20 | 3 | pilot | 0.96 / 1 / 1 / 1 / 1 | 51 / 54 / 58 | 20 | 20 | 124 |
| f < 0.20 | 3 | ×0.5 | 0.03 / 1 / 1 / 1 / 1 | 65 / 74 / 85 | 28 | 27 | 233 |
| f < 0.20, D0-host thr | 3 | pilot | 0.68 / 1 / 1 / 1 / 1 | 52 / 58 / 68 | 21 | 21 | 123 |
| **f < 0.05 (recommended)** | 3 | pilot | 0.01 / 0.67 / **1.00** / 1 / 1 | 72 / **92** / 116 | 34 | 34 | 209 |
| **f < 0.05** | 3 | ×0.5 | 0 / 0.01 / 0.56 / **0.97** / 1 | 118 / 146 / 179 | 54 | 53 | 406 |
| f < 0.05 | 2 | pilot | 0.06 / 0.88 / 1 / 1 / 1 | 64 / 80 / 102 | 44 | 43 | 205 |
| f < 0.05 | 2 | ×0.5 | 0 / 0.02 / 0.70 / 0.99 / 1 | 111 / 135 / 168 | 74 | 70 | 409 |
| f < 0.20 + single MS1 | 3 | pilot | 0 / 0.72 / 1 / 1 / 1 | 80 / 94 / 108 | 35 | 34 | 312 |
| f < 0.20 + single MS1 | 3 | ×0.5 | 0 / 0 / 0.12 / 0.96 / 1 | 150 / 171 / 193 | 63 | 61 | 617 |

**Data volume.**
- **v1 minimum** (f < 0.05, R = 3): about 35 projects (54 at half yields), 100–130 files downloaded at p50–p90 (160–200 at half yields), which is **≈ 0.1–0.2 TB**.
- **v2_sep tier-1 full stop:** about 210 runs (≈ 0.23 TB), or 410 at half yields (≈ 0.45 TB).
- **Pool use:** under 5% of the 1,465 eligible projects, so a not-evaluable outcome burns nothing scarce.
- **Binding criterion:** with f < 0.05 the binding criterion is the capped tail proxy (it reaches only 1.03–1.04× its minimum when the stop occurs). With f < 0.20 it is J = 20.

### 4b. Gate behaviour at the v1-minimum states (R = 3, pilot yields, own threshold)
κ = n_tail(v1) / hyperscore proxy.

| scenario | f < 0.20, κ 1 / 0.5 | **f < 0.05, κ 1 / 0.5** | single-MS1, κ 1 / 0.5 |
|---|---|---|---|
| clean: pass / invalid / fail | 0.94 / 0.07 / 0; 0.92 / 0.08 / 0 | **0.90 / 0.10 / 0.005; 0.87 / 0.08 / 0.055** | 0.95 / 0.05 / 0; 0.93 / 0.06 / 0.01 |
| margin 0.60 (τ 0.10): pass | 0.005; 0.015 | **0.01; 0.005** | 0.01; 0.005 |
| leak 0.90 in 25% of labs: pass | 0.015; 0.01 | 0.005; 0.015 | 0.02; 0.01 |
| contamination δ = 0.02 (hyperscore 0.52, learned 0.53): pass / invalid | 0.22 / 0.78; 0.44 / 0.53 | 0.54 / 0.46; 0.53 / 0.33 | 0.38 / 0.62; 0.54 / 0.43 |
| contamination δ = 0.05: invalid | 1.00; 0.995 | 0.985; 0.89 | 1.00; 0.94 |
| median capped n_tail | 2,334; 1,151 | 1,252; 626 | 2,019; 1,031 |
| G1 clean pass (n_top, J) | 0.985 (15,740; 20) | 0.99 (17,975; 34) | 0.94 (12,313; 34) |
| G3 clean pass (all 10 bins) | 1.00 | 1.00 | 1.00 |

MC SE is ≈ 0.02 near 0.9 and ≈ 0.01 near 0.02.

### 4c. D28 quantities (budget 200 runs; R = 3; own threshold)

| | f < 0.20 | **f < 0.05** | A4/DN3 (D28, reference cell) |
|---|---|---|---|
| P(v1 minimum) | 1.00 | **1.00** (0.97 at ×0.5 yields) | 0.08–0.43 |
| P(clean v1 passes) = P(min) × tail × G1 × G3 | 0.91–0.92 | **0.86–0.89** | ≈ 0.06–0.4 |
| P(marginal leak caught) = P(min) × P(fail \| margin 0.60, τ 0.10) | 0.92–0.93 | **0.905–0.91** | ≈ 0.07–0.4 |
| False pass at the margin \| evaluable | 0.005–0.015 | 0.005–0.01 | 0.003–0.057 |
| Volume / time | ≈ 0.07 TB / hours | **≈ 0.1–0.2 TB / hours** | 0.75–1.9 TB / 3–4 days |

**Caveat that matters.** "P(clean v1 passes)" assumes a valid null: hyperscore tail target fraction 0.5 in expectation. If the residual target excess of section 2 is about 0.02, P(pass) falls to about 0.5 and the outcome is mostly "invalid pool". At 0.05 it is almost surely invalid. This **fails safe**: not evaluable, never a false pass. But it is the main unknown, and it can be learned only by running step 5.

---

## 5. Recommendation
**GO** (D28.2: simulated P(v1 minimum) ≥ 0.8 at feasible cost). P(v1 minimum) is 0.97–1.00 within 200 runs (≤ ~0.2 TB) in every simulated cell. The pool is 1,465 projects from 969 labs. Cost is hours of download and under 1 h of compute; the amendment and its independent review are the long pole.

Conditions to write into the amendment before any listing:
1. **Ladder screen f_max = 0.05**, not 0.20 (section 3c). Removal is pair-closed and label-free, as in section 1.2.
2. **Removal must not use any hit of the null spectrum itself** (section 1.3, R6). I comes from rank-1 IDs at 1%. G0-style swap test with certification frozen.
3. **Strengthen certification against R4.** Use project-pooled TDC at PSM q ≤ 1e-4, or agreement with a second fixed engine. Measure its label-free yield cost in the implementation dry run (per-run 1e-4 is not resolvable).
4. **Drop contaminant pairs from the gate DB.** Use the null's own run as host. The A4 statistics, cap, IUT and α/2 apply unchanged. Every n_min floor is taken from A4 as a minimum.
5. **State the validity-control risk in the amendment:** δ ≥ 0.02 makes "invalid pool" likely. An invalid pool is a completed test under A4.6 (next attempt at α/4). Adopt the strict screens now rather than after a run.
6. **Exclude every walked project from future T4 redraws.** Re-run `c9_power.py` with the final rules and the measured certification yield before the tag.

If the PI or reviewer judges condition 3 or the validity risk unacceptable, the D28 else-branch (T1 under A5) applies unchanged.

## Files
- `scripts/c9/pilot_yield.py`: certification, removal, label-free counts and timings. Output: `per_run/*.json` and `pilot_yield.log`.
- `scripts/c9/certify_strict.py`: certified counts at q 1e-3 vs 1e-4. Output: `certify_strict.json`.
- `scripts/c9/pool_count.py`: metadata-only pool and labs. Output: `pool_count.json`.
- `scripts/c9/c9_power.py`: walk and gate simulation, reusing `reviews/checks/a4_power.py` and `vulcan.gates.a2_stats`. Output: `c9_power.json` and `c9_power.out`.
- `scripts/c9/make_metrics.py`: collects everything into `metrics.json`.

## Limitations
- **Representativeness.** There are 12 runs from 11 projects. D8h is a seeded random draw from the same frame, so it is representative but small. The ×0.5-yield cells cover the uncertainty in level, not in shape.
- **The tail proxy is hyperscore-based.** v1's own tail enters only through κ ∈ {1, 0.5}.
- **The host threshold** was the run's full-search 1% hyperscore threshold. The amendment's own-run mixture threshold will be somewhat higher, which gives fewer tail spectra; this is covered by κ and the D0-host cells.
- **v2_sep's V-rows** were not simulated beyond the full-stop size. V12 needs about 40,000 per label per bin, which drives A_cand ≥ 450,000.
- **The f < 0.05 rule** was chosen after seeing label-free tail counts only. It is stricter and symmetric, but it is a design choice informed by pilot data, and the amendment must say so.
