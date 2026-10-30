# Independent review of Amendment A4 (draft revision 0, commit d42062c)

- **Reviewer:** independent statistician (agent). I had no role in A1–A4, the DN2/DNIT walks, E7/E7b, M7, the T1 driver or any gate run.
- **Date:** 2026-10-28.
- **Read:** `reviews/amendment_A4_draft.md`; `reviews/checks/a4_power.{py,json,out}`, `a4_pool_count.{py,json}`; DESIGN §6 (§6.0, §6.0b, §6.0c, §6.3–6.4, §6.6); RESEARCH_TREE D9–D25; `reviews/amendment_A2_draft.md`; `reviews/amendment_A3_draft.md` and `amendment_A3_review.md`; `reviews/checks/a3_power.{py,out}`; `scripts/a2/IMPLEMENTATION_NOTES.md`, `scripts/a3/IMPLEMENTATION_NOTES.md`; `vulcan/gates/a2_stats.py`, `a3_stats.py`; `config/a2.yaml` (keyword lists); `results/2026-10-27_A2_dn2/{NOTES.md, admission_counts.json, DEVIATIONS.md}` and `data/manifests/dn2_select.log.json` (decision reasons only); `results/2026-10-27_3c_itcid/followup/NOTES.md` §3; `results/2026-10-28_A3_preconditions/NOTES.md`; `bench/sealed_run.py`, `bench/t1/a3.py` (headers).
- **Not read:** anything under `data/sealed/`; any sealed T1–T5 file; `/tmp/vulcan_a2`; any DN2 or DNIT per-label quantity (none exists); anything under `results/2026-10-28_E7b_fair_speed/`. `results/SEALED_LEDGER.md` was checked only to confirm it has no entries (header only). No PRIDE call was made.
- **Local metadata I did use:** the 2026-10-26 crawl table and the cached PRIDE *project records* from the same crawl (`data/cache/pride_meta/projects_p*.json.gz`: title, keywords, description, sample/data protocol). These are the fields the A4 protocol screen will read. I computed only **aggregate** pass rates per tier; I did not list or inspect any project's files, and I compiled no per-project outcomes into this review except the two exclusion leaks in M5.
- **My checks** (scratch, `…/scratchpad/reviewer/`; synthetic plus label-free metadata; nothing in the repo was written or overwritten):
  - `pool.py` and `orgkw.py` rebuild the A4 pool with the same logic as `a4_pool_count.py` (198 / 259 / 135, identical), and compute its genus/lab structure and organism-name keyword self-matches.
  - `proto.py`, `proto2.py` and `proto3.py` compute protocol-screen pass rates on the cached crawl text, under the draft's lists and under variants.
  - `walk_real.py` repeats the S4 walk with the **real** genus/lab structure, a reference-proteome step, and the draft's own yield model (imported from `a4_power`). It runs 200 walks per cell.
  - `tailsim.py`, `tailsim2.py` and `tailsim3.py` test tail-G2 at the draft's v1-minimum walk states (300–400 walks): concentrated and size-correlated leaks, a cluster jackknife-t hull, and cap shares of 20%, 10% and 5%.
  - I imported functions from `a4_power.py` and did not run its `main()`, so `a4_power.json/.out` are untouched.

---

## Verdict: **do not adopt revision 0. Adopt after a revision 1 that fixes B1–B3 and M1–M7, with S4/S5/S6 re-run.**

The architecture is sound:
- label-free admission and stopping;
- the margins are unchanged;
- interval hulls for G1, the TOST rows and G2-original;
- geometric α spending;
- separate agents;
- A1 and A2 are reported.

Most statistic changes really are stricter. But revision 0 has three blocking problems:
1. One added keyword (`IgE`) matches the word "digest". That empties the walk.
2. The evaluability claim rests on a walk model whose key inputs are optimistic by a factor of 2 or more. Re-derived, DN3 is more likely than not to be **not evaluable** under the draft's own rules.
3. The tail-G2 rule is anticonservative under plausible cluster heterogeneity. The cluster cap can make a pass *easier* on given data. Both contradict the draft's "stricter on any given data" principle and its mandatory "≤ 1.5α" wording.

None of these needs a change of architecture. The fixes are local and can be checked by simulation.

---

## BLOCKING

### B1. `IgE` in A4.3.5(a) matches "digest", which removes about 95% of the pool (line 146)
- **The problem.**
  - A2's lists are matched case-insensitively and without word boundaries (`config/a2.yaml` `regex_flags: IGNORECASE`; §6.0b point 2). A4.3.5(a) appends `…|allergen|IgE` to them.
  - A4.1.6 (the stricter reading wins) obliges the implementer to match the additions the same way. `ige` occurs in "d**ige**st", "d**ige**stion" and "ni**ge**r".
- **Measured on the cached crawl text of all 592 candidates** (the same fields as `protocol_fields`), the share whose text contains `ige`:
  - 0.84 (tier I), 0.84 (tier II), 0.76 (tier III);
  - essentially all of these are "digest".
- **Protocol-screen pass rates, unconditional:**

  | Lists | Pass rate |
  |---|---|
  | As drafted | 0.045 (I), 0.046 (II), 0.052 (III) |
  | Same lists without `IgE` | 0.177 (I), 0.166 (II), 0.163 (III) |

- **Consequence.** In my walk re-run with the drafted lists (protocol pass 0.05), P(v1 minimum) = **0.00** in every cell (median 11 admitted projects, J = 10). DN3 would be not evaluable with certainty, and every reached candidate would be burned (A4.3.3 "Downstream").
- **Fix.**
  - Match every A4 *addition* with a leading word boundary, `\b(?:…)`, and match `IgE` case-sensitively as `\bIgE\b`. Say explicitly that §6.0b point 2 ("verbatim, no word boundaries") covers A2's lists only.
  - Keep the dry-run test in A4.4.3, and extend it with a regression test asserting that "digested with trypsin" passes the protocol screen.

### B2. The evaluability claims of A4.7.2 rest on optimistic walk inputs; re-derived, DN3 is likely not evaluable
A4.7.2 claims that "P(v1 min) ≥ 0.95 once tier II is allowed" and P(full) of 0.78–1.00. `simulate_walk` (a4_power.py lines 276–348) has four optimistic inputs.

**(a) Protocol pass rate (A4.7.1, line 325).**
- DN2's realised rate was **19/48 = 0.40**, not 0.45:
  - the genus-slot check precedes the protocol check (IMPLEMENTATION_NOTES 14), so the 10 slot skips are not protocol passes;
  - the 48 were `t4_frame` rows that had already passed BAD_KW on title and keywords at the frame stage.
- A4's frame (A4.3.1; `a4_pool_count.py` lines 79–85) has **no** title/keyword filter. So the per-candidate rate must cover title, keywords and protocol together.
- Measured on the crawl text, the A2-verbatim lists pass:
  - **0.18 (I), 0.17 (II), 0.17 (III)**;
  - conditional on passing title and keywords: 0.30, 0.23, 0.21.
- The dominant hits are search-setting and generic terms in `dataProcessingProtocol` and descriptions: `acetyl` (Protein N-term acetylation as a variable modification), `phospho`, `tissue`, `N-termin`, `standard`, `glyco`, `immunoprecip`, `human`.
- With prefix word boundaries on all lists the rate is still 0.21. With experiment-type terms restricted to the sample-description fields it is 0.24.
- The draft's 0.70 sensitivity scenario corresponds to no rule anyone has adopted.

**(b) Reference-proteome step.** This step is not modelled. In DN2 it passed 16/19.

**(c) Cluster structure.**
- Clustering is modelled as a random link with probability 0.15 (lines 277 and 318), and a project that would overfill a cluster simply opens a new one.
- The real pool is highly concentrated:
  - tier II: Oryza 48, Zea 34, Solanum 32, Triticum 24 (`a4_pool_count.json`);
  - with ≤ 3 projects per cluster, at most Σ_g min(n_g, 3) projects are admissible: **140 (I), 101 (II), 94 (III)**, before lab links merge genera.
- All 38 *Chlamydomonas* candidates (the largest tier-I genus) are excluded by `HLA` matching the organism's own name (`orgkw.py`). The DN2 implementation notes record the same effect at frame level.

**(d) Marker-screen cost and pass rate.**
- A stage-2 exclusion costs R downloads, not min(R, 3) (lines 303–305).
- p_marker should fall with R, because the pooled screen covers up to 24 runs. At 3 runs, 3 of DN2's 6 admitted projects already showed one non-keratin panel member (NOTES table).

**Re-derived** (`walk_real.py`: real genus/lab components with labs normalised by lower-casing and stripping non-alphanumerics, proteome 16/19, p_files 0.93, p_marker 0.5, the draft's yield model and stopping rule, tiers I→II→III, 200 walks per cell):

| p(protocol) | median runs/project | P(full stop) | P(v1 minimum) | median J | median admitted | median cluster-full skips |
|---|---|---|---|---|---|---|
| 0.45 (draft) | 6 / 12 / 24 | **0.44** / 0.99 / 1.00 | 0.84 / 1.00 / 1.00 | 49 / 36 / 28 | 74 / 51 / 38 | 129 / 31 / 15 |
| 0.30 | 6 / 12 / 24 | 0.18 / 0.77 / 0.99 | 0.43 / 0.94 / 1.00 | 39 / 34 / 28 | | |
| **0.18 (A2 lists, measured)** | 6 / 12 / 24 | **0.00 / 0.18 / 0.54** | **0.07 / 0.46 / 0.78** | 27 / 27 / 25 | 38 / 37 / 34 | |
| 0.05 (as drafted, B1) | any | 0.00 | 0.00 | 10 | 11 | |

These are model numbers with a crude per-candidate Bernoulli for the protocol screen, but the direction and size are not in doubt. Under the A2 lists, DN3 reaches the v1 minimum only if projects have many usable runs, and v2_sep tier 1 (which needs the full stop) is unlikely.

**Why it is blocking.**
- A not-evaluable DN3 spends no α, but it permanently burns every reached candidate (A4.3.3).
- If a T1 look happens first, the *next* null-set amendment would be designed by people who have seen T1 (see the sequencing review).

**Fix (all label-free, all before the tag).**
1. Re-run S4 with:
   - the real genus/lab component structure from the frozen frame;
   - the proteome step;
   - protocol-pass rates measured on the frozen crawl text, by tier and per project. The text is local metadata; a deterministic per-project pass flag is better than a Bernoulli;
   - stage-2 costs, and p_marker decreasing in R.
2. Then choose yield levers so that P(v1 minimum) ≥ 0.9 in the pessimistic cell. Options:
   - word-boundary matching of acronym terms (HLA, MHC, UPS, DIA, PRM, SRM, TMT, BHI, FBS, FCS, HEK);
   - restricting experiment-type terms (acetyl, phospho, glyco, N-termin, ubiquitin, crosslink…) to title, keywords, description and sample protocol rather than `dataProcessingProtocol`. Enrichment and labelling are also caught by the frame `enrichment none` filter and the reporter-ion check;
   - raising the per-project file cap;
   - entering tier III earlier.

   Each option changes selection, not any gate. It is legitimate now because no DN3 file has been listed and no label exists (P2). The D16 objection was to changing DN2's list *after* DN2's yield was seen.
3. State the evaluability target and the cell it is met in.

### B3. Tail-G2 is anticonservative under cluster heterogeneity, and the cap is a re-weighting that can ease a pass; the "≤ 1.5α" wording is not supported
**Concentrated leak.**
- The draft's own S5 gives false-pass rates of 0.045–0.060 (leak 0.90 in 25% of clusters) and 0.050–0.098 (leak 1.0 in 20%), against the nominal **0.025** (A4.7.3 table).
- I reproduced 0.060 and 0.090 (κ = 1), and 0.043 and 0.070 (κ = 0.5). The clean pass rate was 0.96.
- Adding a leave-one-cluster-out jackknife-t bound to the hull does **not** help (0.060, 0.0875). The cause is few effective clusters (median J_eff of the capped tail proxy at the v1 minimum is **9.1**), not the interval type.
- The TOST rows show the same thing at J = 10–20: 0.035–0.038 at the margin with τ = 0.005, against 0.025 (S2).

**The cap is not "only removing data" (A4.3.6, line 174).**
- Subsampling clusters at different rates re-weights the estimand. If a leak sits in the *largest* clusters, the cap down-weights exactly those clusters.
- Simulation: the leak is placed in the largest clusters so that the uncapped tail-weighted target fraction is exactly 0.60 (`tailsim3.py`, 300 walks):

  | Gate | Pass rate |
  |---|---|
  | Uncapped | 0.00 |
  | 20%-capped (draft) | **0.067** |
  | 10%-capped | **0.33** |

- So the cap can make a pass easier on given data, which contradicts A4.1.3. A2.1 already lists re-weighting of projects as impermissible *after* the gate. Before the gate it is allowed (label-free), but it is not "stricter".

**Fix (all three parts).**
1. **Intersection–union over capped and uncapped pools.** A row passes only if it passes on both the capped and the uncapped pool. This is never looser than either pool.

   | Scenario | Capped only | Capped and uncapped |
   |---|---|---|
   | Size-correlated leak, 20% cap | 0.067 | 0.00 |
   | Size-correlated leak, 10% cap | 0.33 | 0.00 |
   | Clean scorer at the v1 minimum, 20% cap | 0.93 | 0.88 |
   | Clean scorer at the v1 minimum, 10% cap | 0.92 | 0.88 |

   Apply the same intersection–union to every A-based row, since the cap subsamples A for all of them.
2. **Cap share 10%, not 20%** (`tailsim2.py`):

   | Cap | False pass, leak 0.90 in 25% | False pass, leak 1.0 in 20% | Clean pass at v1 minimum | Capped tail proxy | J_eff |
   |---|---|---|---|---|---|
   | 10% | 0.013 | 0.033 | 0.90 | ≈ 870 (from ≈ 1,256) | 14.5 |
   | 5% | 0 | 0 | 0.57 | | |

   A 5% cap costs too much power. The stopping rule must then be evaluated on 10%-capped counts, and a J_eff criterion should be added (M2).
3. **Wording.** Drop "procedure-level ≤ 1.5α" from the mandatory wording (A4.6 line 300, A4.9 line 429), or qualify it as "nominal; simulated per-attempt false-pass rates under concentrated heterogeneity are …". See also M3, which shows that 1.5α is not correct for G3 even nominally.

---

## MAJOR

### M1. Replacing the project unit by the coarser A4 cluster is not monotone on given data (A4.5.2, line 250)
- A coarser clustering can *narrow* a percentile bootstrap on a specific dataset. For example, two projects in one genus with opposite deviations average out inside the merged cluster.
- So "which can only widen the bound" is true in expectation only. A2's tail-G2 bound used the project bootstrap; A4 drops it.
- **Fix.** UCB = max(CP or delta bound, project-cluster bootstrap, A4-cluster bootstrap). Do the same for every cluster term (G1 jackknife, TOST and G2-original bootstrap): compute both units and take the outer envelope. Then A4.1.3 holds on any data.

### M2. The cap uses only tail-proxy shares, so A_cand and |A| shares stay uncapped (A4.3.6); `cap_factors` silently disables the cap at J < 1/share (a4_power.py line 250)
- In DN2 the tail-proxy rate per A spectrum varied 21-fold between projects (0.0004–0.0085, A4.7.1). A cluster with a low tail rate and large |A| can dominate G1, G3, V12 and V6 while staying under the tail cap.
- The n_min and power simulations assume log-normal weights with log-sd 1 at J = 20. They do not cover that case.
- **Fix.**
  - Set f_j to the minimum of the water-filling factors for the tail proxy, A_cand and |A| (and |B| for B), all label-free.
  - Add J_eff = 1/Σw_j² ≥ 15 on the capped tail proxy and on capped A_cand to both the full-stop and v1-minimum criteria.
  - Where J·share < 1, the rule must return "not evaluable", not "no cap".

### M3. The multiplicity statements need correcting; the geometric scheme itself is sound
- **(a) A1's G3 had no controlled level.**
  - A1's G3 was a point rule. A2.8 shows it passed about 45% of the time at a 1-pp deficit, so "≤ α + α/2" (line 300) is false for G3.
  - A1's tail-G2 had no cluster term, so its α was nominal under the source heterogeneity M3e found.
  - The correct statement is per row: the A4 attempt has nominal level α_k, *conditional* on A1 having failed. The unconditional sequence bound is Σ (actual operating characteristic of each attempt), which for G3 is dominated by A1's point rule.
- **(b) Inconsistency between A4.6 and A4.9.**
  - A4.6 counts "only completed tests", but A4.9 (line 432) requires α/4 for the next attempt after any "not evaluable" outcome.
  - **Fix.** Define a completed test as an attempt in which step 5 ran for that scorer, so that any per-label quantity was computed, whatever the verdict. A label-free stop before step 5 (as in A2) is not a test. A4.9 should then say "α/4 if step 5 ran; α/2 otherwise".
- **(c) Counting A2 as zero is correct.** A2 stopped on label-free counts and no per-label quantity exists (D24, NOTES). The geometric split is a reasonable choice among the honest options:
  - keeping α for A4 and reporting 2α is less protective against test-until-pass;
  - fallback or α-investing schemes do not apply, because A1 rejected nothing.

  Keep it.

### M4. The organism class is broader than A2's and includes taxa usually grown in animal-peptone media (A4.3.2, line 110; Appendix A)
- A4.3.2 says A2's class "is kept as tier I". It is not:
  - A2 used a 40-genus ALLOW list;
  - A4 uses all microbes minus a DENY list.
- The pool contains *Dictyostelium* (19 projects), *Paramecium* (6) and *Tetrahymena* (5). Their standard media are proteose-peptone based (HL5, PPYS/SPP) or bacterised.
- The new media keywords appear in the crawl text of ≤ 2.5% of candidates (`proteose peptone` 0.000, `HL5` 0.002), so the keyword screen will not catch them.
- Animal-peptide carry-over lands in B. A4.10.2 concedes that B contamination can produce a **final fail** of G2-original/G5 for a clean v1.
- **Fix.**
  - Deny Amoebozoa and Ciliophora by lineage, or name the three genera.
  - Deny Mollicutes by lineage. Appendix A lists *Mycoplasma*, *Mycoplasmoides*, *Metamycoplasma* and *Ureaplasma*, but not the 2018 split genera *Mycoplasmopsis*, *Mesomycoplasma* and *Malacoplasma*, which grow in serum media. None is in the current pool, but the taxonomy is rebuilt at the tag.
  - Deny *Leptospira* (BSA medium) and the dermatophytes (*Trichophyton*, *Microsporum*, *Epidermophyton*).
  - Describe tier I honestly as broader than A2's class.

### M5. The exclusions miss projects that prior work reached (A4.3.3)
- **PXD002407** (*Phaseolus vulgaris*, tier II) appears in an alternative T4 draw: `results/2026-10-27_T1driver_dryrun/t4_selector/frame/t4_selection.tsv`, 8 files, tracked. Its PRIDE file listing is cached (`data/cache/pride_meta/files_PXD002407_p000.json.gz`). It is a T4-reachable project whose files have been listed, yet it is in the A4 pool. Item 1 reads only `data/manifests/*.tsv`.
- "Except the two frame files" (line 119) does not name the second file, and `a4_pool_count.py` (line 58) skips only `t4_frame.tsv`.
- **Fix.**
  - Exclude every `PXD` accession appearing in any tracked file under `data/manifests/` or `results/`: selection logs, `*walk_events*`, every `t4_selection*`, fetch logs.
  - Also exclude every accession with a cached PRIDE file listing.
  - Name the frame files explicitly.
  - Before the tag, an agent permitted to open E7b outputs must confirm that the E7b cohort ("2,987 runs/h", untracked `cohort_runs.jsonl`) and the 3i/E9 throughput cohorts are covered. I did not open them.

### M6. The T1 safeguards must be in the tagged text, not in a future "A4.12 (revision 1)" (A4.8, lines 414–418)
- A4.8 incorporates "by reference" safeguards that do not yet exist.
- My separate answer (`reviews/a4_t1_sequencing_review.md`) requires two things:
  - no T1 look before A4 is tagged **and** DN3 steps 2–3 are committed with the v1 minimum met;
  - DN3 step-4b/5 agents and every A4 deviation reviewer must not have seen T1.
- Those requirements must be normative A4 text at the tag. A4.4's "disclose which T1 quantities the adopter has seen" is necessary but not sufficient.

### M7. Options and undefined items remain in normative text (contrary to A3 Q6)
Each of these must be closed in the tagged text:
- **A4.11.3.** The PI "may decide" to reuse DN2. Close it as **no reuse**, consistent with P4 and D24 ("DN2 is spent").
- **A4.7.4.** Admitted RAWs are kept "if quota allows". Fix the policy (check the quota before the tag), and record re-fetch sha256 either way.
- **"Normalised PRIDE lab string"** (A4.3.4). Specify the algorithm in `config/a4.yaml`. Consider also linking by submitter or lab-head name, because the `lab` field mixes PI names and facility names.
- **The 2,000-run hard cap.** Say that reaching it ends the walk and A4.3.8 applies.
- **A4.3.8.** "will then probably be not evaluable": delete "probably", since the row rule decides.
- **The ce_s 24-h rule.** Say who measures, on what, and when, as a dated label-free note before step 5.

---

## MINOR

1. **Admitted vs cluster counts (A4.7.2, line 353).** "37–50 admitted projects in 32–67 clusters" mixes cells. In `a4_power.json` the median admitted count is 37–80 (79–80 at 6 runs per project), and a cell's J can never exceed its admitted count.
2. **Wrong cross-reference (A4.3.6, line 181).** The residual of the concentrated leak is said to be "stated in A4.10.7", but A4.10.7 is about cluster definition. Add the residual as its own item, with the B3 numbers.
3. **Noisy n_min estimates (A4.5.3).** The values come from 1,000 (TOST) and 400 (G1) replicates, with MC SE of about 0.01–0.015 near 0.90. The "smallest n with p̂ ≥ 0.90" search is biased low. Use at least 4,000 replicates, or require p̂ ≥ 0.92. This affects power only.
4. **Clean TOST power assumes one base rate across clusters.** Heterogeneous base recall or keep rates across organisms, with independent T and D allocations, widen the cluster interval. Add a logit-sd of 0.5 in the clean scenario before fixing n_min for V6 and V12.
5. **S6 uses G2-original power without the A4 cluster term** (`G2.g2_power`, a4_power.py line 570). This is immaterial at B ≈ 33 M, but say so.
6. **G1 between-cluster composition.** A clean scorer can show a pooled Spearman ρ ≠ 0 if organisms differ in both mean score and mean classifier logit, for example through precursor mass or peptide length. The jackknife does not remove this; the hyperscore control catches it only if hyperscore shares it. Report the within-cluster pooled ρ (`a2_stats.fisher_pooled_rho`) as reported-only.
7. **Seed predictability.** The walk order is a deterministic function of the frozen frame and seed 47, so anyone can compute it before the tag. The adoption note should state that nobody computed the DN3 order while drafting, the DENY list or the keyword choices. Alternatively, take the walk seed from a value fixed after the tag. The risk is low because no label exists.
8. **Download time (A4.7.4).** "1.2–2.4 h" needs 0.1–0.2 GB/s. DN2's per-stream rate was about 3 MB/s (1.04 GB per 5–7 min), so "≥ 8 ordered streams" gives about 10 h for 0.85 TB. That is still fine for a 0.5–1 day walk, but correct the figure. The G3 exhaustive cost over about 1,800–2,000 mixtures is unmeasured, and B (33–38 M spectra) dominates compute. If the dry run shows G3 is the bottleneck, a seeded label-free subsample of B to a fixed multiple of |A|, fixed at the tag, is a legitimate pre-tag change. It requires an explicit amendment of A2's "all retained B".
9. **v1 closure guard (A4.4.3).** A4 deliberately changes v1's *statistics* (levels, hulls, n_min). The guard should freeze v1's scorer, feature and search closure against `prereg-v0.4.2`, and its statistics against `vulcan/gates/a4_stats.py` at the A4 tag. Say which is which.
10. **G1 clean pass at J = 10.** The text says 0.73; S3 gives 0.713.
11. **Optional: exclude same-lab projects.** Exclude projects whose lab or submitter matches an A1 or DN2 null project, all of them rather than only the contaminated one, so the rule does not depend on outcomes. This blocks a shared-instrument carry-over route.
12. **Tier-III spectra.** They are scored inside QE (D0) hosts. Keep the per-instrument-class reporting (A4.5.5). Without a gate requirement, a QE-specific leak can be diluted (A4.10.1). The strongest cheap remedy is to require that the tail-G2 and G1 point estimates on the A1-class subset do not exceed the margin when that subset has n_tail ≥ 200. This can only make passing harder.

---

## What should not change
- **P1/P2 discipline.** Every admission, cap and stopping quantity is label-free; membership hashes are committed before any gate statistic; the walk agent and the step-5 agent are different; step 5 is one session, v1 committed first.
- **Every margin of A1 and A2 is unchanged.** Not evaluable blocks; a fail is final (P4); DN2 is excluded (no reuse).
- **The hull construction for G1, the TOST rows and G2-original.** The outer envelope of the A2/A3 interval at the halved level and a cluster interval is monotone, so it is stricter on any data (M1 extends this to tail-G2).
- **The validity controls are not hulled** (Wilson 95%, G1 control 99%). Widening them would make contamination *harder* to detect, which would be a relaxation.
- **Geometric α spending, A1 counted as a completed test, A2 counted as zero**, and one level per row shared by v1, v2_sep and ce_s.
- **Two-stage marker screen** (stage 2 pooled over all runs), the file-level SHA-1 and (name, size) guard, cluster-full skips made before any listing, permanent exclusion of every reached candidate, the byte-for-byte host-threshold check, and P1 swap tests (including the cap).
- **Reporting.** Per-tier, per-instrument and per-cluster tables; both interval components; the cap factors; A1 and A2 quoted in every report; ce_s licenses nothing.

## Numbers I re-derived (summary)
- **Pool size.** 198 / 259 / 135 (identical to the draft). Σ_g min(n_g, 3) = 140 / 101 / 94. *Chlamydomonas* (28 + 10 candidates) self-matches `HLA`.
- **Protocol pass, unconditional, per tier I / II / III:**

  | Lists | I | II | III |
  |---|---|---|---|
  | A2 verbatim | 0.18 | 0.17 | 0.17 |
  | As drafted (B1) | 0.045 | 0.046 | 0.052 |
  | Prefix word boundaries | 0.22 | 0.19 | 0.26 |

- **DN2 realised protocol pass:** 19/48 = 0.40. **Proteome step:** 16/19.
- **Walk:** see the B2 table.
- **Tail-G2 at the draft's v1-minimum states:**
  - concentrated-leak false pass 0.060 / 0.090 (κ = 1) and 0.043 / 0.070 (κ = 0.5) under the A4 rule;
  - unchanged with a jackknife-t hull;
  - 0.013 / 0.033 with a 10% cap;
  - size-correlated leak: 0.067 (20% cap) and 0.33 (10% cap) capped-only, 0.00 with the capped-and-uncapped intersection–union;
  - clean pass 0.88–0.96 across these variants.
