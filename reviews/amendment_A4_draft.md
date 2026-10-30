# Amendment A4 (DRAFT, revision 3): one fresh null set (DN3) gating v1 and v2_sep; T1 escrow

- **Status:** draft for adoption as DESIGN.md §6.0d. It is not binding. DESIGN.md, the tags and every A2/A3 file are unchanged.
- **Author:** statistician (agent). I did not build the A2/A3 implementations, the DN2 walk, E7/E7b or M7. I will not run any DN3 step, and I have seen no T1 output.
- **Date:** 2026-10-28.
- **Revisions:**
  - Revision 0 is d42062c.
  - Revision 1 answers the independent review `reviews/amendment_A4_review.md` (verdict on rev 0: do not adopt; 3 blocking, 7 major, 12 minor). Every change is listed in A4.12.
  - The T1-sequencing opinion of the same reviewer, `reviews/a4_t1_sequencing_review.md`, is made normative in A4.8.
  - **Revision 2 implements the PI decisions of RESEARCH_TREE D26:**
    - (1) DN2's 6 admitted projects enter DN3 as walk tier 0, with all their files, under every A4 screen and cap (A4.3.10);
    - (2) the word-bounded keyword rule is kept;
    - (3) T1 runs into escrow (A4.8, normative);
    - (4) the thresholds are unchanged, and DN3 stays sized for v1 and v2_sep.

    Revision 2 is reviewed by a second, new independent reviewer (`reviews/amendment_A4_review2.md`). The changes are listed in A4.13.
  - **Revision 3** answers the second review (verdict on rev 2: do not adopt; 2 blocking, 6 major, 12 minor). Every resolution is in A4.14.
    - The walk and tail figures now come from `reviews/checks/a4_r3_walk.{py,json,out}`. It implements the full selection spec: item 5 on the pool and on tier 0, order-insensitive person keys, the three-factor cap, the fixed tier-0 shuffle, tier-0 marker evidence, the empirical DN2 raw-file distribution as the reference cell, and IUT precedence.
    - Four points need a PI decision before adoption (A4.11).
- **Basis:** DESIGN §6 including §6.0 (A1), §6.0b (A2) and §6.0c (A3); RESEARCH_TREE D9–D26; `reviews/amendment_A2_draft.md`; `reviews/amendment_A3_draft.md` and `reviews/amendment_A3_review.md`; `reviews/checks/a3_power.*`; `scripts/a2/IMPLEMENTATION_NOTES.md`; `scripts/a3/IMPLEMENTATION_NOTES.md`; `results/2026-10-27_A2_dn2/{NOTES.md, admission_counts.json, DEVIATIONS.md}` (label-free records only); `results/2026-10-27_3c_itcid/followup/NOTES.md` (DNIT); `results/2026-10-28_A3_preconditions/NOTES.md`.
- **Checks (synthetic plus label-free/metadata only; seed 131):**
  - `reviews/checks/a4_power.{py,json,out}` (rev 0): DN2 label-free yield model, n_min under the A4 rule, TOST/G1 cluster robustness, idealised walk.
  - **`reviews/checks/a4_r1_walk.{py,json,out}` (rev 1):** the realistic walk on the actual candidate pool and the revision-1 tail-G2 rule. The walk uses:
    - per-project protocol-screen flags, computed on the locally cached 2026-10-26 PRIDE project records;
    - the real genus/lab cluster structure;
    - the reference-proteome step;
    - stage-2 marker costs.
  - **Revision 3:** `reviews/checks/a4_r3_walk.{py,json,out}` implements the full selection spec. **It is the source of every evaluability and tail figure in this revision.**
  - **Revision 2:** `a4_r1_walk.py` adds tier 0 (DN2's admitted projects), with DN2's label-free per-project yields and the raw-file counts from DN2's own listing events (`a4_power.dn2_tier0_projects`).
  - `reviews/checks/a4_pool_count.{py,json}` and `a4_genus_lineage.json`: a metadata-only pool count.
  - No PRIDE listing or call was made for any candidate, no file was downloaded, no spectrum was read, and no label was used.
- **Not read by the author:** any sealed file; any T1 output (none exists; the ledger is empty); any DN2 admission array beyond `len(codes)` and `len(A_spec)`; `/tmp/vulcan_a2`; any E7b output.

---

## A4.0 Why A4, and what it is

**Attempt history.** It is quoted in every report (A4.9).

| Attempt | Null set | Scorers | Outcome | Per-label quantities computed? |
|---|---|---|---|---|
| A1 (`prereg-v0.4.1`) | M3e nulls | v1 | tail-G2, G1, G3 **FAIL**; final for those nulls (D12). The contamination explanation came after the failure and is never used | yes |
| A2 (`prereg-v0.4.2`) | DN2 | v1 | **Not evaluable** on label-free counts: 6 projects / 18 runs, \|A\| 24,581 vs 120k, tail proxy 22 vs 400 (D24). Step 5 never ran | **no** |
| A3 (D22; freeze pending) | DN2 | v2_sep (gated), ce_s (licenses nothing) | Moot: DN2 not evaluable. A3.3 and the A3.1.3 freeze stand | no |
| 3c DNIT | DNIT (IT-CID) | `itcid_linear_pred` | Futile and spent (D20). **Outside A4** | no |
| Dev diagnostics, not verdicts | M3 mixtures / D12 nulls | v2_sep G2-original and G5 pass (E7). ce_s tail/G1 on contaminated nulls (M7) | Diagnostics only | yes (dev) |

**Why A2 failed** (from DN2's label-free records):
- **The frame was too small for the selection rules.** It held about 82 projects from genera on A2's ALLOW list, and A2 took 3 files per project and one project per genus.
- **Per-run yields were low:**

  | Per run | DN2 actual | A2/A3 assumption |
  |---|---|---|
  | A spectra | 1,366 | — |
  | A spectra with an exhaustive candidate | 39% | 60% |
  | Tail proxy (share of A) | 0.09% | 0.33% |
  | B / \|A\| | 40× | ≤ 3× (A3) |

**What A4 is.**
- One new fresh set, **DN3**, selected under rules fixed here, before any DN3 candidate is listed.
- It gates v1 and v2_sep with A2/A3 statistics. Every change to a gate statistic is either stricter on any given data or needed because multiple files per project make spectra correlated within clusters. Each change is simulated.
- ce_s is gated and licenses nothing.

**What A4 cannot promise.** Revision 1's realistic walk model (A4.7.2) shows DN3's evaluability depends on two quantities nobody can know before listing:
- the number of usable RAW files per project;
- the stage-2 marker-screen attrition.

**Revision 3, full specification (`a4_r3_walk`).** In the **reference cell**, each project's raw-file count is drawn from DN2's own 21 listing events (median 12):

| Stage-2 attrition | P(v1 minimum) | P(full stop, v2_sep tier 1) |
|---|---|---|
| e = 0.15 | **0.43** | 0.18 |
| e = 0.35 | **0.08** | 0.01 |

**Sensitivity grid** (assumed median files per project):

| Median files per project | P(v1 minimum), e = 0.15 / 0.35 |
|---|---|
| 6 | 0.01–0.08 |
| 12 | 0.22–0.64 |
| 24 | 0.76–0.98 |

**Tier 0 contributes ≈ 0.** Under A4's own same-lab rule only PXD038243 remains eligible (A4.3.10).

**Evaluability: DN3 is more likely than not to be not evaluable.** This is stated plainly, because a not-evaluable DN3 consumes most of the remaining candidate pool (A4.3.3). The decision to proceed is the PI's (A4.11).

---

## A4.1 Principles (inherited; nothing relaxed)
1. **Inherited rules.** A2 P1–P5 and A3 Q1–Q6 apply, reading DN2 as DN3.
2. **Margins unchanged.**
   - tail-G2: n_tail ≥ 200, 0.60, 0.10;
   - G1: ±0.05;
   - TOST rows: ±0.01 two-sided per bin;
   - G2-original: 0.01 / 0.02;
   - G5: 0.55;
   - every validity control (not hulled; widening it would relax).
3. **Gate-statistic changes are stricter on any given data.**
   - The α split (A4.6).
   - The outer envelope over project-unit and cluster-unit intervals (A4.5.2).
   - The intersection–union rule over the capped and uncapped pools (A4.3.6).
   - n_min re-derived under these rules (A4.5.3).
4. **Selection-frame changes** (A4.3) are fixed before any listing and are label-free. Where one could weaken null validity, a compensating control is added: the two-stage marker screen, lineage denies, the cluster cap with the intersection–union rule, and tier-split reporting.
5. **Reuse is limited to DN2's admitted projects** (D26 option (b); A4.3.10). They enter DN3 as tier 0 under every A4 rule.
   - **Why this is compatible with P4:** P4 forbids re-drawing, extending or pruning a null set *after its gate run*. DN2 never reached step 4 or 5, and no per-label DN2 quantity exists (D24). What is known about these projects is label-free: yields, marker counts and file counts. A4 already uses that kind of information for design.
   - **Everything else stays excluded:** DNIT, A1, the T-sets, the other 52 candidates DN2 reached, and everything listed in A4.3.3.
6. **Stricter reading wins.** Any deviation is a dated note written before the step it affects, and, after any T1 look, is reviewed T1-blind (A4.8).

---

## A4.2 Scorers and roles (closed)

| Scorer | Frozen identity | Role on DN3 | Licenses |
|---|---|---|---|
| v1 `static_gbm_pred` | sha256 `31bfe685…` (A1.1) | A2 analysis with the A4.5/A4.6 changes | tail-G2/G1/G3 for v1 (A2.9.6) |
| v2_sep | model, pipeline and code at the A3 freeze commit (`config/a3.yaml` at `prereg-v0.4.3`) | A3 rows V1–V14 and tiers with the A4.5/A4.6 changes | as A3.2.6 |
| ce_s | sha256 `aa19c8f3…` | A3.2.7 | nothing |
| all others | — | not gated | — |

**Precondition.** `prereg-v0.4.3` exists before A4 is tagged. DN3 is valid only for v1 and for the scorers hashed in `config/a3.yaml` at that tag.

**ce_s cost rule (closed).**
- Before step 5, the step-4b agent measures label-free the ce_s scoring wall time on the first 20 mixtures and extrapolates it to all mixtures. This is written as a dated note.
- If the projection exceeds 24 h, ce_s is dropped in full.
- Since ce_s licenses nothing, this cannot steer any claim.

---

## A4.3 DN3 selection

### A4.3.1 Frame (metadata only)

**Source.** A committed snapshot of the 2026-10-26 crawl table:
- `data/cache/pride_meta/all_projects_classified.tsv`, sha256 `30b35340…`, copied to `data/manifests/a4_frame_source.tsv` at the tag;
- the crawl's project records `projects_p*.json.gz`, which supply the text that the protocol screen reads;
- their sha256 are frozen in `config/a4.yaml`;
- nothing is re-crawled.

**Row filter.**
- `exclude_reasons` ∈ {"", "no_quant_or_phospho_evidence"};
- Thermo, DDA, enrichment none, label ∈ {LFQ, UNK};
- one organism.

**Instruments.** Every listed instrument must match one of:
- **A1 class:** A1's regex, unchanged;
- **reserve:** Orbitrap Fusion / Fusion Lumos / Eclipse / Ascend / Exploris 120. Reserve instruments are used only in tier III, and each run must pass A4.3.5(c).

**Why unlabelled "no quant evidence" projects are admitted.** A null spectrum does not need quantification. Labelling and enrichment are excluded by the frame filter, by the keyword screen and per run by the reporter-ion check.

### A4.3.2 Organism class (a priori, by lineage)
Lineages come from UniProt taxonomy per genus. The table is frozen at the tag as `config/a4_taxonomy.tsv` with its sha256. An unresolved genus is excluded.

**Tier I: microbes.**
- Bacteria, Archaea, Fungi, green algae and other eukaryotic microbes.
- **This is broader than A2's 40-genus ALLOW list** (a deny-list design rather than an allow-list), and is stated as such.

**Tier II: land plants (Embryophyta), except *Arabidopsis*.**

**Always excluded:**
- **lineages:** Metazoa, viruses, Apicomplexa, Kinetoplastea, Metamonada, Parabasalia, Entamoebidae, Microsporidia, **Amoebozoa, Ciliophora, Mollicutes**. Amoebozoa and Ciliophora are usually grown in proteose-peptone (animal-derived) or bacterised media. Mollicutes need serum media.
- **genera:** the DENY list (Appendix A), now including *Leptospira*, the dermatophytes and the Mollicutes split genera;
- **species:** *P. aeruginosa*, *A. fumigatus*, the *B. cereus* group, *C. diphtheriae*;
- every A1 null species, the training species (human, mouse, *E. coli*), and *Arabidopsis*.

**Why plants.** Tier I alone cannot reach the v1 minimum in most scenarios (A4.7.2). Land plants share no host material, are grown without serum, and are about as distant from human as fungi. Their risks are:
- handling contamination, covered by certification, the screens and the marker screen;
- multi-organism samples, covered by the A4 contamination keywords.

### A4.3.3 Exclusions (everything prior work touched or could reach)
The exclusion list is frozen as `config/a4_excluded_accessions.tsv` with its sha256. It is the union of:
1. every `PXD\d{6}` string in any **tracked** file under `data/manifests/` or `results/` at the tag, except `data/manifests/t4_frame.tsv` and `data/manifests/a4_frame_source.tsv`. This covers:
   - all manifests;
   - every selection log (DN, DN2, DNIT, T4, the alternative T4 draw in `results/2026-10-27_T1driver_dryrun/t4_selector/`), walk-event files and fetch logs;
   - E9/D8.
2. every accession with a cached PRIDE file listing (`data/cache/pride_meta/files_PXD*`);
3. `config/a2_excluded_accessions.tsv`, and every candidate reached by the A2, DNIT and A1 walks, **except DN2's 6 admitted projects**: PXD006818, PXD018953, PXD044078, PXD011591, PXD062500, PXD038243. They are carved out of items 1–3 and of the file-level guard, and enter only as tier 0 (A4.3.10);
4. `results/2026-10-28_E7_chimeric_v2/train_v2_folds.tsv`;
5. **same-lab projects:** every project (tier 0 included) that shares any person/lab key (A4.3.4) with any A1 null project (`g2null.tsv` and the g2null "tried" list), any DNIT walk candidate, or any DN2 walk candidate other than the six tier-0 accessions. This applies to all such projects, not only the contaminated ones, and closes the shared-instrument carry-over route.
   - It removes 33 pool projects (26 / 5 / 2 in tiers I / II / III; 17 of them protocol passers) and 5 of the 6 tier-0 projects (A4.3.10).
   - Projects sharing a key with an *admitted* tier-0 project are linked into its cluster.

**Before the tag,** an agent permitted to open E7b outputs confirms that the accessions of the E7b cohort and of the 3i/E9 throughput cohorts are in the list, and adds any that are missing. This is a dated note.

**File-level guard.** A RAW is unusable if its PRIDE SHA-1, or its (name, size), matches any file in any committed `*.sha256` or `cache_raw.sha256`.

**After DN3,** its accessions join `select_t4.py`'s exclusions, and every candidate the walk *reaches* is excluded from any future null set.

**Pool at drafting:** 576 exclusion accessions remove 351 frame rows (`a4_r1_walk.json`). Item 5 (person/lab keys, `a4_r3_walk.json`) then removes 33 more, leaving tiers I / II / III = **145 / 253 / 129 projects (60 / 73 / 54 protocol passers)**. The table below is before item 5.

| Tier | Projects | Genera | Pass the rev-1 protocol screen (A2 verbatim lists for comparison) |
|---|---|---|---|
| I | 171 | 92 | 73 (0.43; verbatim 0.19) |
| II | 258 | 48 | 77 (0.30; verbatim 0.17) |
| III | 131 | 55 | 54 (0.41; verbatim 0.17) |

### A4.3.4 Seeded walk order and clusters
**Order.**
- Within a tier: sort by accession, then shuffle with `random.Random(47)`.
- Tiers are walked in the order **0** → I → II → III. A later tier is entered only when the earlier ones are exhausted without a full stop.
- Tier 0 (DN2's admitted projects) is ordered the same way: sorted by accession, then shuffled with `random.Random(47)`.

**Person/lab keys (order-insensitive).**
- Take the crawl `lab` field and every `labPIs` and `submitters` name in the cached project record.
- Lower-case each, split into alphanumeric tokens, drop single-letter tokens, sort the tokens and join them. For example, "Zhang Yongqian" = "Yongqian Zhang", and "Brook L. Nunn" = "Brook Nunn".
- A project's key set is the set of these strings.

**Cluster.** A connected component of *admitted* projects. Two projects are linked if they share a genus or any key.

**Cluster slots.**
- A candidate whose cluster would exceed **5** admitted projects is skipped as "cluster full", before any listing.
- Only admitted projects take slots.

### A4.3.5 Per-project procedure (in order; each step logged)

**(a) Protocol screen (revision 1; replaces "A2 lists verbatim").** Matching is case-insensitive with a **leading word boundary** `\b(?:…)` on every term. The terms split into two groups:

| Group | Terms | Fields searched |
|---|---|---|
| Contamination and multi-organism | `serum`, `FBS`, `FCS`, `fetal`, `blood`, `plasma`, `patient`, `clinical`, `skin`, `wound`, `swab`, `sputum`, `biopsy`, `macrophage`, `cell line`, `HeLa`, `HEK`, `K562`, `brain heart`, `BHI`, `blood agar`, `peptone from (meat|casein)`, `Homo sapiens`, `human (database\|proteome\|sample\|cell\|tissue\|plasma\|serum)`, `co-?culture`, `host`, `infect`, `symbio`, `nodul`, `mycorrhiz`, `endophyt`, `aphid`, `herbivor`, `nematode`, `infest`, `microbiome`, `metagenom`, `allergen`, `proteose peptone`, `meat extract`, `HL-?5`; and `\bIgE\b` **case-sensitive** | **every** field: title, keywords, description, sample protocol, data protocol |
| A1 experiment-type and generic | A1's list (HLA, MHC, glyco, SILAC, DIA, phospho, TMT, iTRAQ, acetyl, N-termin, crosslink, …), plus `tissue`, `spike`, `UPS`, `standard`, `benchmark`, `mixed proteome`, `entrapment`, `human`, `patient` | the **sample-description** fields only: title, keywords, description, sample protocol |

- **SDRF rule:** A2.3 step 2, unchanged.
- **Why this differs from A2's verbatim lists.**
  - The verbatim lists remove about 82% of candidates. Most of the removals come from search-settings boilerplate in `dataProcessingProtocol` ("acetylation (Protein N-term)", "phospho" as a variable modification, "standard parameters") or from organism names ("HLA" inside "C**hla**mydomonas").
  - None of these indicates contamination. Contamination terms are still searched in every field, and contamination is caught again by certification, the identity and family screens and the marker screen.
  - The change is made before any DN3 listing and is label-free (P2). D16's objection applied to changing DN2's list *after* DN2's yield was known.
- **Tests:**
  - "digested with trypsin" passes;
  - "Chlamydomonas" does not trigger HLA;
  - an "acetylation (Protein N-term)" data protocol passes;
  - a "fetal bovine serum" sample protocol is excluded.

**(b) Files.**
- List the project once. A RAW is usable if it is `.raw`, 0.3–1.5 GB, does not match BAD_FILE, and passes the A4.3.3 file guard.
- Order the usable files with `random.Random(47 + accession digits)`. The selected set is the first min(n_usable, **48**).
- A project with ≥ 1 usable file continues.

**(c) Run checks** (label-free, after download, before any human-DB search). A run is unusable, and is replaced by the next file in order, if any of these holds:
- fewer than 99% of its MS2 scans are FTMS HCD;
- TMT/iTRAQ reporter ions are detected (the production pilot detector);
- it has no MS1 scans.

**(d) Two-stage marker screen.** The A2.3 step-3 rule, panels and threshold are unchanged.
- **Stage 1:** the first ≤ 3 selected runs.
- **Stage 2:** if stage 1 passes, all selected runs, pooled.
- The project is excluded if either stage triggers.

**(e) Admission.** A2.2 (a)–(e) per run, unchanged. The P1 swap test runs on the first admitted project and on the cluster cap (A4.3.6).

### A4.3.6 Cluster cap and the capped/uncapped intersection–union rule (label-free)
**Cap factors.**
- For each cluster j, compute water-filling retention factors f_j so that no cluster holds more than **10%** of the capped total.
- Compute them separately for the **tail proxy**, for **A_cand** and for **|A|**, and take the **minimum** of the three.
- These are label-collapsed counts, as in A2.9.3.

**Applying the cap.**
- Cluster j's A spectra are subsampled uniformly: keep a spectrum iff `u(sha256(source run | store index | 47)) < f_j`. FS and A-syn units follow their A donor.
- B is capped the same way, using its own |B| shares.
- If J·0.10 < 1, the cap is undefined and DN3 is **not evaluable**. It never silently means "no cap".

**Two pools.** The *capped* membership and the *uncapped* membership are both committed (A2.9.4 step 3).

**Intersection–union rule.** Every A-based row (tail-G2, G1, G3(A), V1–V3, V12–V14, FS rows, and ce_s's A rows) must **pass on both the capped and the uncapped pool**. B-based rows use the same rule, with B capped and uncapped. The verdict takes the first of these that applies:
1. **invalid** if either pool is invalid;
2. **fail** if either pool fails;
3. **not evaluable** if either pool is not evaluable;
4. **pass** otherwise.

This keeps A2's within-pool precedence, where a validity control can only turn a verdict into "not evaluable".

**Why both pools.**
- The cap alone re-weights clusters. If a leak sits in the largest clusters, capping hides it: the capped-only rule passes 0.25–0.31. Under the intersection–union rule the pass rate is **0.00** (`a4_r3_walk` tail_r3).
- The uncapped pool alone is anticonservative when a leak is concentrated in a few clusters: uncapped-only false pass 0.06–0.13.
- The intersection–union rule is never looser than either pool.

### A4.3.7 Stopping rule (label-free; evaluated on capped counts after each admitted project)

**Full stop.** All of the following hold:
- J ≥ **20** clusters;
- J_eff = (Σx)² / Σx² ≥ **10**, both for the capped tail proxy and for capped A_cand;
- capped tail proxy ≥ **1,200**;
- capped A_cand ≥ **450,000**.

**v1 minimum.** The same criteria, but with capped A_cand ≥ **120,000**.

**Hard caps.** 2,000 admitted runs, or all tiers exhausted. Reaching either ends the walk, and A4.3.8 applies.

| Criterion | Serves | Basis |
|---|---|---|
| J ≥ 20; J_eff ≥ 10 | every clustered row | G1 clean pass is 0.71 at J = 10 vs 0.90 at J = 20 (S3). **J_eff ≥ 10 is essentially implied by the 10% cap** (Σw² ≤ 0.1 when every share is ≤ 0.1; the minimum observed was 10.4). It is kept only as a guard against the min-of-factors exceeding the share. The protection is the cap plus the intersection–union rule; the residual size is stated in A4.10.6. A binding J_eff ≥ 15 was simulated in rev 1 (P(v1 minimum) 0.01–0.61) and is not adopted |
| proxy ≥ 1,200 | V1, V13, v1 tail-G2 | At κ = n_tail / proxy = 0.5, the capped n_tail is ≈ 630, the clean pass is 0.85 and the clean *fail* is 0.06 (`a4_r3_walk`). At proxy 800 the clean fail was 0.19 (rev-1 exploration). A fail is final, so 1,200 is kept |
| A_cand ≥ 450,000 | V12 (v2_sep tier 1) | ≈ 40,000 per label per bin for a whole-row clean pass of 0.86 (S1b) |
| A_cand ≥ 120,000 | G1 and G3 for v1 | G1 n_min (A4.5.3) |

**Thresholds as functions of the final n_min.** At the tag, the A_cand thresholds are set to max(the values above, 10 × final G1 n_min, 10 × n_row(V12) / 0.9), where n_row(V12) is the per-bin n at which the whole 10-bin V12 row passes with probability ≥ 0.85 under the final rule. This is label-free, and it can only raise them.

B rows never bind (B ≈ 40|A|).

### A4.3.8 End of walk without a full stop
- **v1 minimum met:** steps 4, 4b and 5 run. Rows below their own n_min are not evaluable under their own rules.
- **v1 minimum not met:** DN3 is **not evaluable** for every scorer. Steps 4 and 5 do not run. The result is final under P4.
- A project still being processed when the stop condition is met is recorded as "reached, not used". It is not counted, and it is excluded from future sets.

### A4.3.9 Downstream
- **Mixtures:** A2.9.3 unchanged. D0 hosts, 20,000 per mixture, permutation `default_rng(47)`, all *retained* A/B/FS spectra. Both pools are built from the same mixtures: the uncapped pool searched, the capped pool as a membership subset.
- **Step 4b:** A3.2.2 unchanged.
- **Host threshold:** recomputed, and it must equal DN2's 33.875 byte for byte. If not, the run stops, and the cause is investigated T1-blind.

### A4.3.10 Tier 0: DN2's admitted projects (D26 option (b))

**Projects:**

| Accession | Organism | DN2 \|A\| per run | DN2 raw files listed |
|---|---|---|---|
| PXD006818 | *Streptomyces lividans* | 117 | 76 |
| PXD018953 | *Caballeronia mineralivorans* | 1,854 | 12 |
| PXD044078 | *Bacillus pumilus* | 2,165 | 4 |
| PXD011591 | *Rhodobacter capsulatus* | 917 | 24 |
| PXD062500 | *Deinococcus radiodurans* | 3,012 | 8 |
| PXD038243 | *Emiliania huxleyi* | 129 | 9 |

**Every A4 rule applies, with nothing grandfathered from DN2:**
- **Frame and organism.** The frame, organism-class, instrument and DENY rules (A4.3.1–A4.3.2). All six pass on metadata.
- **Protocol screen.** The revision-1 rule (A4.3.5(a)) is evaluated afresh. On the cached crawl text, **PXD011591 fails**: "human" appears in its project description ("… transporters ABCB7 in humans …"). The other five pass.
- **Same-lab rule (A4.3.3 item 5) applies to tier 0.** This is the stricter reading, as A4.1.6 requires. Using order-insensitive person keys (A4.3.4) on the crawl lab plus the PRIDE `labPIs` and `submitters` names:

  | Project | Shares a person/lab key with | Result |
  |---|---|---|
  | PXD018953, PXD044078 | "Jean Armengaud", whose facility also produced five DN2 candidates excluded by the marker screen for human/bovine serum markers (PXD054267/69/72/85, PXD018045) and other DN2 candidates | **excluded** |
  | PXD006818 | "Tassos Economou" (PXD006819, a DN2 candidate) | **excluded** |
  | PXD062500 | "Yongqian Zhang" (PXD027969, a DN2 *D. radiodurans* candidate, listed as "Zhang Yongqian") | **excluded** |
  | PXD011591 | a person key, and the protocol screen | **excluded** |
  | PXD038243 *Emiliania huxleyi* | none | **eligible** |

  **These outcomes are recorded now, before the tag, so that they cannot be revisited.** Tier 0 is therefore {PXD038243}. D26 option (b) is implemented as written, under every A4 rule, but it yields one project. Keeping any of the others would require an explicit, reasoned carve-out of item 5 in a dated PI decision before the tag (A4.11). The draft does not propose one: item 5 exists precisely for a facility with documented serum contamination.
- **Files.** The files are re-listed and re-picked: the first min(n_usable, 48) usable RAWs in `random.Random(47 + digits)` order.
  - The file-level guard is carved out **only** for the project's own DN2 files, i.e. matches in `dn2.sha256` for that accession. Every other sha256 match stays unusable.
  - DN2's downloaded RAWs are reused only if they are picked, and only after a sha256 match.
- **Run checks, two-stage marker screen and admission (A2.2) are recomputed under A4 for every selected run.** That includes the DN2 runs. DN2's admission arrays are **never opened**, because the stores hold T∪D candidate identities. The only comparison allowed is between the sha256 of the recomputed npz and `results/2026-10-27_A2_dn2/cache_admission_stores.sha256`/`membership.sha256`, and it is logged.
- **Caps.** The cluster slots, the cluster cap and the stopping rule apply as for any other project.

**Clusters.** Later *Emiliania* candidates, and projects sharing a person key with PXD038243, link into its cluster.

**Disclosure (required in the adoption note).**
- The design of tier 0 is informed by DN2's label-free per-project yields. They are reproduced in the table above and used in `a4_r1_walk`/`a4_r3_walk`.
- No DN2 label, tail or gate quantity exists or was used.
- Tier 0 is walked first, whatever its yield, so its inclusion cannot be tuned to any outcome.
- Tier 0's order is `random.Random(47).shuffle` of the accession-sorted list. The revision-2 simulation had a constant sort key, which is fixed in `a4_r3_walk`.

**Multiplicity.** DN2 was not a completed test (A4.6), so using its projects spends no α.

**Expected contribution** (`a4_r3_walk`, 300 walks): median 0 admitted (p90 1); P(v1 minimum) changes by ≤ 0.02 (within MC error). Its stage 2 uses an extra exclusion rate of 0.5, because PXD038243 already showed bovine CSN1S1 at 3 runs.

---

## A4.4 Order of operations and roles

**Adoption preconditions, in this order (normative; second review M3).** The A3 freeze hashes `vulcan/`, `bench/`, `scripts/a3/`, `scripts/e7/` and `scripts/e7b/` (`a3_config.freeze_status`). Any A4 code in those trees must therefore land **before** the A3 freeze commit, and nothing may land in them afterwards until DN3 step 5 is committed and the T1 escrow run has finished. There is no fallback route.
1. A4 (this revision or later) is adopted. Its review(s) and A4.12–A4.14 are committed.
2. The A4 code in frozen trees is committed and tested: `vulcan/gates/a4_stats.py`, the `bench/` escrow code, `bench/t1/claims_a4.py`, and the stratum-1 outputs of `bench/t1/a3.py`. **The A3.1.6 dry run** (v1 columns bitwise identical with and without the A3 additions) **is re-done** after these changes.
3. **The A3 freeze commit**, then `prereg-v0.4.3`.
4. The A4 files outside frozen trees (`scripts/a4/`, `config/a4*`), then the tag `prereg-v0.4.4`. Steps 3–8 below are part of this.

Any other planned change to frozen trees (e.g. the D27 GPU MS1 merge) waits until DN3 step 5 and the escrow run are committed, or is done on a branch that the freeze does not see.

Each step of that order, in detail:
1. `prereg-v0.4.3` exists (step 3 of the order).
2. This text, its review and A4.12 are committed. The PI adopts the text as DESIGN §6.0d with a dated note that states the adopter's exposure to T1 (expected: none) and to DN2 label-free records.
3. **Implementation** is committed and tested on synthetic data:
   - `config/a4.yaml` (every constant in Appendix B);
   - the frozen frame, records, taxonomy and exclusion files with sha256;
   - `scripts/a4/select_dn3.py`;
   - `vulcan/gates/a4_stats.py` (outer envelopes, intersection–union; `a2_stats` and `a3_stats` unchanged);
   - `scripts/a4/gates_a4.py`.
   - **Tests:**
     - a dry-run walk on a synthetic frame;
     - the A4.3.5(a) regressions;
     - P1 for the cap;
     - no per-label counts in `scripts/a4/`;
     - every A4 rule on synthetic pools, including both pools;
     - the closure guard. It freezes v1's scorer, feature and search closure against `prereg-v0.4.2`, and v1's *statistics* against `vulcan/gates/a4_stats.py` and `scripts/a4/gates_a4.py` at the A4 tag.
4. **n_min final values** (A4.5.3) are re-derived with ≥ 4,000 replicates and the rule p̂ ≥ 0.92, with between-cluster heterogeneity of the base rate (logit-sd 0.5) in the clean scenario. The **G3 exhaustive-search cost** is measured on dev mixtures. Both are committed in `config/a4.yaml`.
5. **Escrow implementation (A4.8)** is committed and tested:
   - the `sealed_run.py` one-look rule and `--escrow` mode;
   - `bench/escrow.py open`;
   - `bench/t1/claims_a4.py`;
   - the stratum-1-only outputs in `bench/t1/a3.py`;
   - the D0/D1 escrow dry run with its no-metric console check.

   It lands before the A3 freeze commit (ordering above).
6. **Disk.** The Lustre quota is checked, and the RAW policy is fixed: admitted RAWs are **kept** under `data/cache/a4/raw` (≈ 0.9–2 TB), plus our own sha256. If the quota is insufficient, the tag waits.
7. Tag `prereg-v0.4.4`.
8. **After the tag:** the T1 escrow run (A4.8), by an operator who takes no DN3 role. It runs on an **otherwise idle node**, not overlapping the walk's admission compute, so that T1 speed stays claims-grade. If it does overlap, T1 speed is labelled "not claims-grade".

**Steps (each committed before the next):**
- **2–3:** walk, run checks, marker screens, admission, cap, stopping rule. Commit the selection log, `dn3.tsv`, sha256 and PRIDE SHA-1, both membership hashes, label-free counts and swap tests.
- **4:** mixtures.
- **4b:** by an agent that does not run step 5.
- **A3.6 step 3:** an end-to-end label-free run plus a synthetic-label dry run.
- **5:** one session, by an agent that ran none of steps 2–4b, wrote none of A2–A4, and is T1-blind. v1 first (commit), then v2_sep and ce_s (commit).
- **6:** NOTES and the RESEARCH_TREE decision.

---

## A4.5 Gate statistics

### A4.5.1 Reused unchanged
- **v1:** A2.4–A2.8 (pools, units, (c-final), validity controls, FS rule, G2-original/G4/G5).
- **v2_sep:** A3.2.3–A3.2.6 (τ_s, tail membership, merge maps, A-syn, V1–V14, tiers, production rule).
- **ce_s:** A3.2.7.

### A4.5.2 Change 1: the outer envelope over project-unit and cluster-unit intervals, for every bound-based row
A coarser clustering is stricter only in expectation (review M1). A4 therefore computes each row's interval with **both** units and takes the outer envelope with the A2/A3 analytic interval. The cluster terms per row type:

| Row type | Rule |
|---|---|
| **tail-G2** (V1, V13, v1, ce_s) | UCB = max(CP or delta bound, project bootstrap, A4-cluster bootstrap). B = 10,000, `default_rng(83)` |
| **G1** (V2, V14) | Fisher-z CI (Bonett–Wright) ∪ project jackknife-t CI ∪ cluster jackknife-t CI (leave-one-unit-out Fisher z, t_{J−1}) ⊂ (−0.05, 0.05). The hyperscore control is unchanged (99% Fisher CI ∋ 0) |
| **G3, V4, V5, V6, V12** | Newcombe ∪ project bootstrap ∪ cluster bootstrap ⊂ [−0.01, +0.01] |
| **G2-original** (V7, V8) | max(§6.4 UCB, project bootstrap, cluster bootstrap) |
| **G5** | point rule, unchanged |

**Simulated** in `a4_power` S2/S3, with heterogeneity at the margin:
- the A2-form rules pass up to 0.11 (TOST) and 0.05 (G1);
- the A4 rule passes ≤ 0.038 and ≤ 0.007.

### A4.5.3 Change 2: n_min under the A4 rule
These are provisional values from rev 0 (J = 20). The final values come from adoption precondition 4 and can only be equal or larger.

| Row | A2/A3 | A4 (provisional) |
|---|---|---|
| G1 n_top | 7,200 | 11,100 |
| V6 per label per bin | 22,700 | 30,600 |
| V12 per label per bin | 19,200 | 26,800 |
| G3, V3–V5 | 2,000 | 2,000 (one-bin clean pass 0.999 at recall 0.998) |
| tail-G2 n_tail | 200 | 200 |
| V8 | 107,570 | 107,570 |

Also: **fewer than 10 clusters ⇒ every clustered row is not evaluable.**

### A4.5.4 Change 3: the α split (A4.6)

### A4.5.5 A1-class check (gating, stricter)
If the A1-instrument subset has n_tail ≥ 200, its tail-G2 **point estimates** must not exceed the margins (target fraction ≤ 0.60; difference ≤ 0.10). This is evaluated on **both** pools (intersection–union), and an exceedance is a fail. It guards against a QE-specific leak diluted by tier III.

### A4.5.6 Reported, never used to pass
- Everything A2/A3 report.
- Per tier, per instrument class, per cluster.
- Leave-one-cluster-out for the tail rows.
- Each interval component.
- Both pools.
- The cap factors.
- The within-cluster pooled ρ (`a2_stats.fisher_pooled_rho`) next to G1, which shows organism-composition effects.

---

## A4.6 Multiplicity across attempts

**What counts as a completed test.** An attempt for a scorer is a *completed test* iff step 5 ran for that scorer, so that some per-label quantity was computed, whatever the verdict. A label-free stop before step 5, as in A2, is not a test.

**Ledger** (rows A1 failed or never validly tested: tail-G2, G1, G3, and the A3 V-rows built from them):

| Scorer | Completed tests before DN3 | DN3 is |
|---|---|---|
| v1 | 1 (A1: fail). A2 was not a test | test 2 |
| v2_sep | none | test 1 |
| ce_s | none | test 1 |

**Rule (geometric spending, fixed for all future attempts).** Completed test k uses level α/2^(k−1) on every re-tested row. For DN3:
- tail-G2: one-sided 0.975;
- G1: two-sided 0.995;
- TOST rows: two-sided 0.95.

The same levels apply to all scorers, because V3 *is* v1's G3. G2-original and G5 keep their levels (v1 passed G2-original; re-testing a passed gate can only fail).

**What can be claimed (corrected from rev 0).**
- **Conditional on A1's failure,** DN3's nominal per-row level is α/2. That is the statement used in reports.
- **The unconditional bound over attempts is not 1.5α.** For G3, A1's rule was a point rule with no controlled level: it passed about 45% at a 1-pp deficit (A2.8). A1's tail-G2 had no cluster term.
- **Simulated per-attempt false-pass rates of DN3's tail-G2** at the margin, under the realistic walk (A4.7.3):
  - ≤ 0.01 with heterogeneity τ = 0.10;
  - 0.023–0.057 when the leak is concentrated in 20–25% of clusters (`a4_r3_walk` tail_r3; MC SE ≈ 0.012). The upper end exceeds the nominal 0.025 (A4.10.6);
  - against a nominal 0.025, with MC SE ≈ 0.01.
- These figures are reported. No "procedure-level 1.5α" sentence is used.

**Why not α for DN3, or α = 0?**
- Keeping α makes test-until-pass cheaper.
- α = 0 is impossible, because A1 has already spent α.
- Ignoring A1 (because its nulls were contaminated) would be the post-hoc argument D12 forbids.

**Next attempt.** After a DN3 that completed step 5 without a pass: α/4. After a label-free not-evaluable DN3: α/2 again.

**Across scorers and within scorers:** unchanged (A3.2.5).

---

## A4.7 Power and feasibility

### A4.7.1 Inputs from DN2 (label-free; `a4_power` S0)

| Quantity | Value | Use |
|---|---|---|
| \|A\| per run | mean 1,366, median 1,276 (18 runs, one with 0) | log-normal; between-project log-sd 1.35, within 0.32, zero-run fraction 0.056 |
| A_cand / A | 0.391 (per project 0.16–0.74) | logit-normal |
| tail proxy / A | 0.000895 (per project 0.0004–0.0085) | log-normal, sd ≈ 1 |
| B / A | 39.9 | B never binds |
| effective projects (A, proxy, A_cand) | 3.7, 3.6, 4.0 of 6 | motivates the cap and J_eff |
| reference proteome found | 16 / 19 | walk |
| marker screen, stage 1 | 6 / 12 excluded | p = 0.5; extra stage-2 exclusion e ∈ {0.15, 0.35} for R > 3, costing R downloads |
| RAW size | 1.04 GB mean | volume |

The **label ICC** cannot be estimated without labels. It is assumed through the τ grid. Under a clean scorer on a valid null, every cluster's tail target fraction is 0.5 by exchangeability, so τ matters only under the alternative.

### A4.7.2 Evaluability (full specification, `a4_r3_walk`, 300 walks per cell)
**What the walk implements:**
- the real pool after item 5;
- deterministic protocol flags;
- order-insensitive person/lab keys and genus links (≤ 5 projects per cluster);
- the proteome step (16/19);
- ≤ 48 files per project;
- the 10% cap on min(proxy, A_cand, |A|) factors and the A4.3.7 criteria;
- stage-1 marker pass 0.5, with an extra stage-2 exclusion e (cost R downloads);
- tier 0 = {PXD038243}, with stage-2 exclusion 0.5.

**Reference cell.** Each project's raw-file count is drawn from DN2's 21 label-free listing events (p10/p25/p50/p75/p90 = 4/6/12/16/36), with usable files ~ Binomial(n, 0.8).

| files per project | e | P(v1 minimum) | P(full stop) | runs at end (p10–p90) | J | admitted | files downloaded (p90) |
|---|---|---|---|---|---|---|---|
| **empirical (reference)** | 0.15 | **0.43** | **0.18** | 713 (598–852) | 35 | 59 | 1,042 (1,168) |
| **empirical (reference)** | 0.35 | **0.08** | **0.01** | 573 (448–691) | 31 | 49 | 1,087 (1,234) |
| median 6 | 0.15 / 0.35 | 0.08 / 0.01 | 0.01 / 0.01 | 477 / 390 | 36 / 32 | 60 / 52 | 741 / 766 |
| median 12 | 0.15 / 0.35 | 0.64 / 0.22 | 0.37 / 0.07 | 864 / 713 | 34 / 30 | 56 / 47 | 1,200 / 1,293 |
| median 24 | 0.15 / 0.35 | 0.98 / 0.76 | 0.93 / 0.59 | 1,186 / 1,102 | 29 / 27 | 45 / 42 | 1,536 / 1,890 (2,209) |

**Without tier 0:** P(v1 minimum) is 0.41 / 0.08 in the reference cell, and within ±0.02 in every cell. Tier 0 contributes ≈ 0.

**What e should be.** 3 of DN2's 6 admitted projects already showed one non-keratin panel member at 3 runs, and stage 2 pools up to 48 runs. So e = 0.35 is at least as plausible as 0.15.

**Earlier figures are superseded:**
- The revision-1/2 figures (`a4_r1_walk`: 0.08–1.00) omitted item 5, the person links and the |A| factor, used a non-shuffled tier 0, and assumed the file distribution. They are kept for the record only.
- Alternatives simulated in revision 1, not adopted: proxy ≥ 800 (higher clean final-fail risk, 0.19 at κ = 0.5); J_eff ≥ 15; ≤ 8 per cluster or ≤ 96 files (≤ +0.06).

**Reading.**
- In the cell anchored on DN2's own file counts, DN3 reaches the v1 minimum with probability ≈ 0.1–0.45, and v2_sep's full stop with ≈ 0.01–0.2.
- The binding constraints are the pool's cluster structure under the same-lab rule and the 10% cap that validity requires. File caps and organisms are not binding.

### A4.7.3 Gate power and size (`a4_r3_walk` tail_r3)
Tail-G2 is evaluated at the v1-minimum states of the reference cell (e = 0.35, tier 0 included): median J = 31, capped J_eff ≥ 10.4, maximum capped share 0.14. Results use the IUT precedence of A4.3.6. There are 300 states per column; MC SE ≈ 0.012 near 0.04 and ≈ 0.02 near 0.9.

| scenario: verdict rate (intersection–union, A4 levels) | κ = 1 (n_tail ≈ 1,260) | κ = 0.5 (n_tail ≈ 630) |
|---|---|---|
| clean: pass | **0.91** | **0.85** |
| clean: invalid pool / fail | 0.09 / 0.00 | 0.08 / 0.06 |
| margin 0.60, τ = 0.10: pass | 0.003 | 0.010 |
| mean 0.60, leak 0.90 in 25% of clusters: pass (uncapped alone) | 0.033 (0.083) | 0.023 (0.057) |
| mean 0.60, leak 1.0 in 20% of clusters: pass (uncapped alone) | **0.057** (0.13) | 0.040 (0.11) |
| uncapped mean 0.60, leak in the largest clusters: pass (capped alone) | 0.000 (0.31) | 0.000 (0.25) |

The clean "invalid pool" route (≈ 0.08–0.09) is an additional way to end up not evaluable (A2.11.4). It is not a fail.

**Joint clean pass, conditional on DN3 being evaluable** (rows from `a4_power` S6 with the tail figures above; approximate):

| Scorer / tier | At | Joint clean pass |
|---|---|---|
| v1 | v1 minimum | ≈ 0.75–0.85 |
| v2_sep tier 1 | full stop | ≈ 0.75–0.80 (V12 binding) |
| v2_sep tier 2 | full stop | ≈ 0.65–0.75 |

These are recomputed with the final n_min values (adoption precondition 4). Multiplied by the evaluability above, the **unconditional** chance that a clean v1 passes on DN3 is about 0.06–0.4 in the reference cell.

### A4.7.4 Feasibility (re-measured in the implementation dry run)

| Item | Estimate |
|---|---|
| Files downloaded | ≈ 740–1,890 (median by scenario; p90 ≤ 2,210); reference cell ≈ 1,040–1,090; including stage-1 and stage-2 exclusions |
| Volume | ≈ 0.85–2.0 TB |
| Download time | DN2's per-stream rate was ≈ 3 MB/s, so 8 ordered streams take ≈ 10 h per 0.85 TB (≈ 23 h at 2 TB). If E9's 0.1–0.2 GB/s aggregate is reachable, 1.2–5.5 h |
| Admission compute | ≈ 11,300 certified spectra per run × ≈ 25 ms ≈ 70–90 process-hours. With 16 deterministic workers: ≈ 5 h |
| Disk | RAWs stream through `/tmp` (≈ 85 GB free, shared). Stores ≈ 0.13 GB per run on Lustre (≈ 0.1–0.15 TB). Admitted RAWs kept (≈ 0.9–1.5 TB, quota checked before the tag) |
| Mixtures | (\|A\| + \|B\| + \|FS\|) ≈ 45–65 M / 20,000 ≈ 2,300–3,300, each with a D0 host |
| Step 4 search and v1 features | ≤ 1–2 GPU-h |
| Step 4b | MS1 ≈ 6–9 core-h; chimeric v2 ≈ 13 s per mixture ≈ 8–12 h single stream (parallelisable) |
| G3 exhaustive | unmeasured; measured before the tag (adoption precondition 4). If it is the bottleneck, any B subsample must be a separate dated amendment before the tag |
| ce_s | ≈ 1 GPU-h |
| T1 escrow run (A4.8) | one `sealed_run.py` invocation of the frozen T1 driver; hours, and it can overlap the walk |
| Total wall time | walk ≈ 1–2 days, steps 4/4b ≈ 1 day, step 5 hours: **≈ 3–4 days** |

**Engineering.** The walk prefetches and admits in parallel, but decisions are taken strictly in walk order. A crash restart is a resumption (DN2 DEVIATIONS 3).

---

## A4.8 T1 in escrow (normative; D26 decision 3, built on `reviews/a4_t1_sequencing_review.md`)
D12's block on the T1 look is amended **only in this escrow sense**:
- T1 may be *computed* before DN3 step 5;
- **no one reads any T1 output** (the PI, any agent or any person) until an unseal trigger (A4.8.4).

Every claim condition of D12, A2.9.6, A3.3 and A4.9 is unchanged.

### A4.8.0 Replacement text (adopted with A4; dated notes, not reinterpretations)
- **D12:** "FDR claims and the T1 look stay blocked until gates pass" → "FDR claims stay blocked until v1 passes the null gates. T1 may be computed only under the escrow of DESIGN §6.0d A4.8, and is read only at an A4.8.4 trigger."
- **A3.3.1, first bullet, and A3.3.7 row 1:** "If v1 fails or is not evaluable, there is no T1 look under A3" → "If v1 fails or is not evaluable on DN3, there is no v1 or v2_sep claim. T1 is reported in full under A4.8.5."
- **A3.3.2 and A3.3.8:** the single look also computes the stratum-1-list inputs of A4.8.1.4.
- **Consequence, stated once:** unsealing spends T1 as a confirmatory set for every scorer, present or future (e.g. a v3), whatever DN3 shows. This reverses A3.3.1's rationale of not spending T1 on a scorer without validated FDR, and the PI accepts this with D26.

### A4.8.1 Preconditions for the escrow run (enforced in code)
1. **Tags.** `prereg-v0.4.3` and `prereg-v0.4.4` exist, and HEAD descends from both.
2. **Clean tree.** The tree is clean, and the A3 freeze commit and tree hashes match (A3.1.3). The A4 order of A4.4 guarantees the escrow code is inside the freeze.
3. **Committed and hashed before the run** (hashes recorded in `config/a4.yaml`): the T1 driver, the escrow tools, and the **claim script** `bench/t1/claims_a4.py`.
   - The claim script is a pure function of (DN3 verdicts per scorer and tier; the A3.1.8 dev-precondition verdicts of D23, i.e. G4, arm M and the §6.7 MC; the escrowed T1 outputs) to the A3.3.7/A4.9 claim table and report.
   - It is reviewed T1-blind before the escrow run.
4. **Every branch in one run.** The driver computes everything any DN3 outcome could need:
   - the v1 primary and every A1/A3 comparator;
   - v2_sep's family-S inputs for the full list and for the **stratum-1-only list**;
   - the leave-one-run-out backstops;
   - every entrapment verdict (arms P/C, per merged stratum, protein groups);
   - speed and quant.

   **Stratum-1-only list (defined now; second review M6):**
   - The list is the v2_sep accepted units of merged stratum 1 at a common `q_sep` cutoff.
   - Matched arm-P FDP is computed on that list's own paired counts (target, entrapment and their decoys restricted to stratum-1 pair winners), with `evaluate.py`'s first-crossing rule.
   - Ratios, bootstrap and the leave-one-run-out minimum are as in A3.3.4, against every A3.3.3 comparator.

   A D0/D1 escrow-mode dry run shows that every output is produced and that no file appears outside the escrow root (A4.8.2.1). It also shows that every console line matches the **whitelist of line templates** (stage names, job counters, the exit code). This dry run is committed.
5. **Code location.** See A4.4: all of this lands before the A3 freeze commit. There is no other route.

### A4.8.2 One look, and the escrow mechanics
**One look.** `bench/sealed_run.py` gains a mandatory one-look rule for T1: it refuses to start if `results/SEALED_LEDGER.md` already holds a T1 line, except an escrow retry under item 5 below. A test covers this.

**Tag rule in escrow mode.** Ledger-only commits and escrow-manifest-only commits (`results/T1_escrow/`) are skipped when `sealed_run.py` looks for the tag, so retries keep the same base commit.

**Escrow mode** (`--escrow`):
1. **Escrow root.** The root is a node-local directory **outside the repository**, `/tmp/vulcan_escrow_T1_<ledger-id>/`, created mode 0700. The driver is called with `--out <root>/out` **and** `--work <root>/work`. Every cache it writes (predcache, predcache_a3, indexes, PSM tables, mokapot, chimeric outputs) lives under the root. Stdout and stderr go to `<root>/logs`. The terminal and the agent see only whitelisted heartbeat lines and the exit code.
2. **Sealing.** At exit, whatever the exit code:
   - the whole root is tarred;
   - the tar is **encrypted to a public key whose private key is held by the PI off the cluster** (`age` or `gpg`; the public key's fingerprint is recorded in `config/a4.yaml` at the tag);
   - the plaintext root is deleted, and a check confirms it is gone;
   - the ciphertext goes to `data/escrow/T1/<ledger-id>.tar.age`;
   - the ciphertext's sha256 and its total size only (no file names or sizes) go to `results/T1_escrow/MANIFEST.md`, which is committed;
   - **the T1 wall time and the exit code's detail go inside the escrow.** The ledger "result" line records only "sealed" or "failed".
3. **Leak test.** After the escrow-mode dry run on D0/D1, a test asserts that no file newer than the run start exists outside the escrow root, under `/tmp/vulcan_t1_*`, `results/`, `data/cache/` or `work/`. The only exceptions are the ledger and the manifest.
4. **Failure.** A non-zero exit produces a committed failure summary that contains only the stage name, the exception type and the traceback's file:line:function frames, with no messages and no locals. Debugging is done **only** by reproducing the failure on D0/D1.
5. **Retries.**
   - At most two, each a new ledger line "T1-escrow-retry n".
   - A retry may follow only an engineering fix that touches no T1 analysis parameter. The fix is reviewed by a T1-blind agent. It must re-pass the A3.1.6 bit-identity dry run, and it changes the hashes in `config/a4.yaml` by a dated deviation.
   - If the fix touches a frozen tree, the A3 freeze must be re-taken (new freeze commit and new `prereg-v0.4.3.x` tag, with a dated note) **before** the retry.
   - Every earlier ciphertext stays sealed and is opened with the final one. **The last successful run is the result of record.**
6. **Opening.** `bench/escrow.py open` is the only route. It:
   - verifies the ciphertext sha256;
   - refuses unless an A4.8.4 trigger commit exists;
   - appends an "unsealed" ledger line naming the trigger commit;
   - only then asks the PI to decrypt, and extracts to `results/<date>_T1/`.

### A4.8.3 Access control and log
- **What is logged.** Every action on the escrow (create, seal, retry, open) is a ledger line. Any other read of the root during the run, or any decryption outside `escrow.py open`, is a protocol violation and is reported with T1.
- **What the cryptography covers.** Encryption makes the sealed ciphertext unreadable by any agent on the cluster. During the run (hours), the plaintext root is readable by same-uid processes. The run therefore happens on a node and at a time when no other agent of this project is active (A4.4 step 8), and no other agent may access `/tmp/vulcan_escrow_T1_*`.
- **Who runs it.** The escrow run's operator runs no DN3 step and reviews no A4 deviation.
- **Who must stay T1-blind until the unseal:** every A4 artefact author (notes, deviations, reviews, adoption notes, step-5 NOTES), and the **author, reviewer, adopter, walk agent and step-4b/5 agents of any later null-set amendment**. Each such artefact carries a "T1 exposure: none" statement.
- **During escrow,** nobody computes anything from escrowed data, and no re-runs are done to "check" results.

### A4.8.4 Unseal trigger (whichever comes first)
- **(a) DN3 step 5 complete.** The commit of DN3 step-5 results for all gated scorers (A4.4 step 5).
  - A DN3 that ends **not evaluable does not unseal T1.** The escrow then stays sealed until the step-5 commit of the next null-set attempt for v1, or until trigger (b).
  - This follows D26's wording ("until A4 step 5") and the sequencing opinion's safeguards 1 and 3.
- **(b) Time limit.** It fires on a **fixed calendar date** for the project's final-report phase. The PI states this date in the A4 adoption note and may not move it later. Trigger (b) may not fire **before DN3 steps 2–3 are committed**.
  - When it fires, the PI writes a dated RESEARCH_TREE decision before unsealing, recording the status of A4 (and of any later null-set attempt) at that moment.
  - The trigger cites no T1 information; none exists.

If the unseal happens under (b) while A4 is unfinished:
- **A4 runs to completion.** Stopping it, or a later null-set attempt, other than by its own label-free rules is a protocol violation, reported with T1.
- The A4 agents still to act (4b, 5, deviation reviewers) must be T1-blind.
- The A4 rules cannot change.
- A4 is reported whatever T1 showed.

### A4.8.5 Reporting after the unseal
- **T1 is reported in full**, computed by the pre-hashed claim script. No output is withheld.
- **Labels while DN3 step 5 is incomplete:** "not gate-validated; A4 pending; v1 failed tail-G2/G1/G3 under A1 (D12); A2 not evaluable (D24)".
- **After DN3, the A3.3.7/A4.9 claim table applies:**
  - if v1 fails or is not evaluable, the label is "FDR not gate-validated", with no FDR, ID-at-FDR or superiority claim for any scorer;
  - v2_sep claims require its tier.
- **No production default and no "validated" wording** before DN3 step 5.
- **No external release** of T1 numbers before the unseal. After the unseal and before DN3 step 5, any release carries the labels above.
- A4 is reported whatever T1 showed. T1 is reported whatever A4 shows.

### A4.8.6 Why escrow is sound
- Under escrow, nobody's information set contains T1 before DN3 step 5, or before the fixed date. So no A4 design, execution, deviation or persistence decision, and no later null-set design, can depend on T1.
- The T1 numbers are a fixed function of frozen code and sealed data.
- **Residual visible information** is only:
  - success or failure;
  - the failure stage and code frames;
  - the ciphertext size, which is compressed and encrypted (the residual leak is small, accepted, and stated).

---

## A4.9 Reporting commitment
`results/<date>_A4_gates/NOTES.md` is written whatever the outcome. It contains:
- one section per scorer, v1 first;
- every row on both pools, on A, B, FS, A-syn and the A1 offsets;
- per-cluster, per-tier and per-instrument tables;
- both interval components;
- the cap factors and the validity controls.

**Required wording:**
- **v1 pass:** "passed tail-G2/G1/G3 on fresh nulls DN3 under A4, at its second completed test (A1: failed on contaminated nulls, D12; A2: not evaluable before any test, D24); per-row nominal level α/2 conditional on the A1 failure; simulated per-attempt false-pass rates at the margins ≤ 0.01 under diffuse heterogeneity and up to ≈ 0.06 under leaks concentrated in a fifth to a quarter of clusters (nominal 0.025)".
- **v2_sep:** A3.3.7's wording, with "A3/A4 (DN3), first completed test".
- **A fail** is final under P4.
- **Not evaluable** keeps claims blocked. The next attempt follows A4.6.

---

## A4.10 Residual risks
1. **Evaluability is poor in the reference cell** (A4.7.2): P(v1 minimum) ≈ 0.08–0.43. If DN3 is not evaluable, the remaining PRIDE pool for this organism class is largely consumed, and T1 stays sealed until the next attempt's step 5 or the fixed date (A4.8.4).
2. **Organism and instrument heterogeneity.** A QE- or microbe-specific leak can be diluted by tiers II and III. This is mitigated only partly, by the A1-class point check and by reporting.
3. **B contamination can produce a final *fail*** of G2-original or G5 (A ∪ B is about 97% B). This is A2's rule, and it is kept.
4. **The validity controls at large n** can flag a residual hyperscore bias of about 0.03 as an invalid pool.
5. **The label ICC, κ and the tail rate of new organisms are assumed.** At κ = 0.5, a clean scorer fails tail-G2 (finally) about 6% of the time, and is "invalid pool" about 8% (A4.7.3). This is under the pessimistic independent-tails model. The real learned and hyperscore tails overlap, so the true rate should be lower.
6. **Concentrated leaks.** A margin-level leak concentrated in 20–25% of clusters passes the intersection–union rule 0.023–0.057 of the time (MC SE ≈ 0.012), against the nominal 0.025. The cause is few effective clusters (capped J_eff ≈ 10–14).
   - A binding J_eff ≥ 15 would fix this but makes DN3 almost certainly not evaluable (rev 1).
   - **The PI accepts this residual in the adoption note** (option (a) of second review M5), and every v1/v2_sep report quotes it (A4.9).
7. **v2_sep tier 1's joint clean pass is ≈ 0.75–0.80 at the full stop.** A chance not-evaluable or fail outcome is possible.
8. **The cluster definition from genus and person/lab keys** can miss shared facilities, and person-name matching can give false links. The order-insensitive keys reduce the second problem.
11. **Escrow residuals** (A4.8.6). The plaintext is readable by same-uid processes during the run; this is mitigated by an idle node and the access rule. The ciphertext size and the success/failure status are visible.
12. **The A3 freeze ordering** (A4.4) blocks other changes to `vulcan/` and `bench/` (e.g. D27) until DN3 step 5 and the escrow run are complete, unless they are done on a branch.
9. **New parallel-walk code.** It adds engineering risk.
10. **Seed predictability.** The walk order is computable from the frozen frame and seed 47. The adoption note states that nobody computed it while the rules were drafted. The author did not; the walk simulations used other order seeds.

---

## A4.11 Points for the PI (revision 3): decisions needed before adoption
1. **Tier 0 vs the same-lab rule (second review B1).**
   - Revision 3 applies item 5 to tier 0, the strict reading. Only **PXD038243** remains; PXD018953 and PXD044078 (the Armengaud facility, whose five other DN2 candidates were excluded for serum markers), PXD006818 and PXD062500 drop out. D26 option (b) thus yields one project and ≈ 0 evaluability gain.
   - Keeping any of the others needs an **explicit, reasoned carve-out** in a dated PI decision before the tag. The draft recommends against it for the Armengaud pair.
2. **Not-evaluable DN3 keeps T1 sealed (B2).** Revision 3 follows D26's wording: unseal at DN3 step 5 or the fixed date, not at "not evaluable". If the PI wants unsealing on not-evaluable, that is a new decision, and it accepts that the next null-set amendment is designed by T1-exposed people. That amendment would then carry the disclosure.
3. **A fixed calendar date for the time-limit trigger (M2).** It goes in the adoption note, and (b) cannot fire before DN3 steps 2–3 are committed.
4. **An escrow key (M1).** The PI holds a private key off the cluster. Its public-key fingerprint goes in `config/a4.yaml`.
5. **Accept the concentrated-leak residual (M5):** ≈ 0.06 against a nominal 0.025 (A4.10.6).
6. **Proceed or not, given the evaluability (M4).**
   - In the reference cell, P(v1 minimum) is 0.08–0.43 and P(full stop) 0.01–0.18. A clean v1's unconditional chance of passing on DN3 is ≈ 0.06–0.4.
   - The cost is about 0.75–1.9 TB of downloads and 3–4 days, and a not-evaluable result consumes most of the remaining pool.
   - The levers that would raise evaluability all weaken validity or were rejected (proxy 800; no same-lab rule; larger cap share). The draft makes no recommendation to change them.
7. **Freeze ordering (M3).** A4's code in `vulcan/` and `bench/` lands before the A3 freeze. D27's GPU MS1 merge must wait until DN3 step 5 and the escrow run are complete, or be done on a branch.

## A4.12 Changes after the independent review (`reviews/amendment_A4_review.md`)

| Review item | Revision 1 |
|---|---|
| **B1** `IgE` matched "digest" | Leading `\b` on every term; `\bIgE\b` case-sensitive; regression tests (A4.3.5(a)) |
| **B2** optimistic walk | New `a4_r1_walk.py`: real pool, deterministic protocol flags, real clusters, proteome step, stage-2 costs. Protocol screen changed (field split, `\b`); ≤ 5 projects per cluster; ≤ 48 files. Evaluability stated by cell, including the bad cells (A4.0, A4.7.2, A4.11.1) |
| **B3** concentrated leaks; cap re-weights; 1.5α | 10% cap on the min over proxy, A_cand and \|A\| factors; intersection–union over the capped and uncapped pools for every row; simulated (A4.3.6, A4.7.3). The 1.5α wording is withdrawn (A4.6, A4.9) |
| **M1** coarser unit not monotone | Outer envelope over the project and cluster units for every bound row (A4.5.2) |
| **M2** cap only on the proxy; silent no-cap | Cap on the min of the three factors (and \|B\| for B); J_eff ≥ 10 on both capped quantities in the stopping rule; J·share < 1 ⇒ not evaluable (A4.3.6, A4.3.7). The reviewer's J_eff ≥ 15 was simulated and costs too much evaluability (A4.7.2); 10 is used, with the measured false-pass rates reported |
| **M3** multiplicity | "Completed test" = step 5 ran; conditional nominal statement; the G3 caveat; next-attempt rule made consistent (A4.6) |
| **M4** broad tier I; peptone media | Amoebozoa, Ciliophora and Mollicutes denied by lineage; *Leptospira*, dermatophytes and Mollicutes split genera denied; tier I described as broader than A2's class (A4.3.2) |
| **M5** exclusion leaks (PXD002407) | Every PXD in tracked `data/manifests/` and `results/`, plus cached listings, plus same-lab projects; frame files named; E7b and 3i cohort check before the tag (A4.3.3) |
| **M6** T1 safeguards by reference | Normative in A4.8 |
| **M7** open options | DN2 reuse closed (no); RAW policy fixed (keep; quota checked before the tag); lab normalisation specified, including PI and submitter names; hard cap → A4.3.8; "probably" removed; ce_s rule specified (A4.2) |
| minor 1 counts | Tables rebuilt from `a4_r1_walk` |
| minor 2 cross-reference | Residual 6 added |
| minor 3/4 n_min noise and heterogeneity | Adoption precondition 4 (≥ 4,000 replicates, p̂ ≥ 0.92, logit-sd 0.5) |
| minor 5 G2 power without the cluster term | Noted: immaterial at B ≈ 33–60 M |
| minor 6 G1 composition | Within-cluster pooled ρ reported (A4.5.5) |
| minor 7 seed predictability | Adoption-note statement (A4.10.10) |
| minor 8 download time; B compute | Corrected (≈ 10–23 h at 8 streams); G3 cost measured before the tag; a B subsample only by a separate pre-tag amendment |
| minor 9 closure guard | Split: v1 scorer/search closure vs A4 statistics (A4.4.3) |
| minor 10 G1 at J = 10 | 0.71 |
| minor 11 same-lab exclusion | Adopted (A4.3.3 item 5) |
| minor 12 tier-III dilution | A1-class point-estimate check, gating (A4.5.5) |
| **Not adopted** | J_eff ≥ 15 (replaced by J_eff ≥ 10 with simulated size; see M2). Nothing else |

---

## A4.13 Changes in revision 2 (PI decisions D26)

| D26 item | Revision 2 |
|---|---|
| (1) DN2's 6 admitted projects enter DN3 with all files | Tier 0, walked first, under every A4 rule (A4.1.5, A4.3.3 carve-out, A4.3.4, A4.3.10). PXD011591 fails the protocol screen. The walk check includes tier 0 (`a4_power.dn2_tier0_projects`, `a4_r1_walk`). Evaluability +0.02 to +0.06 at rev 2. **Superseded in rev 3:** under item 5, tier 0 = {PXD038243}, ≈ 0 gain |
| (2) word-bounded keyword rule | Kept, unchanged from revision 1 (A4.3.5(a)) |
| (3) T1 escrow | A4.8 rewritten as a normative escrow protocol: preconditions, one-look and escrow mechanics in `sealed_run.py`, access log, unseal triggers (A4 complete, or the PI's dated time-limit decision), reporting |
| (1, cont.) thresholds unchanged, sized for v1 and v2_sep | Unchanged (A4.3.7) |
| (D26 item 4) second review | `reviews/amendment_A4_review2.md`; resolutions in A4.14 (revision 3) |

## A4.14 Changes after the second review (`reviews/amendment_A4_review2.md`; revision 3)

| Review item | Revision 3 |
|---|---|
| **B1** tier 0 vs same-lab item 5 | Item 5 applies to tier 0 (stricter reading). Person/lab keys are order-insensitive and include `labPIs` and `submitters`. Tier 0 = {PXD038243}. The outcomes are recorded before the tag (A4.3.3, A4.3.4, A4.3.10). A carve-out is left to an explicit PI decision (A4.11.1) |
| **B2** unseal on not-evaluable | Removed. T1 stays sealed until DN3 step 5, the next attempt's step 5, or the fixed date. Later null-set authors, reviewers and agents must be T1-blind (A4.8.3, A4.8.4) |
| **M1** escrow leaks | Escrow root outside the repo, with `--out` and `--work` inside it; whole-root tar encrypted to the PI's off-cluster key; manifest gives sha256 and total size only; wall time inside the escrow; leak test; whitelist console check; idle node (A4.8.2, A4.8.3, A4.4.8) |
| **M2** discretionary time limit | Fixed calendar date in the adoption note, not movable; not before DN3 steps 2–3; "if time allows" removed; abandonment is a violation (A4.8.4) |
| **M3** A3-freeze collision | Normative order: A4 code in frozen trees → A3 freeze → `prereg-v0.4.3` → A4 files → `prereg-v0.4.4`. The A3.6 route is removed. D27 waits or uses a branch (A4.4, A4.8.1.5) |
| **M4** optimistic evaluability | New `a4_r3_walk.py`: item 5 on the pool and tier 0, person links, three-factor cap, tier-0 shuffle fixed, tier-0 marker evidence (e = 0.5), empirical DN2 file-count reference cell. The figures replace all earlier ones (A4.0, A4.7.2, A4.11) |
| **M5** J_eff ≥ 10 non-binding; stale sizes | Stated as implied by the cap. The residual concentrated-leak size (0.023–0.057) is stated, and the PI accepts it in the adoption note. All size figures come from `a4_r3_walk` (A4.3.6, A4.3.7, A4.6, A4.10.6). J_eff ≥ 12 was not adopted, for the evaluability cost |
| **M6** amendment by "reading as"; stratum-1 estimator | Explicit replacement text for D12, A3.3.1, A3.3.7 row 1 and A3.3.2/3.3.8. "Unsealing spends T1" stated. The stratum-1-list estimator is defined (A4.8.0, A4.8.1.4) |
| minor 1 IUT precedence | invalid > fail > NE > pass; the simulation was re-run with it (A4.3.6, `a4_r3_walk`) |
| minor 2 clean invalid rate | Reported (A4.7.3, A4.10.5) |
| minor 3 `peptone from (meat\|casein)` dropped | Restored to the contamination group. `tryptone` was **not** added: lysogeny broth (LB) contains it and is near-universal for bacteria, so it would remove most of tier I. Bovine casein is caught by the marker screen (Panel B: CSN1S1/CSN2/CSN3) |
| minor 4 thresholds vs final n_min | A_cand thresholds defined as the max with functions of the final n_min (A4.3.7) |
| minor 5 file-guard carve-out | Limited to the project's own `dn2.sha256` files (A4.3.10) |
| minor 6 DN2 npz | Never opened; sha256 comparison only (A4.3.10) |
| minor 7 A1-class check | Moved to A4.5.5 as gating, on both pools |
| minor 8 retries | Hash and freeze handling, T1-blind review, traceback frames, D0/D1 reproduction, last success is the result of record (A4.8.2.5) |
| minor 9 `sealed_run.py` tag rule and wall time | Escrow-manifest commits skipped for the tag check; wall time inside the escrow (A4.8.2) |
| minor 10 speed under load | Idle node, or T1 speed labelled "not claims-grade" (A4.4.8) |
| minor 11 stale references | Fixed (A4.13, A4.10.5, A4.3.10; the `a4_r1_walk` docstring is superseded by `a4_r3_walk`) |
| minor 12 claim-script inputs | The A3.1.8 / D23 verdicts are named (A4.8.1.3) |
| **Not adopted** | `tryptone` (minor 3, reason above). J_eff ≥ 12 (M5 option b; evaluability cost). Nothing else |

## Appendix A: DENY list (frozen at tag)

**Lineages:** Metazoa, Viruses, Apicomplexa, Kinetoplastea, Metamonada, Parabasalia, Entamoebidae, Microsporidia, Amoebozoa, Ciliophora, Mollicutes.

**Genera:**
- *Homo, Mus, Escherichia, Arabidopsis, Saccharomyces, Pyrococcus, Trypanosoma, Staphylococcus, Synechococcus*;
- *Shigella, Salmonella, Klebsiella, Enterobacter, Yersinia, Acinetobacter, Mycobacterium, Mycobacteroides, Neisseria, Haemophilus, Bordetella, Brucella, Francisella, Legionella, Listeria, Clostridioides, Clostridium, Bacteroides, Phocaeicola, Parabacteroides, Bifidobacterium*;
- *Lactobacillus, Lacticaseibacillus, Lactiplantibacillus, Limosilactobacillus, Apilactobacillus, Lactococcus, Leuconostoc, Pediococcus, Fructobacillus*;
- *Helicobacter, Campylobacter, Vibrio, Aliivibrio, Chlamydia, Chlamydophila, Borrelia, Borreliella, Treponema, Leptospira*;
- *Mycoplasma, Mycoplasmoides, Metamycoplasma, Mycoplasmopsis, Mesomycoplasma, Malacoplasma, Ureaplasma*;
- *Streptococcus, Enterococcus, Porphyromonas, Prevotella, Tannerella, Fusobacterium, Aggregatibacter, Moraxella, Pasteurella, Actinobacillus, Coxiella, Rickettsia, Ehrlichia, Anaplasma, Wolbachia, Buchnera, Sodalis, Photorhabdus, Xenorhabdus, Serratia, Proteus, Citrobacter, Cronobacter, Stenotrophomonas, Achromobacter, Burkholderia, Aeromonas, Cutibacterium, Micrococcus, Kocuria, Schaalia, Gardnerella, Akkermansia, Faecalibacterium, Roseburia, Ruminococcus, Eubacterium, Blautia, Veillonella, Endozoicomonas*;
- *Candida, Cryptococcus, Malassezia, Histoplasma, Paracoccidioides, Coccidioides, Pneumocystis, Blastomyces, Talaromyces, Emergomyces, Madurella, Trichophyton, Microsporum, Epidermophyton*;
- *Vairimorpha, Anncaliia, Nosema, Encephalitozoon, Durusdinium, Cladocopium, Symbiodinium, Breviolum, Naegleria, Acanthamoeba, Candidatus*.

**Species:** *P. aeruginosa*, *A. fumigatus*, the *B. cereus* group, *C. diphtheriae*.

## Appendix B: constants for `config/a4.yaml`

**Seeds**

| Use | Seed |
|---|---|
| Walk order | `random.Random(47)` per tier |
| File pick | `random.Random(47 + digits)` |
| Mixtures | `default_rng(47)` |
| Cap hash | 47 |
| A-syn | 47 |
| Bootstrap | `default_rng(83)` |
| Checks | 131 |

**Levels**

| Rule | Level |
|---|---|
| tail-G2 | 0.975 one-sided |
| G1 | 0.995 two-sided |
| TOST | 0.95 two-sided |
| G2-original | 0.95, plus the cluster max |
| Validity controls | unchanged |

**Intervals:** outer envelope over the analytic, project and A4-cluster intervals.

**n_min** (provisional; final values from adoption precondition 4)

| Row | n_min |
|---|---|
| G1 | 11,100 |
| V6 | 30,600 |
| V12 | 26,800 |
| G3 | 2,000 |
| tail | 200 |
| V8 | 107,570 |
| Clustered rows | J ≥ 10 |

**Walk**

| Item | Value |
|---|---|
| Files per project | ≤ 48 |
| Projects per cluster | ≤ 5 |
| Cap share | 10%, min over the proxy, A_cand and \|A\| factors (\|B\| for B) |
| Rows | intersection–union over the capped and uncapped pools |
| Full stop | J ≥ 20; J_eff ≥ 10 (capped proxy and capped A_cand); proxy ≥ 1,200; A_cand ≥ 450,000 |
| v1 minimum | as the full stop, with A_cand ≥ 120,000 |
| Hard cap | 2,000 runs |
| Tiers | 0 (DN2's admitted projects) → I → II → III |
| Tier-0 accessions | PXD038243 eligible. PXD006818, PXD018953, PXD044078 and PXD062500 are excluded by item 5; PXD011591 by the protocol screen and item 5 |
| Person/lab keys | order-insensitive tokens of the crawl `lab`, `labPIs` and `submitters`; single letters dropped |
| IUT precedence | invalid > fail > not evaluable > pass |

**Escrow:** the root is `/tmp/vulcan_escrow_T1_<id>/` (outside the repo); the ciphertext is `data/escrow/T1/<id>.tar.age`, encrypted to the PI's public key (fingerprint in `config/a4.yaml`); `results/T1_escrow/MANIFEST.md` holds the sha256 and total size only; at most 2 retries; the unseal triggers are A4.8.4 (a) DN3 step 5 and (b) the fixed date.
