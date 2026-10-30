# Amendment A6 (DRAFT, revision 1): one null-spectrum gate test on human "removed-peptide" nulls (HN1), gating v1

- **Status:** draft for adoption as DESIGN.md §6.0e under A5.1.3 (GO branch). It is not binding. DESIGN.md, the tags and every A1–A5 file are unchanged.
- **Author:** statistician (agent), fresh context, T1-blind. I did not write A1–A5, their reviews, the D28 review or the C9 pilot, and I ran no gate step. I will not run any A6 step.
- **Date:** 2026-10-28.
- **Revisions:** revision 0 was reviewed by a fresh T1-blind agent (`reviews/amendment_A6_review.md`; verdict "adopt with changes": 3 blocking, 5 major, 5 minor). Revision 1 resolves every item (A6.17). The main changes: A6 gates **v1 only**; every label read moves to step 5; the swap test recomputes I from raw PSMs on every project; a throughput-based tag criterion; a status file for the A5 look.
- **Mandate:** RESEARCH_TREE D30 (C9 GO; draft A6 with the pilot's six conditions), A5.1.3 (GO branch), A5.10 (S6 floor), and the PI time limit D31 (node time limit: selection/admission and the single gate run finish by **18:00 CET today**).
- **Read:** RESEARCH_TREE D9–D30; DESIGN §6 including §6.0–§6.0d; `reviews/pi_proposal_D28_review.md`; `reviews/amendment_A2_draft.md`, `amendment_A3_draft.md`, `amendment_A4_draft.md` (rev 3), `amendment_A5_draft.md` (rev 1); `results/2026-10-28_C9_pilot/{NOTES.md, metrics.json, c9_power.json, certify_strict.json}`; `scripts/c9/*`, `scripts/a2/*` (headers, `gates_a2.py`), `scripts/a3/*` (headers, `gates_a3.py` V6 path), `vulcan/gates/*` (`a2_stats`, `a2_admission`), `vulcan/index/build.py`, `vulcan/gpu/engine.py` (signatures).
- **Checks run for this draft (DEV data, label-free on every null candidate; seed 131):**
  - `reviews/checks/a6_cert_yield.{py,json,out}`: the final admission rules on the pilot's 12 DEV runs (two-engine certification, f_max = 0.05, contaminant pairs dropped). Full-index labels were read only for certification, the run's hyperscore 1% threshold and the identified set, exactly as in the pilot.
  - `reviews/checks/a6_power.{py,json,out}`: `c9_power` re-run with the final rules (pilot condition 6).
- **Not read:** any sealed file; `results/SEALED_LEDGER.md` beyond its header; any T1/T5 output (none exists); any per-label quantity on any null candidate (none was computed by me or, per its NOTES, by the pilot).
- **Exposure:** T1 none; T5 none; C9 pilot: its label-free outputs (counts, simulations), which A6 uses for design (disclosed in A6.13).

---

## A6.0 What A6 is (first paragraph; quoted in every report of its result)

**A6 is one null-spectrum gate test on a fresh set of human "removed-peptide" null spectra (HN1). For v1 (`static_gbm_pred`) it is the second completed test, run at level α/2, after v1 FAILED tail-G2, G1 and G3 under A1 (D12) and after A2 (DN2) was not evaluable before any test and A4 (DN3) was not run. It is designed and run T1-blind, under the A5.10 floor, before the T1 look. Its result enters the T1 labels and void rules through A5.4, A5.7 and A5.9.3, whatever it is.**

**Attempt history (quoted in every report, A6.12):**

| Attempt | Null set | Scorers | Outcome | Per-label quantities? |
|---|---|---|---|---|
| A1 (`prereg-v0.4.1`) | M3e species and offset nulls | v1 | tail-G2, G1, G3 **FAILED**; final (D12). The contamination explanation is post hoc, unconfirmed, and not used | yes |
| A2 (`prereg-v0.4.2`) | DN2 | v1 | not evaluable on label-free counts (D24); step 5 never ran | no |
| A3 (D22) | DN2 | v2_sep (gated), ce_s (licenses nothing) | moot (DN2 not evaluable) | no |
| A4 (draft rev 3) | DN3 | v1, v2_sep, ce_s | drafted, reviewed twice, **not adopted, not run** (D28.1) | no |
| C9 pilot (D28.2) | 12 DEV human runs | none (label-free) | GO: simulated P(v1 minimum) ≥ 0.8 (D30) | no |
| **A6 (this)** | **HN1** | **v1 (test 2, α/2)**; v2_sep and ce_s **not run** (A6.2) | — | — |

**What A6 changes in earlier text (dated, before any HN1 listing):**
1. **A2.4, last sentence of rule 2** ("Human spectra are never used as donors …") does not apply to HN1. A6.3–A6.4 answer each of A2.4's three routes by construction (removal of the true peptide, its relatives and every pair a fixed score identified in the project), and A6.11 states the residual risks and their direction. This is a new null construction, not a change of any margin or statistic.
2. **A5.1.3 time box** is shortened by PI decision D31: selection, admission and the single gate run end by **18:00 CET on 2026-10-28** (A6.10). Otherwise A6 is recorded **"stopped: node time limit"** (D31.3).
3. **Scope.** A6 gates **v1 only**. D30 asked for v2_sep and ce_s as well; A6.6.2 shows that v2_sep's tier 1 cannot be evaluated on single-hypothesis human nulls under A3 and the A5.10.3 floor, and ce_s licenses nothing and does not fit D31. Neither is run, and A6 is not a test for either (A5.10.1).
4. Nothing else. Every A5.10.3 floor item is met or tightened (A6.14, compliance table).

---

## A6.1 Principles (inherited; nothing relaxed)
1. **A2 P1–P5 and A3 Q1–Q6 apply**, reading DN2 as HN1. In particular:
   - **P1, label-swap invariance:** every admission and removal decision is a function of the spectrum, the frozen certification outputs and the *unordered* T∪D with its pair map;
   - **P2:** memberships and label-free counts are hashed and committed before any per-label quantity exists;
   - **P4, one shot:** a fail is final under this protocol; nobody re-draws, extends or prunes HN1 after step 5.
2. **No null-spectrum-driven removal (pilot condition 2).** No hit of a null spectrum itself (its rank ≥ 2 PSMs, chimeric secondary hits, its winner in any null search, or any score it gives a candidate) is used to remove a candidate, a pair or the spectrum. Removal is fixed before any null search is run (A6.4.4).
3. **Certification is frozen first.** The certified list, the identified set I and the per-run thresholds are hashed and committed before the reduced database is built or searched (A6.4.6).
4. **Stricter reading wins.** Any ambiguity found during implementation is resolved by the stricter reading and recorded under A6.9.
5. **Report regardless of outcome** (A6.12).

---

## A6.2 Scorers, roles and levels (closed)

| Scorer | Frozen identity | Role in A6 | Licenses |
|---|---|---|---|
| v1 `static_gbm_pred` | sha256 `31bfe685…` (A1.1) | **Gated: completed test 2** | tail-G2/G1/G3 for v1 (A2.9.6 wording, A6.12) |
| v2_sep | — | **not run** (A6.6.2); no per-label quantity; not a test | nothing; designation stays "candidate"; next test is still test 1 |
| ce_s | — | **not run** (licenses nothing; does not fit D31); not a test | nothing |
| every scorer | — | **HN1 is spent for every scorer once v1's step 5 starts** (any later amendment needs a new null set) | — |

**Levels (A5.10.2 and A4.6).** Every re-tested v1 row runs at α/2:
- tail-G2: one-sided 0.975;
- G1: two-sided 0.995;
- TOST rows (G3, V3, V4, V12): two-sided 0.95 (i.e. TOST at 2.5% per side);
- G2-original and G5: their unchanged levels (A4.6: not re-tested gates keep their levels).


**Tag and hook into the A5 look (closed).**
- **Tag.** A6 is tagged `prereg-v0.4.5` (the name A5.1.3 gives it) on the commit that holds the adopted text and the A6 code. Under D31, `prereg-v0.4.3` and `prereg-v0.4.4` may not exist yet; A6 does not need them, because it gates only v1, whose identity is fixed by A1.1 and guarded against `prereg-v0.4.2` (A6.8.4). The out-of-order tag number is confirmed by the PI at adoption (A6.16).
- **Code placement.** `scripts/a6/`, `config/a6.yaml` and `tests/a6/` are outside every A3/A5 freeze tree. No A6 file changes after the tag.
- **Status file (the only channel to the look).** A6 writes `results/2026-10-28_A6_hn1/a6_status.json` with a closed enum `status` ∈ {`PASS`, `FAILED`, `not evaluable`, `invalid pool`, `not completed`, `label-free stop`, `not tagged`, `stopped: node time limit`}, the step-5 commit sha (or null) and the sha256 of the gate output. The A5 claim script reads this fixed path; the path and the enum are committed before the A5 freeze. If the file is absent, unverifiable, or not committed when the look starts, the A5.9.3 state is **"stopped: node time limit"**. No one sets the state by hand.

---

## A6.3 HN1 selection (metadata, then label-free per project)

### A6.3.1 Frame (metadata only)
- **Source:** `data/manifests/t4_frame.tsv` (sha256 recorded in `config/a6.yaml`) and the cached crawl records `data/cache/pride_meta/projects_p*.json.gz` (sha256 recorded). Nothing is re-crawled.
- **Row filter (the C9/E9 frame, unchanged):** `eligible == 1`; organisms == "Homo sapiens (human)"; DDA; label LFQ; enrichment none; vendor Thermo; instruments match `scripts/e9/select_e9.py INSTR` (QE, QE Plus, QE HF, HF-X, Exploris 240/480: the A1 OT-HCD class); title + keywords + experiment types do not match `select_e9.BAD_KW`.

### A6.3.2 Exclusions (union; frozen as `config/a6_excluded_accessions.tsv` with sha256)
1. every accession in any `data/manifests/*.tsv` except `t4_frame.tsv` (dev; sealed T1–T5; the T4 selection; tier S; train, train_m7, train_itcid; predictor-training accessions; DN, DN2, DNIT; D8/E9);
2. every project in any T4 draw under `results/**/t4_selector/**/t4_selection.tsv` (T4-reachable);
3. every project E9 tried (`data/manifests/e9_d8.log.json`, 166), which includes the C9 pilot's runs; and PXD001468 (D0);
4. every `PXD\d{6}` string in any **tracked** file under `data/manifests/`, `results/` or `config/` at the tag (A4.3.3 item 1; includes `results/2026-10-28_E7_chimeric_v2/train_v2_folds.tsv` and `config/a2_excluded_accessions.tsv`);
5. every accession with a cached PRIDE file listing `data/cache/pride_meta/files_PXD*` (A4.3.3 item 2);
6. **same-lab rule (A4.3.3 item 5):** every project sharing a person/lab key (A6.3.3) with an A1 null project (`g2null.tsv` and its "tried" list), a DNIT walk candidate or a DN2 walk candidate.

**Pool at drafting** (`a6_power.json`, metadata only): the pilot's 1,465 projects minus 50 (15 by item 4, 35 by item 6) = **1,415 projects**.

### A6.3.3 Walk order and clusters
- **Order:** sort the pool by accession, then `random.Random(S).shuffle`, with **S = int(first 8 hex digits of the `prereg-v0.4.5` tag's commit sha, 16)**. No one can compute the order before the tag. (The pilot's simulation used seeds 47+i on the same pool, so seed 47 is not used.)
- **Person/lab keys (A4.3.4, unchanged):** the crawl `lab` field plus every `labPIs` and `submitters` name in the cached record; lower-cased, split into alphanumeric tokens, single letters dropped, tokens sorted and joined. A project with no key gets the key "?<accession>".
- **Cluster:** a connected component of **admitted** projects linked by any shared key. Re-uploads are handled by the file guard (A6.3.4(b)).
- **Cluster slots:** a candidate whose cluster would exceed **5** admitted projects is skipped as "cluster full", before listing.

### A6.3.4 Per-project procedure (in walk order; every step logged in `selection_log.tsv`)
- **(a) Listing.** List the project once. A RAW is usable if it is `.raw`, 0.1–2.0 GB, does not match `BAD_FILE` (as E9), and passes the file guard (b).
- **(b) File guard.** A RAW is unusable if its PRIDE SHA-1, or its (name, size), matches any file in any committed `*.sha256`/`cache_raw.sha256`, or any RAW already selected in HN1.
- **(c) Files.** Order the usable RAWs with `random.Random(S + accession digits)`. Download in that order until **3** runs pass (d), or the list ends.
- **(d) Run checks (label-free, before any search):** ≥ 99% of MS2 scans FTMS HCD; no TMT/iTRAQ reporter ions (the production pilot detector); MS1 present; decodable. A failing run is replaced by the next file.
- **(e) Admission.** A project is **admitted** iff ≥ 1 run passes (d). Admitted projects take a cluster slot whatever their yield. Each admitted run goes through A6.4.

### A6.3.5 Cluster cap and the two pools (A4.3.6, adopted unchanged)
- Water-filling factors f_j per cluster so that no cluster holds more than **10%** of the capped total; computed separately for the tail proxy, A_cand and |A|; the **minimum** of the three is used. If J·0.10 < 1, HN1 is not evaluable.
- Cluster j's A spectra are kept iff `u(sha256(run | store index | 47)) < f_j`.
- The **capped** and **uncapped** memberships are both committed. Every row must pass on both (intersection–union), with precedence **invalid > fail > not evaluable > pass**.

### A6.3.6 Stopping rule (label-free; evaluated on capped counts after each admitted project)
The walk **stops** at the first of:
1. **v1 minimum reached** (A4.3.7 values): J ≥ 20; J_eff ≥ 10 for both the capped tail proxy and capped A_cand; capped tail proxy ≥ 1,200; capped A_cand ≥ max(120,000, 10 × final G1 n_min);
2. **hard cap:** 400 admitted runs;
3. **walk deadline: 16:30 CET, 2026-10-28** (the last project being processed is "reached, not used" and excluded from future sets);
4. pool exhausted.

**Why the walk stops at the v1 minimum** (and not at a v2_sep full stop): A6 gates only v1 (A6.6.2), the time limit (D31), and stopping there keeps the validity controls at the size the GO decision was simulated for (A6.11, risk R7).

**End of walk without the v1 minimum** (stops 2–4): HN1 is **not evaluable for every scorer**. Steps 4–5 do not run. This is a label-free stop, not a completed test (A5.10.2: v1's next test stays at α/2).

### A6.3.7 Downstream
- Every **listed** project (admitted, or rejected after listing, or "reached, not used") gets role `dev-null` in `data/manifests/hn1.tsv`, is excluded from every future null set, and is added to `scripts/data/select_t4.py`'s exclusions. Projects skipped as "cluster full" before listing were seen as metadata only and are **not** burnt. **Check:** with the new exclusions, the metadata-only T4 selector reproduces the committed T4 selection byte for byte; if not, stop and investigate T1-blind before step 4.

---

## A6.4 Spectrum admission (per admitted run; label-free on every null candidate)

### A6.4.1 Certification (two engines; pilot condition 3)
A spectrum s of an admitted run is **certified** with peptide p iff all hold:
1. **Hyperscore (as the pilot, `scripts/c9/pilot_yield.certify`):** production closed search of the run against the **full** human + contaminant prod index (E2 SearchParams; report 5); rank-1 PSM is a human target; per-run paired double competition (hash tie-break seed 11); PSM-level q ≤ 0.001 and peptidoform q ≤ 0.01; p not A1.3-near the contaminant FASTA (`vulcan.gates.a2_admission.certify_run`).
2. **Second engine (Sage, the frozen comparator binary `tools/bin/sage`, sha256 in `config/a3.yaml`):** Sage searches the project's admitted runs (mzML from ThermoRawFileParser, `conda_envs/trfp`) with the T1 driver's plain `sage` comparator configuration, database = the prod index's human + contaminant target FASTA with that configuration's decoys. Its rank-1 PSM for s is a target with `spectrum_q ≤ 0.001`, the same I/L-collapsed stripped sequence as p, and the same charge.
3. A run whose mzML conversion or Sage search fails is unusable and is replaced (A6.3.4(c)).

**Why two engines, not project-pooled q ≤ 1e-4.** Per-run q ≤ 1e-4 is not resolvable by TDC (pilot: 11 of 12 runs give 0). Project pooling does not fix it: 9 of the 11 pilot projects certify 0 spectra at pooled q ≤ 1e-4 (`a6_cert_yield.json`). Two-engine agreement keeps **98.7%** of certified spectra and **98.5%** of A_cand, but only **86%** of the tail proxy (118 of 137). The 1.5% of A_cand where Sage disagrees carries 14% of the tail: a 10× enrichment (tail rate ≈ 5.3% vs 0.51%). That is the label-free signature of false certifications (risk R4): their true peptide is still in the database. Disclosure: this tail-retention figure was seen before the rule was fixed (A6.13).

### A6.4.2 Single precursor hypothesis (A2.2(e))
Each certified spectrum carries exactly one hypothesis in every null search: certified mono m/z, charge z, isotope error 0. The MS1 co-isolation count is reported, never used as a filter.

### A6.4.3 Identified set I (per project; symmetric)
I is the set of **unordered T∪D pairs** whose paired double-competition winner, in any admitted run of the project, has hyperscore ≥ that run's peptidoform-level 1% threshold, **whichever label won**. Only rank-1 PSMs of the full-index certification search enter. Each pair is represented by both member sequences (I/L-collapsed).

*Change from the pilot, and why:* the pilot built I from accepted **targets** only. A falsely accepted target (≈ 1% of I) was selected because it scored well on some spectrum s′. If s′ is correlated with a null spectrum s (a repeated MS2 of the same precursor, or a chimeric partner), removing that pair removes a target that would also have scored well on s by chance, which is pro-decoy. Building I from winning *pairs* of either label removes this route symmetrically. The yield effect is ≈ 1% more pairs in I.

### A6.4.4 Removal set R_P (per project; fixed before any null search)
A T∪D pair is in R_P iff **either member** meets any of:
- **(i) near I:** I/L-identical to, Hamming ≤ 2 from any window of, or sharing a 6-mer with any sequence in I or any isobaric variant (N↔GG, Q↔AG), all modforms (`vulcan.gates.homolog.HomologFilter`, as the pilot);
- **(ii) ladder (pilot condition 1):** it lies in the single-hypothesis precursor window (±10 ppm, isotope 0) of **any** certified p of the project and shares **≥ 0.05** of p's b2…b(L−1), y2…y(L−1) ions (z = 1, ±20 ppm, cross-type, all modforms; `pilot_yield.f_scores`). This is applied project-wide, so a pair removed for one certified spectrum is removed for every null spectrum of the project (stricter than the pilot's per-spectrum rule, and needed so that one reduced index serves the project);
- **(iii) contaminants (pilot condition 4):** its sequence occurs in a protein of the contaminant FASTA (`data/fasta/contaminants_frankenfield2022.fasta`) or of Panel B (A2 appendix), or is A1.3-near (6-mer or Hamming ≤ 2) to a **non-human** protein of either.

R_P is closed under pairing by construction (the pair is the unit). It is a function of (the frozen certification outputs, I, the unordered T∪D with its pair map, the fixed FASTAs), and of nothing else.

*(iii) replaces A2.2(c)'s winner filter.* The A2 filter dropped a null spectrum when a pipeline's winner was near a contaminant. That is a decision driven by hits on the null spectrum (forbidden by A6.1.2). Removing the same neighbourhoods from the database beforehand is at least as strict: no such pair can win. The A2.2(c) check is still computed at step 4 as a diagnostic; it is vacuous by construction and must find nothing. Human Panel H proteins are part of the human database and are not removed.

*A2.2(d) (family screen) is not used.* It is vacuous: the reduced database is a subset of the full one, in which p won.

### A6.4.5 Strata and label-free counts
- **A** = certified spectra of admitted runs. **A_cand** = A spectra with ≥ 1 candidate (≥ 4 matched peaks) in the production search of the reduced index (A6.5.1).
- **Tail proxy** = A_cand spectra whose label-collapsed max hyperscore over their reduced-index candidates is ≥ the run's own hyperscore 1% threshold (from the full-index certification search).
- There is **no stratum B** (uncertified human spectra are not nulls) and no offset nulls. **FS** (A2.4 rule 2: every A spectrum's fragments shifted by +2/7 Th, certified precursor kept) is searched against the **reduced** index (with the full index, the true peptide's precursor and RT would remain).

### A6.4.6 Freeze, swap test and label guard (pilot condition 2)
1. **Freeze.** Per project: certified list (store index, p, z, mono m/z), I, the per-run thresholds and R_P are written, sha256-hashed and committed **before** the reduced index is built.
2. **Swap test (P1), on every admitted project** (set operations; seconds per project): rebuild the gate index from the same FASTA with target and decoy labels swapped. Keep frozen only what is legitimately label-dependent: the certified list (p, z, mono m/z) and the per-run threshold scalars. **Recompute everything else from the raw full-index rank-1 PSM table and the swapped index**: the per-run double competitions and their pair winners, I (pairs whose winner is ≥ the frozen threshold), R_P, the reduced index, A and A_cand membership and the tail proxy. The I pair set, the R_P pair set and both memberships must be **bit-identical** to the unswapped run. A mismatch stops A6 (implementation stop; not a test).
   - **Unit test (synthetic labels, `tests/a6/`):** a decoy-won pair above the threshold is in I; a target-won pair just below it is not; swapping labels leaves I unchanged.
   - **Count assertion per project:** |I| equals the number of distinct pairs whose rank-1 pair winner is ≥ the run threshold, counted by the certification module without any split by label.
3. **Label guard.** In steps 2–4, only the certification module opens `labels.arrow`/`proteins.arrow`, and only to index rows of the **host** full-index search (rank-1 PSMs and their pair partners: certification, the per-run thresholds and I). A test asserts that no label array is ever indexed at a row that came from a reduced-index search or a null copy, and that no step 2–4 output carries a per-label key (as `tests/a2/test_no_label_counts.py`). Every other label read happens in `gates_hn1.py` at step 5.

---

## A6.5 Searches and mixtures (step 4; label-free)

### A6.5.1 Reduced index (per project)
The reduced index is the prod index with the rows of R_P removed, built by `scripts/a6/reduced_index.py` as a row subset of the frozen `L1Index` arrays (`pf_id, mass32, seq_off, length, mask, n_var, missed, partner`; `residues` shared; `partner` remapped; a `row_map` back to full-index rows). `vulcan.gpu.engine.GpuIndex` and `engine.search`/`search_exhaustive` are used **unchanged**.
- **Exactness test (before the tag, dev D0):** with R = ∅ the PSM tables are bit-identical to the full index; with the pilot's R, every returned row maps outside R and has the same hyperscore as for the same (spectrum, row) on the full index.
- If this cannot be done outside the frozen trees, A5.1.4's fallback applies (A6 code lands before the A3 freeze); if that cannot happen by the tag deadline (15:15 CET), A6 is "not tagged" (a label-free stop).

### A6.5.2 Self-host mixture (one per admitted run)
- **Host:** the null's own run, **unmodified**: all its MS2 spectra, production hypotheses, searched against the full prod index (this is the certification search, reused).
- **Nulls:** provenance-flagged **copies** of the run's A spectra (prov 20) and of their FS versions (prov 22), each with the single certified hypothesis, searched against the project's reduced index.
- The two PSM tables are concatenated (`spec_idx` offset; `pf_row` mapped back to full-index rows) **before** any feature, calibration or scoring step. Step 4 then computes, label-free: features (`scripts/m3/features.py`, as for the A2 mixtures), `StaticModel.apply` scores for hyperscore, v1 and v1 noRT, and `psm_competition` winners.
- **Each pipeline's own per-mixture 1% threshold** (`double_competition` on the mixture, as `scripts/a2/gates_a2.collect`), tail membership and every target/decoy split are computed **only at step 5**, by `gates_hn1.py`. They use labels of null-copy rows, so they are per-label quantities (P2).
- **RT** is the run's real RT; no quantile mapping (the host is the source run).
- **Why not the pilot's "own run with certified spectra replaced":** the certified spectra are a large share of the run's confident IDs (≈ 25–60%). Replacing them removes those targets from the run's competition, inflates every per-run threshold and degrades run-level calibration, which is not the deployment regime. Copies leave the host exactly as deployed; the nulls add ≈ 10% spectra, mostly low-scoring, as the A2 mixtures added 20,000 to a D0 host.
- **Why not a D0 host:** it needs RT quantile mapping and a mixture with two databases per host; the self-host is simpler and closer to deployment.

### A6.5.3 What is computed in step 4 (label-free; committed with sha256)
Per mixture: the PSM tables, the feature tables, the scores and winners (row ids only) of every pipeline (hyperscore, v1, v1 noRT) for every null copy, the production and exhaustive searches of every A copy on the reduced index (for G3; row ids and scores only), the capped and uncapped memberships, the per-cluster, per-project and per-instrument label-free counts, and the A2.2(c) diagnostic (expected empty; it reads sequences, not labels).

---

## A6.6 Gate rows

All rows use A2/A3 statistics with the A4 changes adopted by A5.10.3 (outer envelope, cluster cap, intersection–union, final n_min) at the A6.2 levels. "Project" = PRIDE accession; "cluster" = A6.3.3 component.

### A6.6.1 v1 (on HN1)

| Row | Pool | Rule |
|---|---|---|
| tail-G2 | A | A2.5 (n_tail ≥ 200; UCB_abs ≤ 0.60; UCB_diff ≤ 0.10; hyperscore validity control: two-sided 95% Wilson CI ∋ 0.5) with **UCB = max(CP or delta, project bootstrap, cluster bootstrap)**, B = 10,000, `default_rng(83)`, one-sided 0.975. **Difference clause, stricter:** max(UCB_diff, UCB_abs − min(p̂_hyperscore, 0.5)) ≤ 0.10, so a hyperscore tail inflated by contamination never loosens it |
| G1 | A | A2.6 rule at two-sided 0.995; pass iff Fisher-z (Bonett–Wright) ∪ project jackknife-t ∪ cluster jackknife-t ⊂ (−0.05, 0.05); n_top ≥ final G1 n_min (≥ 11,100); hyperscore control (99% CI ∋ 0) |
| G3 | A (A ∪ B ≡ A) | A2.8 TOST per bin at two-sided 0.95: Newcombe ∪ project bootstrap ∪ cluster bootstrap ⊂ [−0.01, +0.01]; ≥ 2,000 per label in each top-half decile; production = `engine.search` on the reduced index, exhaustive = `search_exhaustive` on the reduced index |
| G2-original, G5 | A | §6.4 margins at their unchanged levels; G2-original UCB = max(§6.4 UCB, project bootstrap, cluster bootstrap) |
| FS rows | FS | A2.4 rule 3: tail-G2, G1, G5 if evaluable; an evaluable fail blocks; not evaluable does not block |
| Subset check (A4.5.5, generalised by A5.10.3 item 7) | A | for each instrument class (QE/QE Plus; QE HF/HF-X; Exploris 240/480), each cluster and **each project** with n_tail ≥ 200: tail-G2 point estimates ≤ 0.60 and ≤ 0.10, on both pools; an exceedance is a **fail** |

**v1 verdict:** pass iff tail-G2, G1, G3, G2-original and G5 all pass, no evaluable FS row fails and no subset check fails, each on both pools. Reported, never used to pass: noRT, per project, per cluster, per instrument class, leave-one-cluster-out, each interval component, both pools, the cap factors.

### A6.6.2 v2_sep (tier mapping per A3.2.4, and why it is not run)

| A3 row | Tier | HN1 pool | Status in A6 |
|---|---|---|---|
| V1 tail-G2 | 1 | A, rank 1 | as A3, A6 levels and envelope |
| V2 G1 | 1 | A | as A3 |
| V3 G3 | 1 | A | = v1's G3(A), asserted byte-identical |
| V4 G3, A ∪ B all units | 1 | A (B empty) | = V3 |
| **V6 class 1** (neighbour-filter symmetry, rank-1 winners of null spectra with ≥ 2 groups) | 1 | **none**: A spectra carry one hypothesis and B does not exist | **not evaluable by construction** |
| V7 G2-original | 1 | A | as A3 |
| V9 (a), (b) G5 | 1 | A | as A3 |
| V10 FS | 1 | FS | as A3 |
| V12 top-5 truncation symmetry | 1 | A | as A3, n_min ≥ 26,800 per label per bin |
| V5, V6 class 2, V8, V9 (c) | 2 | none (no null secondary units) | not evaluable by construction |
| V13, V14 (A-syn) | 2 | A-syn is not built in A6 | not evaluable |

**Consequence (stated plainly):** under A3's rule "not evaluable blocks the tier" and the A5.10.3 floor, **v2_sep's tier 1 is not evaluable on HN1**, because V6 class 1 needs null spectra with ≥ 2 precursor groups. A single-hypothesis human null has none, and real secondary groups of human spectra are not nulls. Building synthetic multi-group nulls (A-syn groups attached to A spectra) is a new construction that cannot be implemented, dry-run and reviewed within D31.

**Rule (closed).** v2_sep is **not run** in A6: no step 4b, no per-label quantity. A6 is not a test for v2_sep (A5.10.1), its next test is still test 1, and HN1 is spent for it once v1's step 5 starts (A6.2). Its A5.9.3 state is "A6: not run (V6 class 1 not evaluable by construction)", and the T1 title stays "Secondary arm; not gate-validated". A later v2_sep test needs null spectra with ≥ 2 precursor groups (for example synthetic secondary groups attached to certified spectra), designed, simulated and reviewed in its own amendment.

### A6.6.3 ce_s
**Not run.** It licenses nothing (A3.2.7), its scoring path is not implemented for self-host mixtures, and D31 leaves no time. A6 is not a test for ce_s.

### A6.6.4 n_min (A5.10.3 item 3)
Final values are derived before the tag by `reviews/checks/a4_power.py`'s A4.5.3 rule with ≥ 4,000 replicates, p̂ ≥ 0.92 and logit-sd 0.5, and are never below the A4 provisional table (G1 11,100; G3 2,000 per label per bin; tail 200). If that run cannot finish before the tag, the provisional table is used; it is a floor either way, and A5.10.3 item 3's re-derivation is then recorded as not done (a disclosed gap, A6.16). Fewer than 10 clusters makes every clustered row not evaluable.

---

## A6.7 Multiplicity and the attempt ledger (A5.10.1–A5.10.2)

| Scorer | Completed tests before A6 | A6 is | Level | After A6 |
|---|---|---|---|---|
| v1 | 1 (A1: FAILED) | test 2 | α/2 | completed without a pass → next test α/4; label-free stop → next test α/2 |
| v2_sep | none | not a test (not run) | — | test 1 remains, at α, on a new null set |
| ce_s | none | not a test (not run) | — | — |

- **Completed test (A4.6/A5.10.1), applied literally:** as soon as any per-label quantity of v1 on HN1 exists, A6 is v1's completed test 2, whatever happens next: pass, FAILED, not evaluable at step 5, **invalid pool**, or any stop during step 5 (time limit, guard failure, crash), recorded as "not completed". Each of these except PASS moves v1's next test to α/4.
- **Label-free stop:** any stop before the first per-label quantity (not tagged, walk without the v1 minimum, implementation stop, time limit before step 5). v1's next test stays at α/2.
- Reports state: "per-row nominal level α/2, conditional on A1's failure". No "procedure-level" level over attempts is claimed (A4.6). Simulated per-attempt false-pass rates at the margins are quoted (A6.11).

---

## A6.8 Order of operations, roles and time (D31)

### A6.8.1 Roles (each a separate agent context; all T1-blind)
| Role | Who | Must not |
|---|---|---|
| Author | this agent | run any step |
| Reviewer | a fresh agent (A6.15) | have written any A-amendment or run any step |
| Implementer | a fresh agent, commissioned as soon as this revision is committed: `scripts/a6/`, `config/a6.yaml`, `tests/a6/`, the dry run, the throughput measurement and the a6_power re-run | run steps 2–5 |
| **Selection/admission agent (S)** | a fresh agent: steps 2–4 (walk, downloads, certification, I, R_P, freeze, swap tests, reduced indexes, mixtures, features, label-free counts) | run step 5; open any `labels.arrow` except through the certification module |
| **Gate agent (G)** | a fresh agent: step 5 only | have written A6, implemented it, or run steps 2–4 |
| Deviation reviewer | a fresh agent per deviation (A6.9) | have run the step the deviation affects |

The commissioning prompt of every role is committed verbatim with its artefact.

### A6.8.2 Steps and clock (CET, 2026-10-28)
| Step | Content | Deadline |
|---|---|---|
| 0 | This revision and its review committed | 14:00 |
| 1 | PI adoption (verbatim, or with T1-blind-reviewed changes): DESIGN §6.0e and a RESEARCH_TREE note with the adopter's exposure; implementer commissioned at once | 14:15 |
| 2 | Implementation and tests (synthetic), including the `a6_stats` acceptance test (A6.8.3); dev dry run on D0 and two D8h runs (label-free counts; synthetic-label gate path); **end-to-end throughput r** (admitted runs per hour, including listing, download, mzML, Sage, certification, I, R_P, reduced-index search); cert-yield and **a6_power re-run with the dry-run yields** (pilot condition 6) | 15:15 |
| 3 | **Tag criterion (label-free):** let N = r × (16:30 − tag time). Tag `prereg-v0.4.5` iff the re-run gives P(v1 minimum within N runs) ≥ 0.8 at ×1 yields **and** ≥ 0.5 at ×0.5 yields. Otherwise record **"not tagged: infeasible"** (a label-free stop; no project is listed, so no pool is burnt) | 15:15 |
| 4 | Agent S: walk to the A6.3.6 stop with admission inside the walk (prefetching downloads, decisions strictly in walk order), freeze, swap tests, reduced indexes, mixtures, features and label-free scores; commit memberships, counts and hashes | walk ends by 16:30; step 4 committed by 17:10 |
| 5 | Agent G: one session, from the tag commit plus the committed step-4 records, clean tree: v1 rows, `a6_status.json`, commit | 17:50 |
| 6 | NOTES and the RESEARCH_TREE decision | 18:00 |

From the current figures (A6.10), the tag criterion needs N ≳ 170 runs (×0.5 yields: P = 0.33 at 150 and 0.86 at 200 runs), i.e. r ≳ 135 admitted runs/h for a 15:15 tag. Whether this is reachable is measured, not assumed; E9's ≈ 580 runs/h was download-only.

**Hard stop and resources.**
- A6 uses the GPU and a lock file `/tmp/vulcan_a6.lock` that are not shared with the A5 freeze/G0 run, the speed track's bench lock or the look. If they cannot be separated, A6 work yields to them and its clock is unchanged.
- **All A6 work stops at 17:55 CET** whatever its state. A6.7's completed-test rule then applies. Partial outputs are committed and reported, never re-run.
- **Time limit (D31.3).** If step 5 is not committed by 17:55, the status is **"stopped: node time limit"**; it is a completed test iff any per-label quantity of v1 exists (A6.7).

### A6.8.3 Implementation (minimal; reuse first)
All A6 code lives outside the frozen trees (`scripts/a6/`, `config/a6.yaml`, `tests/a6/`) and imports frozen code unchanged (A5.1.4).

| File | Reuses (exact targets) | New code |
|---|---|---|
| `scripts/a6/select_hn1.py` | `reviews/checks/a6_power.pool` (frame, exclusions, person/lab keys); `scripts/e9/stream.py` (`download`: size and PRIDE SHA-1 verified, `process_batch`); `scripts/e9/select_e9.py` (`INSTR`, `BAD_KW`, `BAD_FILE`) | seed from the tag; cluster slots; stopping rule on capped counts; selection log |
| `scripts/a6/admit_hn1.py` | `scripts/c9/pilot_yield.py` (`certify`, `ladders`, `f_scores`, `null_stage` logic, `single_feature_flags`); `reviews/checks/a6_cert_yield.py` (Sage agreement, contaminant drop); `bench/t1/comparators.py` (`sage_config`, `run_sage`); `conda_envs/trfp` for mzML | symmetric I (A6.4.3); project-wide R_P (A6.4.4); freeze files; swap test on every project |
| `scripts/a6/reduced_index.py` | `vulcan.index.build.L1Index` arrays; `vulcan.gpu.engine.GpuIndex`, `search`, `search_exhaustive` | row-subset index + `row_map`; exactness test |
| `scripts/a6/mixtures_hn1.py` | `scripts/a2/build_mixtures.py` (`build`/`search` pattern; prov codes 20/22); `vulcan.gates.a2_fs` (single hypothesis, FS shift); `scripts/m3/features.py` (feature tables); `vulcan.rescore.static.StaticModel`; `scripts/m3e/gates_a1.neutral_rt` | self-host concatenation (A6.5.2); label-free winners |
| `scripts/a6/a6_stats.py` | `vulcan.gates.a2_stats` (`tail_g2_a2`, `cluster_bootstrap_tail`, `pool_validity`, `g1_gate_a2`, `g3_gate_a2`, `tost_bin`); `reviews/checks/a4_power.cap_factors` (copied verbatim) | envelope over project and cluster units for every bound-based row; G1 jackknife-t hulls; TOST project/cluster bootstrap; cap subsample; IUT precedence; the stricter difference clause; subset check |
| `scripts/a6/gates_hn1.py` | `scripts/a2/gates_a2.py` (`evaluate`, `collect`, `g3` structure; `double_competition` thresholds) | one-session runner; refuses unless the step-4 record is committed and unmodified; writes `a6_status.json` |
| `config/a6.yaml` | — | every constant of Appendix A, file hashes, the clock |

**`a6_stats` acceptance test (before the tag; part of step 2).** On fixed synthetic inputs with known structure (clean, margin, concentrated leak, invalid pool), `a6_stats` must (i) reproduce `vulcan.gates.a2_stats` exactly where it calls it, (ii) reproduce, on the same synthetic states and seeds, the cap factors of `reviews/checks/a4_power.cap_factors` and the per-pool verdicts and IUT precedence of `scripts/c9/c9_power.tail_iut` (bounds to 1e-12), and (iii) for the new interval components (G1 jackknife-t, TOST and G2-original bootstraps) show one-sided coverage ≥ 0.95 at nominal in 1,000 synthetic replicates. A failure means **"not tagged"**.

### A6.8.4 Guards
- **v1 closure guard:** `scripts/a3/v1_closure_guard.py` (unchanged) passes for v1's scorer, feature and search closure against `prereg-v0.4.2`, plus `scripts/a6/` at the A6 tag.
- **Clean tree** at step 5; `gates_hn1.py` verifies `config/a6.yaml`, the scorer and classifier hashes and the committed step-4 record.
- **No per-label output** from steps 2–4 (test).

---

## A6.9 Deviations
- Any departure from this text is a **dated note committed before the step it affects**, with the stricter reading chosen where the text is open.
- A deviation that **changes a rule** (a threshold, a pool, a level, a stopping criterion) is not allowed; it stops A6 (label-free stop if before step 5).
- An **implementation correction** (a bug, as D16) is allowed only before any per-label quantity exists, only if label-free, and is reviewed by a fresh T1-blind deviation reviewer before the affected step resumes. A step affected by a corrected bug is re-run **once**, with the same seeds.
- Every deviation, its review and its effect on label-free counts are reported in NOTES.

---

## A6.10 Feasibility and power (`a6_power.json`; final rules; 400 walks per cell; seed 131)

**Inputs.** Final-rule yields on the 12 DEV runs (`a6_cert_yield.json`): certified 63,692; A_cand 23,049; tail proxy 118. With one-engine certification and the same removal (contaminants dropped) the counts are 64,507 / 23,406 / 137; the retention figures of A6.4.1 are relative to these. (The pilot's f < 0.05 cell, without the contaminant drop, was 64,507 / 23,527 / 139.) Pool 1,415 projects; clusters by person/lab keys; R = 3; one-engine certification as a sensitivity cell. The simulation does not model symmetric I or the project-wide ladder rule (both remove slightly more); the step-2 re-run uses dry-run yields that include them.

**Evaluability (walk stops at the v1 minimum):**

| Cell | P(v1 minimum) within 100 / 150 / 200 / 300 runs | runs at stop p10 / 50 / 90 | files downloaded p50 (p90) | projects p50 | J p50 |
|---|---|---|---|---|---|
| **final rules, pilot yields** | 0.47 / 0.99 / **1.00** / 1.00 | 80 / **102** / 130 | 113 (144) ≈ 0.12 (0.16) TB | 38 | 37 |
| final rules, yields ×0.5 | 0 / 0.33 / 0.86 / **1.00** | 134 / 163 / 207 | 180 (230) ≈ 0.19 (0.25) TB | 61 | 56 |
| one-engine certification (sensitivity), pilot yields | 0.66 / 1.00 / 1.00 / 1.00 | 74 / 92 / 116 | 103 (130) | 35 | 33 |

- The binding criterion is the capped tail proxy (1.03× its minimum at the stop); J is 1.8–2.8× its minimum. No project was skipped as "cluster full" in the median walk.
- P(v1 minimum) within 300 runs is 1.00 in every cell (the simulation's cap was 2,000 runs; the 400-run hard cap is not separately simulated). **The binding constraint is the walk deadline (16:30 CET), not the run cap:** within 100 runs P = 0.47 (×1) and 0 (×0.5). Hence the throughput-based tag criterion (A6.8.2 step 3).
- A v2_sep full stop (A_cand ≥ 450,000, A4) would need p50 213 runs (×0.5: 418), but it is not used (A6.6.2).

**Gate behaviour at the stop** (A4 levels, intersection–union over both pools; 200 states per cell; MC SE ≈ 0.02 near 0.9, ≈ 0.01 near 0.02). κ = n_tail(v1) / hyperscore proxy. Capped n_tail p50 ≈ 1,250 (κ 1) and ≈ 630 (κ 0.5).

| Scenario (pilot yields) | κ = 1 | κ = 0.5 |
|---|---|---|
| clean: pass / invalid pool / fail | **0.935 / 0.065 / 0** | **0.87 / 0.06 / 0.07** |
| margin 0.60 (τ = 0.10 between clusters): pass / fail | 0.015 / 0.91 | 0.02 / 0.91 |
| margin 0.60, homogeneous: pass | 0.00 | 0.00 |
| leak 0.90 in 25% of clusters: pass | 0.01 | 0.02 |
| residual excess δ = 0.01 (hyperscore 0.51, learned 0.515): pass / invalid | 0.805 / 0.19 | 0.825 / 0.10 |
| **δ = 0.02: pass / invalid** | **0.53 / 0.455** | **0.59 / 0.27** |
| δ = 0.05: invalid | 0.975 | 0.855 |
| G1 clean pass (n_top ≈ 19,800, J ≈ 36); at ρ = 0.05 | 0.99; 0.00 | |
| G3 clean pass, all 10 bins | 1.00 | |

At ×0.5 yields the tail rows are within MC error of these (clean pass 0.925 / 0.88; margin pass 0.00 / 0.01; concentrated-leak pass 0.03 / 0.03).

**D28 quantities (final rules):**

| | Pilot yields | Yields ×0.5 |
|---|---|---|
| P(v1 minimum) within 300 runs (within 100 / 150 / 200) | 1.00 (0.47 / 0.99 / 1.00) | 1.00 (0 / 0.33 / 0.86) |
| P(clean v1 passes), **upper bound** (P(min) × tail × G1 × G3, treated as independent; G2-original, G5, FS and the subset check omitted; deadline not modelled) | **≤ 0.86–0.93** | ≤ 0.87–0.92 |
| P(margin leak caught) | **0.91** | 0.92–0.93 |
| False pass at the margin \| evaluable | 0.015–0.02 (nominal 0.025) | 0.00–0.01 |
| False pass, leak in 25% of clusters \| evaluable | 0.01–0.02 | 0.03 |
| P(clean v1 passes) if δ = 0.02 | 0.53–0.59 (the rest mostly "invalid pool") | 0.58–0.64 |

These assume a valid null (δ = 0). δ cannot be measured without labels; it is the main unknown (R7). The contamination cells model learned = 0.5 + 1.5δ (the learned scorer exploits partial evidence more). The reverse case (hyperscore inflated more than v1, e.g. residual homologs with the wrong RT) would loosen A2's difference clause by δ_h − δ_l; A6's stricter difference clause (A6.6.1) removes that, and the absolute clause UCB ≤ 0.60 is unaffected by it. A leak of +0.10 in v1 on top of any δ_h then faces the margin cells above (pass ≤ 0.02).

---

## A6.11 Residual risks and their direction

| # | Residual | Direction | Mitigation / what happens |
|---|---|---|---|
| R1 | Unremoved homologs (Hamming ≥ 3, no shared 6-mer, f < 0.05) carrying partial true evidence | pro-target, every scorer, but **not necessarily equally** | f_max 0.05; I per project. Leads to **invalid pool** (hyperscore) or **fail** (learned). If hyperscore gains more than v1 (δ_h > δ_l), A2's difference clause would be loosened; A6's stricter difference clause prevents it, and the absolute clause is unaffected |
| R2 | Unidentified co-isolated peptides | pro-target | I removes identified co-isolates' neighbourhoods; single hypothesis |
| R3 | Contaminants (keratins, trypsin, BSA; bovine serum in human cell culture) | pro-target | A6.4.4(iii); certification excludes contaminant-near p |
| R4 | False certification (true peptide still in the DB) | pro-target, strong; possibly more for hyperscore than for v1 (wrong RT) | two-engine agreement (A6.4.1), which removed a 10× tail-enriched 1.5% of A_cand on dev; residual unknown; validity control plus the stricter difference clause |
| R5 | Population shift to well-identified spectra | changes the estimand | reported (as A2.11.1) |
| R6 | Label-dependent removal | **pro-decoy (false pass)** | forbidden by A6.1.2; symmetric I (A6.4.3); swap test (A6.4.6) |
| R7 | **Invalid pool** from a residual target excess δ in the hyperscore tail (pilot condition 5) | not evaluable (completed test) | **At δ ≈ 0.02, "invalid pool" is likely (A6.10); at δ ≈ 0.05 almost certain.** An invalid pool is a completed test: v1's next attempt runs at α/4. The strict screens are adopted now for this reason; nothing is loosened after a run |
| R8 | Concentrated leaks (a fifth to a quarter of clusters) | false pass | simulated rates in A6.10, against nominal 0.025 |
| R9 | The +0.05 tail blind spot of the A1 margins (A2.11.2) | — | unchanged; tail control beyond it rests on entrapment |
| R10 | v2_sep untestable on single-hypothesis nulls (A6.6.2) | — | v2_sep stays "never tested"; a later amendment needs multi-group nulls |
| R11 | Time limit (D31) | label-free stop or "not completed" | A6.8.2 clock; partial outputs reported |
| R12 | Self-host: null copies duplicate their host spectra; run-level features see the twin and the host's sub-threshold PSMs | pro-target (real but unidentified peptides have host PSMs; decoys rarely do) → fail or invalid, not a false pass | **closed rule:** every v1 feature that counts host PSMs, proteins or charges of a candidate peptide is listed by the implementer before the tag and reported with this direction; it is never a reason to drop spectra or change the pipeline |

---

## A6.12 Reporting (whatever the outcome)
`results/2026-10-28_A6_hn1/NOTES.md`, written whatever happens (pass, fail, not evaluable, invalid pool, label-free stop, not tagged, node time limit). It contains:
- the A6.0 first paragraph and attempt history;
- every v1 row on both pools, with each interval component, the validity controls, the cap factors, per project, per cluster, per instrument class, leave-one-cluster-out and noRT;
- v2_sep and ce_s: "not run", with the A6.6.2/A6.6.3 reasons;
- the selection log, label-free counts, swap tests, deviations, timings and data volume;
- the pilot disclosures (A6.13).

**Required wording:**
- **v1 pass:** "passed tail-G2/G1/G3 on fresh human removed-peptide nulls (HN1) under A6, at its second completed test (A1: FAILED, D12; A2: not evaluable before any test, D24; A4: not run, D28); per-row nominal level α/2 conditional on the A1 failure; simulated per-attempt false-pass rates at the margins as in A6.10". In T1 labels this is A5.9.3's "GO, A6: v1 PASS" row.
- **v1 fail:** "FAILED" (upper case); final under P4; all v1 FDR-based claims void (A5.7.1.3).
- **Not evaluable / invalid pool / stopped:** A5.9.3's last row, naming the status.

---

## A6.13 Disclosures
1. **f_max = 0.05 was chosen after seeing label-free tail counts** (pilot §3c: the tail is 2.6× enriched in winners sharing 15–20% of p's ions). It is stricter and symmetric, and no label was seen.
2. **Two-engine certification was chosen after measuring its label-free yield and tail retention** on the 12 DEV runs (A6.4.1). The alternative (project-pooled q ≤ 1e-4) was infeasible on the same data. No label of any null candidate was seen.
3. **The walk seed is derived from the tag commit** because the pilot's simulation used seeds 47+i on the same pool.
4. **The author has seen the pilot's label-free outputs** and the a6 checks. The PI's exposure to them is stated in the adoption note.
5. **The DEV runs used for yields** (D0 and the 10 D8h projects) are excluded from HN1.

---

## A6.14 A5.10.3 floor compliance

| Floor item | A6 |
|---|---|
| 1. A1 margins and n-minima; A2 P1–P5; validity controls; A3 V-rows with tiers; IUT within a scorer | A6.1, A6.6 (all margins unchanged; V-rows mapped, unevaluable rows stay not evaluable) |
| 2. Outer envelope (analytic, project, cluster), B = 10,000, rng 83 | A6.6.1 |
| 3. n_min under the A4 rule, never below provisional; < 10 clusters → not evaluable | A6.6.4 |
| 4. Cluster at least as coarse as A4.3.4 | A6.3.3 (person/lab keys; the source project is the donor, so no extra key) |
| 5. Cap and capped/uncapped IUT | A6.3.5 |
| 6. Label-free stopping rule committed before listing | A6.3.6 |
| 7. Subset check | A6.6.1 (instrument classes, clusters and projects) |
| 8. Reporting | A6.12 |
| 9. Exclusions at least A4.3.3 + T-sets, train/dev, predictor training, earlier null sets, C9 runs | A6.3.2 |
| 10. Label-free evaluability simulation committed before listing | A6.10 (`a6_power`), re-run at step 2 |
| 11. Closed text | yes; PI points (A6.16) are decided at adoption |

Also met: A5.10.4 (T1-blind roles, one open null-set amendment, commissioning prompts committed), A5.1.4 (code placement).

---

## A6.15 Review
Revision 0 is reviewed by a fresh, T1-blind agent (`reviews/amendment_A6_review.md`). The resolutions are listed in A6.17.

## A6.16 Points for the PI (decide at adoption)
Each point has a default written into the text; adoption without comment accepts the default.
1. **Scope: v1 only (A6.2, A6.6.2–A6.6.3).** Default: v2_sep and ce_s are not run; A6 is not a test for either; HN1 is spent for every scorer once v1's step 5 starts. This departs from D30's "gates v1 and v2_sep": under A3's tier rule and the A5.10.3 floor, V6 class 1 makes v2_sep's tier 1 not evaluable on single-hypothesis nulls, so running it could only spend its test 1. v2_sep stays "candidate" with test 1 still at α.
2. **Invalid-pool risk (R7, pilot condition 5).** Default: accept. At δ = 0.02, P(invalid pool) ≈ 0.27–0.46 and P(clean v1 passes) ≈ 0.53–0.59; an invalid pool is a completed test and moves v1's next test to α/4.
3. **Two-engine certification with Sage (A6.4.1).** Default: accept. It costs mzML conversion and one Sage search per project inside the clock, which the throughput criterion measures. Using the comparator binary as a certifier does not touch any comparator result.
4. **Self-host mixture (A6.5.2)** instead of the pilot's "own run with certified spectra replaced". Default: accept.
5. **Stricter changes to the pilot design:** symmetric I (A6.4.3), project-wide ladder removal (A6.4.4(ii)), contaminant neighbourhoods removed from the database instead of A2.2(c)'s winner filter, and the stricter difference clause (A6.6.1). Default: accept. Their yield cost enters the step-2 re-run.
6. **Tag number and the A5 hook (A6.2).** Default: A6 is tagged `prereg-v0.4.5` although `prereg-v0.4.3`/`v0.4.4` may not exist yet; the A5 look learns A6's state only from `a6_status.json` (closed enum; absent = "stopped: node time limit"). The PI confirms the out-of-order tag and asks the A5 implementer to commit the claim-script hook before the A5 freeze.
7. **Throughput-based tag criterion (A6.8.2 step 3).** Default: accept; "not tagged: infeasible" is a label-free stop that lists no project. The PI should expect this outcome to be likely unless admission reaches ≳ 135 runs/h.
8. **n_min re-derivation (A6.6.4).** Default: if the ≥ 4,000-replicate re-derivation cannot finish before the tag, the A4 provisional table is used and the gap is disclosed. The alternative is to treat a missing re-derivation as "not tagged".
9. **Clock and hard stop (A6.8.2).** Default: as written; all A6 work stops at 17:55 CET.

## A6.17 Changes after the review (`reviews/amendment_A6_review.md`)

| Item | Resolution |
|---|---|
| **B1** swap test cannot fail | Accepted. The swap test now runs on every admitted project and recomputes I, R_P, the reduced index and memberships from the raw rank-1 PSMs and the swapped index, freezing only the certified list and the threshold scalars. Added a synthetic-label unit test (decoy-won pair above threshold ∈ I) and a per-project count assertion on \|I\| (A6.4.6) |
| **B2** step 4 reads labels on null candidates | Accepted. Per-mixture thresholds, tail membership and every target/decoy split move to step 5 (`gates_hn1.py`). Step 4 commits label-free winners, scores and features only. The certification module may index labels only at host full-index rows; a test asserts no label read at reduced-index or null-copy rows (A6.4.6.3, A6.5.2–A6.5.3) |
| **B3** clock infeasible; pre-tag criterion ignores the clock; A4 statistics untested | Accepted. The implementer is commissioned at adoption. The tag criterion now uses the measured end-to-end throughput r and the time to the walk deadline (P(v1 minimum within r × window) ≥ 0.8 at ×1 and ≥ 0.5 at ×0.5), otherwise "not tagged: infeasible". Added the `a6_stats` acceptance test (reproduce `a2_stats`, `a4_power.cap_factors`, `c9_power.tail_iut`; coverage of the new interval components), a separate GPU/lock and a 17:55 hard stop. v1 only (M5). The clock is re-timed (tag 15:15, walk 16:30, step 5 17:50) |
| **M1** no hook into the A5 look; tag order | Accepted. Closed-enum `a6_status.json` at a fixed path, read by the A5 claim script; absent = "stopped: node time limit". The out-of-order `prereg-v0.4.5` is a PI confirmation point (A6.2, A6.16.6) |
| **M2** difference clause can pass falsely under differential contamination | Accepted in substance. The difference clause is now max(UCB_diff, UCB_abs − min(p̂_hyperscore, 0.5)) ≤ 0.10 (never more lenient than A2). R1/R4 now state that contamination need not be equal across scorers. The reverse-direction case is argued in A6.10: with the stricter clause, a +0.10 learned leak faces the absolute clause as in the margin cells; it was not re-simulated, because of the time limit |
| **M3** power summary optimistic; 400-run figure inferred | Accepted. "P(clean v1 passes)" is labelled an upper bound with the omitted rows listed; the 300-run figure is quoted; the deadline, not the run cap, is named as binding, with P(min) at 100/150/200 runs |
| **M4** subset check narrower than floor item 7 | Accepted. Each project with n_tail ≥ 200 is added; A6.14 updated |
| **M5** close v2_sep and ce_s; declare HN1 spent | Accepted. A6 gates v1 only; v2_sep and ce_s are not run (not tests); HN1 is spent for every scorer once v1's step 5 starts (A6.0, A6.2, A6.6.2–A6.6.3, A6.7) |
| minor 1 (counts) | Accepted: A6.10 now gives both one-engine and pilot counts |
| minor 2 (completed-test rule) | Accepted: generalised to any stop during step 5 (A6.7) |
| minor 3 (R12) | Accepted: closed rule in R12 |
| minor 4 (reuse targets) | Accepted: exact files and functions named in A6.8.3 |
| minor 5 (cluster-full not burnt) | Accepted (A6.3.7) |

Not changed: the swap-test idea of testing I's construction was also raised by the author during the review (an interim edit), and is superseded by B1's fuller version.

---

## Appendix A: constants for `config/a6.yaml`

| Item | Value |
|---|---|
| Walk seed S | int(sha[:8], 16) of the `prereg-v0.4.5` commit |
| File order | `random.Random(S + accession digits)` |
| Runs per project | 3 |
| RAW size | 0.1–2.0 GB |
| Projects per cluster | ≤ 5 |
| Cap share | 10%, min of the proxy, A_cand and \|A\| factors; hash seed 47 |
| Stop | v1 minimum: J ≥ 20, J_eff ≥ 10 (proxy and A_cand), capped proxy ≥ 1,200, capped A_cand ≥ max(120,000, 10 × G1 n_min); hard cap 400 runs; walk deadline 16:30 CET |
| Certification | hyperscore PSM q ≤ 0.001, peptidoform q ≤ 0.01, tie-break seed 11; Sage rank-1 target, spectrum_q ≤ 0.001, same I/L sequence and charge |
| Removal | HomologFilter over I (pairs of either label at the run's 1% peptidoform threshold) + isobaric variants; f ≥ 0.05 project-wide; contaminant FASTA ∪ Panel B pairs, A1.3-near to their non-human proteins |
| Levels | tail 0.975 one-sided; G1 0.995 two-sided; TOST 0.95 two-sided; validity controls unchanged |
| Bootstrap | B = 10,000, `default_rng(83)` |
| n_min (floor) | G1 11,100; G3 2,000; tail 200; clustered rows J ≥ 10 |
| Tail-G2 difference clause | max(UCB_diff, UCB_abs − min(p̂_hyperscore, 0.5)) ≤ 0.10 |
| Clock (CET) | tag by 15:15 (throughput criterion); walk ends by 16:30; step 4 by 17:10; step 5 by 17:50; hard stop 17:55; NOTES by 18:00 |
| Status file | `results/2026-10-28_A6_hn1/a6_status.json`, closed enum (A6.2) |
