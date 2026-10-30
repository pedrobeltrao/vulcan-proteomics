# Amendment A2: implementation notes (adoption precondition 2)

Text: `reviews/amendment_A2_draft.md` (to become DESIGN §6.0b). Config: `config/a2.yaml`. Written 2026-10-27, before any
DN2 file was listed, downloaded or searched. Nothing below was tuned on any null data. The rule used throughout is
A2.9.5's: **where the draft is open, the stricter reading wins.** Each choice is also fixed in the docstring of the
module that implements it.

## What exists

| A2 item | Code | Status |
|---|---|---|
| A2.9.1–A2.9.3 DN2 selector | `scripts/a2/select_dn2.py` | Code only, never run. Refuses to run without tag `prereg-v0.4.2` and an unmodified A2 tree. Dry-run on a synthetic table: `tests/a2/test_select_dryrun.py` |
| Panels (Appendix) | `vulcan/gates/a2_panels.py`, `scripts/a2/resolve_panels.py` → `config/a2_panels.tsv` | Resolved 2026-10-27 (UniProt REST and GPM FTP only): 177 entries. FASTA sha256 in config |
| A2.2(a) certification and A2.3 step 3 marker screen | `scripts/a2/own_search.py`, `vulcan/gates/a2_admission.py` | Pure logic tested. GPU path not run |
| A2.2(b)–(e) admission, label-free counts, P1 swap test | `scripts/a2/admit.py`, `vulcan/gates/a2_identity.py` | P1 tested on synthetic T/D lists (`tests/a2/test_admit_swap.py`). The real-index swap test (`admit.py --swap-test`) runs at step 3 |
| Tail proxy threshold | `scripts/a2/host_threshold.py` | Reads labels of the two D0 *host* runs only and stores 2 scalars |
| A2.4 FS nulls and A2.2(e) hypotheses | `vulcan/gates/a2_fs.py`, `scripts/a2/build_mixtures.py` | Pure parts tested (`tests/a2/test_fs.py`) |
| A2.5 / A2.6 / A2.8 gates | `vulcan/gates/a2_stats.py` (new functions; `vulcan/gates/stats.py` is unchanged) | Unit tests and MC checks |
| Step-5 driver | `scripts/a2/gates_a2.py` | Evaluation layer tested on synthetic pools. Its finalize→gates order is enforced by a committed-hash guard |
| MC checks against `reviews/checks/a2_power.*` | `scripts/a2/power_checks.py` → `results/2026-10-27_A2_implementation/power_checks.{json,out}` | see below |
| A2.9.4 "no per-label counts before step 5" | `tests/a2/test_no_label_counts.py` | Static (AST/grep) and runtime checks |

**Not executed on any real data:** own search, admit, build_mixtures, gates_a2, host_threshold. The GPU wiring follows
`scripts/m3e/{g2mix,gates_a1,g3}.py`, but it has not been exercised. The first step-3 run must start with
`host_threshold.py` and must include `admit.py --swap-test` on the first admitted project.

## Ambiguity resolutions (stricter reading) and fixed degrees of freedom (A2.11 item 5)

### Panels and contaminants
1. **Panel resolution.**
   - Query: `gene_exact:<symbol> AND organism_id:<taxid> AND reviewed:true`.
   - *Every* reviewed hit is kept, including all isoforms. A member with no hit falls back to the entry name, e.g. `OVAL_CHICK` (chicken OVAL has the gene name SERPINB14).
   - An unresolved member is a hard error.
   - Several members resolve to more than one entry (e.g. FLG 22, DMKN 16, ALB 9, bovine CSN3 2). They are kept: more panel sequence means more removal and more detection.
2. **"cRAP"** = GPM cRAP ∪ Frankenfield 2022 (A1's contaminant list). The union is the stricter reading.
3. **Panel K** (keratins, which never trigger) is the Appendix list of 10 (KRT1/2/5/6A/9/10/14/16/17/77), not the shorter list in A2.3. Keratins *do* count as contaminants for certification and for the A1.3 filters.

### Certification and marker screen
4. **Certification (A2.2a).**
   - The estimator runs per run, not pooled; the per-run +1 is more conservative.
   - The peptide-level unit is the peptidoform, as in production. PSM-level q is TDC over the run's rank-1 PSMs.
   - "Own-proteome target" means the peptide maps to at least one `OWN|` entry and is not A1.3-near any contaminant or panel entry, keratins included.
5. **Marker screen (A2.3 step 3).**
   - Pooling is the union of accepted peptides over the project's runs.
   - A peptide shared by several panel members counts for each of them.
   - Peptides A1.3-near the own proteome are discarded.
   - A member is detected with ≥ 2 distinct I/L-collapsed peptides. The project is excluded at ≥ 2 detected H/B members.

### Identity screen (b), family screen (d), A1.3 filter (c), single precursor (e)
6. **f(c) details (A2.2b).**
   - **Definition.** f(c) is computed exactly in numpy. The draft's "K2 kernel" route is not used, because K0 (min_peaks 15, 150–2000 m/z, top-150 peaks) would drop ions of short peptides and *lower* f.
   - **Ions of p:** b2…b(L−1) and y2…y(L−1), z = 1. **Ions of c:** *all* b1…b(Lc−1) and y1…y(Lc−1).
   - **Matching:** a b ion may match a y ion. No m/z range is applied. The tolerance is ±20 ppm of the p ion.
   - **Masses:** Sage residue masses, static CAM, Met-ox from the mask bit. Index windows use the index's own Sage f32 `mass32`.
   - **W:** the certified mono m/z, ±10 ppm, isotope errors {0, 1, 2}.
7. **Hamming clause (A2.2b).**
   - "T∪D peptide" means the gate index's peptide entries (stripped sequences), as written, not arbitrary protein windows.
   - p is near if any isobaric variant of p is within Hamming ≤ 2 (I = L) of an equal-length T∪D peptide. Variants are every subset of non-overlapping rewrites N→GG, Q→AG, GG→N and AG→Q, so mixed cases are covered.
   - Above 4,096 variants, p is treated as near.
8. **A1.3 filter (c) is applied twice.**
   - **(c-pre), in admit:** on the exhaustive hyperscore top candidate over T∪D and its pair partner (label-free).
   - **(c-final), in `gates_a2.py --finalize`:** on the matched pair of *every* compared pipeline, before any statistic.
   - (c-final) is also applied to FS spectra, which the draft does not require; this is stricter.
   - The A1.3 reference set is source proteome ∪ cRAP ∪ Frankenfield ∪ Panel H (incl. keratins) ∪ Panel B.
9. **(d) family screen.** Ties are removed (≥). The max is taken over every candidate in W, label-free.
10. **(e) single precursor.**
    - In the gate search, A and FS spectra carry exactly one hypothesis: certified mono m/z, certified z, isotope error 0.
    - The screens (b) and (d) use the wider W (isotope set 0–2). The screened candidate set therefore contains the searched one.
    - The other options were rejected: keeping the isotope set adds lighter-mass hypotheses that are not "the certified precursor".

### FS
11. **FS (A2.4).**
    - Donors are *all* A spectra.
    - Each peak is shifted as float32(m/z) + float32(2/7).
    - FS spectra share the null permutation with A and B.
    - **G5 on FS counts as evaluable when there is ≥ 1 target and ≥ 1 decoy winner.** This is the strict reading of "an alarm on any of these blocks". At tiny n a chance alarm (AUC > 0.55) blocks. **The PI may want to replace this with a minimum n before the tag.** Any such change must be dated and made before the tag.
    - Tail-G2 and G1 on FS are evaluable only when their own n_min and validity conditions hold.

### Selection (A2.9)
12. **Exclusions.**
    - The run-time union of `data/manifests/*.tsv` (except `t4_frame.tsv`) is combined with `config/a2_excluded_accessions.tsv`, the 300 accessions at adoption, sha256 recorded.
    - The adoption list is needed because three of those manifests (`d7b`, `d10cf`, `e9_d8`) are untracked today and are absent from the tag.
    - Non-PXD frame rows are dropped.
13. **ALLOW details.**
    - *Bacillus* excludes the whole *B. cereus* group (cereus, anthracis, thuringiensis, weihenstephanensis, mycoides, pseudomycoides, cytotoxicus, toyonensis, wiedmannii, paranthracis, tropicus).
    - *Sulfolobus*/*Saccharolobus* share one genus slot.
    - Species-restricted genera (*Corynebacterium glutamicum*, *Propionibacterium freudenreichii*, *Aspergillus niger/nidulans*, *Pseudomonas putida*) are matched on the second word.
14. **Genus slot.** Only an *admitted* project takes its genus slot, as in A1 (the slot was taken on selection). A project excluded by the protocol check, file count, proteome or marker screen does not block its genus. This is A1's semantics, not a strictness choice.
15. **Order inside the walk.**
    - The metadata checks (genus slot, protocol text, SDRF) run *before* any file listing, so a project that would be excluded is never listed.
    - Then: list files → pick → reference proteome → download → own search + marker screen → admit → stopping rule.
16. **Protocol check (A2.3 step 2).**
    - The regex is applied to `description`, `sampleProcessingProtocol` and `dataProcessingProtocol`, plus title and keywords (stricter).
    - SDRF organism part must be exactly "whole organism" or "cell culture"; any other value, "not applicable" included, excludes.
    - A missing record or SDRF never excludes.
17. **Reference proteome.** Only a UniProt *reference* proteome (`proteome_type:1`) is accepted. The first hit matching the full frame organism name wins, else genus + species. If there is none, the project is skipped (logged).
18. **BAD_KW quirks, kept verbatim.**
    - A1's list has no word boundaries and is case-insensitive, so `HLA` matches "C**hla**mydomonas", `DIA\b` matches "me**dia**" and `UPS` matches "gro**ups**". A2 says "A1's list plus …", so it is reused unchanged.
    - Consequence: a *Chlamydomonas* project that names its genus in the title or keywords is dropped at the frame filter.
    - This only removes projects. Pinned by a test so that nobody "fixes" it silently.
19. **Steps 2 and 3 interleave.**
    - The stopping rule needs admit's label-free counts after each project, so the selector calls own search, marker screen and admit inside the walk.
    - A2.9.4's "each step committed before the next" therefore applies as follows:
      - the selection log, the manifest, the RAW sha256 values, the admission membership hashes, the counts and the swap-test result are committed together at the end of the walk;
      - that commit happens before step 4.
    - No gate statistic exists at that point.
20. **Tail proxy threshold.**
    - Computed from a production hyperscore search of each D0 host run alone, using that run's own 1% threshold.
    - The proxy uses the *max* over the two hosts. This counts fewer spectra, so the walk can only run longer.
21. **Cap at 16.** If 16 admitted projects are reached first, the walk stops and the gates run (P4).

### Gates
22. **Tail-G2 (A2.5).**
    - **Cluster bootstrap:**
      - cluster = project;
      - one (B, J) index matrix with B = 10,000, `default_rng(83)`, used for both the absolute and the paired difference statistic;
      - type-7 quantile.
    - **Degenerate resamples:** an empty learned tail counts as 1.0, and an empty learned or hyperscore tail counts as a difference of +1.0. J < 2 gives a bound of +∞.
    - **Validity control:** two-sided 95% Wilson CI as a closed interval. An empty hyperscore tail is invalid.
    - **Precedence:** n_tail < 200 → not evaluable; invalid pool → not evaluable (invalid pool) *even if a clause fails*, because A2.1 says the control "can only turn a verdict into not evaluable". Clause outcomes are still reported.
23. **G1 (A2.6).**
    - n_top is the size of the learned top decile (score ≥ 90% quantile, as A1).
    - The control is hyperscore's ρ on its own top decile of its own winners, on the same pool. It must contain 0 (closed interval). If it is undefined (< 10 spectra), the pool is invalid.
24. **G3 (A2.8).**
    - Both labels must reach 2,000 in each top-half decile.
    - Bottom-half deciles are merged upward from decile 0 until *both* labels reach 2,000.
    - A leftover group containing decile 4 cannot merge upward. It is **not evaluable**: no downward merging and no testing of an under-sized group.
    - The pass condition is lo ≥ −0.01 and hi ≤ +0.01 (closed).
25. **Pools.** Gated A is the admitted A after (c-final). B is the B candidates after (c-final). The A1 offsets come from the D0 hosts in the first two mixtures and are reported only.
26. **A2.9.1 "`bench/sealed_run.py`'s T4 assertion".** No such assertion exists today. Instead `dn2.tsv` was added as an optional manifest to `scripts/data/select_t4.py` `DEV_MANIFESTS`, so `assert_disjoint` (and `--check`) cover DN2 once it exists. The committed T4 selection still passes the check.
27. **Frozen classifier.** `/tmp/vulcan_m3e/seqclf.pkl` was hashed today, before any DN2 listing. Its sha256 (4ec134eb…) is in the config, and the selector and the gate run both check it. The copy is node-local; if it is lost, A2.6 item 4 applies.

## MC checks (`results/2026-10-27_A2_implementation/power_checks.{json,out}`; seed 131)
- **Agreement with the draft:** the implemented code agrees with `reviews/checks/a2_power.json` in all 55 comparisons, within 4 combined MC SE:
  - S3 coverage of max(CP, cluster): 0.958 / 0.979 / 0.994. These are identical to the draft because the same RNG stream gives the same resamples, which also confirms that `cluster_bootstrap_tail` reproduces the draft's estimator exactly.
  - S5 G1 pass probability (e.g. n_top 7,200: 0.913 at ρ = 0 vs 0.905 analytic; 0.006 at ρ = 0.05).
  - S6 G3 TOST, every cell with n ≥ 1,000 (e.g. 2,000/bin, recall 0.995: 0.879 vs 0.8945).
  - S7 FS δ gaps.
- **A2 is never looser than A1:** the implemented tail-G2 never passes where the A1 rule fails on the same data.
- **Full A2 tail-G2 rule** (independent tails, 10 homogeneous projects; 500 replicates per cell):

  | n_tail | true fraction 0.50 | 0.55 | 0.60 |
  |---|---|---|---|
  | 400 | 0.83 | 0.31 | 0.004 |
  | 750 | 0.93 | 0.48 | 0.012 |
  | 1500 | 0.95 | 0.77 | 0.010 |

  - Each cell is the A2 pass rate.
  - For comparison, the A1 rule on the same data at true fraction 0.50 passes 0.88 / 0.99 / 1.00.
  - **Validity control:** it alone sends a clean pool to "not evaluable (invalid pool)" in 4–7% of replicates. That is the 95% Wilson CI by construction, i.e. A2.11 item 4 quantified. The rest of the gap to A1 is the cluster max rule and the independent-tails design.
- **Fires at +0.10:** at every n_tail ≥ 400 a true fraction of 0.60 passes ≤ 1.2%; the rest are fail or invalid.
