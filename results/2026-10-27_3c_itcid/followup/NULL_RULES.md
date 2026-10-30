# IT-CID fresh null set "DNIT": admission rules (A2 adapted to the IT-CID instrument class)

- **Written and committed 2026-10-27, before any DNIT file is listed (PRIDE file API), downloaded or searched.**
  While drafting, only the local metadata frame `data/manifests/t4_frame.tsv` was read (organism / instrument / label
  columns, to enumerate the organism class below and count candidates), as A2 did for its ALLOW list. No PRIDE file
  listing, no spectrum, no search.
- **Purpose:** the A1 gates still pending for the frozen IT-CID scorer `itcid_linear_pred`
  (sha256 667ce39cb653026baed40bb18b25f5cbe196b192e2d5a7f2f8e7aa3618c23b89): **tail-G2, G1, G3**. G2-original and G5 on
  A ∪ B are also computed (A2.7) and reported. Arm M is out of scope (no IT-CID yeast dev run).
- **Normative basis:** DESIGN §6.0 (A1) and §6.0b (A2, `reviews/amendment_A2_draft.md`, implementation `config/a2.yaml`,
  `scripts/a2/IMPLEMENTATION_NOTES.md`). Everything below that is not listed as an IT-CID adaptation is A2 **verbatim**,
  including every A1 margin and every A2 minimum-n / validity rule. Where this text is open, the stricter reading wins.
- **Machine-readable:** `config/a2_itcid.yaml` (same commit). Code under `scripts/a2_itcid/` reads it; tagged A2 code
  (`scripts/a2`, `vulcan/gates/a2_*`, `config/a2.yaml`) is imported, never edited.

## P. Principles (A2.0, unchanged)
P1 label-swap invariance of every admission decision (tested on the real gate index, as A2); P2 admission from metadata,
own-species search and label-collapsed quantities only, membership hashed and committed before any gate statistic;
P3 no relaxation; **P4 one shot**: DNIT is drawn once by the seeded rule below; a fail is final for `itcid_linear_pred`
under this protocol; nothing is re-drawn, extended or pruned after the gate run; P5 reported whatever the outcome.

## 1. Exclusions (mechanical, checked by code before the walk)
1. Every accession in the `accession` column of every `data/manifests/*.tsv` except `t4_frame.tsv` at walk time
   (dev, sealed, tier S, T4 selection, all training sets incl. `train_itcid*` and its 30+30 walked candidates,
   `null_itcid*` (DN3c and its walked candidates), predictor training, DN, D8, ...).
2. **DN2 isolation** (PI rule from 3c, applied unchanged through `scripts/m3c/select_train_itcid.dn2_guard`): every
   accession the A2 DN2 walk has seen (`data/manifests/dn2*`, `results/2026-10-27_A2_dn2/{selection_log.tsv,
   walk_events.jsonl}`), **and** every frame row that passes the DN2 frame filter (A2.3 step 1 as coded in
   `scripts/a2/select_dn2.py`; the set of projects the DN2 walk *could* reach), **and** any organism whose genus is on
   the DN2 ALLOW list (`config/a2.yaml` `allow`). So DNIT shares neither a project nor a genus with DN2.
3. Organisms of any target, entrapment or previous-null database: *Homo sapiens* (gate DB), *Mus musculus*
   (3c training), *Arabidopsis* (arm C), *Pyrococcus furiosus* (G0 / D1), *Saccharomyces cerevisiae* (DN3c, D7, arm M),
   *Escherichia coli*, *Trypanosoma brucei*, *Staphylococcus aureus*, *Synechococcus*, *Bacillus cereus* (A2.9.1 list).
4. ProteoBench / non-PXD rows.

## 2. Organism class and genus list (fixed now)
Class (A2.9.2, unchanged): free-living, non-host-associated, axenically culturable unicellular organisms (incl.
filamentous fungi, as A2's *Neurospora*, *Trichoderma*, *Abortiporus*) grown in animal-free media. Excluded: human or
animal commensals and pathogens, parasites, symbionts, plant pathogens, plants, animals, enrichment / mixed /
environmental cultures. Matched on genus (species restriction where the genus contains pathogens):

*Acidithiobacillus, Alteromonas, Blastobotrys, Blastococcus, Dehalococcoides, Desulfovibrio, Desulfurella,
Dinoroseobacter, Fragilariopsis, Geodermatophilus, Kluyveromyces, Leptospirillum, Magnetospirillum, Modestobacter,
Nannochloropsis, Paenarthrobacter, Podospora, Prochlorococcus, Prymnesium, Pseudoalteromonas, Psychroflexus, Rhodococcus,
Saccharopolyspora, Schizophyllum, Sulfobacillus, Sulfurimonas, Thalassolituus, Trichococcus, Vibrio (natriegens only)*.

Metadata-only count while drafting (frame rows passing §3 with these genera): ~33 projects in 29 genera (before
file listing, protocol text and every later step).

## 3. Frame filter (metadata)
Single organism; DDA; Thermo; label LFQ; enrichment none; `eligible` = 1; instrument field exactly one of
`LTQ Orbitrap`, `LTQ Orbitrap XL`, `LTQ Orbitrap Velos`, `LTQ Orbitrap Velos Pro`, `Orbitrap Velos Pro`,
`LTQ Orbitrap Elite`, `Orbitrap Fusion`, `Orbitrap Fusion Lumos` (IT-CID-capable class; the analyzer is verified on the
spectra in §6). Title + keywords + experiment types must not match A2's BAD_KW (`bad_kw_a1` + `bad_kw_a2_add`, verbatim).

## 4. Protocol text / SDRF (A2.3 step 2, verbatim)
PRIDE project record (`description`, `sampleProcessingProtocol`, `dataProcessingProtocol`, title, keywords) and SDRF if
present; exclude on BAD_KW ∪ `protocol_kw_add`; SDRF with > 1 organism or an organism part other than whole organism /
cell culture excludes. Missing text or SDRF is not a reason to exclude.

## 5. Seeded walk, files, download
- Order: candidates sorted by accession, shuffled with `random.Random(47)` (SEEDS[2]; DN used 11, DN2 and DN3c 23).
- One project per genus. Per project: ≥ 3 usable RAW files (`.raw`, **0.15–1.5 GB** (IT-CID adaptation: the 3c
  IT-CID size range; ion-trap MS2 files are small), name not matching A2's `bad_file_regex`); 3 files picked with
  `random.Random(47 + int(accession digits))`.
- Reference proteome: UniProt `(organism_name:"<genus species>") AND (reference:true)`, first hit whose organism starts
  with the frame organism's full name, else with genus + species (the D16-corrected A2 rule); none -> project skipped.
- **Download verified by PRIDE SHA-1** (the `checksum` of the PRIDE v3 file record) and size; a file without a PRIDE
  checksum, or failing SHA-1 after 3 full attempts, makes the project skipped (logged); sha256 recorded in
  `data/manifests/dnit.sha256` and DATA.md (role dev, set DNIT).
- Steps run in A2's order: genus slot -> frame -> protocol/SDRF -> list files -> pick -> proteome -> download ->
  §6 class check -> own search + marker screen -> admission -> stopping rule.

## 6. IT-CID class check on the spectra (new, label-free; LESSONS 7)
Each of the 3 runs must have ≥ 95% of MS2 with analyzer IT and activation CID, and no MS3 (3c rule). A failing run
rejects the project (replaced by the next in the walk).

## 7. Own-species search, certification, marker screen (A2.2a / A2.3 step 3)
Database and certification as A2 (own reference proteome ∪ GPM cRAP ∪ Frankenfield ∪ Panel H ∪ Panel B, pseudo-reverse
decoys, hyperscore, per-run paired double competition, q_PSM ≤ 0.001 and q_peptide ≤ 0.01, own-target only, not
A1.3-near any contaminant / panel entry). **IT-CID adaptation:** the search uses the 3c production IT-CID parameters
(`scripts/m3c/common.py` SEARCH: fragment ±0.5 Da (`fragment_unit` "da"), precursor ±10 ppm, isotope errors 0..2,
report_psms 5, K1 bins 200 ppm, otherwise E2 SearchParams defaults; index `default_params()`). Marker screen verbatim
(≥ 2 non-keratin Panel H/B members with ≥ 2 distinct peptides -> project excluded).

## 8. Spectrum admission (A2.2 b–e)
- Gate DB = the production DB of the IT-CID hosts: human_contam prod (`/tmp/vulcan_e2/index/human_contam_prod`, the
  index 3c searched D10cf with).
- (b) identity screen: f(c) as A2 (`vulcan.gates.a2_identity`), **IT-CID adaptation: ions coincide within ±0.5 Da**
  (the IT-CID fragment tolerance) instead of ±20 ppm; threshold 0.20 and the Hamming clause unchanged. (A wider
  tolerance only removes more spectra.)
- (d) family screen: exhaustive hyperscore over (T∪D) ∩ W with the IT-CID search parameters; removed iff max ≥
  hyperscore(p) (ties removed).
- (c) A1.3 filter (pre and final) and (e) single precursor hypothesis: verbatim.
- P1 swap test: on the first admitted project, gate index rebuilt with `swap_labels=True`; membership sha256 must be
  identical, else the procedure stops.

## 9. Mixtures and nulls
- Hosts: the 29 D10cf fractions (IT-CID label-free human; the hosts of 3c's G2), host of mixture m = fraction m mod 29.
  TMT-class nulls are not covered (stated limitation; D4 TMT is not a host).
- Mixture recipe as A2 (20,000 nulls per mixture, A1 cut rule, permutation `default_rng(47)`, RT quantile-mapped onto
  the host, each null used once, single hypothesis for A spectra, production hypotheses otherwise).
- A1-recipe offsets of the host (2 copies, first two mixtures): generated, **reported only** (A2.4 rule 1).
- **FS nulls: not generated (IT-CID adaptation, declared before data).** A fragment shift cannot make IT-CID fragment
  evidence chance-only: at ±0.5 Da matching every shift δ lands within tolerance of an integer-spaced (isotope / ±1 Da)
  relation of a real fragment for some ion. FS is therefore "not evaluable", which under A2.4 does not block; the
  prior-dominated regime is gated on stratum A, as A2.4 states.
- Tail proxy threshold (stopping rule only): max over the 29 D10cf hosts of the hosts' own hyperscore 1% threshold
  (production IT-CID search, host alone), as A2's "max over hosts".

## 10. Gates (step 5; A2.5–A2.8 unchanged)
- Scorer under test: frozen `itcid_linear_pred`; reference: hyperscore (IT-CID production search). Also reported: the
  same scorer with RT features neutralized (as A2's `gbm_noRT`).
- Features: `vulcan.rescore.itcid` with every Prosit prediction read from the **pinned prediction store**
  (`data/cache/koina_pinned`, item 1 of this follow-up); the store is filled (label-free) before the gate run and
  frozen during it.
- tail-G2 on A (n_tail ≥ 200; UCB_abs ≤ 0.60; UCB_diff ≤ 0.10; cluster bootstrap 10,000 / seed 83; hyperscore Wilson
  validity); G1 on A (n_top ≥ 7,200; ±0.05 at 99%; hyperscore control); G3 TOST on A and A ∪ B (≥ 2,000 per label in each
  top-half decile; exhaustive vs production top-K, all mixtures); G2-original and G5 on A ∪ B, G5 also on A.
- G1 sequence classifier: the frozen A2 classifier (`/tmp/vulcan_m3e/seqclf.pkl`, sha256 4ec134eb…). Its `rt_pred`
  input is, by its training definition, **AlphaPeptDeep's** RT prediction (not Prosit iRT); it is computed for the
  winners with the M3 AlphaPeptDeep predictor (`vulcan.predict.apd`, label-blind sequence property). Never retrained.
- Overall verdict: as A2 (`overall()`): tail-G2(A), G1(A), G3(A), G3(A∪B), G2-original, every G5 must pass; not
  evaluable blocks.

## 11. Power / stopping rule (A2.9.3 counts; label-free from admission)
Walk until **all** hold: ≥ 10 admitted projects; |A| ≥ 120,000; A spectra with an exhaustive candidate ≥ 72,000;
tail proxy ≥ 400; cap **16** admitted projects; or the candidate list is exhausted.
**If the walk ends (cap or exhaustion) with any of these four unmet, power is declared insufficient: no mixture is
built and no gate is run; tail-G2/G1/G3 for `itcid_linear_pred` stay "not evaluable" and IT-CID FDR claims stay
blocked.** (Stricter than A2, which would run the gates anyway; this avoids spending the one-shot set on an
underpowered run.)

## 12. Order of operations (each step committed before the next)
1. This file + `config/a2_itcid.yaml` committed (before any listing).
2. Walk + downloads + class check + own search + marker screen + admission; commit `dnit.tsv`, `dnit.sha256`, DATA.md,
   per-project label-free counts, membership sha256, P1 swap-test result.
3. Power check (§11). Stop here if insufficient.
4. Mixtures, search, features (pinned store fill), final (c); commit the final membership record.
5. Gates once; results reported in `results/2026-10-27_3c_itcid/followup/` whatever the outcome.
No code under `scripts/a2_itcid/` prints or stores a per-label count before step 5. Steps 2–4 and step 5 are
separate commits; the step-5 code is committed with step 4 and not edited afterwards.

## Dated notes (each written before the step it affects)
1. **2026-10-27 23:2x, before the walk (no DNIT file listed yet).** Tail-proxy threshold (§9, stopping rule only):
   5 of the 29 D10cf hosts accept no target at 1% when searched alone (sparse SEC fractions), so their 1% threshold
   is undefined (`tail_threshold` returns +inf). Taken literally, "max over hosts" would be +inf and the tail proxy
   identically 0, so the stopping rule could never be met and power would always be declared insufficient. The max is
   therefore taken over the 24 hosts with a defined threshold: **33.142** (range 28.35–33.14). Only host-run labels
   were read (as `scripts/a2/host_threshold.py` does for D0); this affects only the label-free stopping proxy.
2. **2026-10-27 22:52, walk run 1 stopped after 6 decisions (no RAW file completed, nothing ingested or searched).**
   §5 said "a file without a PRIDE checksum (in the PRIDE v3 file record) makes the project skipped". Run 1 skipped
   PXD001518 because its 2014 file records carry no `checksum` field. PRIDE also publishes its SHA-1 sums in the
   project's `checksum.txt` (the source the A2 DN2 operator wrapper already uses as fallback). Both are PRIDE SHA-1, so
   the narrower rule was an implementation slip, not a protocol choice. Correction: expected SHA-1 = API `checksum`,
   else the project's `checksum.txt`; a file with neither is still unverifiable (project skipped). The walk is
   **re-run once from the start** (same seed 47, same order, `DNIT_WALK_RUN=2`); run 1's events stay in
   `walk_events.jsonl`. Only metadata had been seen (6 projects: 4 protocol/SDRF rejections, 1 checksum skip,
   1 partial download); this cannot depend on any gate quantity.
3. **2026-10-27 23:25, after the first admitted project (label-free counts only; no mixture, no gate quantity).**
   PXD014351 (*Dinoroseobacter shibae*): 76,002 spectra, 26,431 certified, of which the (b) identity screen at
   ±0.5 Da removed 23,603 (89%) -> n_A = 461, tail proxy 0. At unit fragment tolerance a ≥ 20% b/y coincidence with
   some T∪D candidate in the precursor window is common, so stratum A is ~50x smaller per project than A2 assumed.
   Futility rule added now (it can only end the walk earlier with "power insufficient"; it cannot admit, drop or
   re-weight data, and no gate is run either way): **after ≥ 3 admitted projects, stop and declare power insufficient
   if the projection at the cap, (cumulative n_A / n_admitted) × 16, is below 60,000 (half of the 120,000 minimum)**.
   Otherwise the walk continues under §11.
4. **2026-10-27 23:56, outcome of the walk (run 2): power insufficient, walk stopped by the note-3 futility rule.**
   28 decisions: 15 protocol/SDRF rejections, 7 with < 3 usable RAW, 3 without a PRIDE SHA-1 (API or checksum.txt),
   3 admitted (all runs 100% IT-CID, marker screen passed): PXD014351 *Dinoroseobacter shibae* n_A 461 (26,431
   certified), PXD008499 *Nannochloropsis gaditana* n_A 24 (7,475), PXD011268 *Trichococcus patagoniensis* n_A 96
   (2,991). Totals: |A| = 581, A with candidate 580, **tail proxy 0**; projection at the cap 3,099 < 60,000. The
   (b) identity screen at ±0.5 Da removed 86–89% of certified spectra in every project. **No mixture was built and no
   gate was run; tail-G2 / G1 / G3 for `itcid_linear_pred` remain "not evaluable"** and IT-CID FDR claims stay
   blocked. DNIT is spent under P4 (a later attempt needs a new amendment and must report this one). The P1 swap test
   was not run (no stratum A is used). Manifests: `data/manifests/dnit.tsv`, `dnit.sha256`, `dnit_select.log.json`.
