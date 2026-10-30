# Domain review of DESIGN.md v0.1

Reviewer: domain (search engines, DIA, XL/LiP/CF-MS, phospho). Inputs: PROMPT.md, DESIGN.md, lit/02, lit/04 (lit/01, lit/03 skimmed). 2026-10-26.

**Verdict.** The statistical scaffolding is unusually good: pre-registration, paired TDC, Group-walk, entrapment, RESET and the G1–G3 gates. The front end that decides sensitivity is close to absent. Against Sage and MSFragger on real PRIDE data, sensitivity comes as much from five front-end steps as from the scorer:
- mass calibration;
- charge and isotope handling;
- MS1 features;
- ion models per activation type;
- parameter inference.

The proposed store schema makes several of these *impossible*. The scope is also several tool-years wide for one GPU and four cores.

Severity: **[B]** blocker, **[M]** major, **[m]** minor.

## 1. Blockers

**[B1] The store discards what calibration, chimeric search, quant and XL need.** L0 keeps MS2 only and applies top-150, precursor removal and deisotoping irreversibly at ingest. Without MS1 you cannot do any of these:
- recalibrate precursors;
- re-determine the monoisotopic peak or charge;
- detect co-isolated precursors (DDA+ works from MS1 envelopes);
- compute TMT isolation purity;
- run LFQ, MBR or DIA MS1 scoring.

The schema also lacks:
- activation type (HCD/CID/ETD/EThcD) and supplemental energy;
- analyzer (OT/IT/TOF) and resolution;
- master-scan linkage, needed for DSSO MS2→MS3 and SPS-MS3 TMT;
- SPS ions, FAIMS CV, and the scan filter string.

**Fix:**
- Store all MS1 and MS2 centroids together with these fields.
- Make top-N, deisotoping and precursor removal search-time *views* that depend on the analyzer and activation, and can be ablated.

**[B2] There is no calibration or parameter-inference stage, and that is the PRIDE problem.** Repository reprocessing has to infer, per run:
- precursor and fragment mass error and its drift;
- tolerances;
- analyzer: ion-trap CID needs ~0.4 Da bins, and a high-resolution fragment index fails silently on it;
- activation type;
- labels (TMT/iTRAQ/SILAC/dimethyl), enrichment and enzyme;
- common artefacts (protein N-terminal acetylation, pyro-Glu, carbamylation).

PRIDE metadata and SDRF are sparse or wrong for most projects. MSFragger's calibration and parameter optimisation measurably contribute to its sensitivity.

**Fix:** add an "L1.5 pilot pass" per run. It runs a fast, wide-tolerance search of a 5–10% subsample and produces:
- recalibration curves;
- empirical tolerances;
- the analyzer and activation class;
- a reporter-ion check;
- a mass-offset histogram for choosing mods.

Add **E0** (IDs with/without calibration vs Sage/MSFragger).

**[B3] The regression gate leaks the test set into development.** Every merge is tuned on Tier S, and the same datasets or close relatives are then reported. §5 protects *training* data but not heuristic tuning (top-N, bins, K, thresholds). **Fix:**
- Freeze a dev/test split by dataset in DATA.md now.
- Tier S is dev, and nothing measured on it is reported.
- Test sets are opened only at milestones and are logged before anyone looks at the results.
- The localization calibration sets must be disjoint from the sets used for FLR claims (§3).

**[B4] This scope cannot deliver on this hardware.** The plan has eight fronts:
- the DDA core;
- the cross-encoder;
- bi-encoder retrieval and clustering;
- open/PTM search;
- DIA including diaPASEF 4D tracing;
- XL, LiP and CF;
- three quant modes;
- repository scale-out.

Each of DIA-NN, xiSEARCH and FragPipe took person-years. See §6.

## 2. Sensitivity gaps against Sage, MSFragger and DIA-NN

**[M1] Charge, isotope error and MS1 features.** The design does not handle any of these:
- precursors with no assigned charge, which are common for Thermo low-m/z and ion-trap MS2. Try z=2–4, and let all charge hypotheses of a spectrum compete in *one* TDC;
- monoisotopic misassignment (isotope error −1…+3). This interacts with deamidation (0.984 vs 1.003 Da) and pollutes open-search offset histograms;
- MS1 feature finding.

**Fix:** add a batched, GPU-friendly MS1 isotope-envelope feature finder in the biosaur2/Dinosaur style. One feature table then serves:
- charge and isotope correction;
- DDA+ chimeric candidates. Without MS1, enumerating every mass in the isolation window inflates the search space and weakens TDC;
- LFQ, MBR and TMT purity;
- DIA MS1 scores.

**[M2] The fragment model is implicitly b/y only, which the design's own benchmarks rule out.** S1 Pfu and PXD000138 contain ion-trap CID. S2 Ferries and PXD014337 contain EThcD. **Fix:** an ion-type table per activation class, with fragment charge up to z−1:
- HCD: b/y, a2, immonium ions, H₃PO₄ loss;
- CID-IT: low-resolution bins, −98 losses;
- ETD/EThcD: c, z•, z+1, plus b/y.

Index according to the run's class, and put neutral losses in scoring rather than in the index.

**[M3] Peak preprocessing.** A global top-150 truncates three kinds of spectra:
- crosslink spectra, which hold two peptides;
- chimeric and DIA spectra;
- dense Astral and EThcD spectra.

Deisotoping is valid only at high resolution. **Fix:** use top-k per 100 Th as a view, deisotope only at high resolution, and ablate both on S0/S2/S3.

**[M4] Predictors on heterogeneous data.** Stage-B features need three things, none of which PRIDE metadata supplies:
- a per-run iRT→RT mapping;
- an NCE mapping. NCE is not comparable across Orbitrap generations, and timsTOF uses CE ramps;
- coverage of the modification. Crosslinks and open offsets have no predictions, so those strata get structurally different feature vectors.

**Fix:**
Fine-tune per instrument class, pick NCE empirically from pilot IDs, report Stage-B gains per instrument and stratum, and handle prediction-less strata explicitly in Group-walk.

**[M5] The frozen scorer should not be the default.** One frozen scorer will not be calibrated across Velos ion-trap CID, QE, Astral, TMT and EThcD. G2 tests behaviour under the null, not calibration across strata. Training labels taken from Sage reanalyses would import Sage's biases and its ~1% wrong labels. **Fix:**
- Make RESET-mokapot over Stage-A plus predicted features the default; MS2Rescore and Oktoberfest-class gains are proven.
- Run the cross-encoder as a research arm that must win E6 within each instrument class.
- Train on ProteomeTools and MassIVE-KB.
- Add **G4**: compare decoy-score quantiles per (instrument, activation, charge) before pooling.

**[M6] DIA: wrong first target, and interference is under-designed.** S5, PYE plasma on diaPASEF, is the hardest possible entry point:
- 4D tracing;
- narrow IM peaks;
- plasma dynamic range;
- 25 Th windows, where a spectrum-centric tier 1 has poor recall.

DIA-NN's sensitivity comes from interference handling (elution-reference fragments, fragment exclusion), which a peak-group cross-encoder does not replace. **Fix:**
- Start with Astral nDIA (ProteoBench DIA-Astral plus the Astral entrapment module, ~20 GB), where a spectrum-centric tier 1 is defensible.
- Then do QE variable-window DIA.
- Do diaPASEF last.
- Add an interference-aware fragment-selection operator.
- Measure tier-1 recall against DIA-NN IDs before building tier 2.

**[M7] Bruker DDA has no benchmark.** If Bruker is in scope, add ddaPASEF IM features, frame-summed MS2, CCS prediction and a PXD014777 subset.

## 3. PTM localization

**[M8] Retrieval truncates positional isomers.** With top-K=10–50, the isomers of multiply phosphorylated, S/T-rich peptides crowd each other, and the true peptide, out of K. **Fix:** localization should be a separate, exhaustive operator. For the best sequence, it should:
- enumerate *all* site configurations (hundreds, so cheap);
- score site-determining ions, aware of neutral losses and of c/z ions under EThcD;
- add predicted-intensity evidence.

Decoy residues (Ala/Leu) enter only at this step, never in identification (Ramsbottom). Add a localization G1 test: the neural head must not separate decoy-AA configurations from real ones on sequence alone.

**[M9] The validation design is not split, and the synthetic sets have known traps.** **Fix:**
- Calibrate on PXD000138 (Marx) and PXD009449 (ProteomeTools-PTM).
- Hold out Ferries PXD007058 for FLR claims.
- Use the Ferries U2OS TiO2 runs (real background) to check decoy-AA FLR against synthetic FLR.
- Use PXD014525 later for DIA.

Handle known traps: isomers co-occurring in a pool, synthesis impurities, CID phosphate migration.

Report FLR separately for mono vs multi-phospho and for adjacent sites (≤2 residues apart). Separate near-isobaric pairs: acetyl vs trimethyl, citrullination vs deamidation.

**[M10] Open-search localization.** For unknown offsets claim only enrichment (PTM-Shepherd-style); validate by searching PXD009449 *open*. The G2 null (precursor shifted by +7–50 Da) breaks under open or offset search, because the search absorbs the shift as an offset. Use wrong-species spectra instead.

## 4. XL-MS, LiP-MS, CF-MS

**XL [M11].** The core mechanisms are right: cleavable doublets give linear searches, α-first offset search, FDR computed at each level, self vs heteromeric. Missing:
- **One-sided matches.** A good α with a random β is the dominant false-CSM mode. Require fragment evidence for each peptide, and use features aware of min(α, β).
- **Competing explanations.** Linear, mono-link and loop-link explanations must compete in the same search (as in xiSEARCH).
- **MS3 and doublets.** DSSO MS2-MS3 needs the B1 scan linkage. Doublet detection needs fragment charge from high-resolution spectra.
- **Database reduction** from a linear pre-search must shrink targets and decoys symmetrically.
- **Self/heteromeric labels** must handle shared and homologous peptides.
- **Small counts.** When TD<DD, follow xiFDR conventions. Cross-check against xiFDR and IMP-X-FDR.

S3's 95-peptide database does not stress the O(n²) path; add a human background.

**LiP [M12].**
- **DIA dependency.** PXD031627 is DIA, so LiP depends on *semi-specific DIA*, the costliest DIA mode. Check PXD031616 and LiP-Quant for DDA runs, and do DDA LiP first.
- **Half-tryptic calls.** These must exclude protein N-termini, signal/transit peptides and Met removal. They must also flag in-source fragments, which co-elute with their fully tryptic parent.
- **The deliverable is peptide × run matrices** with cut-site coordinates, for MSstatsLiP's TrP correction. Do not reimplement that correction.
- **Truth.** Agreement with Spectronaut/MSstatsLiP shows concordance, not accuracy. Use known targets from PXD015446 (e.g. rapamycin–FKBP1A) as the truth set.

**CF-MS [M13].**
- **The DDA concept is being tested on DIA data.** S6 is CF-*DIA*, where "MBR" means library extraction in every fraction, so "transfer only within the same fraction index" is meaningless there. Define instead:
  - a run-specific detection q-value per fraction;
  - a library-borrowed provenance flag;
  - series-level global FDR.
- **Strict same-index transfer creates missing values that destroy elution profiles,** and transfers between adjacent fractions are often correct. Use PIP-ECHO-style transfer FDR (decoy peptides plus random-RT competition) with a profile-continuity prior.
- **E10's "known non-interacting pairs" is ill-defined:** absence from CORUM does not mean the proteins do not interact. Measure instead:
  - PR-AUC for CORUM intra-complex vs inter-complex pairs;
  - inflation of elution-peak width;
  - an entrapment-in-fractions control.
- Five fractions are a smoke test only.

## 5. Quantification and benchmarks

**[M14] Quant is under-specified and has no MS1 basis.** Specify each mode.
- **DDA LFQ:** MS1 features, alignment, and MBR with PIP-ECHO-style transfer FDR. Use Sage's per-target MS1 decoy as the baseline.
- **TMT:**
  - impurity correction;
  - reporter tolerance per plex (~3 mDa for TMTpro 32/35);
  - SPS-MS3 linkage;
  - isolation purity;
  - no normalisation for TPP.
- **SILAC** is common in legacy PRIDE and is absent. Support it or declare it out of scope.

A QuantUMS-style model is research; start with Top-N, weighted-median MaxLFQ and directLFQ. Benchmark ProteoBench DDA-QExactive first, then DIA-Astral, reporting ε with number quantified and completeness.

**[M15] The benchmarks do not represent PRIDE.** Add:
- **A stratified PRIDE sample** of 30–50 projects across years, instruments, labels and enrichments, with no manual parameters. Measure the fraction of runs that succeed, how well inferred parameters agree with metadata, and IDs against Sage at verified FDP. This is the headline test for the stated goal, and it validates B2.
- **Others:**
  - **TMT in Tier S:** PXD016766 (2.5 GB; Sage's 18 s is published).
  - **HLA-I immunopeptidomics in Tier M:** the non-specific index is 22.6 GB in Sage, about the 3090's 24 GB, and neural rescoring gains most here.
  - **ddaPASEF** (PXD014777).
  - **Low-input DIA entrapment** (PXD049412), where FDR control actually fails.
  - **ProteoBench DDA-Astral.**
- **Primary DDA FDR test:** PXD001468 with FDRBench paired entrapment, as in Wen 2025. Pfu S1 is a single low-resolution run with a few thousand IDs. It is a smoke test with too few expected false IDs at 0.1% to decide anything.

**[m1]** The entrapment decision rule has no tolerance, unlike the regression gate. Use the same binomial tolerance in both, or tools that are valid will come out "inconclusive" from noise alone.

## 6. Scoping, sequencing, minimum credible deliverable

**[M16] Hardware fairness.** Beating 4-core Sage with a 3090 is near-certain and uninformative. Report Sage per-core throughput, a cost-normalised comparison, and Sage at 32–64 threads (measured, or labelled as extrapolated).

TRFP ingest on 4 cores will probably dominate end-to-end time. Report it separately, as per-device rates.

**Minimum credible deliverable:** an open, GPU-accelerated **DDA** engine with:
- calibrated ingest that infers parameters;
- closed, open/offset and semi-specific search;
- RESET rescoring with predicted features;
- exhaustive in-search localization with decoy-AA FLR on held-out synthetic sets;
- TMT and MS1 LFQ;
- picked-protein and **global repository FDR shown flat in N** (PXD001468 plus the PRIDE sample);
- cleavable XL and semi-specific LiP as demonstrations of stratified FDR.

All at entrapment-verified FDP with honest speed.

**Order:**
1. Harness: dev/test freeze, E1, and Sage/Comet with entrapment on PXD001468.
2. DDA core: full-centroid store, pilot calibration, charge/isotope/MS1 features, ion tables, E2, paired TDC, picked protein.
3. RESET plus predicted features, with L6 global FDR and E9 brought forward. Both are cheap and central to the prompt.
4. PTM: open/offset search, localization, FLR on held-out sets, chimeric search (E7).
5. Assays on the shared offset index: XL, LiP-DDA, CF provenance and transfer FDR.
6. TMT/LFQ on ProteoBench DDA. This can run in parallel with steps 4–5.
7. Astral nDIA with entrapment.

Research arms behind go/no-go gates: cross-encoder (E6), bi-encoder (E5), clustering, diaPASEF, QuantUMS clone.

**Cut or defer:**
- bi-encoder-first retrieval;
- clustering/propagation, whose gain is unclear once each member must be verified;
- the cross-encoder as the default scorer;
- diaPASEF and semi-specific DIA;
- proteome-wide non-cleavable XL beyond α-first;
- TPP-specific features.

## 7. Minor

- **[m2]** Semi-specific decoys need a recipe per class. Generate them from decoy proteins so that target and decoy counts match.
- **[m3]** Mass-offset strata and tiny heteromeric XL groups need a minimum-size rule for Group-walk, and the power of each stratum should be reported.
- **[m4]** Split FAIMS multi-CV runs into per-CV pseudo-runs.
- **[m5]** f32 m/z (≈0.06 ppm) is fine. Apply recalibration as a stored offset.
- **[m6]** Exclude run-level tolerances and calibration from the index content hash, so recalibration never forces a rebuild.

## 8. What is good

- Statistics:
  - pre-registration with fixed seeds;
  - paired double competition;
  - Group-walk instead of TDC per stratum;
  - RESET for on-data learning;
  - lower-bound, combined and paired entrapment over several seeds.
- G1–G3, especially target/decoy recall symmetry for pruning stages.
- Picked protein groups (no razor), Rosenberger's global constraint, q×N reporting.
- Arrow batch operators with fixed schemas, a content-addressed amortised index, and a cohort mass sweep. This is the right architecture for scale.
- MBR that is opt-in, design-aware, has a separate transfer q-value, and never counts transfers as IDs.
- XL FDR per level and split self/heteromeric. LiP strata kept separate. No global normalisation for CF or TPP.
- The CPU/GPU split is decided by measurement (E2), and speed is reported per device.
- lit/04's accession verification is careful.
