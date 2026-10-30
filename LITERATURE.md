# LITERATURE.md: synthesis for decision-making

Compiled 2026-10-26 from the four surveys: [01 DDA engines](lit/01_dda_engines.md), [02 DIA and quantification](lit/02_dia_quant.md), [03 statistics, FDR, scoring](lit/03_stats_fdr_scoring.md), [04 assays and benchmarks](lit/04_assays_benchmarks.md). Section numbers like "01 §3" point into those files. **[unverified]** markers are inherited from the sources. Accession checks marked "PRIDE API" were run on 2026-10-26 for this synthesis (section 7). Dev node: 4 CPU cores, ~500 GB RAM, 1x RTX 3090 (04); most published speed numbers below are on 32-64 cores, so compare per core or per GPU.

## 1. Executive summary

The fastest open CPU baseline is a fragment-ion index (MSFragger, Sage). Sage (MIT, Rust) searches 12 TMT files in 18 s on 32 vCPU, and nothing on GPU has yet beaten it for closed tryptic search with matched sensitivity. GPUs pay off for open and semi-specific search, dense similarity kernels, neural reranking, localization and embedding clustering. The real costs at PRIDE scale (2.85 PB) are raw-file decoding and I/O, open or PTM search spaces, neural rescoring and global FDR, not closed-search arithmetic. The statistical state of the art is weaker than the speed state of the art. DIA-NN, Spectronaut and EncyclopeDIA do not reliably control FDR under valid entrapment estimators (worst at protein level and low input). Match-between-runs (MBR) has no validated transfer FDR in DIA. Naively merged per-study 1% lists reach 28% precursor FDR (MassIVE-KB). Learned scorers can break target-decoy "equal chance" through sequence priors and train/test leakage. Our plan: a Sage-style indexed core extended with a partitioned, persistent, cohort-amortised index and a columnar zero-copy store; chimeric retrieve-then-rerank scoring; peptide-level paired competition with RESET or static models; entrapment as a built-in harness; stratified FDR through Group-walk; and master-level picked-group FDR. We beat DIA-NN InfinDIA and Sage on speed per unit hardware and the incumbents on entrapment-validated FDR, PTM FLR and transfer FDR.

## 2. State of the art by area

### 2.1 DDA engines and acceleration (01 §1-3, 7)
- **Fragment index is the standard.** MSFragger (Kong 2017, doi:10.1038/nmeth.4256): open 500 Da search of 41,820 spectra in 5.4 min on 4 cores, ~150x Comet and X!Tandem (815/976 min). Index memory grows with search space: 1.6 GB tryptic, 55.8 GB semi-tryptic, 86.6 GB phospho. 32-bit peptide IDs cap it at ~2x10^9 peptides.
- **Sage** (Lazear, doi:10.1021/acs.jproteome.3c00486, MIT): 12 TMT files in 18 s vs MSFragger 101 s, Comet 345 s. Open search 44 s vs MSFragger 225 s. Its two-level layout (fragment-m/z buckets re-sorted by precursor mass) is cache-friendly and a natural GPU tile. It runs out of memory (OOM) on huge databases. Tide v4.3 (streaming, precursor-sorted sweep) handles 25.2 M spectra x 5.8 B peptides in 16.3 h and 2 GB, where Sage and MSFragger OOM.
- **Chimeric search is worth +15-57% peptides.** More than 2/3 of DDA MS2 spectra are chimeric. MSFragger-DDA+ (doi:10.1038/s41467-025-58728-z) is 4x faster than MaxQuant and at least 20x faster than CHIMERYS+PD. CHIMERYS is commercial (doi:10.1038/s41592-025-02663-w). No open DDA+-class engine with entrapment-validated FDR exists.
- **GPUs.** AlphaPept's GPU PSM step took 77 s vs 9 s multithreaded CPU (memory-bound). GiCOPS (doi:10.1038/s41598-023-43033-w; CC BY-NC-ND) gets 1.2-2x closed and 3-5x open over its CPU predecessor. HOMS-TC (tensor cores, library OMS) is 31x over ANN-SoLo-GPU. Near-memory work finds search is limited by reference-data movement, not compute.
- **Licences.** MSFragger (Fragmatics commercial licensing), GiCOPS, MS-GF+, CHIMERYS, ionbot and pFind are benchmark-only. Reusable: Sage, timsrust, ThermoRawFileParser, Casanovo, PTM-Shepherd, ANN-SoLo/HOMS-TC.

### 2.2 GPU and embedding retrieval (01 §3-5)
- **Learned scores win sensitivity.** Casanovo-DB: +31-40% peptides vs XCorr, +88-102% vs Andromeda at 1% FDR (no entrapment check done; 03 §9).
- **Bi-encoder first is unproven.** The 2026 framework paper (doi:10.1021/acs.jproteome.5c01153) estimates biencoder+ANN at 2-3 h per 10^6 spectra x 10^8 peptides vs ~1-2 h for a fragment index. Its GPU dot-product rate (2x10^7/s) is conservative by 50x or more [FAISS>10^9/s on A100]. Exhaustive cross-encoder: ~2-12 days. Recommended: index top-K, then neural rerank.
- **Where embeddings clearly pay.** GLEAMS embedded 669 M spectra (32-d), clustered 511 M into 60 M clusters and added 132 M PSMs (+71%). HyperSpec clustered 21 M spectra in 24 min on GPU. Embeddings suit cross-run dedup, unknown-mass-shift retrieval and pre-filters.
- **De novo as side channel.** Casanovo v5 runs at ~235 spectra/s per GPU; pi-PrimeNovo is up to 89x faster than autoregressive models. This is 1-2 orders slower than index search, so run it only on unexplained spectra.

### 2.3 DIA (02 §1-4)
- **DIA-NN 2.x is the speed and sensitivity bar.** InfinDIA pre-search scores 210-440 M target+decoy pairs/min on Orbitrap/Astral and 29-51 M/min on timsTOF (64-core 7980X). Its q-values may deviate by 1.5-2x. Peptidoform FDR can reach 5-10% without peptidoform mode. Enterprise-only "Knowledge base" must be off for entrapment.
- **Open alternatives are slower or narrower.** AlphaDIA (Apache-2.0): Astral 21 min HeLa gives >120k precursors, entrapment FDP 0.2% precursor / 1.08% protein. It is slow with predicted libraries on Bruker. DIA-BERT is GPU-only and Orbitrap-only (+51% proteins vs DIA-NN, claim not independently validated). Beta-DIA targets diaPASEF. DIA-CLIP (zero-shot dual encoder) is unreplicated.
- **Spectrum-centric DIA works.** MSFragger-DIA / diaTracer / FragPipe: 144 vs 337 min (Spectronaut 17) on 20 ccRCC runs, entrapment FDP 1.9%. diaTracer enables open/offset PTM search on diaPASEF.
- **Instrument drives design.** Astral 2-Th nDIA is close to "DDA with complete sampling" (spectrum-centric fits). diaPASEF needs 4D tracing or IM-aware extraction and is everyone's slow path (4-8x).
- **Prediction cost is amortisable for closed human tryptic, not for open spaces.** AlphaPeptDeep: 7.9 M precursors in ~10 min on one GPU. Open/multi-PTM spaces (10^8-10^9 precursors) cost 1-30 GPU-hours. Hence lazy, cached prediction.

### 2.4 Quantification and MBR (02 §7, 9)
- **Roll-up scaling.** directLFQ is linear: 10,000 DIA samples in ~16 min vs ~2.5 days for iq (doi:10.1016/j.mcpro.2023.100581). QuantUMS (doi:10.1038/s41587-026-03131-2) is linear, trained on 6-100 runs, with >1.4x lower low-intensity log2-SD, and ratio compression remains at low input.
- **MBR is the weakest statistical link.** At nominal 1%, measured false-transfer FDP on single-cell data: PIP-ECHO 1.0%, MaxQuant 2.0%, IonQuant 3.4%, FlashLFQ 6.1% (doi:10.1021/acs.jproteome.5c00065). Lim 2019: 44% of yeast proteins transferred into yeast-free runs (PXD014415). IonQuant: 19-38x faster than MaxQuant, with m/z-shifted decoy transfers. Sage controls MS1 integration with a decoy per target. No DIA tool reports a validated transfer FDR.
- **Cohort recipe.** Empirical library from 20-100 runs, per-run search, cross-run steps on cached .quant files (two raw passes). quantms: 29,354 RAW files, ~2 min per file; quantmsdiann ~17 node-min per DIA run.
- **Incremental stability.** Quantities should not change when runs are added (Top-N against references, or frozen QuantUMS parameters).

### 2.5 FDR theory, entrapment, rescoring and neural-scoring pitfalls (03 §1-5, 9)
- **TDC needs equal chance and independence.** With the +1 correction, peptide-level TDC controls FDR (Barber-Candes; He 2015). PSM-level TDC is liberal because of multiplicity. Paired double competition gives +17% peptides at controlled FDR (Lin 2022, doi:10.1021/acs.jproteome.2c00282). Multi-decoy competition (Emery 2020) gives up to +50% at small thresholds. Keep an explicit target-decoy pairing id everywhere.
- **Percolator and mokapot can fail** (Freestone 2024, doi:10.1021/acs.jproteome.3c00902), notably when runs are combined. RESET (decoy split into train and estimate halves; factor 2 single-decoy, 1 two-decoy) is provably valid and slated for Percolator 4.0. Static pretrained models avoid leakage too.
- **Entrapment estimators (Wen 2025, doi:10.1038/s41592-025-02719-x; FDRBench).** Combined N_E(1+1/r)/(N_T+N_E) is an upper bound. Paired (r=1) is tighter. N_E/(N_T+N_E) is only a lower bound and can show failure, never control. "Sample" is invalid. With valid estimators, DIA-NN 1.8.1, Spectronaut 18.7 and EncyclopeDIA 2.12 all failed somewhere (single cell, protein level); DDA Tide, Sage, MS-GF+ and MSFragger+MSBooster looked valid at peptide level. DIA-NN 2.x was not tested.
- **Neural-scoring leaks.** (i) Sequence priors violate equal chance (sequence-only classifiers separate decoys from targets at AUC 0.64 for reversed, up to 0.81 for DIA-NN-style; Reznikov 2026). (ii) A near-perfect target/decoy split from one feature such as -log MS-GF+ E-value in MS2Rescore/Oktoberfest is a red flag (Winkelhardt 2026). (iii) Training on PRIDE and reprocessing PRIDE is a leakage default. (iv) ANN/top-k retrieval must treat targets and decoys symmetrically (03 §9.2).
- **Groups.** Per-group TDC is not rigorous. Group-walk (+26% by modification mass difference, +4% by charge) is. CONGA gives rigorous open search with narrow-search power.
- **FDR vs FDP.** FDR is an expectation. Report TDC-UB or FDP-SD bounds and variability across >= 5 decoy seeds (decoy-induced variability 4.7% at 1% FDR on 15k spectra).

### 2.6 PTM localization and FLR (01 §2; 03 §6)
- PSM FDR, FLR and site FDR are different quantities. Decoy amino acids (Ala, Leu best) give an in-situ FLR (Ramsbottom 2022, doi:10.1021/acs.jproteome.1c00827). Synthetic-calibrated scores are often miscalibrated on real data.
- Single-evidence sites in PhosphoSitePlus/PeptideAtlas have 82-98% estimated FDR (Kalyuzhnyy 2022). Camacho 2023 aggregates PSM to site to multi-study with tracked FLR.
- Locard-Paulet 2020: across 22 pipelines the same score threshold yields very different FLRs.
- Every engine does localization as a separate post-pass. Fusing site-determining-ion evidence into top-K scoring is open (dense, regular, GPU-friendly). Tools to reuse or compare: PTMProphet, ptmRS, Ascore, LuciPHOr2, DeepFLR, PTM-Shepherd (offset annotation, no per-offset FDR).
- Ground truth: PXD000138, PXD007058, PXD009449, PXD014525 (section 5).

### 2.7 XL-MS (03 §8; 04 A1)
- Non-cleavable search is O(n^2) in linkable peptides. Remedies: MS-cleavable doublets (linear searches), alpha-first offset search with an open-search index (pLink 2: 40x pLink 1), sparse-matrix formulation (MS Annika 3.0, our nearest GPU precedent), database reduction (reduce decoys the same way), ML rescoring (Scout).
- FDR must be set at the claim level (CSM, residue pair, PPI) and separately for self and heteromeric links. At 5% CSM-FDR, actual PPI error was 35-60% (Lenz 2021; heteromeric decoy matches 10.6x more likely). Merge to PPI first, then estimate FDR on heteromeric links only (gives 4.8-4.9% at 5%).
- At nominal 5%, actual false rates range 2.4-32% by pipeline (Beveridge 2020; pLink 11.3%, xi 2.5%, XlinkX 0%). All six engines exceeded nominal 1% in Matzinger 2022. Scout is >200x faster than xiSEARCH on a small test.

### 2.8 LiP-MS (04 A2)
- Half-tryptic peptides must be searched (semi-specific, ~10-20x candidate space). Use stratified FDR for half vs fully tryptic. Quantify at peptide/cut-site level, normalised to the paired TrP control (MSstatsLiP, FLiPPR). Never transfer by MBR between LiP and TrP runs. DIA libraries must include semi-tryptic precursors.

### 2.9 CF-MS (04 A3)
- Elution-profile correlation is the signal, so cross-fraction MBR can manufacture false co-elution (mechanistic argument; no paper quantifies it, 04). Needs per-fraction provenance, same-fraction-index transfers only (or +-k window with its own transfer FDR), no per-fraction total-intensity normalisation, peptide-by-fraction matrices out, and series-level (not per-run) FDR. PPI-level FDR needs decoy complexes (CCprofiler) or a reference set (PrInCE, EPIC). Skinnider & Foster reanalysed 206 CF-MS experiments (CFdb, PXD042664).

### 2.10 Repository-scale reanalysis and global FDR (01 §7; 03 §4)
- Scale: PRIDE has 42,036 datasets / 2.85 PB (Aug 2024). ionbot reanalysed ~924 M human spectra; MassIVE reanalysed 2.2 B spectra / 392 datasets; GLEAMS embedded 669 M.
- Merged per-dataset 1% PSM lists give 28% library precursor FDR (MassIVE-KB, doi:10.1016/j.cels.2018.08.004). Classic protein FDR over-counts decoys with data size (MAYU; Savitski 2015). Picked Protein Group FDR (The 2022, doi:10.1016/j.mcpro.2022.100437): razor/Occam grouping becomes anticonservative, giving no protein groups at 1% entrapment FDR beyond ~300 experiments. ProteomicsDB gives 18,000 groups at 1% group FDR from ~410 M PSMs.
- Required: global competition per claim level, plus a global analyte constraint (pass run-specific and global q <= 1%; Rosenberger 2017, doi:10.1038/nmeth.4398). Novel proteins need HPP evidence (>= 2 unique non-nested 9+ aa peptides).

### 2.11 File formats and I/O (01 §6)
- mzML is slow and large. mzMLb/HDF5 and Aird are faster. The emerging standard is Parquet-based mzPeak (doi:10.1021/acs.jproteome.5c00435; spec in development, no hard benchmarks yet). Sage writes Parquet, ThermoRawFileParser 2.0 emits Parquet (schema v0.2), Casanovo v5 moved to Lance. Bruker: timsrust (MIT, pure Rust). Thermo needs the .NET RawFileReader (redistributable). No native SCIEX reader on Linux (04 B2).
- Decision: decode and centroid once into a memory-mappable columnar store aligned with mzPeak's point layout. Treat conversion throughput as a first-class benchmark. Raw-read cost alone for 2.85 PB at 1 GB/s per node is ~33 node-days (01 estimate).

## 3. Competitor numbers we must beat

Hardware is as published. Our dev node (4 cores + RTX 3090) is much smaller than most rows, so report per-core and per-GPU-hour and re-run baselines on our node.

| Tool | Dataset / accession | Metric | Value | Hardware |
|---|---|---|---|---|
| Sage | PXD016766, 12 TMT files | Time / PSMs at 1% | 18 s / 100,448 | AWS c5ad.8xlarge, 32 vCPU |
| MSFragger | same | Time / PSMs | 101 s / 97,844 | same |
| Comet | same | Time / PSMs | 345 s / 92,333 | same |
| Sage | PXD001468, open -100..+500 Da | Time / PSMs | 44 s / 24,682 | 32 vCPU |
| Sage +chimera | same | Time / PSMs | 70 s / 28,724 (+16%) | 32 vCPU |
| MSFragger | same | Time / PSMs | 225 s / 24,911 | 32 vCPU |
| Sage open + Met-ox | PXD001468 | Time, RAM | 68 s, 7.3 GB (MSFragger 340 s, 11.55 GB) | 32 vCPU |
| MSFragger | HEK293 41,820 spectra, open 500 Da | Time | 5.4 min (Comet 815) | 4-core Xeon |
| Tide v4.3 | PXD011515 (12.5 M spectra, 7.1 B peptides) | Time, RAM | 179 h, 20 GB (Sage OOM; MSFragger 448 h, 124 GB) | 4 threads |
| Tide v4.3 | 539-file set, 5.8 B peptides | Time, RAM | 16.3 h, 2 GB (derived ~430 spectra/s) | 4 threads |
| Sage / Tide / MSFragger | PXD028806 (221 GB) | Time | 431 s / 673 s / 793 s | 4 threads (Tide paper) |
| MSFragger-DDA+ | 62 mzML | Time | 68.8 min (>= 20x CHIMERYS+PD) | workstation |
| AlphaPept | PSM step, 260k spectra vs 5 M entries | Time | CPU 9 s; GPU 77 s | multithreaded CPU / GPU |
| GiCOPS | open search | Speed-up | 3-5x vs HiCOPS | RTX A6000 + 2x Xeon 5215 |
| HOMS-TC | iPRG2012, HEK293 library OMS | Speed-up | 31x vs ANN-SoLo-GPU | RTX 4090 |
| Casanovo v5 | de novo | Throughput | ~235 spectra/s | one GPU |
| DIA-NN 2.3 InfinDIA | Astral 40 SPD, PXD044991 [instrument tag disputed, section 7] | Pairs/min | 215 M | 64-core 7980X |
| DIA-NN InfinDIA | timsTOF 60 SPD PXD017703 / 70 min PXD034222 | Pairs/min | 51 M / 29 M | same |
| DIA-NN InfinDIA | QE HF-X 65 min, PXD005573 [K-GG label doubtful] | Pairs/min | 440 M | same |
| FragPipe-MSF vs DIA-NN lib-free vs Spectronaut 17 | 20 ccRCC runs | Time | 144 min vs 2x slower vs 337 min | workstation |
| CHIMERYS | LFQbench | Speed | 4.9x faster than DIA-NN, 1.7x slower than Spectronaut | 8-core VM |
| AlphaPeptDeep | human library 7.9 M precursors | Time | ~10 min GPU / ~100 min CPU | one GPU |
| directLFQ | 10,000 DIA samples | Time | ~16 min (iq 2.5 d) | 64 cores |
| IonQuant MBR | HeLa / two-proteome | Time | 37 / 28 min (MaxQuant 699 / 1,056) | desktop |
| quantms | 29,354 RAW | Time | ~2 min per file, ~10 h per dataset | HPC |
| quantmsdiann | 2,300 single-cell runs | Time | 2.2 h on 300 nodes (~17 node-min/run) | HPC |
| HyperSpec | clustering | Time | 21 M spectra in 24 min | GPU |
| Entrapment bars (accuracy) | Wen 2025 | FDP at 1% claimed | DIA-NN single cell 2.3-4.7% (lower bound); AlphaDIA 0.2% / 1.08%; FragPipe-MSF 1.9% | n/a |
| MBR transfer bars | single-cell | False transfer FDP at 1% | PIP-ECHO 1.0%, MaxQuant 2.0% | n/a |
| Quant bars | HYE ProteoBench | Epsilon, CV | Not extracted. DIA-NN Astral nDIA CV <7% (Guzman); QuantUMS 1.3x lower protein SD | n/a |

## 4. Open problems and opportunities

1. **Compressed, partitioned, persistent fragment index** (cohort-amortised; hash-keyed by FASTA, enzyme, mods, version). Target PXD011515-scale (7 B peptides) in <64 GB where Sage OOMs.
2. **Cohort sort-once, sweep-once execution** (Tide idea on a Sage index) plus spectrum clustering and dedup (GLEAMS suggests 511 M to 60 M clusters, a 5-10x scoring reduction [estimate]).
3. **GPU closed search that beats Sage.** Needs batching of thousands of spectra per index slice, SpMM over binned spectra, or restricting the GPU to open search, rerank, localization and clustering. Prototype and compare on identical candidates before committing.
4. **Open, entrapment-validated chimeric search** (DDA+ quality). Rank-stratified FDR.
5. **Retrieve-then-rerank neural scoring** with static or RESET-valid, calibrated, percolator-free scores. Go/no-go for bi-encoder-first: recall@K >= 95% of stage-1 top-1 on entrapment benchmarks.
6. **In-search PTM localization** with decoy-residue FLR, site-level FDR, mass-offset search with a curated catalogue and per-group FDR (Group-walk).
7. **A validated transfer FDR for DIA and cross-run borrowing**, and CF-MS-aware (fraction-restricted) MBR.
8. **Fast open (MIT/Apache) DIA engine**, especially on timsTOF where CPU engines are 4-8x slower and sparse 4D extraction on GPU is open ground.
9. **Single raw pass** with per-candidate summaries stored columnar for library refinement, MBR and quant; incremental addition of runs with frozen parameters.
10. **Lazy cached prediction** keyed by (peptidoform, charge, NCE, instrument) for open/PTM spaces.
11. **Global FDR master lists** (picked peptide / protein group / PTM site) with absolute false-count reporting, and tested at N in {1, 10, 100, 1000} runs.
12. **Assay-specific engines**: XL SpMM on GPU with level-aware FDR; LiP semi-specific index with TrP pairing; the CF-MS MBR harm quantified (small publishable experiment, 04).

## 5. Benchmark datasets (consolidated)

"Verified" = PRIDE API in survey 04 (live) or in this synthesis (re-check marked "re-checked 2026-10-26"). Sizes are total project size. See 04 Part B for the full list and the S/M/L tiering.

| Accession | Purpose | Size | Instrument | Verified? |
|---|---|---|---|---|
| PXD001468 | Chick 2015 HEK293: DDA speed, open search, DDA entrapment (Sage/MSFragger reference) | 116.9 GB, 24 raw | Q Exactive | Yes (re-checked) |
| PXD016766 | Sage TMT benchmark (12 files used) | 16.7 GB, 48 raw | Orbitrap Eclipse | Yes (re-checked) |
| PXD011515 | Metaproteome, huge search space (Tide paper) | not in 04 | Q Exactive | Exists (re-checked); size not checked |
| PXD028806 | Normal-size speed test (Tide paper) | not in 04 | Q Exactive Plus | Exists (re-checked) |
| PXD001077 | Tiny DDA entrapment (P. furiosus) | 0.5 GB | LTQ Orbitrap Velos | Yes |
| PXD000138 | Marx 2013 synthetic phospho library, FLR gold standard | 69.4 GB, 193 raw | LTQ Orbitrap Velos | Yes (re-checked) |
| PXD007058 | Ferries 2017 phospho pools + U2OS, 8 fragmentation schemes | 32.6 GB, 96 raw | Orbitrap Fusion | Yes (re-checked) |
| PXD009449 | ProteomeTools 21 PTMs | 150.1 GB | Fusion Lumos | Yes (re-checked) |
| PXD014525 | DIA phospho with synthetic spike-ins | 699 GB | Q Exactive HF(-X) | Yes (re-checked) |
| PXD004732 | ProteomeTools part I, synthetic unmodified truth | 1.09 TB | Orbitrap Fusion | Yes |
| ProteoBench DIA-entrapment-Astral | DIA FDR, leaderboard | 5.65 GB tar | Orbitrap Astral | Yes (server) |
| PXD041421 | A549/K562 diaPASEF batch effects, Bruker entrapment | 532 GB | timsTOF Pro | Yes (re-checked) |
| PXD012988 | DIA reproducibility (entrapment) | 112.7 GB | QE HF | Yes (re-checked) |
| PXD028735 | Van Puyvelde HYE multi-vendor, ProteoBench DDA-QE | 3.43 TB | TTOF, Synapt, QE HF-X, timsTOF | Yes |
| PXD070049 | HYE "Generation Beta", modern instruments | 3.11 TB | Astral, timsTOF Ultra, ZenoTOF | Yes (re-checked) |
| PXD062685 | ProteoBench diaPASEF HYE | 44 GB | timsTOF SCP | Yes (re-checked) |
| PXD002952 | LFQbench HYE (SCIEX format, low priority) | 470 GB | TripleTOF | Yes |
| PXD017703 | diaPASEF (Meier), InfinDIA speed set | 826 GB | timsTOF Pro | Yes |
| PXD014777 | HYE ddaPASEF; IonQuant MBR | 204 GB | timsTOF Pro | Yes |
| PXD014415 | Lim 2019, MBR false-transfer two-proteome | not in 04 | Fusion Lumos | Exists (re-checked) |
| PXD046283 / PXD046453 | Astral nDIA (Guzman) | 346 GB | Orbitrap Astral | Yes (PXD046453 re-checked) |
| PXD049412 | Astral single-cell, two-proteome | 527 GB | Astral, Exploris 480 | Yes (re-checked) |
| PXD003881, PXD001819 | DDA LFQ spike-ins (IonStar, UPS1) | 40 / 61 GB | Fusion / Velos | Yes |
| PXD040205, PXD016647, PXD029009, PXD015087 | QuantUMS/directLFQ/CHIMERYS quant sets | n/a | timsTOF / QE / TripleTOF | Exist (re-checked) |
| PXD056598 | PYE plasma, 1,116 runs, 12 sites (per 02) | n/a | multi-vendor | **No: PRIDE API 404** |
| JPST003358 | ProteoBench DIA-plasma PYE subset | 14.1 GB tar | timsTOF | Yes (server, 04) |
| PXD080610 | ProteoBench low-input (per 02) | n/a | Astral, timsTOF Ultra 2 | **No: PRIDE API 404** |
| PXD042704, PXD034525, PXD023325 | Wen 2025 DIA entrapment sets | n/a | Astral, Lumos, Eclipse | **Not on PRIDE** (404; Panorama/jPOST). PXD023325 maps to JPST000971 |
| MSV000084000 | Wen 2025 entrapment | n/a | n/a | Not verified |
| PXD036062 | Messner yeast KO cohort (per 02, "verify") | n/a | SWATH | **No: PRIDE API 404** |
| PXD005573 | InfinDIA speed set | n/a | QE HF | Exists, but is the Bruderer DIA-optimisation project (re-checked; 02's "K-GG" label doubtful) |
| PXD044991 | InfinDIA speed set | n/a | PRIDE instrument tag: QE HF; keyword: Astral | Exists (One-Tip, re-checked), instrument ambiguous |
| PXD014337 | XL: Cas9 95-peptide library | 9.9 GB | QE HF-X, Lumos | Yes (re-checked) |
| PXD029252 | XL: ribosomal library, 6 linkers | 112.7 GB | QE HF-X, Exploris, Eclipse | Yes (re-checked) |
| PXD042173 | XL: Scout 256 recombinant proteins | 233.8 GB | Fusion Lumos | Yes (re-checked) |
| PXD019120, PXD024148 | XL: Lenz E. coli SEC (BS3/DSSO; sSDA) | 331 / 38 GB | QE HF / Lumos | Yes (PXD019120 re-checked) |
| PXD027161 | XL on Bruker (MaxLynx) | 3.4 GB | timsTOF Pro | Yes |
| PXD031627, PXD031616 | LiP/TrP (MSstatsLiP protocol) | 110 / 260 GB | Fusion Lumos | Yes (run-to-condition map to be read from SDRF) |
| PXD015446, PXD034606 | LiP-Quant; PELSA | 1.95 TB / 1.35 TB | QE / Exploris | Yes |
| PXD059608, PXD044083 | CF-DIA (small Thermo fractions) | 241 / 125 GB | QE HF / Exploris | Yes (PXD059608 re-checked) |
| PXD007038, PXD010288, PXD036711, PXD049340, PXD012009, PXD042664 | SEC-SWATH, mini-complexome, PCP-SILAM, CFdb | 0.17-2.6 TB | various | Yes (PXD042664 re-checked) |
| PXD013774, PXD009495 | TPP | 56 / 248 GB | Q Exactive | Yes |
| Glyco (Kawahara HUPO, GlyConnect 2943) | Glycopeptide consensus | n/a | n/a | **Not verified** (04) |

## 6. Known pitfalls we commit to avoiding

1. Claiming FDR control from N_E/(N_T+N_E), the "sample" estimator, or agreement with another tool's list. Anything but upper bound <= alpha at every tested alpha on every dataset is "inconclusive", never a pass (03 §10.3 E5).
2. Training on decoys and estimating on the same decoys; PSM-level cross-validation folds; pooling folds with different calibration. Use RESET or a static model, with folds grouped by peptide/stem.
3. Train/test overlap: project-disjoint model folds; never score a project with a model trained on it (including PRIDE-trained models on PRIDE).
4. Learned sequence priors and leaking features. Gate with the sequence-only decoy-AUC test (threshold > 0.6), the null-spectrum target-win-fraction test (0.5 +- 99% CI), and a near-perfect separation alarm. Use foreign-species entrapment as well as shuffled.
5. Asymmetric retrieval or pruning of decoys vs targets in ANN/GPU stages; re-run entrapment after every retrieval change.
6. PSM-level FDR on multiplicity-prone data; merged per-run q-values at repository scale; classic protein FDR; razor/Occam grouping for FDR.
7. Independent per-group TDC (use Group-walk), independent per-offset TDC in open search, and pooling half- and fully-tryptic classes.
8. Tuning features, hyperparameters or cutoffs on the entrapment results of test data. Changing decoy generation, entrapment construction, estimator formulas or seeds after seeing results (fix them, with the seed list, in DESIGN.md).
9. XL: CSM-level FDR claimed as PPI-level; pooled self and heteromeric; reduced search database with unreduced decoys.
10. MBR transfers counted as identifications; MBR across LiP/TrP or non-matching CF-MS fractions; global per-fraction normalisation of CF-MS or TPP series.
11. Benchmark hygiene: CV on complete cases only, run-specific q-values presented as global, imputation-like defaults, comparing against DIA-NN with Knowledge base on, single decoy seed, small test lists without FDP bounds, speed numbers without hardware and peak RAM.
12. Treating published ID gains from unreleased or commercial neural scorers (Tesorai, pUniFind, DIA-BERT, DIA-CLIP) as established until re-tested with valid entrapment.
13. Using SCIEX .wiff or controlled-access data where no Linux reader or open access exists; relying on PXD000647 (Legionella, not TPP).

## 7. Cross-file inconsistencies and how they resolve

Checked against the PRIDE API (https://www.ebi.ac.uk/pride/ws/archive/v3/projects/ACC) on 2026-10-26.

| # | Issue | Resolution |
|---|---|---|
| 1 | **PTM synthetic sets: PXD000138 vs PXD007058.** 01 refers to "the synthetic PXD000138 set" and "the Kuster synthetic set" interchangeably (pyLucXor 51,468 sites, 98.55% accuracy attributed to PXD000138). 04 lists both. | **Resolved.** Different datasets. PXD000138 = Marx 2013 (Kuster, LTQ Orbitrap Velos; doi:10.1038/nbt.2585), >100k synthetic peptides. PXD007058 = Ferries 2017 (Eyers lab, Orbitrap Fusion; doi:10.1021/acs.jproteome.7b00337), 8 fragmentation schemes. 01's R6 "PXD000138 and the Kuster synthetic set" is the same dataset. The pyLucXor numbers stay as stated in 01, not re-checked. |
| 2 | **Sage benchmarks: PXD016766 vs PXD001468.** | **Resolved.** PXD016766 (Eclipse TMT benchmark) is Sage's TMT-MS2 benchmark; PXD001468 (Chick 2015 HEK293) is the open-search benchmark. 01 and 04 agree. 04 is right that PXD001468 does not apply to Sage's TMT or LFQ. The "HEK293 41,820 spectra" MSFragger timing is from the same Chick data family (single-run subset); not re-checked. |
| 3 | **Wen 2025 DIA entrapment accessions.** 02 lists PXD042704, PXD034525, PXD023325 as PRIDE; 04 says they are Panorama/jPOST. | **Resolved in 04's favour**: all three return 404 on PRIDE. PXD023325 corresponds to JPST000971. 02 labels PXD023325 "1-cell/100-cell Eclipse"; 04 lists the instrument as "LTQ Orbitrap (metadata)". Instrument unresolved. |
| 4 | **PYE plasma accession.** 02: PXD056598, 1,116 runs. 04: ProteoBench subset JPST003358. | **Unresolved.** PXD056598 returns 404 on PRIDE (may be private, in another repo, or mis-cited). Use JPST003358 / ProteoBench only. |
| 5 | **PXD080610 (02, ProteoBench low-input) and PXD036062 (02, Messner yeast KO).** | **Unresolved.** Both 404 on PRIDE. PXD080610 may be unreleased. Do not plan around them. |
| 6 | **InfinDIA speed-table accessions in 02.** PXD044991 is listed as "Astral 40 SPD"; PRIDE tags it Q Exactive HF (title "One-Tip", keyword Astral). PXD005573 is listed as "QE HF-X 65 min, K-GG"; PRIDE shows the Bruderer DIA-parameter-optimisation project on QE HF with no K-GG mention. | **Partly unresolved.** Treat 02's accession-to-condition mapping as unverified. Re-run DIA-NN on our node on named files and report our own numbers, not the table. |
| 7 | **Wen 2025 numbers: FDP vs lower bound.** 02 reports DIA-NN "FDP 2.3-4.7%" on single cell; 03 states they are lower-bound estimates at 1%. | **Resolved in 03's favour.** They are lower bounds, so true FDP is at least that. Ranges like "0.7-1.3% Astral" in 02 are the same type of value. |
| 8 | **Sage activity.** 02: last release v0.14.7 (Mar 2024), development slowed. 01: repo active again with commits in Oct 2026 (named PTMs, radical z-ions, CI harness), last tags v0.14.7 / v0.15.0-beta.2. | **Resolved.** 01 is more recent. Releases are old; HEAD is active. Benchmark against HEAD and record the commit hash. |
| 9 | **Tesorai Search.** 01: "no decoy-based validation, FDR red flag". 03: authors report entrapment checks (estimator unverified; JMB 2026). | **Unresolved**; keep as unverified in either direction until estimators are inspected. |
| 10 | **Zero-shot / pretrained DIA claims.** 02 and 03 both flag DIA-BERT and DIA-CLIP gains as unvalidated; 02's "+51% proteins" and "FDP 0.36% vs 1.1%" are paper claims. | Consistent; no change, but 02 table row "Conservative" is a claim, not a result. |
| 11 | **Glyco and ISB18 / MSV000084000.** | **Unresolved** (04): no verified accessions. |
| 12 | **CHIMERYS completeness (02) vs speed numbers (01).** 02: 4.9x faster than DIA-NN; 01: DDA+ at least 20x faster than CHIMERYS+PD. | Not contradictory (different data, DIA vs DDA, cloud vs PD). Do not combine them. |
| 13 | **Speed hardware mismatch.** Published rows use 32-64 cores; 04 says our node has 4 cores and one RTX 3090. Derived throughput comparisons (spectra/s) are derived, not published. | Always rerun baselines on our node; report per-core and absolute. |

Remaining unresolved items: 4, 5, 6 (partly), 9, 11, and the instrument/mapping details in 3.
