# 02 — DIA analysis and quantification: state of the art (as of Oct 2026)

Scope: DIA search engines (peptide-centric, spectrum-centric, feature-free, deconvolution), predictors and their inference cost, FDR in DIA, label-free/isobaric quantification incl. MBR and alignment, benchmark datasets and metrics, and cohort-scale processing. Numbers come from primary papers or official docs. Items marked *(verify)* were not confirmed against the primary text.

---

## 0. TL;DR

- **DIA-NN (2.x, closed source, now Academia-limited vs paid Enterprise) is the speed and sensitivity reference.** Its new **InfinDIA** pre-search scores **210–440 million target+decoy pairs per minute** on Orbitrap/Astral data on one 64-core AMD 7980X. It falls to 29–51 M/min on timsTOF. That sets a hard CPU bar for any GPU engine. **QuantUMS** (Nature Biotechnology 2026) is the reference quant module. It trains on 6–100 runs, has linear complexity, and reduces ratio compression.
- **No widely used DIA tool reliably controls peptide-level FDR on every dataset.** Wen et al. (Nat Methods 2025) found DIA-NN at a claimed 1% FDR reaching up to 2.3–4.7% precursor FDP on single-cell data, and Spectronaut 3.8–7.6%. CHIMERYS reports that DIA-NN and Spectronaut underestimate run-specific FDR. Entrapment-validated FDR is therefore something a new engine can win on outright, not just a box to tick.
- **MBR/identity propagation is the weakest statistical link.** On a single-cell dataset, at a nominal 1% PIP FDR, the measured false-transfer FDP was 3.4% for IonQuant, 6.1% for FlashLFQ and 2.0% for MaxQuant. PIP-ECHO (JPR 2025) was the only method that controlled it. DIA "MBR" (a two-pass empirical library) has the same problem in a different form.
- **Open, GPU-capable alternatives exist but are slower or narrower.** AlphaDIA is Apache-licensed and uses numba plus PyTorch transfer learning, but predicted-library search on Bruker data is slow. DIA-BERT needs a GPU and is Orbitrap-only. Beta-DIA uses a GPU and targets diaPASEF. CHIMERYS (LASSO deconvolution) is commercial. Sage is MIT-licensed, has LFQ with decoy-controlled MS1 integration and a wide-window mode, but is not a real DIA engine.
- **Benchmarks:** we should standardise on ProteoBench's HYE modules (Astral, diaPASEF, ZenoTOF, low-input, plasma-PYE), the ProteoBench Astral entrapment module, LFQbench (PXD002952), Van Puyvelde (PXD028735), PYE plasma (PXD056598, 1,116 runs across 12 sites) and two-proteome MBR sets (PXD014415).

---

## 1. Computational anatomy of DIA search

DIA analysis has to (i) generate candidate precursors (from a library or predicted from FASTA), (ii) extract fragment and precursor XICs (plus ion mobility on timsTOF), (iii) find and score peak groups, (iv) estimate FDR with decoys, then (v) do cross-run steps: library refinement, MBR, protein inference and quantification. The paradigms are:

| Paradigm | Representative tools | Core idea | Cost scaling |
|---|---|---|---|
| Peptide-centric, targeted extraction | DIA-NN, Spectronaut, OpenSWATH, EncyclopeDIA, MaxDIA, Skyline | For each library precursor, extract XICs in the RT/IM window and score co-elution and library similarity | ∝ (#precursors × RT-window width × fragments); RAM ∝ library (DIA-NN: ~0.5 GB per 1M precursors) |
| Spectrum-centric (pseudo-MS/MS or direct) | MSFragger-DIA, diaTracer, DIA-Umpire | Search each DIA MS2 spectrum (or a deconvolved pseudo-spectrum) with a fragment index, then rescore peptide-centrically | ∝ #spectra × candidates per index lookup; insensitive to library size, so suited to open/PTM searches |
| Feature-free / convolutional | AlphaDIA | Aggregate signal over RT × IM × fragments with learned convolution kernels; no centroiding or feature finding | Dense tensor ops; GPU-friendly, but slow with large predicted libraries on Bruker |
| Linear deconvolution of chimeric spectra | CHIMERYS | Treat each MS2 as a non-negative LASSO mixture of predicted spectra for candidates | ∝ #spectra × candidate set; needs accurate RT and intensity prediction |
| Learned end-to-end scorers | DIA-BERT, Beta-DIA, DIA-CLIP, DreamDIA | A transformer/CNN scores the extracted peak-group tensor; pretrained across many runs | GPU required; pretrained scorers aim to remove per-run semi-supervised training |
| Index-based massive pre-search | DIA-NN InfinDIA | A "novel spectral representation" scores peptides against all spectra of a run at very high rate; the output is an empirical library for a normal search | 2×10⁸–4×10⁸ targets/min on 64 cores (Orbitrap/Astral) |

Data volume for scale: Astral narrow-window DIA (nDIA) runs at ~200 Hz MS/MS, so a 30-min gradient yields ~360k MS2 scans (2 Th windows, 300 windows over m/z 380–980). On timsTOF diaPASEF every frame adds a dense IM dimension, which is why DIA-NN's InfinDIA rate drops 4–8× on PASEF data (51 M/min at 60 SPD and 29 M/min at 70-min gradients, versus 215 M/min on Astral at 40 SPD).

---

## 2. Tools

### 2.1 DIA-NN (Demichev et al., Nat Methods 2020, doi:10.1038/s41592-019-0638-x)
**Algorithm.** Peptide-centric extraction. For each precursor it picks the best co-eluting fragment as an "elution reference", which removes interference. It computes ~70+ peak-group scores and uses an ensemble of fully connected neural networks as the target-decoy classifier. In the default "NNs (cross-validated)" mode, each network is trained for 4 epochs on a PSM subset and scores only PSMs dissimilar to its training set; this takes 15–30 s per ~100k-precursor run. "NNs (fast)" trains for 1 epoch on everything and takes a few seconds. **Library-free mode** predicts spectra, RT and IM from FASTA with built-in deep-learning predictors (<2 min per million precursors on a 16-core desktop). **MBR** is a two-pass procedure: pass 1 builds an empirical library from all runs, and pass 2 re-searches every run against it. For large studies the docs recommend building the empirical library from **20–100 high-quality runs**, then searching all runs with MBR off.

**Version history (GitHub releases).** 1.8 (Jun 2022): improved library-free mode, Linux build. 1.8.1 (Apr 2023): plexDIA and peak-height quant. 1.9 (Jun 2024): peptidoform scoring, QuantUMS, Slice/midia/Synchro-PASEF. 1.9.2 (Oct 2024): redesigned NN classifier, Astral calibration, "Conservative ML", `--nn-fold 4`. **2.0 (Jan 29 2025):** proteoform confidence mode (decoys made by mutating one residue, so q-values reflect correct precursor m/z and PTM set), fine-tuning of DL predictors for PTMs, Scanning-PASEF, non-specific digests. **2.3:** InfinDIA and beta DDA support. **Current: 2.7.0.**

**Licensing changed in 2.x.** "DIA-NN Enterprise" (full functionality, paid via Aptila Biotech / Thermo Fisher) versus "DIA-NN Academia" (limited functionality, non-profit only). Enterprise-only features include the **Knowledge base**, which scores PSMs against prior data and is "a hybrid library + library-free approach". The docs say it must be **turned off for entrapment validation**. CPU affinity control and some QuantUMS metrics (`MaxLFQ.Empirical.Quality`) are also Enterprise-only. DIA-NN is closed source and Windows/Linux only.

**Speed claims.** Up to "1000 runs per hour" (README headline; this presumably means an empirical library, short gradients and a large server). InfinDIA benchmark (DIA-NN 2.3.0 Enterprise, AMD 7980X with 64 cores and 128 threads, fixed mass accuracy, calibration library):

| Data | Speed (M target+decoy pairs/min) |
|---|---|
| QE HF-X 65 min, K-GG, PXD005573 | 440 |
| QE HF 30 min, PXD005573 | 315 |
| Astral 30 min phospho, max 5 var mods, MSV000093613 | 275 |
| Astral 40 SPD, PXD044991 | 215 |
| Astral 40 SPD, **unspecific digest** | 210 |
| timsTOF 60 SPD, PXD017703 | 51 |
| timsTOF 70 min, PXD034222 | 29 |

The docs also say InfinDIA q-values "may deviate from externally controlled FDR estimates by a factor of up to approximately 1.5–2x in either direction". Speed is bought with calibration accuracy, and InfinDIA is meant to be followed by a regular two-pass search.

**Peptidoform FDR.** The docs state that without dedicated peptidoform scoring, "the effective peptidoform FDR may reach 5–10% for library-free analyses". That number matters for our PTM goals.

**Large-cohort workflow (official guidance).** (1) Build an empirical library from a subset of runs. (2) Search each batch with fixed mass accuracies and scan window. (3) Train QuantUMS on 10–100 runs (`--quant-train-runs`, `--quant-params` for reuse). (4) Run cross-run steps only with "Reuse .quant files". In MBR mode the pass-2 .quant data are held in RAM. The authors call MBR "just a convenience feature" for large studies.

### 2.2 QuantUMS (Kistner, Grossmann, Sinn, Demichev; Nat Biotechnol 2026, doi:10.1038/s41587-026-03131-2; bioRxiv 10.1101/2023.06.20.545604)
- Models the log-variance and log-bias of each MS1 and fragment signal as a function of intensity and a correlation-based quality score. Precursor quantities are variance-weighted aggregates of those signals. Protein quantities come from a MaxLFQ variant that uses a **weighted median of precursor log-ratios**, with weights per run pair.
- Hyperparameters are fitted by autodiff with Armijo backtracking. The fit reflects the LC-MS setup and matrix rather than the experimental design, so parameters transfer between experiments with the same setup.
- **Linear complexity in number of runs**, versus quadratic/cubic for classical MaxLFQ.
- Reported gains: >1.4× lower log2-SD for low-intensity precursors, and 1.3× (human) and 1.2× (E. coli) for proteins. A "high-accuracy" mode removes most ratio compression. Benchmarks: PXD040205 (K562 + E. coli mixes) and PXD022216 (50 CLL samples).
- Outputs per-quantity quality scores (`Quantity.Quality`, `PG.MaxLFQ.Quality`) for filtering. **Limitation (docs):** it cannot remove ratio compression for low-abundance analytes at 200–500 SPD nanogram loads. With mixed bulk and single-cell runs, training favours the majority type.

### 2.3 Spectronaut / directDIA (Biognosys; commercial, closed)
- Industry standard alongside DIA-NN. directDIA (library-free via the Pulsar search) and "directDIA+". Version 19 (2024) added IM-aware timsTOF processing.
- **Spectronaut 20 (May 2025):** new **Kuiper** search engine with proprietary predictors, "70% faster runtimes vs 19" and "~8 min per raw file on average" for HLA-I unspecific searches. Also Linux and cloud support.
- Algorithm details are sparsely published, so it is hard to reproduce. Missing-value handling (q-value "sparse" vs "percentile" imputation-like reporting) affects benchmark CVs and completeness. CHIMERYS's LFQbench comparison found Spectronaut's completeness fell to 61% only after an extra empirical-FDR filter. Wen 2025 found Spectronaut's FDP on single-cell data reached 3.8–7.6% (precursor) at a nominal 1%.
- In the FragPipe paper's runtime comparison, Spectronaut 17 took **337 min** versus FragPipe's 144 min for 20 ccRCC runs.

### 2.4 AlphaDIA (Wallmann, Skowronek et al., Mann lab; Nat Biotechnol 2025, doi:10.1038/s41587-025-02791-w; **Apache-2.0**)
- **Feature-free search.** Signals are aggregated over RT × IM × fragments with learned convolution kernels, with no centroiding or feature detection. By default 5 candidates are kept per precursor. Up to 47 features feed a small MLP (layers 100-50-20-5, ~10.8k parameters). FDR is count-based target-decoy at precursor and protein level.
- **End-to-end transfer learning.** Starting from fully predicted libraries (AlphaPeptDeep), it iteratively fine-tunes the RT, IM, MS2 and charge models on confident IDs from the data (max 50 epochs, 5 warm-up). On dimethyl-labelled data, RT R² rose from 0.69 to 0.99, median RT error fell from 317 s to 11 s, and precursors rose from 65k to 96k (+48%). This enables generic PTM DIA.
- **Numbers (preprint):** Astral HeLa 21 min (60 SPD): >120k precursors and ~9.5k protein groups. timsTOF empirical-library benchmark: 81.5k peptides and 7,649 PGs. Entrapment (Arabidopsis library): **0.2% precursor and 1.08% protein** FDP at nominal 1%; repeated transfer-learning rounds gave <0.5% false Arabidopsis IDs. Median CV 7.7% on HeLa dia-PASEF. Runtime stayed "below acquisition time" on a 32-core machine for rapid gradients.
- Quant via **directLFQ**. Supports MBR as a second search step and HPC distribution.
- **Weakness:** the docs warn that "feature-free search takes a long time with fully predicted libraries" on Bruker .d. Python/numba code is hard to compare with the C++ engines.

### 2.5 MSFragger-DIA, diaTracer, FragPipe-DIA (Nesvizhskii lab)
- **MSFragger-DIA** (Yu et al., Nat Commun 2023, doi:10.1038/s41467-023-39869-5) searches every DIA MS2 spectrum with the fragment index **before** feature detection. It enumerates precursor charges, runs a two-stage calibration search, then does peptide-centric peak tracing and rescoring with a greedy step that blocks peak re-use. Libraries are built with EasyPQP, rescored with MSBooster and Percolator, and **quantified with DIA-NN**.
- Speed: on 20 ccRCC runs FragPipe-MSF was fastest, with DIA-NN library-free 2× slower, MaxDIA 10× slower, and Spectronaut 17 at 337 vs 144 min. On phospho data it was ≥6× faster than DIA-NN library-free. Entrapment FDP over 20 runs: DIA-NN lib-free 1.2%, FP-MSF 1.9%, hybrid 1.8%, so the spectrum-centric route is slightly more liberal.
- Weaknesses: no diaPASEF in MSF-DIA itself, PTM localisation "an open question", and the speed advantage shrinks with many variable mods.
- **diaTracer** (Li et al., Nat Commun 2025, doi:10.1038/s41467-024-55448-8) does 3D (m/z, RT, IM) peak tracing on diaPASEF and builds precursor-resolved pseudo-MS/MS spectra. This allows library-free, **open/mass-offset PTM searches on diaPASEF**. It is integrated in FragPipe.
- Licensing: MSFragger, IonQuant and diaTracer are free for academic use but closed; commercial use is licensed (Fragmatics). FragPipe calls DIA-NN for quant, which inherits DIA-NN's licence.

### 2.6 EncyclopeDIA / Scribe / Carafe (Searle; MacCoss/Noble labs)
- **EncyclopeDIA** (Searle et al., Nat Commun 2018, doi:10.1038/s41467-018-07454-w) searches DDA spectrum libraries or DIA **chromatogram libraries** built by gas-phase fractionation (GPF). The Walnut/PECAN mode is library-free. Searle 2020 (Nat Commun, doi:10.1038/s41467-020-15346-1) introduced Prosit-predicted libraries empirically corrected with 6×GPF runs. **Scribe** (a DDA library search) has been in the main build since v2.12.30 (Dec 2022). Apache-2.0.
- In Wen 2025, EncyclopeDIA's peptide-level FDR looked fine but **protein-level FDP was 4.3–9.1% at a claimed 1%**, which shows protein inference/FDR is a separate failure point.
- **Carafe** (Wen et al., Nat Commun 2025, doi:10.1038/s41467-025-64928-4) fine-tunes AlphaPeptDeep-style models **directly on DIA data** to produce experiment-specific in-silico libraries. It is integrated into Skyline and avoids the DDA-trained predictor mismatch.

### 2.7 OpenSWATH + PyProphet + TRIC/DIAlignR (Röst/Aebersold lineage; OpenMS BSD-3, PyProphet BSD)
- OpenSWATH (Röst et al., Nat Biotechnol 2014, doi:10.1038/nbt.2841) does targeted extraction. PyProphet is the semi-supervised LDA/XGBoost scorer. Rosenberger et al. (Nat Methods 2017, doi:10.1038/nmeth.4398) defined **run-specific, experiment-wide and global** peptide/protein FDR contexts and a "global analyte constraint" to stop false positives accumulating across large datasets. This is the conceptual basis for repository-scale FDR.
- **TRIC** (Röst 2016, doi:10.1038/nmeth.3954) and **DIAlignR** (Gupta et al., MCP 2019) align chromatograms across runs. DIAlignR's hybrid global/local alignment matched plasma peaks with 98% accuracy versus 76% for the previous method. Multi-run MST/progressive alignment (2023, PMC10616189) gave a 3× lower quant error rate and 50% lower cross-site CV for abundant proteins over 11 LC-MS setups.
- Recent: native Bruker .d (TDF) input via OpenTIMS, and Parquet (.xipm/.oswpq) outputs. Mobi-DIK handles diaPASEF.
- Fröhlich 2022 found OpenSWATH needed TRIC to avoid false E. coli calls. Throughput is well below DIA-NN.

### 2.8 MaxDIA (Sinitcyn et al., Nat Biotechnol 2021, doi:10.1038/s41587-021-00968-7)
Part of MaxQuant: free, closed source. It supports library-based and library-free "discovery" modes with in-silico libraries, BoxCar-DIA and bootstrap DIA. In FragPipe's benchmark it was **~10× slower than FragPipe** and less sensitive than DIA-NN and FragPipe. It is the slowest mainstream option and rarely used for large cohorts.

### 2.9 Sage (Lazear, JPR 2023, doi:10.1021/acs.jproteome.3c00486; **MIT**)
A Rust fragment-index search engine. It does not have a dedicated peptide-centric DIA engine, but it has a **wide-window mode** ("enables WWA/PRM/DIA searches") and **chimeric search**. Its **LFQ** uses FlashLFQ-style MS1 extraction over all charge states and isotopologues with **peptide-specific time warping** across files. It scores normalised spectral angle, intensity and ΔRT, and **controls MS1-integration FDR with a decoy peptide for every quantified target**. That is a clean model for MBR FDR. Sage keeps all data in RAM and reads each file once. The last release is v0.14.7 (Mar 2024), so development has slowed. Wen 2025 found Sage's DDA FDR control valid with both shuffled and foreign entrapment.

### 2.10 CHIMERYS (Frejno et al., MSAID/TUM; Nat Methods 2025, doi:10.1038/s41592-025-02663-w; commercial)
- **Spectrum-centric, acquisition-agnostic deconvolution.** Each MS2 is modelled as a non-negative linear combination of Prosit-predicted spectra, solved with an L1-regularised LASSO ("explain the most intensity with the fewest precursors"). Targets and decoys get predicted RT and intensity. The same engine handles DDA, DIA and PRM.
- DDA HeLa 2 h: 238,795 PSMs; >85% of MS2 spectra got ≥1 PSM and >66% of identified spectra were chimeric. LFQbench DIA (per file): 529,993 PSMs, 66,888 peptide groups, 7,331 PGs; >82% of identified spectra contained >1 precursor. Run-specific FDR matched entrapment eFDR across 1.4–20.4 Th windows, while **DIA-NN and Spectronaut underestimated it**. Data completeness after eFDR correction: CHIMERYS 41%, DIA-NN 30%, Spectronaut 61%. Median CVs: CHIMERYS 26.9%, DIA-NN 29.1%, Spectronaut 29.2%. Speed: **4.9× faster than DIA-NN and 1.7× slower than Spectronaut** on LFQbench (8-core VM, 64 GB).
- Limitations: Thermo raw only at publication; modifications outside the predictor's training set are missed; isobaric peptides are reported as the top-scoring representative only. It is available via Proteome Discoverer ≥3.1 and the MSAID cloud.

### 2.11 Deep-learning-native newcomers
- **DIA-BERT** (Liu et al., Westlake/Guo lab; Nat Commun 2025, doi:10.1038/s41467-025-58866-4). An encoder-only transformer on peak-group matrices (330×16, MS1+MS2), pretrained on **276M precursors from 952 DIA files**. A quant head was trained on 34M synthetic precursors (Synthedia). Reports +51% proteins and +22% precursors over DIA-NN on 5 cancer cohorts, and is more conservative (precursor FDP 0.36% vs 1.1%). **Needs a GPU, and v1.0 is Orbitrap-only** (not recommended for timsTOF/TripleTOF). Open source.
- **Beta-DIA** (bioRxiv 10.1101/2024.11.19.624419; PyPI). Targets diaPASEF. Two DL scorers (DeepProfile for 2D RT×IM co-elution, DeepMall for spectrum similarity) feed an ensemble. Runs on an RTX 4090 GPU with 12 CPU cores versus DIA-NN on 24 cores, with similar runtime. Reports +24–33% peptides and +3–10% PGs at 1% **external** (Arabidopsis entrapment) FDR.
- **DIA-CLIP** (arXiv 2602.01772, Feb 2026). A **dual-encoder contrastive** peptide–spectrum representation for **zero-shot** PSM scoring with no per-run semi-supervised training. Claims up to +45% proteins and −12% entrapment hits. This is directly relevant to our bi-encoder plan, but the claims are unreplicated.

---

## 3. Tool comparison table

| Tool (latest) | Paradigm | Library-free | Accel. | Licence / source | Quant | Notable speed number | FDR caveats |
|---|---|---|---|---|---|---|---|
| DIA-NN 2.7 | Peptide-centric + NN ensemble; InfinDIA index pre-search | Yes (built-in DL predictors, fine-tunable) | CPU only (C++) | Closed; Academia (limited, free) vs Enterprise (paid) | QuantUMS, MaxLFQ, legacy Top-N; plexDIA | InfinDIA 210–440 M targets/min (Orbitrap/Astral), 29–51 M/min (timsTOF) on 64 cores; "≤1000 runs/h" | Single-cell FDP 2.3–4.7% at 1% (Wen 2025); peptidoform FDR 5–10% without peptidoform mode; Knowledge base incompatible with entrapment |
| Spectronaut 20 | Peptide-centric; Pulsar/Kuiper | directDIA+ | CPU (Win/Linux) | Commercial | Built-in MaxLFQ-like, imputation options | ~8 min/file (HLA, v20); 337 min / 20 runs (v17) | Single-cell FDP 3.8–7.6% (Wen 2025) |
| AlphaDIA (2.x) | Feature-free convolution + MLP; transfer learning | Yes (AlphaPeptDeep) | numba CPU; PyTorch GPU for TL | **Apache-2.0** | directLFQ; DIA multiplexing | < acquisition time (32 cores, short gradients); slow on Bruker with predicted libs | Entrapment 0.2% prec / 1.08% prot at 1% |
| FragPipe (MSFragger-DIA + diaTracer) | Spectrum-centric index search then peptide-centric rescoring | Yes (DIA → library → DIA-NN quant) | CPU (Java) | Free academic, closed core | DIA-NN quant / IonQuant | 144 min / 20 runs; 2× faster than DIA-NN lib-free | 1.9% FDP (entrapment, 20 runs) |
| CHIMERYS 4 | LASSO deconvolution with Prosit predictions | Yes (prediction-only) | Cloud/AWS | Commercial (PD/MSAID) | Fragment-level (PD) | 4.9× faster than DIA-NN on LFQbench | Matches eFDR across window widths |
| EncyclopeDIA / Carafe | Library search; GPF chromatogram libraries; DIA-trained predictors | Walnut / Prosit / Carafe | CPU (Java) | Apache-2.0 | Fragment XIC | — | Protein FDP 4.3–9.1% at 1% (Wen 2025) |
| OpenSWATH + PyProphet + DIAlignR | Targeted extraction; semi-supervised scorer | Needs library (predicted OK) | CPU (C++/Python) | BSD (open) | Fragment XIC + alignment | Slowest open pipeline | Global/experiment-wide FDR contexts defined here |
| MaxDIA | Library / discovery | Yes | CPU (.NET) | Free, closed | MaxLFQ | ~10× slower than FragPipe | — |
| Sage 0.14.7 | Fragment index; wide-window / chimeric | DDA-like only | CPU (Rust) | **MIT** | LFQ with decoy-controlled MS1 integration; TMT | Order of magnitude faster than DDA engines | DDA FDR valid (Wen 2025) |
| DIA-BERT 1.0 | Transformer peak-group scorer | Yes | **GPU required** | Open | Transformer peak area | — | Conservative (0.36% vs DIA-NN 1.1%) |
| Beta-DIA | DL co-elution + spectrum scorers (diaPASEF) | Yes | GPU (RTX 4090) | Open (PyPI) | — | ≈ DIA-NN with half the cores | Evaluated with external entrapment |

---

## 4. Instrument specifics that drive engine design

- **timsTOF diaPASEF** (Meier et al., Nat Methods 2020; DIA-NN IM module: Demichev et al., Nat Commun 2022, doi:10.1038/s41467-022-31492-0). There is a 4th dimension (1/K0), and windows are placed in m/z × IM space. ProteoBench's dataset uses 8 TIMS ramps × 3 windows of 25 Th, a 0.96 s cycle and 1/K0 0.64–1.37. Variants: **Slice-PASEF** (all ions, maximal sensitivity, from 1 frame), **Synchro-PASEF**, **midia-PASEF** and **Scanning-PASEF** (DIA-NN 2.0). Data are sparse 2D (scan × TOF) frames in TDF/SQLite with compressed binary. Extraction must handle IM peaks of a few scans. This is where generic engines lose 4–8× throughput (InfinDIA) and where AlphaDIA's predicted-library search is slowest. IM prediction (IM2Deep, AlphaPeptDeep CCS) is mandatory for library-free work.
- **Orbitrap Astral** (nDIA: Guzman et al., Nat Biotechnol 2024, doi:10.1038/s41587-023-02099-7). ~200 Hz MS/MS with 2-Th windows gives near-DDA-clean spectra. The paper reports ~10k proteins in 30 min, >170k precursors, and **14k+ protein groups in a 3-species mix with median CV <7% (nDIA) vs <19% (DDA)**, at throughputs from 180 SPD (5-min gradients) to 40 SPD. Data were processed with Spectronaut 17/18 and DIA-NN 1.8.1. MS1 spectra at 240k resolution are excellent, which is why QuantUMS gains most here. Low input: Bubis et al. (Nat Methods 2025) report up to 7,400 PGs from 250 pg HeLa at 50 SPD and 1,801–5,300 PGs per single A549 cell. **Astral Zoom** (2025) adds multi-pass resolving power for TMTpro 35-plex.
- **ZenoTOF 8600 / Zeno SWATH** and **Scanning SWATH** (Messner et al., Nat Biotechnol 2021) make the quadrupole position a scoring dimension in DIA-NN.

**Consequence for us:** narrow-window Astral DIA is close to "DDA with complete sampling", so a spectrum-centric, fragment-index GPU engine (MSFragger-DIA/CHIMERYS-like) becomes competitive there. diaPASEF needs 4D tracing (diaTracer-like) or IM-aware peptide-centric extraction.

---

## 5. Predictors (spectrum / RT / IM) and inference cost

| Model | Predicts | Size / speed (published) | Notes |
|---|---|---|---|
| Prosit (Gessulat et al., Nat Methods 2019, doi:10.1038/s41592-019-0426-7) | MS2 (HCD/CID, TMT), iRT | Prosit-Transformer: 24 min for 1.4M peptides (per AlphaPeptDeep paper) | Basis of CHIMERYS and Oktoberfest; mostly served via Koina |
| AlphaPeptDeep (Zeng et al., Nat Commun 2022, doi:10.1038/s41467-022-34904-3) | MS2 (transformer, ~4M params), RT (~0.7M), CCS (~0.7M) | **Human library: 2.6M peptides / 7.9M precursors in ~10 min on one GPU, ~100 min CPU**; MS2 for 1.4M peptides in 35 s (40× Prosit-Transformer) | Few-shot TL: 10 peptides → PCC90 48%→87%; 500 RT points → R² 0.927→0.986; TL "within minutes, even without GPU" |
| MS2PIP (Gabriels 2019; Declercq et al., NAR 2023) | MS2 (XGBoost / DL) | Fast CPU | Many fragmentation/instrument models |
| DeepLC (Bouwmeester et al., Nat Methods 2021, doi:10.1038/s41592-021-01301-5) | RT | CPU-fast; atomic-composition encoding generalises to unseen mods | Calibrated per run |
| IM2Deep (Devreese, Declercq et al., JPR 2025) | CCS, incl. multi-conformer | DeepLC-style encoding | Handles unseen mods and multiple conformers |
| Chronologer (Wilburn/Searle) | RT | — | Used in EncyclopeDIA/Skyline ecosystem |
| DIA-NN built-in | MS2, RT, IM | **<2 min per 1M precursors on 16 CPU cores**; fine-tunable in 2.x | Closed |
| Koina (Lautenbacher et al., Nat Commun 2025, doi:10.1038/s41467-025-64870-5) | Serves 30+ models (Prosit, MS2PIP, AlphaPeptDeep, DeepLC, IM2Deep, UniSpec, Chronologer) over Triton/KServe | **>200k predictions/s peak, median >30k/s; ~1.5 min per 2M spectra** for an average model (4×GTX1080 servers) | Network-bound; DIA-NN's local predictor was faster on a 36-core server |

**Cost arithmetic.** A human tryptic library (~5M precursors) costs ~10 min on CPU (DIA-NN) or ~5 min on one GPU (AlphaPeptDeep), so it is negligible when amortised across a cohort. **An open, semi-specific or multi-PTM space (10⁸–10⁹ precursors) costs 1–30 GPU-hours just to predict**. DIA-NN's InfinDIA avoids this by pre-searching without full predictions (it predicts only for the survivors). So prediction should be lazy: predict only for candidates surviving a coarse index or embedding filter, and cache predictions keyed by peptidoform + charge + NCE + instrument.

---

## 6. FDR in DIA: the evidence

- **Wen, Freestone, Riffle, MacCoss, Noble, Keich (Nat Methods 2025, 22:1454; doi:10.1038/s41592-025-02719-x; FDRBench, Apache).** Of the common entrapment estimators, one is invalid, one is only a lower bound, and the valid "combined" estimator is FDP̂ = N_E(1+1/r)/(N_T+N_E). They propose a tighter **paired** estimator for r=1. At a claimed 1%:
  - DIA-NN precursor/protein FDP: Astral 0.7–1.3% / 1.5–2.1%; Lumos 0.9–1.7% / 2.1–2.6%; **single-cell Eclipse 2.3–4.7% / 2.0–3.5%** (inflation up to ~48% precursor and ~44% protein reported for single-cell).
  - Spectronaut: Astral 0.8–1.5% / 1.6–2.3%; **single-cell 3.8–7.6% / 3.0–5.6%**.
  - EncyclopeDIA: peptide ~0.4–1.2%, **protein 4.3–9.1%**.
  - DDA tools (Tide+Percolator-RESET, Sage, MS-GF+) appeared valid.
  - DIA datasets used: PXD042704, PXD034525, PXD028735, PXD041421, PXD017703, PXD023325 (1-cell/100-cell Eclipse), PXD012988 (timsTOF), MSV000084000. The exact per-dataset mapping is in their Suppl. Table 1 *(verify mapping)*.
- **ProteoBench now runs this as a live module** ("DIA ion entrapment (Astral)", global precursor FDR, lower-bound/combined/paired FDP, 1:1 shuffled entrapment FASTA). We can submit directly to it.
- **The usual failure modes are protein-level and low-input.** Run-specific vs global q-values, protein-group FDR (picked-group FDR is preferable), MBR second-pass libraries (pass 1 is filtered at a global library q-value, pass 2 at run-specific q-values), and "knowledge base" or prior-data scoring all break entrapment assumptions.
- **Proteoform/peptidoform FDR:** DIA-NN 2.0's single-residue-mutation decoys target this explicitly. Localisation needs a separate module (see the PTM report).

---

## 7. Quantification

### 7.1 Protein roll-up
- **MaxLFQ** (Cox et al., MCP 2014, doi:10.1074/mcp.M113.031591) solves a least-squares problem over pairwise median peptide log-ratios. It is O(n²) in runs (O(n³) with a dense solve). Implementations: MaxQuant, DIA-NN, **iq** (Pham et al., Bioinformatics 2020, doi:10.1093/bioinformatics/btz961; a fast C++ MaxLFQ in R), and IonQuant.
- **directLFQ** (Ammar et al., MCP 2023, doi:10.1016/j.mcpro.2023.100581; Apache) normalises by "intensity-trace shifting" in log space, merging the most similar traces first by single-linkage clustering. It is **linear in samples: 10,000 DIA samples in ~16 min vs ~2.5 days for iq; 10,000 DDA samples in ~10 min vs >2 weeks for MaxLFQ; 100k samples in ~2–3 h** (64 cores, 128 GB). Accuracy: yeast SD 0.21 vs iq 0.26 vs Spectronaut 0.31 on PXD016647. CVs on the ~900-sample clinical set PXD029009 were similar to iq and DIA-NN. Weakness: per-sample *linear* shifts cannot correct nonlinear, RT-dependent distortions.
- **QuantUMS** (§2.2): variance-aware, linear complexity, with per-value quality scores.
- **Top-N / Top-1** with RT-binned normalisation is DIA-NN's recommended path when quantities of existing runs must **not change as runs are added**. That property matters for incremental repository-scale processing.

### 7.2 Normalisation
Global median/quantile normalisation, RT-dependent normalisation (DIA-NN's "RT-dependent" mode; the legacy recipe uses ~100 RT bins with median or half-sample-mode of the 40% least-variable precursors, then progressively coarser bins), directLFQ trace shifting, and VSN. Batch correction uses reference-run precursor factors. In multi-species benchmarks, normalisation must be anchored on the constant-ratio species (human in HYE). Otherwise it inflates errors.

### 7.3 Match-between-runs / identity propagation and its FDR
- **Lim, Paulo, Gygi (JPR 2019, doi:10.1021/acs.jproteome.9b00492; PXD014415).** 20 human-only versus 20 human+yeast runs. MaxQuant MBR increased IDs by ~40%, but **44% of identified yeast proteins were transferred into at least one yeast-free sample**. Only 2.7% survived MaxLFQ filtering.
- **IonQuant** (Yu, Haynes, Nesvizhskii, MCP 2021, doi:10.1016/j.mcpro.2021.100077). MBR FDR from LDA on 5 features (log intensity, isotope KL divergence, ppm error, ΔIM, ΔRT). **Decoy transfers are m/z-shifted by +11×1.0005 Th at the same RT/IM**, and an EM mixture model gives ion/peptide/protein-level MBR FDR. It was **19–38× faster than MaxQuant** (37 vs 699 min HeLa; 28 vs 1,056 min two-proteome), with peptide coverage up from 57% to 79% at 1% MBR FDR. Datasets: PXD014415, PXD010012, PXD014777 (3-organism timsTOF), MSV000085230.
- **PIP-ECHO** (Solivais et al., JPR 2025, doi:10.1021/acs.jproteome.5c00065; FlashLFQ). Two competitions: target vs decoy *peptide* transfers (identification errors) and predicted-RT vs **randomised-RT** peaks (peak-matching errors). Measured FDP at a nominal 1% on single-cell data: **PIP-ECHO 1.0%, MaxQuant 2.0%, IonQuant 3.4%, FlashLFQ v1 6.1%**. At 5% FDP it found 47% more DA proteins than MaxQuant. It also found that existing tools sometimes give *fewer* discoveries with MBR on, because transfer errors accumulate. Data: PXD014415, PXD003881, PXD037527, PXD057758.
- **DIA analogue.** DIA-NN's two-pass MBR re-searches every run against a cohort-wide empirical library, so "transfer" errors appear as library precursors that pass run-level q-values in runs where they are absent. Two-proteome DIA designs (Fröhlich 2022, Bubis HeLa/yeast, PYE) measure this directly. Wen 2025 shows the worst inflation at low input, which is exactly where MBR-style borrowing is used most. **No DIA tool reports a dedicated, validated "transfer FDR".**

### 7.4 Cross-run RT alignment
Approaches: global LOWESS/spline to a reference (MaxQuant, DIA-NN library RT recalibration); peptide-specific warping (Sage); chromatogram-level hybrid alignment with MST/progressive multi-run trees (DIAlignR; 3× lower quant error); optimal-transport alignment (Alignstein, PMC9633278). For thousands of runs, a star alignment to a library iRT scale is O(n), whereas pairwise schemes are O(n²). Alignment errors turn directly into false transfers.

### 7.5 Isobaric (TMT) and non-isobaric multiplexing
- TMTpro 18-plex is standard. **TMTpro 32/35-plex** adds ²H-isotopologues 3 mDa apart, which need >90k resolving power at m/z 130 (Orbitrap MS3 at 90k or Astral Zoom "TMT HR mode"). A 32-plex with RTS-SPS-MS3 gave the same protein count as 18-plex at 2× throughput. Thermo app notes report 4,299 PGs at 440 samples/day on Astral Zoom *(vendor data)*.
- Main quant issue: **ratio compression from co-isolation** (SPS-MS3, RTS and narrow isolation mitigate it). Software: TMT-Integrator/FragPipe (bioRxiv 2025.05.27.656447), Sage (built-in TMT), MSstatsTMT. DIA-NN does not yet support isobaric reporter quant (planned).
- **plexDIA** (Derks et al., Nat Biotechnol 2023, doi:10.1038/s41587-022-01389-w): mTRAQ/dimethyl/SILAC with DIA. DIA-NN uses a **decoy channel** to estimate channel-assignment FDR (`Channel.Q.Value`). AlphaDIA also supports DIA multiplexing.

---

## 8. Benchmark designs and metrics

### 8.1 Metrics
- **Precision:** CV per precursor/protein over replicates (report median and the fraction <20%), always together with the number of quantified analytes and completeness. A CV computed only on complete cases rewards tools that drop hard analytes.
- **Accuracy:** log2 fold-change error against the expected ratio per species. ProteoBench's **epsilon** = observed − expected log2FC per precursor, summarised as mean/median |ε|, plotted against #precursors quantified (with a min-runs filter). LFQbench reports accuracy (median deviation), precision (SD/IQR of log-ratios) and species separation.
- **Differential-call ROC / pAUC:** treat spike-in species as positives and the constant species as negatives. Report TPR at an empirical FDP (e.g. 1%/5%), not at nominal q. Fröhlich 2022 found non-parametric tests (SAM, ROTS) beat parametric ones.
- **False-transfer / false-detection rate:** count foreign-species IDs in samples lacking that species (Lim; PIP-ECHO; Fröhlich "human-only" samples).
- **Entrapment FDP** (combined + paired; FDRBench / ProteoBench), at precursor, peptidoform, protein and global levels.
- **Runtime:** wall-clock per run plus cross-run steps, with hardware, threads and RAM, separating library prediction, search and quant (as in DIA-NN's InfinDIA table).

### 8.2 Benchmark datasets

| Dataset | Design | Instrument(s) | Runs | Accession | Use |
|---|---|---|---|---|---|
| LFQbench (Navarro et al., Nat Biotechnol 2016, doi:10.1038/nbt.3685) | HYE: A = 65% H / 30% Y / 5% E; B = 65/15/20 | TripleTOF 5600/6600, 32- and 64-window SWATH | 3+3 per setup | PXD002952 | Classic DIA quant accuracy |
| Van Puyvelde et al. (Sci Data 2022, doi:10.1038/s41597-022-01216-6) | HYE A/B; log2FC H 0, Y −1, E +2 (ProteoBench convention) | TripleTOF 5600/6600+, QE HF-X, Synapt G2-Si/XS, timsTOF Pro; DDA and DIA | >700 | PXD028735 | Cross-platform; ProteoBench AIF module |
| ProteoBench Astral (2 Th nDIA) | Same HYE design, 50 ng, ~16 min μPAC | Orbitrap Astral | 3+3 | "not-yet-released" per docs *(check)* | Epsilon vs depth |
| ProteoBench diaPASEF | HYE, 25 ng, 30 min | timsTOF SCP | 3+3 | via ProteoBench | IM-aware quant |
| ProteoBench ZenoTOF | HYE | ZenoTOF 8600 | — | PXD070049 | Zeno SWATH |
| ProteoBench low-input Astral / timsTOF | HYE, 200 pg | Astral; timsTOF Ultra 2 | — | PXD080610 | Low-input quant |
| Bubis et al. (Nat Methods 2025) | HeLa:yeast 240:10 vs 200:50 pg | Astral, FAIMS | 3+3 | PXD049412 | Single-cell-like two-proteome |
| Guzman nDIA (Nat Biotechnol 2024) | 3-species + HeLa + gradients | Orbitrap Astral | many | PXD046453, PXD046372, PXD046444, PXD046386, PXD046417, PXD046357, PXD046283, PXD046285 | Astral speed/depth reference |
| PYE plasma (Distler et al., Nat Commun 2025, doi:10.1038/s41467-025-64501-z) | Human plasma + low yeast/E. coli spikes (ProteoBench: Y log2FC −1.585, E +1) | 12 sites, DDA and DIA, multiple vendors | **1,116** | PXD056598 | High-dynamic-range, multi-site; protein CV 3.3–9.8% (DIA) |
| Fröhlich et al. (Nat Commun 2022, doi:10.1038/s41467-022-30094-0) | FFPE lymph nodes from many patients + E. coli 1:6/1:12/1:25/none | Orbitrap Eclipse | 92 runs, 1,428 workflows | EGAD00010002223 (controlled access); Zenodo 6379087 | Real inter-patient heterogeneity; DE-call benchmarking |
| Lou et al. (Nat Commun 2023, doi:10.1038/s41467-022-35740-1) | Benchmark mimicking biological complexity, proteome + phospho | Orbitrap + timsTOF | — | see paper | Library × software matrix |
| QuantUMS sets | K562 + E. coli ratios; 50 CLL samples | timsTOF / Orbitrap | — | PXD040205; PXD022216 | Quant precision vs ratio compression |
| directLFQ sets | Mixed-species DIA; ~900-sample clinical; 200 HeLa replicates (DDA) | — | — | PXD016647; PXD029009; PXD015087 | Roll-up scaling and accuracy |
| Two-proteome MBR (Lim et al. 2019) | 20 human vs 20 human+yeast | Orbitrap | 40 | PXD014415 | False-transfer rate |
| IonQuant 3-organism | HYE timsTOF | timsTOF Pro | 6 | PXD014777 | MBR FDR on PASEF |
| Wen 2025 entrapment sets | Human / yeast / mouse / single-cell DIA | Astral, Lumos, QE, Eclipse, timsTOF, TripleTOF | — | PXD042704, PXD034525, PXD028735, PXD041421, PXD017703, PXD023325, PXD012988, MSV000084000 | FDR validation |
| DIA-NN InfinDIA speed sets | K-GG, phospho, Astral 40 SPD, timsTOF 60 SPD / 70 min | various | — | PXD005573, MSV000093613, PXD044991, PXD017703, PXD034222 | Speed comparison against DIA-NN's own numbers |
| CHIMERYS sets | LFQbench, HeLa DDA, etc. | various | — | PXD015087, PXD018171, PXD014877, PXD013615, PXD023218, PXD028735, PXD046453 | Chimeric deconvolution |
| Large cohort: yeast KO proteomes (Messner et al., Cell 2023, doi:10.1016/j.cell.2023.03.026) | 4,699 KO strains, ~2,520 proteins each, >100M peptide quantities | SWATH / DIA-NN | ~thousands | PXD036062 *(verify)* | Cohort-scale scaling/batch test |

---

## 9. Cohort-scale DIA processing

- **quantms** (Dai et al., Nat Methods 2024, doi:10.1038/s41592-024-02343-1; nf-core, open). Reanalysed **83 PRIDE datasets, 29,354 raw files, 13,132 human samples** (DDA-LFQ, TMT, DIA-LFQ via DIA-NN), giving 16,599 proteins and 1.03M peptides. **quantmsdiann** (Research Square preprint rs-10319687) runs DIA-NN per file in parallel under SDRF control. It reanalysed **2,300 single-cell runs in 2.2 h on 300 HPC nodes**, about 660 node-hours or **~17 node-minutes per run**, with results identical to single-machine DIA-NN. Upgrading DIA-NN 1.8.1 to current gave up to +17% PGs.
- **Standard scalable recipe (DIA-NN docs):** empirical library from 20–100 runs → per-run search (embarrassingly parallel) → cross-run steps on .quant files → QuantUMS trained on 10–100 runs. The cost is **two passes over raw data** plus library generation. MBR mode keeps pass-2 data in RAM, which limits scale.
- **The global FDR problem at repository scale:** per-run 1% FDR across 10⁴ runs gives an ever-growing union of false precursors. Rosenberger's global analyte constraint, picked-protein FDR and FDRBench-style entrapment are needed. Library-from-subset strategies also bias against analytes that occur only in rare conditions.
- **Roll-up scale:** directLFQ (linear, 10⁴ samples in ~16 min) and QuantUMS (linear) replace O(n²) MaxLFQ. "Quantities must not change when runs are added" requires fixed-reference designs (Top-1 against reference runs, or frozen QuantUMS parameters).

---

## 10. Weaknesses and open problems we can exploit

1. **FDR validity is not guaranteed** for DIA-NN, Spectronaut or EncyclopeDIA, especially at low input and at protein level. InfinDIA's q-values are off by up to 1.5–2×. An engine that reports **entrapment-validated (paired-FDP) FDR at precursor, peptidoform, protein and global levels by default** is a real differentiator.
2. **MBR / transfer FDR has no principled DIA treatment.** The DDA solutions (IonQuant m/z-shift decoys, Sage per-target decoys, PIP-ECHO decoy peptides plus random-RT competition) carry over to DIA's second pass and to cross-run borrowing.
3. **Closed, license-gated reference tools.** DIA-NN 2.x has moved advanced features (Knowledge base, affinity control, some QuantUMS metrics) to a paid Enterprise edition. Spectronaut and CHIMERYS are commercial. MSFragger is closed. Only AlphaDIA, EncyclopeDIA, OpenSWATH, Sage and the DL newcomers are open, and none matches DIA-NN's speed. **A fast, open (Apache/MIT) DIA engine is an unfilled niche.**
4. **No GPU in mainstream tools.** DIA-NN, Spectronaut and FragPipe are CPU-only. DIA-BERT and Beta-DIA use GPUs for scoring, but are not faster end to end, and DIA-BERT is Orbitrap-only.
5. **timsTOF is everyone's slow path.** It is 4–8× slower in InfinDIA, slowest in AlphaDIA predicted-library mode, and needed a separate tool (diaTracer) for spectrum-centric FragPipe. Sparse 4D extraction on GPU is open ground.
6. **Prediction cost explodes for open/PTM spaces.** Eager prediction of 10⁸–10⁹ precursors costs GPU-hours to GPU-days. Lazy, filter-first prediction with caching is not standard (InfinDIA's internals are unpublished).
7. **Two passes over raw data** (empirical library, then re-search) dominate I/O at repository scale. A single pass that stores compact per-run candidate summaries (DIA-NN's .quant is a proprietary version of this) would halve I/O.
8. **Per-run semi-supervised training** (NN ensembles, Percolator-like) is slow, overfits at low input, and gives run-dependent calibration. Pretrained scorers (DIA-BERT, DIA-CLIP) point to globally calibrated, percolator-free scoring, which suits repository-wide FDR.
9. **Quant accuracy at low input:** ratio compression persists (QuantUMS docs). Interference-aware and deconvolution-based (CHIMERYS-like) quant can still improve.
10. **Benchmark hygiene:** CV-on-complete-cases, run-specific q-value reporting and imputation-like defaults make published comparisons hard to reproduce. ProteoBench now gives a neutral, public leaderboard.

---

## 11. Implications for our design

1. **Two-tier DIA search: index pre-search, then full extraction for survivors.**
   - Tier 1 (GPU): fragment-index or embedding retrieval over all MS2 spectra of a run (Astral/Orbitrap/Zeno). For PASEF, run it over diaTracer-style 4D-traced pseudo-spectra. It emits candidate (precursor, RT, IM) tuples.
   - Tier 2 (CPU or GPU): DIA-NN-style co-elution scoring only for those candidates.
   - **Target to beat:** DIA-NN InfinDIA at 2–4×10⁸ targets/min on 64 cores (Orbitrap/Astral) and 3–5×10⁷ (timsTOF). Our GPU tier must be clearly faster per dollar on timsTOF, where the CPU bar is lowest.
2. **Lazy, cached prediction.** Predict MS2, RT and IM only for Tier-1 survivors. Cache by (peptidoform, charge, NCE, instrument class). Run few-shot fine-tuning per instrument/cohort, not per run (AlphaPeptDeep: 10–500 examples suffice; Carafe: train on DIA, not DDA). Budget: human library ≤5 GPU-min; amortise over all runs of a cohort.
3. **Globally calibrated scoring.** Use a pretrained cross-encoder over peak-group tensors (DIA-BERT-like, but multi-instrument) with a light per-cohort calibration layer rather than per-run NN training. Validate against DIA-CLIP's zero-shot claim. This supports repository-level global FDR (Rosenberger's global constraint plus picked protein/PTM-site FDR) without per-run semi-supervised drift.
4. **FDR validation built in.** Ship FDRBench-style paired entrapment (r=1) as a first-class mode. Report combined and paired FDP at precursor, peptidoform, protein-group and global levels, and submit to ProteoBench's Astral entrapment module. Do **not** use prior-data "knowledge base" scoring in validated runs.
5. **Transfer-FDR for cross-run borrowing.** For any MBR/second-pass/library-from-subset step, generate (a) decoy-peptide transfers, (b) m/z-shifted or random-RT decoy peaks (IonQuant / PIP-ECHO), and (c) per-target MS1 decoys (Sage). Report a separate transfer q-value. Evaluate on PXD014415, the Bubis HeLa/yeast set and PYE "species-absent" samples.
6. **Single pass over raw data.** Convert each raw file once into a columnar, zero-copy store holding centroided fragments, an IM-aware index and per-candidate XIC summaries. Cross-run steps (library refinement, MBR, quant) then run on these summaries, never on raw data again, which generalises DIA-NN's .quant reuse. This also enables incremental addition of runs with frozen quant parameters.
7. **Quant stack.** (a) Fragment/MS1 signal selection with a QuantUMS-style learned variance model, trained once per LC-MS setup and frozen. (b) Linear-time roll-up: directLFQ-style trace shifting or weighted-median MaxLFQ with linear-cost reference anchoring, plus RT-dependent normalisation. (c) Per-value quality scores. (d) An option where existing quantities are invariant to added runs. Isobaric: SPS/RTS-aware TMT reporter extraction with co-isolation (interference) correction. plexDIA with decoy channels.
8. **Benchmark gate before any speed optimisation lands:**
   - ProteoBench HYE modules (Astral, diaPASEF, ZenoTOF, low-input): |ε| vs #precursors.
   - LFQbench PXD002952 and Van Puyvelde PXD028735 (cross-vendor).
   - PYE PXD056598 (plasma dynamic range, multi-site).
   - Entrapment (Wen sets and the ProteoBench Astral entrapment module).
   - MBR false-transfer sets (PXD014415, PXD049412).
   - Speed against DIA-NN's published InfinDIA table (same accessions: PXD044991, PXD017703, PXD034222, PXD005573) on stated hardware.
   - Metrics: median CV with the number quantified, |ε|, DE-call pAUC at empirical FDP, entrapment FDP, wall-clock and cost per run.
9. **Licensing strategy.** We cannot embed DIA-NN, MSFragger, Spectronaut or CHIMERYS. Use DIA-NN only as an external comparator (Academia edition, documented version). Build on Apache/MIT components (AlphaPeptDeep/AlphaBase, Sage ideas, directLFQ, FDRBench, OpenMS file I/O, timsrust/OpenTIMS for TDF).
