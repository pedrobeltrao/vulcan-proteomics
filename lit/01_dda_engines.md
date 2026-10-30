# 01. DDA database search engines and how they are accelerated

*Literature review for the ultra-fast repository-scale proteomics engine. Compiled 2026-10-26. Numbers come from primary sources (DOIs/URLs inline and in the reference list). Items marked **[unverified]** come from background knowledge and were not re-checked against the source during this review. Derived throughputs (spectra/s) are my own arithmetic from the published runtimes and are labelled as **derived**.*

---

## 0. Bottom line (TL;DR)

1. **The fragment-ion index is the dominant CPU design.** MSFragger introduced it (2017), and Sage, Comet-FI, AlphaPept and MSFragger-DDA/DIA/Glyco/Labile all use it. On CPU it gives 5–150× speed-ups over classical precursor-window scoring, and the gain grows with search-space width (open search ≈150× vs Comet/X!Tandem). Its weakness is **memory**. The index grows with peptidoforms (MSFragger: 1.6 GB tryptic human → 55.8 GB semi-tryptic → 86.6 GB phospho). Sage and MSFragger both **run out of memory (OOM)** on metaproteomics- or proteogenomics-scale databases, where the streaming Tide v4.3 succeeds with 2–20 GB.
2. **GPUs have not yet beaten a well-engineered CPU fragment index for closed search.** In AlphaPept the GPU PSM step (77 s) was slower than multithreaded CPU (9 s) because the step is memory-bound. GiCOPS gets only 1.2–2× over its CPU predecessor for closed search, against 3–5× for open search. GPUs win in three places: **wide/open searches**, **dense similarity kernels** (HOMS-TC on tensor cores is 31× faster than ANN-SoLo-GPU), and **neural encoders/rerankers**.
3. **Learned scoring is where the sensitivity gains are.** Casanovo-DB detects 31–40% more peptides than XCorr and 88–102% more than Andromeda at 1% FDR before Percolator. Chimeric/full-isolation-window search (MSFragger-DDA+, CHIMERYS) adds 15–57% more peptides, and >2/3 of DDA MS2 spectra are chimeric. Exhaustive cross-encoder scoring is intractable (estimated ~2–12 days per 10⁶ spectra against 10⁸ peptides). The practical architecture is **index retrieval → top-K neural rerank**.
4. **Bi-encoder retrieval is not yet proven as a primary candidate generator at full-proteome sensitivity with rigorous FDR.** The 2026 JPR framework paper (Kalogeropoulos et al.) estimates biencoder+ANN at ~2–3 h query for 10⁶ spectra × 10⁸ peptides. That is **no faster than a fragment index (~1–2 h)** and uses 40–80 GB. Embeddings are clearly valuable for **clustering/deduplication across repositories** (GLEAMS: 669 M spectra) and for **spectral-library open search**.
5. **At PRIDE scale, I/O dominates.** PRIDE held 2.85 PB in 42,036 datasets in Aug 2024. Converting Thermo RAW and Bruker TDF to mzML is a large share of wall time. Columnar formats are converging on **Parquet/Arrow**: mzPeak (JPR 2025, Rust reference, reads Thermo/Bruker directly), Sage's parquet output, ThermoRawFileParser 2.0 parquet, and Casanovo v5 moving HDF5→Lance. We should convert once into a centroided, columnar, memory-mappable store and never touch mzML in the hot path.
6. **Licensing.** Sage (MIT, Rust; active again with commits in Oct 2026), Comet (Apache-2.0 per authors), Crux/Tide (Apache-2.0), Casanovo/InstaNovo (Apache-2.0), timsrust (MIT), ThermoRawFileParser (Apache-2.0), ANN-SoLo/HOMS-TC (Apache-2.0) and PTM-Shepherd (Apache-2.0) are reusable. MSFragger is **academic-only**, with commercial licensing via Fragmatics LLC since May 2025, so we can benchmark against it but not embed it. GiCOPS is CC BY-NC-ND. MS-GF+ is UC academic/non-profit only. MaxQuant, Mascot, CHIMERYS, ionbot and pFind are closed.

---

## 1. Classical and indexed DDA search engines

### 1.1 MSFragger family (Nesvizhskii lab)
- **Core paper:** Kong et al., *Nat Methods* 2017, doi:10.1038/nmeth.4256.
  - **Index:** an inverted fragment-ion index. Theoretical b/y fragments are binned by m/z; each bin lists peptide IDs sorted by precursor mass, so a precursor-tolerance window becomes a contiguous slice within every fragment bin. The bin width is chosen dynamically from the precursor tolerance. Peptide IDs are 32-bit, giving a **~2×10⁹ peptide limit**. Large searches are split into index "slices".
  - **Efficiency:** "~1.5 comparisons per candidate peptide", and >80% of fragment comparisons contribute to the score.
  - **Memory (human UniProt, tryptic, 1 missed cleavage):**
    - base index: 1.6 GB
    - +Met-ox, N-term acetyl: 2.9 GB
    - HLA non-specific 9–11-mers: 22.6 GB
    - semi-tryptic: 55.8 GB
    - phospho: 86.6 GB
  - **Speed (HEK293, 41,820 spectra, quad-core Xeon E3-1230v2):**

    | Search | MSFragger | Comet | X!Tandem |
    |---|---|---|---|
    | Narrow (100 ppm) | 0.4 min | 1.7 min | 1.7 min |
    | Open (500 Da) | 5.4 min | 815 min | 976 min |

    The open search is ≈150× faster than Comet. *Derived:* ≈130 spectra/s on 4 cores for a 500 Da open search, and ≈1,700 spectra/s narrow.
  - **IDs on the same data:** narrow search gave 456,548 PSMs and 115,755 peptides. Open search gave 609,897 PSMs (+34%). The full 24-run dataset ran in <30 min on one workstation.
  - **Score:** hyperscore = log(N_b!·N_y!·ΣI_b·ΣI_y), the X!Tandem-style score.
- **Mass-offset ("multinotch") search:** the search allows a list of discrete precursor offsets instead of a continuous window. It is cheap with a fragment index because the unmodified fragments are shared across offsets.
- **MSFragger-Glyco** (Polasky et al., *Nat Methods* 2020, doi:10.1038/s41592-020-0967-9) and **MSFragger-Labile** (Polasky et al., *MCP* 2023, PMC10182319) handle modifications whose fragment ions lose the mass. They search the peptide backbone with a precursor delta mass and add diagnostic/Y-ion and remainder-mass logic. Version 4.3 (Jun 2025) sped up Glyco/Labile for HCD/CID.
- **MSFragger-DIA** (*Nat Commun* 2023, PMC10338508) searches DIA MS2 spectra directly before feature detection. It is relevant as proof that the same index handles wide-window spectra.
- **MSFragger-DDA+** (Yu, Deng, Nesvizhskii, *Nat Commun* 2025, doi:10.1038/s41467-025-58728-z).
  - **Method:** it searches every MS2 against **all peptides in the full isolation window**, then runs MS1 XIC detection, shared-fragment removal and rescoring.
  - **Gains:** +57% peptide sequences vs MSFragger's standard DDA mode (Searle NCE data). Gains are +15% (AspN/GluC CID) to +45% (trypsin HCD) on Richards et al. data. It is similar to CHIMERYS on wide-window data.
  - **Speed:** 4× faster than MaxQuant/MetaMorpheus, and ≥20× faster than CHIMERYS+Proteome Discoverer (62 mzML files in 68.8 min).
  - **FDR:** checked with a shuffled entrapment database.
- **Recent versions:** 4.1 (May 2024) loads Thermo RAW much faster and adds DDA+ for ddaPASEF. 4.2 (May 2025) reorders search and post-processing to "significantly reduce memory" and adds Stellar low-res DIA. 4.4/4.4.1 (Dec 2025/Jan 2026). [MSFragger changelog](https://msfragger.nesvilab.org/CHANGELOG.html)
- **Licence:** academic/non-commercial only. Commercial licensing has been handled by Fragmatics LLC since 1 May 2025. Closed source (Java jar).

### 1.2 Sage (Lazear)
- **Paper:** Lazear, *J Proteome Res* 2023, 22:3652, doi:10.1021/acs.jproteome.3c00486. **Repo:** [github.com/lazear/sage](https://github.com/lazear/sage), **MIT, Rust**.
- **Status:** the last tagged releases are v0.14.7 (Mar 2024) and v0.15.0-beta.2 (2024). The repo was **active again in Oct 2026**, with commits adding named PTMs, radical z-ions for ETD/ECD and a CI benchmark harness ("sage-infra").
- **Index design** (from the Sage design blog, [lazear.github.io/sage](https://lazear.github.io/sage/)):
  1. Digest peptides, deduplicate them and sort by mass.
  2. Generate all b/y fragments and sort the whole fragment list by fragment m/z.
  3. Cut the list into fixed-size buckets (e.g. 16 fragments) and **re-sort each bucket by precursor mass**.
  4. This gives a two-level, B-tree-like array: the outer level does a binary search on fragment m/z, and the inner level does a binary search on precursor mass.
  - In the worked example, only ~1,500 of ~4,000 precursor-matched candidates are scored at all. This layout is very cache-friendly and a natural template for a GPU port.
- **Features:**
  - LDA rescoring, target-decoy FDR (picked peptide/protein), RT prediction, TMT quantification and LFQ.
  - "Chimera" mode (a second-pass search after subtracting the matched peaks) and wide-window/dynamic-tolerance mode for WWA/PRM/DIA.
  - Input: mzML(.gz), MGF (v0.14.6), Bruker .d via timsrust (ddaPASEF in v0.14.4; diaPASEF beta in v0.15).
  - Parquet output (v0.14.0), and streaming from S3/GCS/Azure.
  - Python bindings: **sagepy** (Rustims/imspy; Teschner et al., *JPR* 2025, PMC12053931; MIT).
- **Benchmarks** (AWS c5ad.8xlarge, 32 vCPU, 64 GB):

  | Search | Engine | Runtime | PSMs @1% FDR | Notes |
  |---|---|---|---|---|
  | TMT-MS2, PXD016766, 12 files | Sage | 18 s | 100,448 | |
  | | MSFragger | 101 s | 97,844 | |
  | | Comet | 345 s | 92,333 | |
  | | IdentiPy | 3,128 s | 96,465 | |
  | Open (−100…+500 Da), PXD001468 | Sage | 44 s | 24,682 | |
  | | Sage + chimera | 70 s | 28,724 | +16% PSMs |
  | | MSFragger | 225 s | 24,911 | |
  | Open + Met-ox, 2 missed cleavages | Sage | 68 s | 27,500 | 7.3 GB RAM |
  | | MSFragger | 340 s | 27,589 | 11.55 GB RAM |

- **Weaknesses:**
  - The whole index lives in RAM, so Sage **OOMs on huge databases** (fails on PXD011515 metaproteome, 7.1 B peptides; Tide 2025 paper).
  - Localization is limited to what the variable-mod enumeration gives; there is no explicit site-probability or FLR output **[unverified for latest commits]**.
  - The LDA rescorer is weaker than mokapot/Percolator plus predicted features.

### 1.3 Comet and Comet-FI (Eng; Schweppe lab)
- **Classic Comet** is a SEQUEST-XCorr reimplementation with E-values from the score histogram (Eng et al., *Proteomics* 2013, doi:10.1002/pmic.201200439). Apache-2.0 per the authors; [github.com/UWPR/Comet](https://github.com/UWPR/Comet).
- **Comet-FI** (McGann … Schweppe, bioRxiv 2024, doi:10.1101/2024.10.11.617953; Comet ≥ v2024.02.0) adds a fragment-ion index.
  - **Method:** it prefilters candidates by matched-fragment count (≥2 by default) and then computes full XCorr and E-values. The index is persisted to `.idx` files.
  - **Speed-ups vs classic Comet** (i9-11900, 12 threads, 64 GB):

    | Search type | Speed-up |
    |---|---|
    | TMT fractionated | 1.7–2.3× |
    | HLA-I | 6.5× |
    | Non-enzymatic | 18× |
    | Open search | ~50× |
    | Astral narrow-window (45.3 vs 155.6 min) | 3.4× |
    | Real-time search (tryptic / phospho / HLA-I) | 4.9× / 6× / 14.7× |

  - The open search is **still ~3.8× slower than MSFragger**.
  - **Lesson:** XCorr plus an index is viable, and the index mainly helps when the candidate count is large.

### 1.4 Tide / Crux (Noble lab)
- **Tide v4.3** (Kertesz-Farkas et al., bioRxiv 2025, doi:10.1101/2025.04.01.646675; Apache-2.0, [crux-toolkit](https://github.com/crux-toolkit/crux-toolkit)) **does not use an inverted fragment index**. It streams a precursor-sorted peptide index against spectra sorted by an **external merge sort**, keeping only N spectra in memory for N files. Removing redundant fragment calculations made it 2–7× faster than v4.2.
- **Benchmarks** (Xeon E5-2640v4, 4 threads):

  | Dataset | Tide v4.3 | MSFragger | Sage |
  |---|---|---|---|
  | PXD011515 (12.5 M spectra; 31 M proteins, 7.1 B peptides) | 179 h, 20 GB, 351 k PSMs | 448 h in 32 batches, 124 GB, 210 k PSMs | OOM |
  | 539-file PEAKS set (25.2 M spectra, 5.8 B peptides) | 16.3 h, **2 GB** | OOM | OOM |
  | PXD028806 (221 GB spectra, normal size) | 673 s, 9 GB | 793 s, 64 GB | **431 s**, 15 GB |

  *Derived* for the 539-file set: ≈430 spectra/s on 4 threads.
- **Takeaway:** sorting spectra by precursor mass and sweeping the peptide index once is the scalable pattern for huge search spaces. It maps directly onto **cohort-level amortisation**: sort spectra from thousands of runs once, then sweep once.
- **GPU-Tide** (BMC Bioinf 2018, doi:10.1186/s12859-018-2559-6) is 2.7× (Tide) and 5.8× (Tide-Hi) faster at 50 Da tolerance. GiCOPS beat it by >10× (§3).

### 1.5 MaxQuant/Andromeda, X!Tandem, MS-GF+, Mascot, others
- **Andromeda** (Cox et al., *JPR* 2011, doi:10.1021/pr101065j) uses a binomial probability score over the top-q peaks per 100 Th. Closed source, free.
  - MaxQuant 2.6.x runs on .NET 8 on Linux. It is slow on large cohorts: AlphaPept is ~3× faster on Orbitrap, and DDA+/FragPipe is 4× faster.
  - Andromeda is the weakest base score in the Casanovo-DB comparison (§5).
- **X!Tandem** (Craig & Beavis, *Bioinformatics* 2004, doi:10.1093/bioinformatics/bth092) introduced the hyperscore and an E-value from a survival-function fit. Uses a multi-pass "refinement" for unanticipated mods. Artistic licence. Its open search is two orders of magnitude slower than an index (976 min in the MSFragger benchmark).
- **MS-GF+** (Kim & Pevzner, *Nat Commun* 2014, doi:10.1038/ncomms6277) computes spectral E-values from a generating function. These are well calibrated and robust across fragmentation types, but slow and memory-hungry. Licence: UC Regents, **educational/research/non-profit use only** ([LICENSE](https://github.com/MSGFPlus/msgfplus/blob/master/LICENSE.txt)). quantms still uses it as one of its engines.
- **Mascot** (Perkins et al., *Electrophoresis* 1999) uses a probability-based MOWSE score (−10·log10 P). Commercial and closed; it is a de facto reference in many PRIDE submissions, so reprocessing must reproduce Mascot-era results.
- **AlphaPept** (Strauss et al., *Nat Commun* 2024, doi:10.1038/s41467-024-46485-4; Apache-2.0) is written in Python+Numba with HDF5 storage. A 1.6 GB RAW becomes 200 MB of HDF5.
  - **Speed:** 2.07 min/file on a 12-core workstation and 1.4 min/file on a cluster; 200 HeLa runs took 280–482 min.
  - **PSM step** (260 k spectra vs 5 M entries):

    | Implementation | Time |
    |---|---|
    | Pure Python | 10 h |
    | Numba | 80 s |
    | Multithreaded Numba | **9 s** |
    | GPU (CuPy) | **77 s** |

    The GPU is slower because the step is memory-bound. Feature finding got 150× from the GPU.
- **Open-pFind** (Chi et al., *Nat Biotechnol* 2018, doi:10.1038/nbt.4236) uses a **sequence-tag index** of ~5-aa tags to generate candidates for open search.
  - It confidently identifies 70–85% of spectra.
  - It is ">40× faster than three other open-search engines and 2–3× faster than restricted searches".
  - It reanalysed ~25 M human spectra, giving 14,064 proteins.
  - Closed source (free).
- **ionbot** (Degroeve et al., bioRxiv 2021, doi:10.1101/2021.07.02.450686) is an ML-integrated open-modification search with MS²PIP/DeepLC features. It adds up to 12% more peptides through co-eluting lower-rank matches. It ran a **PRIDE-wide reanalysis of ~924 M human spectra (25,000 RAW files) and ~600 M mouse spectra (14,500 files)**. Closed and cloud-based.
- **Dear-PSM** (He et al., *Smart Med* 2024, doi:10.1002/SMMD.20240014) inverts the index: it **indexes the experimental spectra**, then streams peptides. It claims 7× speed and 240× less memory than mainstream engines and allows 3–20 variable mods per peptide. This is a useful idea for cohort-level search where spectra outnumber peptides.
- **MetaMorpheus** (G-PTM-D) is open-source (MIT) and does mod discovery by augmenting the database. It is slower than DDA+ (above) **[licence unverified]**.

---

## 2. Open/mass-tolerant search and PTM localization

### 2.1 Open-search strategies
- **Precursor-window open search (±150–500 Da)** with a fragment index (MSFragger, Sage, Comet-FI). Only unmodified fragments match, so the modified-side ions are lost. MSFragger adds "localization-aware" rescoring with the delta mass placed on each residue. Open search raises PSMs by 34% on HEK293 (§1.1), but FDR must be stratified by mass-shift group. Open search inflates candidate counts and lets mass-shift-specific false positives through.
- **Mass-offset / multinotch:** a curated list of offsets. It is much tighter than an open window and is the right default for repository reprocessing with a known PTM catalogue (Unimod plus empirical PTM-Shepherd output).
- **Tag-based** (Open-pFind, TagGraph, ionbot-like) generates candidates from de novo tags. It is fast for open search and naturally handles mass shifts at any position.
- **Spectral-library open search** (ANN-SoLo, HOMS-TC; §4) uses shifted dot products against libraries, which can be predicted (ANN-SoLo 2025 with Prosit).

### 2.2 Post-search characterisation and localization tools

| Tool | Approach | Reference | Licence |
|---|---|---|---|
| **PTM-Shepherd** | Clusters open-search mass shifts; rescores each delta at each position → localization enrichment per residue; RT shift; spectral similarity; multi-experiment comparison | Geiszler et al., *MCP* 2021, PMC7950090 | Apache-2.0 |
| **PTMProphet** | Bayesian mixture/EM over site-determining ions; outputs site probabilities; part of TPP | Shteynberg et al., *JPR* 2019, doi:10.1021/acs.jproteome.9b00205 | TPP (LGPL) **[unverified]** |
| **phosphoRS / ptmRS** | Cumulative binomial probability of matched site-determining ions, tolerance-aware; integrated in Proteome Discoverer | Taus et al., *JPR* 2011, doi:10.1021/pr200611n | Closed/free binary |
| **Ascore** | Peak-depth-optimised probability of site-determining ions | Beausoleil et al., *Nat Biotechnol* 2006, doi:10.1038/nbt1240 | Open implementations (PNNL) |
| **LuciPHOr2** / **pyLucXor (2026)** | Target-decoy for localization (decoy residues) → explicit FLR | Fermin et al., *Bioinformatics* 2015; pyLucXor *BIO Web Conf* 2026 | Open |
| **DeepFLR** | Predicted phosphopeptide spectra + target-decoy FLR | *Nat Commun* 2023, doi:10.1038/s41467-023-38035-1 | Open |
| **DeepRescore2** | phosphoRS + predicted spectrum similarity + RT features | *MCP* 2024, PMC10831110 | Open |

- **FLR validation:** Ramsbottom et al. (*JPR* 2022, PMC9251759) propose **decoy amino acids**, i.e. scoring phospho on residues that cannot carry it. **Alanine and leucine** are the most reliable decoys. They show that synthetic-calibrated localization scores are often mis-calibrated on real data. pyLucXor on the synthetic PXD000138 set localizes 51,468 sites at 1% FLR with 98.55% true accuracy. The synthetic phosphopeptide library (>100,000 peptides with known sites; Kuster lab, TUM) is the gold standard. Use both, since this feeds our benchmarking work package.
- **Key gap:** in every major engine, localization is a **separate post-processing pass** that re-enumerates isoforms per PSM. With a fragment index this can be fused: compute site-determining-ion evidence for all positional isomers during the final scoring of the top-K candidates. The work is cheap and regular, a dense small-matrix job well suited to the GPU.

---

## 3. GPU and other hardware-accelerated search engines

| System | Year | Hardware | Reported speed-up | Notes |
|---|---|---|---|---|
| **GPU-Tide** (BMC Bioinf, doi:10.1186/s12859-018-2559-6) | 2018 | NVIDIA GPU | 2.7× vs Tide, 5.8× vs Tide-Hi (50 Da) | XCorr on GPU; superseded |
| **IP2-GPU ProLuCID / PaSER** (Bruker, Yates/Park) | 2020 | GPU server attached to timsTOF | Keeps up with 120 Hz PASEF acquisition (real-time) | Commercial; shows GPU closed search can meet real-time rates |
| **HiCOPS** (Haseeb & Saeed, arXiv 2102.02286) | 2021 | 72 nodes / 1,728 cores | "500×" vs existing tools on TB-scale databases; minutes vs weeks | Distributed-memory partitioning of the index across nodes |
| **GiCOPS** (Haseeb & Saeed, *Sci Rep* 2023, doi:10.1038/s41598-023-43033-w) | 2023 | RTX A6000 48 GB + 2× Xeon Gold 5215 | 1.2–2× (closed), 3–5× (open, ±100–500 Da) vs HiCOPS; >10× (up to 50×) vs GPU-Tide | Hyperscore; transfer limits closed search (kernel only 50% of time vs 85% in open); CC BY-NC-ND |
| **AlphaPept GPU** | 2024 | Numba/CuPy | **GPU slower than multithreaded CPU** for PSM scoring (77 s vs 9 s) | Memory-bound; GPU good for feature finding (150×) |
| **ANN-SoLo GPU** (JPR 2019, doi:10.1021/acs.jproteome.9b00291) | 2019 | Faiss GPU | ~10× vs CPU ANN-SoLo | Feature-hashed low-dim vectors + ANN for library OMS |
| **HyperOMS** (arXiv 2211.16422) | 2022 | GPU, HDC | "up to 100×", 1,337× energy efficiency vs OMS GPU tool | Binary hypervectors |
| **HOMS-TC** (*Bioinformatics* 2023, PMC10323168) | 2023 | RTX 4090 tensor cores (INT8 / XOR-popcount) | **31× (geo-mean) vs ANN-SoLo-GPU**; 51× on search stage; iPRG2012 (1.16 M lib spectra), HEK293 (2.99 M) | IDs similar to SpectraST/MSFragger, slightly fewer than ANN-SoLo; Apache-2.0 |
| **HyperSpec** (*JPR* 2023, PMC10243109) | 2023 | GPU HDC clustering | 21 M spectra clustered in 24 min (vs 4 h) | Repository-scale dedup |
| **FiCOPS** (bioRxiv 2026, doi:10.64898/2026.02.15.706012) | 2026 | Intel Stratix-10 FPGA | 3.5× vs CPU tools; 3–5× lower power than CPU/GPU; open search 4× vs MSFragger | On-instrument target |
| **Cross-domain OMS** (Pinge, Song et al., IEEE JETCAS 2026, arXiv 2607.17569) | 2026 | Near-memory/in-storage | >100× speed, >40,000× energy | Finds "performance is dominantly limited by reference data movement rather than computation" |
| **cMGVB** (*Proteomes* 2026, doi:10.3390/proteomes14030034) | 2025/26 | GPU multi-node | Not extracted | GPU PTM-combination search **[details unverified]** |

**What the GPU literature shows:**
- No open-source GPU fragment-index engine **beats Sage/MSFragger on closed tryptic search while matching sensitivity on public benchmarks**. GiCOPS's comparisons are against HiCOPS and GPU-Tide, not Sage. This gap is our opportunity.
- The bottleneck is **data movement**: host↔device transfer of spectra and index slices, plus random gathers in the index. Kernels that convert scoring into dense, regular work win, e.g. binary hypervectors with popcount, INT8 tensor-core dot products, or batched dense matrix products over binned spectra.
- Open search is the regime where the GPU helps most, because there are many candidates per spectrum.

---

## 4. Spectral library search

- **SpectraST** and **MSPepSearch** are the classical CPU reference tools for dot-product library search.
- **ANN-SoLo** (Bittremieux et al., *JPR* 2018, doi:10.1021/acs.jproteome.8b00359; GPU 2019; Apache-2.0) uses a cascade: a narrow search first, then open search. It feature-hashes high-res spectra into low-dimensional vectors and uses ANN (Faiss) candidate selection, then a shifted dot product.
  - **2025 update** (bioRxiv doi:10.1101/2025.05.20.655174): **unified library and database search via Prosit-predicted libraries**, neutral-loss-aware vectors, and spectrum- and peptide-level decoys.
- **HOMS-TC** and HyperOMS: see §3.
- **Flash entropy search** (Li & Fiehn, *Nat Methods* 2023, doi:10.1038/s41592-023-02012-9; metabolomics) queries **1 billion spectra in <2 s** with an inverted index on fragment peaks. It is essentially the fragment-index idea applied to library spectra. It is directly portable to repository-scale "find all spectra of this peptide" queries.
- **MassIVE-KB** (Wang et al., *Cell Syst* 2018) is a library built from >1 B spectra with 5.95 M peptides. **Spectroscape** (PMC10560257) does real-time query of a spectral archive. **FASST** (MassIVE) does repository-scale spectrum search (ProteomeXchange 2026 NAR update).
- **Implication:** with Prosit/AlphaPeptDeep-class predictors, "library search" and "database search" converge. Predicted intensities can serve as the reranker feature, or as a dense vector for GPU scoring, on top of fragment-index candidates.

---

## 5. Learned scoring, de novo sequencing as candidate generator, and embedding retrieval

### 5.1 Learned PSM scores
- **Casanovo-DB** (Bioinformatics 2024, PMC11211853) uses the Casanovo transformer's likelihood as the PSM score. Peptides detected at 1% FDR on E. coli, yeast and human, compared with each base score:

  | Compared with | Extra peptides from Casanovo-DB |
  |---|---|
  | XCorr | **+31–40%** |
  | X!Tandem hyperscore | +52–69% |
  | Andromeda | +88–102% |

  Percolator helps Casanovo-DB more than the other scores. Casanovo v5 `db_search` reports +61–82% vs XCorr.
- **CHIMERYS** (Frejno et al., *Nat Methods* 2025, PMC12074992) explains each spectrum as a **non-negative LASSO mixture of predicted spectra** (Prosit intensity + RT).
  - **>2/3 of DDA MS2 spectra contain multiple precursors.**
  - Peptide groups +10% (plasma) to +35% (FFPE) vs Sequest HT.
  - Precursor gains vs MSFragger: phospho +36%, ubiquitin +45%.
  - Entrapment validation is consistent for 1.4–20.4 Th windows.
  - Commercial (Thermo Proteome Discoverer 3.1, cloud/Kubernetes).
- **Tesorai Search** (bioRxiv 2024, doi:10.1101/2024.08.19.606805) is a pretrained spectrum–peptide model that scores "without using decoys".
  - It reports +12% / +9% / +21% peptides vs FragPipe / PEAKS / Proteome Discoverer.
  - It ran 250 immunopeptidomics searches in 45 min.
  - Free for academics via a webserver; the company is commercial.
  - The lack of decoy-based validation is a **red flag for FDR**.
- **pUniFind** (arXiv 2507.00087, 2025) is pretrained on >100 M open-search spectra and covers >1,300 modifications. It uses a cross-modal joint scorer plus zero-shot de novo. It reports +42.6% immunopeptides; speed and code are not reported.
- **DeepSearch** (arXiv 2405.06511) is a contrastive encoder–decoder for end-to-end search with zero-shot variable PTMs; no speed numbers.

### 5.2 De novo sequencing (as a candidate generator / open-search front end)

| Model | Key facts |
|---|---|
| **Casanovo** (*Nat Commun* 2024, doi:10.1038/s41467-024-49731-x; Apache-2.0) | Transformer trained on 30 M MassIVE-KB PSMs. **v5** (bioRxiv 2025, doi:10.1101/2025.07.25.666826): 3.5× faster inference, **~235 spectra/s on a GPU**, vectorised beam search 3.1–3.5× faster, HDF5→**Lance** data backend, peptide score = product of amino-acid scores (better calibrated) |
| **InstaNovo / InstaNovo+** (*Nat Mach Intell* 2025; Apache-2.0) | Transformer + multinomial-diffusion refinement; 28 M training spectra; InstaNovo-P for phospho |
| **π-HelixNovo** (*Brief Bioinform* 2024, doi:10.1093/bib/bbae021) | Encodes the **complementary spectrum** to recover missing ions; aa precision/recall 0.765/0.758 on the nine-species set |
| **π-PrimeNovo** (*Nat Commun* 2025, doi:10.1038/s41467-024-55021-3) | **Non-autoregressive** with CUDA precise-mass-control decoding; **up to 89× faster inference** than autoregressive models; the fastest option if de novo is in the hot path |
| **Others** | ContraNovo (contrastive), PLMNovo (2025, PLM-aligned embeddings, CC-BY-NC-ND), Winnow (de novo FDR calibration, arXiv 2509.24952), NovoBench (benchmark), MassNet (billion-scale training corpus, bioRxiv 2025) |

- **Use for us:** at ~10²–10³ spectra/s per GPU, de novo is **1–2 orders of magnitude slower than fragment-index search per spectrum**. It is worth running only on spectra left unexplained after the fast search, e.g. 50–70% of a typical run. There it can (a) generate tags or candidates for open/error-tolerant search and (b) flag novel peptides.

### 5.3 Embedding / bi-encoder retrieval

| Method | Key facts | Licence/code |
|---|---|---|
| **GLEAMS** (Bittremieux et al., *Nat Methods* 2022, doi:10.1038/s41592-022-01496-1) | **32-d** spectrum embedding trained on 30 M MassIVE-KB PSMs; embedded **669 M spectra** from 227 human datasets on 4× RTX 2080 Ti; clustered 511 M spectra into 60 M clusters; ID propagation added **132 M PSMs (+71%)** | BSD |
| **yHydra** (bioRxiv 2021, doi:10.1101/2021.12.01.470818) | **Joint** spectrum/peptide Euclidean embedding trained on ~20 M PSMs; closed, open and error-tolerant search modes | — |
| **SpecEncoder** (Liu, Tao, Ye, Tang, *Bioinformatics* 2024, PMC11211836) | Metric learning on spectra + predicted spectra; +1–2% peptides vs SpectraST; **+6–15% vs MS-GF+/Percolator** | [lkytal/SpecEncoder](https://github.com/lkytal/SpecEncoder) |
| **Foundation spectrum encoder** (Sanders … Noble, arXiv 2505.10848) | De novo pretraining → transferable representations for quality, chimericity, phospho and glyco prediction | — |

- **Framework paper** (Kalogeropoulos, Van Goey, Jenkins, Eloff, *JPR* 2026, 25:2234, doi:10.1021/acs.jproteome.5c01153, PMC13142651).
  - **Assumptions:** S = 10⁶ spectra, P = 10⁸ peptides, A100-80GB plus a 32–64-core CPU.

    | Approach | Complexity | Build | Query | Memory |
    |---|---|---|---|---|
    | Classical | O(S·ρP) | — | ~4–6 h | 10–20 GB |
    | Fragment index | O(P+S) | 2–4 h | ~1–2 h | **100–300 GB** |
    | Biencoder + ANN | O(S(log P + K)) | 4–6 h | **2–3 h** | 40–80 GB |
    | Exhaustive cross-encoder | — | — | **~2–12 days** | — |
    | De novo + rescoring | — | — | 10–15 h | — |

  - **Assumed throughputs:** cross-encoder ~10⁴ pairs/s; GPU dot products ~2×10⁷/s. That GPU figure is very conservative; FAISS on an A100 does >10⁹ dot products/s, so the bi-encoder query time is likely overestimated by ≥10×.
  - **Recommendation:** "classical candidate generation + neural rescoring of top-K."
- **Assessment:**
  - No bi-encoder has yet shown **entrapment-validated** parity with a fragment index on a full human proteome with variable mods.
  - The known failure modes are poor recall for rare PTMs and isobaric/positional isomers, which are near-identical in embedding space. Bi-encoders also do not encode the exact precursor-mass constraint, so it must be applied as a filter.
  - Embeddings are strong for: (i) **cross-run clustering and deduplication**, where the same peptide recurs across thousands of runs and could be searched once per cluster; (ii) open/unknown-mass-shift retrieval, where the fragment index degrades; (iii) cheap **pre-filters** for which spectra deserve expensive treatment.

---

## 6. File formats and I/O

| Format / reader | Key facts | Licence |
|---|---|---|
| **mzML** (PSI) | XML with base64(+zlib/numpress) arrays; slow to parse, large; random access only via the indexedmzML offsets | Open standard |
| **mz5** (Wilhelm et al., *MCP* 2012, doi:10.1074/mcp.O111.011379) | HDF5 mapping of mzML; smaller and faster than mzML; low adoption | Open |
| **mzMLb** (Bhamber et al., *JPR* 2021, doi:10.1021/acs.jproteome.0c00192) | HDF5 with chunked compressed arrays and the XML metadata kept verbatim; "faster than all existing approaches in virtually all cases"; compressed size comparable to vendor files; random access to compressed data orders of magnitude faster than compressed mzML; readers in ProteoWizard and pyteomics | Open |
| **MGF / MS2** | Text peak lists; MS2-only, centroided; trivial to parse but lose MS1 and ion mobility | — |
| **Aird** (PMC8756627; StackZDPD PMC8967824) | Computation-oriented; ZDPD/SIMD codecs: ~80% smaller than mzML, decoding 3× faster than zlib | Open (CSi-Studio) |
| **mzPeak** (Kohlbacher, Van Den Bossche, Klein et al., *JPR* 2025, 24:5329, doi:10.1021/acs.jproteome.5c00435; spec at [psidev.info/mzpeak](https://www.psidev.info/mzpeak); reference impl. [mobiusklein/mzpeak_prototyping](https://github.com/mobiusklein/mzpeak_prototyping)) | Uncompressed ZIP of **Parquet** tables (spectrum metadata, signal, chromatograms). "Point" layout = parallel arrays plus a repeated spectrum-index column, filterable through Parquet page indexes. "Chunked" layout = fixed m/z chunks with δ-encoding or Numpress. Zstd (default level 3), byte-shuffle, float32 ion mobility. Rust reader/writer (Thermo RAW via .NET, Bruker TDF, mzML, MGF), Python (pyarrow + DataFusion SQL), R, TS/WASM, C#. "Competes with and sometimes beats vendor formats for size, and uncompressed mzML for complex random access." Spec still "in development"; no published hard benchmarks yet | Open |
| **ThermoRawFileParser** (Hulstaert et al., *JPR* 2020, doi:10.1021/acs.jproteome.9b00328) | .NET wrapper around Thermo's **RawFileReader** (proprietary but freely redistributable DLL); outputs mzML/MGF/**Parquet**; v2.0.0-dev (Oct 2025) on .NET 8 with Parquet schema v0.2. The PRIDE/quantms standard converter. MSFragger 4.1 also reads RAW directly, faster | Apache-2.0 |
| **rawrr** (Kockmann & Panse, *JPR* 2021, doi:10.1021/acs.jproteome.0c00866) | R/Bioconductor wrapper of RawFileReader (.NET 8); has a spectra/s benchmark function | GPL-3 **[unverified]** |
| **Bruker .d / TDF** | SQLite metadata + `analysis.tdf_bin` (zstd-compressed frames of scan/TOF indices). Vendor SDK (timsdata DLL) is closed | — |
| **timsrust** (MannLabs; MIT) | Pure-Rust TDF/TSF/miniTDF reader with centroiding and pseudo-MS2 primitives; used by Sage (ASMS 2026 poster). Fork: **sage-plus-tdf** (MIT/Apache-2.0). Related: alphatims (Apache-2.0), OpenTIMS, Rustims/imspy (MIT) | MIT |
| **Lance** | Columnar ML format with fast random row access; adopted by Casanovo v5 instead of HDF5 | Apache-2.0 |
| **Sage output** | Parquet results since v0.14 | MIT |

**I/O takeaways:**
- **Decode once, centroid once, store columnar.** For DDA identification only centroided MS2 peaks are needed, plus precursor m/z, charge, isolation window, RT, ion mobility and the MS1 needed for quantification.
- A float32 m/z (or u32 fixed-point) plus f16/u16 intensity, top-N filtered, is roughly 10–50× smaller than profile mzML **[estimate]**. Such a store can be memory-mapped directly into the search's spectrum arrays with zero copy (Arrow/Parquet → Rust slices / CUDA pinned buffers).
- Thermo RAW decoding requires the .NET RawFileReader. This is a CPU-bound, per-file job that parallelises across files and should run as a separate, cached conversion tier.

---

## 7. Repository-scale reanalysis: reported compute numbers

| Effort | Scale | Compute / throughput |
|---|---|---|
| **PRIDE Archive** (2025 update, doi:10.1093/nar/gkae1011) | 42,036 datasets, **2.85 PB** (Aug 2024; 1.35 PB in 2021); ~534 datasets/month; >10% of datasets have >100 RAW files; largest has 7,444 RAW files | — |
| **ionbot PRIDE reanalysis** | 924 M human spectra (25 k RAW) + 600 M mouse spectra (14.5 k RAW); open modification search | Not reported |
| **MassIVE reanalysis** (ProteomeXchange 2026 NAR) | **2.2 B spectra from 392 datasets → 1.1 B PSMs**; MassIVE-KB from >1 B spectra | MODa open search used for unexpected mods |
| **GLEAMS** | 669 M spectra embedded; 511 M clustered; +132 M IDs | 4× RTX 2080 Ti workstation |
| **quantms** (Dai et al., *Nat Methods* 2024, doi:10.1038/s41592-024-02343-1) | 83 datasets, **29,354 RAW files**, 13,132 samples, 1.03 M peptides | **~2 min per RAW file average, ~10 h per dataset on an HPC cluster**; engines Comet/MS-GF+/Sage + ConsensusID; open-source Nextflow |
| **Open-pFind** | ~25 M spectra (human draft proteome) → 14,064 proteins | — |
| **Tide v4.3** | 25.2 M spectra × 5.8 B peptides | 16.3 h on 4 threads, 2 GB (*derived* ≈430 spectra/s) |
| **HyperSpec** | 21 M spectra clustered in 24 min on GPU | — |

**Back-of-envelope target** (my estimate):
- Public human DDA in PRIDE is probably 10⁹–10¹⁰ MS2 spectra.
- At a Sage-like ~10⁴ spectra/s per 32-core node for a closed tryptic search (*derived* order of magnitude from the TMT benchmark), 10¹⁰ spectra take ~10⁶ node-seconds (~12 node-days). That is trivial on a cluster.
- **The real costs are:** (i) reading and decoding petabytes of vendor files (2.85 PB at 1 GB/s per node is about 33 node-days just to read); (ii) open/PTM search spaces, which are 10–100× more expensive; (iii) rescoring with neural features; (iv) global FDR across studies.

---

## 8. Master tool table

| Tool | Category | Speed headline (hardware) | Licence | Code |
|---|---|---|---|---|
| Sage | Fragment-index DB search (Rust) | 12 TMT files in 18 s (32 vCPU), 5.6× MSFragger; open search 5× MSFragger | **MIT** | github.com/lazear/sage |
| MSFragger (+DDA+, Glyco, Labile, DIA) | Fragment-index DB search (Java) | Open 500 Da: 5.4 min/41.8 k spectra (4 cores), 150× Comet | Academic only; commercial via Fragmatics | Closed; github.com/Nesvilab/MSFragger (releases) |
| Comet / Comet-FI | XCorr ± fragment index (C++) | FI: 1.7–50× over classic; open search 3.8× slower than MSFragger | Apache-2.0 (per authors) | github.com/UWPR/Comet |
| Tide (Crux) v4.3 | Streaming XCorr/p-value (C++) | 25 M spectra × 5.8 B peptides in 16 h, 2 GB (4 threads) | Apache-2.0 | github.com/crux-toolkit/crux-toolkit |
| MaxQuant/Andromeda | Probability score (.NET) | ~3–4× slower than AlphaPept/FragPipe | Free, closed | maxquant.org |
| X!Tandem | Hyperscore/E-value (C++) | Open 500 Da: 976 min/41.8 k spectra | Artistic | thegpm.org |
| MS-GF+ | Generating-function E-value (Java) | Slow; robust E-values | UC academic/non-profit | github.com/MSGFPlus/msgfplus |
| Mascot | MOWSE probability | — | Commercial | matrixscience.com |
| AlphaPept | Numba fragment-count search | ~1.4–2 min/file | Apache-2.0 | github.com/MannLabs/alphapept |
| Open-pFind | Tag-index open search | >40× other open engines | Free, closed | pfind.ict.ac.cn |
| ionbot | ML open search | 924 M human spectra reprocessed | Closed, cloud | ionbot.cloud |
| MetaMorpheus | G-PTM-D, chimeras | Slower than DDA+ | MIT **[unverified]** | github.com/smith-chem-wisc/MetaMorpheus |
| CHIMERYS | LASSO chimeric deconvolution | DDA+ ≥20× faster | Commercial | — |
| GiCOPS / HiCOPS | GPU / MPI DB search | 3–5× open over HiCOPS (A6000) | CC BY-NC-ND (GiCOPS) | github.com/pcdslab/gicops |
| GPU-Tide | GPU XCorr | 2.7× Tide | — | — |
| ANN-SoLo | ANN library OMS (GPU) | ~10× over CPU version | Apache-2.0 | github.com/bittremieux-lab/ANN-SoLo |
| HOMS-TC | HDC tensor-core library OMS | 31× ANN-SoLo-GPU (RTX 4090) | Apache-2.0 | github.com/tycheyoung/homs-tc |
| HyperSpec | HDC GPU clustering | 21 M spectra in 24 min | Open | — |
| GLEAMS | 32-d embedding + clustering | 669 M spectra embedded | BSD | github.com/bittremieux-lab/GLEAMS |
| yHydra / SpecEncoder | Joint/metric embedding search | SpecEncoder +6–15% peptides vs MS-GF+ | Open | GitHub |
| Casanovo (v5) / Casanovo-DB | Transformer de novo / learned score | 235 spectra/s per GPU; DB score +31–40% vs XCorr | Apache-2.0 | github.com/Noble-Lab/casanovo |
| InstaNovo(+) | Transformer + diffusion de novo | — | Apache-2.0 (code) | github.com/instadeepai/InstaNovo |
| π-PrimeNovo | Non-autoregressive de novo | Up to 89× faster than autoregressive | Open | GitHub |
| π-HelixNovo | Complementary-spectrum de novo | — | Open | GitHub/OpenI |
| Tesorai / pUniFind / DeepSearch | Pretrained PSM scorers | Tesorai: 250 HLA searches in 45 min | Commercial / unreleased | — |
| PTM-Shepherd | Open-search mass-shift characterisation | — | Apache-2.0 | github.com/Nesvilab/PTM-Shepherd |
| PTMProphet / ptmRS / Ascore / LuciPHOr2 / DeepFLR | Site localization / FLR | — | Mixed (see §2.2) | — |
| ThermoRawFileParser | RAW→mzML/MGF/Parquet | — | Apache-2.0 (vendor DLL proprietary) | github.com/compomics/ThermoRawFileParser |
| timsrust | Bruker TDF reader (Rust) | — | MIT | github.com/MannLabs/timsrust |
| mzPeak | Parquet MS format | "≈ vendor size" | Open | github.com/mobiusklein/mzpeak_prototyping |

---

## 9. Weaknesses and open problems we can exploit

1. **Index memory explosion and the 2³¹ peptide cap.** Sage and MSFragger OOM on large or expanded spaces, and MSFragger had to split PXD011515 into 32 batches. Nobody ships a **compressed, partitioned, out-of-core or GPU-sharded fragment index** that keeps Sage-level speed. Options include delta-coded u32 fragment keys, 16-bit precursor-mass offsets within buckets, and peptide IDs local to each slice.
2. **No cross-run amortisation.** Every engine searches run by run, or at best builds the index once per invocation. Comet-FI `.idx` and MSFragger `.pepindex` caches exist. Two cohort-level gains are unused:
   - **Global sorting of spectra by precursor mass** across all runs (Tide-style), followed by a single sweep of each index slice.
   - **Spectrum clustering and deduplication** (GLEAMS/HyperSpec/falcon), searching cluster representatives first. With 10⁹ spectra and heavy redundancy (GLEAMS: 511 M → 60 M clusters), this is potentially a ~5–10× reduction in scoring work.
3. **GPU closed search is memory-bound and unproven against Sage.** The gather-heavy inverted index suits CPUs. A GPU win needs one of: (a) batching thousands of spectra per index slice, so slices stay in shared memory/L2 while spectra stream through; (b) a dense formulation over a binned spectrum matrix (spectra × m/z bins) and a peptide fragment matrix, i.e. SpMM; (c) using the GPU only for the expensive parts (open search, the reranker, localization).
4. **Weak base scores.** Hyperscore, XCorr and Andromeda leave 30–100% of peptides undetected before rescoring (Casanovo-DB). The industry fix is post hoc rescoring with predicted intensities and RT (MSBooster, MS2Rescore, Oktoberfest, inSPIRE). This is a second pass that **re-reads spectra and re-predicts**. Fusing a predicted-intensity reranker into the search loop for the top-K saves I/O and lets the better score drive candidate selection.
5. **Chimeric spectra are the default, but most engines are single-ID.** Over two-thirds of DDA MS2 spectra are chimeric, and full-isolation-window search gives +15–57% peptides. Sage's chimera mode is a cheap approximation. Open-source DDA+-quality chimeric search with entrapment-validated FDR does not exist, because MSFragger is not open.
6. **FDR in open/learned/chimeric settings is fragile.**
   - Tesorai claims a decoy-free search.
   - Learned scores can "learn the decoys".
   - Chimeric secondary IDs need their own FDR.
   - Mass-shift-stratified FDR in open search is ad hoc.
   - Entrapment validation is sporadic (DDA+ and CHIMERYS do it).

   An engine with **built-in entrapment reporting** and per-stratum (charge, mass-shift, rank, chimera-order) FDR would be a differentiator.
7. **Localization is a post-processing afterthought.** No engine computes calibrated site probabilities in-search with a decoy-residue FLR (Ala/Leu decoys). For open search, localization is only "enrichment" (PTM-Shepherd).
8. **I/O is unaddressed.** mzML parsing and RAW conversion dominate small-search wall time (quantms ~2 min/run end to end). There is no adopted columnar standard yet: mzPeak is in development, and TRFP's Parquet is a v0.2 schema. We can define a lean internal Arrow schema aligned with mzPeak so it stays forward-compatible.
9. **Neural candidate retrieval is not established.** The bi-encoder papers (yHydra, SpecEncoder, pUniFind, DeepSearch) report modest or unvalidated gains, and none benchmark against Sage/MSFragger speed. This is a publishable gap if done carefully, but it is risky as a primary path.

---

## 10. Implications for our design (recommendations)

**R1. Build the core on a Sage-like fragment index in Rust, extended rather than reinvented.**
- Fork or borrow from Sage (MIT): the digest, index, mods and LDA code, plus timsrust and sagepy.
- Keep Sage's two-level layout (fragment-m/z buckets sorted internally by precursor mass), because it maps directly to GPU tiles.
- Add **partitioning by precursor-mass slices** with a fixed memory budget per slice, compressed encoding and an on-disk persistent index keyed by a hash of (FASTA, enzyme, mods, version), built **once per cohort**.
- Target: process PXD011515-scale databases (7 B peptides) in <64 GB, where Sage OOMs today.

**R2. Use cohort-level "sort-once, sweep-once" execution.**
- Ingest all spectra from N runs into a columnar store.
- Globally sort, or bucket, MS2 by precursor neutral mass (per charge-state hypothesis).
- Stream index slices through CPU caches or GPU memory with all spectra in that mass range.
- This turns per-run random access into sequential sweeps (the Tide v4.3 insight combined with the Sage index). It also enables the **dedup/clustering pre-pass**: search cluster medoids first, then propagate and verify.

**R3. CPU/GPU split.**
- **CPU (AVX-512/NEON):**
  - closed tryptic candidate generation (memory-latency-bound gather; CPUs win per AlphaPept and GiCOPS);
  - RAW/TDF decoding;
  - I/O.
- **GPU:**
  - (a) open/mass-offset/semi/non-specific searches, where candidates per spectrum are ≥10³ and dense scoring pays off (GiCOPS 3–5×, HOMS-TC 31×);
  - (b) the top-K neural reranker (cross-encoder or predicted-spectrum similarity);
  - (c) in-search positional-isomer localization;
  - (d) spectrum embeddings for clustering.
- First prototype the GPU kernel as a binned-spectrum × fragment-list sparse-dense product and benchmark it against Sage on identical candidates, **before** committing.

**R4. Use retrieve-then-rerank, not bi-encoder-first.**
- Stage 1: fragment-index top-K (K ≈ 10–50) with a cheap hyperscore-like score.
- Stage 2: a neural reranker on GPU (Casanovo-DB-style likelihood or a split-then-joint encoder) plus predicted intensity and RT features. Cache the spectrum-encoder output so it can be reused for clustering and de novo.
- Use a bi-encoder/ANN only as a **side channel**: (i) for spectra with no good stage-1 hit (open/unknown mass shift), and (ii) for cross-run clustering.
- Bi-encoder-first is justified only if it shows ≥95% recall@K of stage-1 top-1 hits on entrapment-validated benchmarks. That is a clear go/no-go experiment.

**R5. Search chimeric spectra by default.** Search the full isolation window, as DDA+ does, for every MS2 in the stage-1 index; it is nearly free with a fragment index. Then do a stage-2 joint explanation: greedy subtraction or NNLS over predicted spectra in the CHIMERYS style, restricted to K candidates. Report rank-stratified FDR. The expected gain is +15–50% peptides.

**R6. Make localization and PTM support first-class.**
- Mass-offset search should be the default for repository reprocessing, using a curated offset list from Unimod plus a PTM-Shepherd-style empirical scan of a pilot subset.
- Run full open search on demand.
- Compute site probabilities in-search for top-K isoforms (an Ascore/ptmRS-style site-determining-ion score plus predicted-intensity evidence).
- **Estimate FLR with decoy residues (Ala/Leu)** and validate on the synthetic phosphopeptide library (PXD000138 and the Kuster synthetic set).

**R7. FDR plumbing.**
- Emit features for mokapot/Percolator, but aim at a calibrated direct score.
- Build **entrapment-database support into the engine** (shuffled or foreign-proteome entrapment) so every benchmark reports estimated FDP alongside q-values. DDA+ and CHIMERYS set this precedent.
- Stratify FDR by mass-shift group, chimera rank and search stage.

**R8. I/O tier.**
- Write a one-time converter to an Arrow/Parquet store aligned with mzPeak's "point layout": centroided MS2, f32 m/z or u32 fixed-point, u16/f16 intensity, top-N per spectrum, zstd, and spectrum-level metadata columns.
- Thermo uses RawFileReader via ThermoRawFileParser or mzPeak's Rust+.NET path; Bruker uses timsrust.
- Memory-map the store straight into the search (zero-copy).
- Keep mzML/MGF readers for compatibility only.
- Measure conversion throughput as a first-class benchmark: it is probably the PRIDE-scale bottleneck (2.85 PB).

**R9. Baselines and benchmarks to run.**
- **Speed and IDs:**
  - Sage (current HEAD), MSFragger 4.4 (academic licence; DDA, DDA+ and open modes) and Comet-FI on PXD016766 (TMT, the Sage benchmark), PXD001468 (open search) and HEK293 (MSFragger paper).
  - Astral narrow-window data (the Comet-FI benchmark).
  - Tide on PXD011515 for huge search spaces.
- **Hardware:** fixed hardware with a documented core count; report spectra/s, peak RSS and IDs at 1% FDR, plus entrapment FDP.
- **GPU comparators:** GiCOPS (non-commercial licence, benchmark only) and HOMS-TC/ANN-SoLo for library OMS.

**R10. Do not depend on closed or non-commercial components.**
- MSFragger, GiCOPS, MS-GF+, CHIMERYS, ionbot and pFind are benchmark-only.
- Reusable stack: Sage (MIT), timsrust (MIT), ThermoRawFileParser (Apache-2.0; RawFileReader is redistributable under Thermo's licence), Casanovo/InstaNovo (Apache-2.0), PTM-Shepherd (Apache-2.0) and ANN-SoLo/HOMS-TC (Apache-2.0).

---

## References (primary sources)

- Kong AT et al. MSFragger. *Nat Methods* 2017. doi:10.1038/nmeth.4256 (PMC5409104)
- Lazear MR. Sage. *J Proteome Res* 2023;22:3652. doi:10.1021/acs.jproteome.3c00486; design blog https://lazear.github.io/sage/
- Yu F, Deng Y, Nesvizhskii AI. MSFragger-DDA+. *Nat Commun* 2025. doi:10.1038/s41467-025-58728-z (PMC11978857)
- Polasky DA et al. MSFragger-Glyco. *Nat Methods* 2020. doi:10.1038/s41592-020-0967-9; MSFragger-Labile. *MCP* 2023 (PMC10182319)
- MSFragger-DIA. *Nat Commun* 2023 (PMC10338508); MSFragger changelog https://msfragger.nesvilab.org/CHANGELOG.html
- McGann CD … Schweppe DK. Comet fragment-ion indexing. bioRxiv 2024. doi:10.1101/2024.10.11.617953
- Kertesz-Farkas A et al. Fast and memory-efficient searching with Tide. bioRxiv 2025. doi:10.1101/2025.04.01.646675
- Haseeb M, Saeed F. GiCOPS. *Sci Rep* 2023. doi:10.1038/s41598-023-43033-w; HiCOPS arXiv:2102.02286
- GPU-Tide/HiXCorr. *BMC Bioinformatics* 2018. doi:10.1186/s12859-018-2559-6
- Kumar S … Saeed F. FiCOPS. bioRxiv 2026. doi:10.64898/2026.02.15.706012
- Strauss MT et al. AlphaPept. *Nat Commun* 2024. doi:10.1038/s41467-024-46485-4 (PMC10924963)
- Chi H et al. Open-pFind. *Nat Biotechnol* 2018. doi:10.1038/nbt.4236
- Degroeve S et al. ionbot. bioRxiv 2021. doi:10.1101/2021.07.02.450686
- He Q et al. Dear-PSM. *Smart Med* 2024. doi:10.1002/SMMD.20240014
- Frejno M et al. CHIMERYS. *Nat Methods* 2025 (PMC12074992)
- Geiszler DJ et al. PTM-Shepherd. *MCP* 2021 (PMC7950090)
- Shteynberg DD et al. PTMProphet. *JPR* 2019. doi:10.1021/acs.jproteome.9b00205
- Taus T et al. phosphoRS. *JPR* 2011. doi:10.1021/pr200611n
- Beausoleil SA et al. Ascore. *Nat Biotechnol* 2006. doi:10.1038/nbt1240
- Ramsbottom KA et al. Decoy-amino-acid FLR. *JPR* 2022 (PMC9251759)
- DeepFLR. *Nat Commun* 2023. doi:10.1038/s41467-023-38035-1; DeepRescore2 (PMC10831110)
- Bittremieux W et al. ANN-SoLo. *JPR* 2018 doi:10.1021/acs.jproteome.8b00359; GPU *JPR* 2019 doi:10.1021/acs.jproteome.9b00291; 2025 update bioRxiv doi:10.1101/2025.05.20.655174
- HOMS-TC. *Bioinformatics* 2023 (PMC10323168); HyperOMS arXiv:2211.16422; HyperSpec *JPR* 2023 (PMC10243109); Cross-domain OMS arXiv:2607.17569
- Li Y, Fiehn O. Flash entropy search. *Nat Methods* 2023. doi:10.1038/s41592-023-02012-9
- Bittremieux W et al. GLEAMS. *Nat Methods* 2022. doi:10.1038/s41592-022-01496-1 (PMC9189069)
- yHydra. bioRxiv doi:10.1101/2021.12.01.470818; SpecEncoder *Bioinformatics* 2024 (PMC11211836)
- Kalogeropoulos K, Van Goey J, Jenkins TP, Eloff K. A framework for database search with AI models. *JPR* 2026;25:2234. doi:10.1021/acs.jproteome.5c01153 (PMC13142651)
- Yilmaz M et al. Casanovo. *Nat Commun* 2024. doi:10.1038/s41467-024-49731-x; Casanovo-DB *Bioinformatics* 2024 (PMC11211853); Casanovo v5 bioRxiv doi:10.1101/2025.07.25.666826
- Eloff K et al. InstaNovo. *Nat Mach Intell* 2025; π-HelixNovo *Brief Bioinform* 2024 doi:10.1093/bib/bbae021; π-PrimeNovo *Nat Commun* 2025 doi:10.1038/s41467-024-55021-3
- Tesorai Search. bioRxiv doi:10.1101/2024.08.19.606805; pUniFind arXiv:2507.00087; DeepSearch arXiv:2405.06511; spectrum foundation model arXiv:2505.10848
- mzPeak. *JPR* 2025;24:5329. doi:10.1021/acs.jproteome.5c00435; https://www.psidev.info/mzpeak
- mzMLb. *JPR* 2021. doi:10.1021/acs.jproteome.0c00192; mz5 *MCP* 2012 doi:10.1074/mcp.O111.011379; Aird (PMC8756627)
- ThermoRawFileParser. *JPR* 2020. doi:10.1021/acs.jproteome.9b00328; rawrr *JPR* 2021 doi:10.1021/acs.jproteome.0c00866
- Teschner D et al. Rustims/sagepy. *JPR* 2025 (PMC12053931); timsrust https://github.com/MannLabs/timsrust
- Dai C et al. quantms. *Nat Methods* 2024. doi:10.1038/s41592-024-02343-1
- Perez-Riverol Y et al. PRIDE 2025 update. *NAR* 2025. doi:10.1093/nar/gkae1011; ProteomeXchange 2026 update, *NAR* 54:D459
- Eng JK et al. Comet. *Proteomics* 2013. doi:10.1002/pmic.201200439; Craig R, Beavis RC. X!Tandem. *Bioinformatics* 2004. doi:10.1093/bioinformatics/bth092; Kim S, Pevzner PA. MS-GF+. *Nat Commun* 2014. doi:10.1038/ncomms6277; Cox J et al. Andromeda. *JPR* 2011. doi:10.1021/pr101065j; Perkins DN et al. Mascot. *Electrophoresis* 1999
