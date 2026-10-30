# 04 — Specialized assays (XL-MS, LiP-MS, CF-MS, TPP, glyco) and benchmark / gold-standard datasets

Compiled 2026-10-26. Every PRIDE accession marked **verified** in the tables was queried live via
`https://www.ebi.ac.uk/pride/ws/archive/v3/projects/<ACC>` and `/files` (paginated, all pages). File counts,
sizes, instruments and licences below come from those API responses. ProteoBench archives were checked by
HTTP HEAD and directory listing on `proteobench.cubimed.rub.de`. Non-PRIDE ProteomeXchange entries
(Panorama, jPOST) were checked on ProteomeCentral only, so their sizes are not known. The verification script
is in the session scratchpad (`pride.py`/`lsf.py`). Re-run it before downloading because file lists change.

Node context: 4 CPU cores, about 500 GB RAM, 1× RTX 3090 (24 GB). `dotnet`, `mono`, `docker`, `cargo` and
container runtimes are **not** on PATH.

---

## Part A — Specialized assays

### A1. Cross-linking MS (XL-MS)

**The core computational problem.** For a non-cleavable linker (DSS/BS3, PhoX, SDA), the precursor mass is
m(α)+m(β)+m(linker). An exhaustive search has to consider every peptide pair whose summed mass matches the
precursor. That is O(n²) in the number of linkable peptides n, before counting modifications and link sites.
At proteome scale (around 10⁶ linkable tryptic peptides) this means about 10¹² pairs. Inter-protein (heteromeric)
pairs far outnumber intra-protein (self) pairs. In the E. coli lysate study of Lenz et al. 2021, a
4,350-protein database gave about 9.5 M protein pairs, almost all heteromeric. Because most of the search space is
heteromeric, random matches fall mainly on heteromeric links. Lenz et al. found that a decoy match was **10.6× more
likely among heteromeric than among self links**. Pooled FDR therefore gives inter-protein links (PPIs) an
FDR far above nominal.

**Known solutions to the quadratic search space** (each is a candidate for our CPU/GPU design):

1. **MS-cleavable linkers** (DSSO, DSBU, DSBSO; Kao 2011 MCP; Müller 2010 Anal Chem). In MS2 the linker
   fragments and leaves signature doublets (e.g., DSSO alkene/thiol +54/+86 Da stubs, Δ≈31.97 Da; DSBU
   +85/+111 Da). The doublets give each peptide's mass, so the search becomes two **linear** searches plus a
   pairing step. This approach is used by XlinkX (Liu 2015 Nat Methods 12:1179; Klykov 2018 Nat Protoc), MeroX
   (Götze 2015 JASMS; Iacobucci 2018 Nat Protoc 13:2864), MS Annika (Pirklbauer 2021 JPR 20:2560) and Scout
   (Clasen et al. 2024 Nat Methods 21:2327, doi:10.1038/s41592-024-02478-1). Acquisition can be stepped-HCD MS2,
   MS2-MS3 or EThcD.
2. **"Open-modification" α-first search.** Score each linear peptide with an unknown mass offset (the partner
   plus linker) using linear fragments only. Keep the top-k α candidates, then look up β by mass. xiSEARCH
   (Mendes 2019 Mol Syst Biol 15:e8994), Kojak (Hoopmann 2015 JPR 14:2190) and pLink 2 work this way. pLink 2
   (Chen 2019 Nat Commun 10:3404, doi:10.1038/s41467-019-11337-z) adds a fragment index and reports 40× the
   speed of pLink 1 and 3–10× the speed of Kojak. pLink 3 was released in January 2025
   (github.com/pFindStudio/pLink3). This is a mass-offset search that is structurally the same as
   MSFragger-style open search, so **our open-search fragment index can be reused directly**.
3. **Sparse-matrix formulation.** MS Annika 3.0 (Birklbauer et al., Commun Chem 2024,
   doi:10.1038/s42004-024-01386-x; bioRxiv 10.1101/2024.09.03.610962) casts non-cleavable candidate
   generation and scoring as **sparse matrix multiplication** (spectra × fragment-ion matrix). This allows
   proteome-wide DSS-type searches on commodity hardware (shown on C. elegans nuclei). It maps naturally onto
   GPU SpMM (cuSPARSE) and is our most direct design precedent.
4. **Database reduction.** Run a linear pre-search, then restrict the XL search to detected proteins, or to the
   top-N proteins or SEC co-fractions. Lenz 2021 crosslinked within 44 SEC fractions. This is effective, but the
   decoy database must be reduced the same way or the FDR becomes biased (Fischer & Rappsilber 2017).
5. **ML rescoring.** Scout uses an ANN, and Percolator-style rescoring is common. Beveridge 2020 showed that
   validation choice alone changes the real FDR: Kojak with Percolator vs PeptideProphet gave **2.4% vs 32%**
   actual false crosslinks at 5% nominal FDR.

**Tools to know.** xiSEARCH + xiFDR (Rappsilber lab, open source), pLink 2/3 (pFind), MeroX/StavroX, XlinkX
(Proteome Discoverer node), Kojak (+Percolator), MaxLynx (MaxQuant; Yılmaz 2022 Anal Chem,
doi:10.1021/acs.analchem.1c03688; handles timsTOF), Scout (ANN, cleavable, very fast; Liu lab FMP Berlin),
OpenPepXL (OpenMS; Netz 2020 MCP), MS Annika 1/2/3 (PD node + standalone; Mechtler/Dorfer).

**Crosslinker chemistry relevant to search:**
- BS3/DSS (non-cleavable, K/S/T/Y/N-term, 138.068 Da).
- DSSO (cleavable, 158.004 Da, doublets).
- DSBU (cleavable, 196.085 Da).
- PhoX (non-cleavable, phosphonic-acid handle for IMAC enrichment; Steigenberger 2019 ACS Cent Sci). Search it
  like a DSS-type linker and include its enrichment-dependent mono-link and loop-link forms.
- DSBSO (cleavable + enrichable).
- CDI (zero-length-ish, cleavable).
- ADH/DHSO (acidic-residue linkers).
- SDA (photo-reactive, heterobifunctional K–any). The second residue is unconstrained, so the search space is
  much larger.

Our engine needs linker definitions as data, covering mass, specificity per end, cleavage stubs and
mono-/loop-link masses.

**FDR in XL-MS.**
- **Target–decoy pairs.** Decoy rates come from TT/TD/DD pairs: FDR ≈ (TD − DD)/TT (Walzthoeni 2012 Nat
  Methods 9:901).
- **Separate FDR at each level.** FDR must be computed separately at CSM, peptide-pair, residue-pair (ResPair)
  and PPI level. False CSMs rarely cluster, so a 1% CSM-FDR can map to a much higher PPI-FDR (Fischer &
  Rappsilber 2017 Anal Chem, doi:10.1021/acs.analchem.6b03745, "Quirks of error estimation…").
- **Self vs heteromeric.** These should be controlled separately, giving a heteromeric PPI-FDR (Lenz 2021).
  xiFDR implements this and adds "boosting" of pre-filters to maximize results at a target level.

**Ground-truth benchmarks**, all from Mechtler's IMP lab apart from Scout's and Lenz's:
- **Beveridge et al. 2020** (Nat Commun 11:742, doi:10.1038/s41467-020-14608-2; **PXD014337**). Design:
  95 synthetic Cas9 (S. pyogenes) peptides, each with a single K. They were split into 12 groups, crosslinked
  within groups, then pooled, giving about 426 possible true crosslinks. Any link between groups is a known
  false positive. Linkers were DSS, DSBU and DSSO.

  Actual FDR at 5% estimated FDR:

  | Engine | Actual FDR |
  |---|---|
  | pLink | 11.3% |
  | StavroX | 5.2% |
  | xi | 2.5% |
  | Kojak | 2.4–32% |
  | MeroX | 4.9% |
  | XlinkX | 0% |

- **Matzinger et al. 2022** (Nat Commun 13:3975, doi:10.1038/s41467-022-31701-w; **PXD029252**). Design: 141
  peptides from 38 E. coli ribosomal proteins, in groups of 6–10, with up to 1,018 possible links. Library
  variants are main, enrichable and acidic. Linkers were DSSO, DSBU, CDI, DSBSO, ADH and DHSO. Some runs were
  spiked 1:5 into a tryptic HEK digest, and some were SEC-fractionated, to add background. All six engines
  (MS Annika, MaxLynx, MeroX, XlinkX, pLink 2, xiSEARCH) had **real FDR above the 1% estimate**. Error was
  higher for inter-protein links. The study also provides the **IMP-X-FDR** tool for computing real FDR from
  the library design.
- **Clasen/Lima et al. 2024 (Scout)** (Nat Methods; **PXD042173**). Design: 256 recombinant human
  proteins/fragments in 32 interaction groups of 8, crosslinked with DSSO. Links between groups are false.
  This gives a ground truth at CSM, ResPair and PPI level at near-proteome scale, about 10× the size of the
  peptide libraries. Scout was >200× faster than xiSEARCH on 4 raw files with a 540-protein database.
- **Lenz et al. 2021** (Nat Commun 12:3564, doi:10.1038/s41467-021-23666-z; **PXD019120** BS3+DSSO,
  **PXD024148** sSDA). Design: E. coli lysate, 44 SEC fractions crosslinked separately, then SCX × hSAX, giving
  90 runs per linker on a QE HF. Result: 590 PPIs at 1% heteromeric PPI-FDR. This is not a strict ground
  truth, but it is the reference realistic proteome-wide dataset for FDR at PPI level.

### A2. Limited proteolysis MS (LiP-MS)

**Principle.** Feng et al. 2014 (Nat Biotechnol 32:1036, doi:10.1038/nbt.2999) and Schopper et al. 2017 (Nat
Protoc 12:2391, doi:10.1038/nprot.2017.100). A native lysate is briefly digested with proteinase K (PK,
broad-specificity), typically 1:100 for 5 min at 25 °C. The sample is then denatured and fully digested with
trypsin. Structural changes alter PK accessibility. This shows up as changes in **half-tryptic peptides** (one
end from PK, one from trypsin) and in fully-tryptic peptides spanning PK sites. A parallel **trypsin-only (TrP)
control** measures protein abundance.

**Applications and analysis tools:**
- Piazza 2018 Cell 172:358 (E. coli metabolite–protein interactome).
- **LiP-Quant** (Piazza et al. 2020 Nat Commun 11:4200, doi:10.1038/s41467-020-18071-x; **PXD015446** +
  PXD018204/PXD019902): drug-target deconvolution from dose–response LiP plus an ML score.
- **LiP-MS in situ** (Cappelletti 2021 Cell 184:545).
- **MSstatsLiP**: Malinovska, Cappelletti et al. 2023 Nat Protoc (doi:10.1038/s41596-022-00771-x;
  **PXD031627**, **PXD031616**), with the MSstatsLiP Bioconductor package by Kohler/Vitek. It fits peptide-level
  models and **adjusts LiP peptide fold changes for TrP protein-level changes**, so that abundance changes are
  not called structural changes.
- **FLiPPR** (Manriquez-Sandoval et al. 2024 JPR 23:2332, doi:10.1021/acs.jproteome.3c00887): FragPipe post-processor. It merges
  ions into cut-sites, uses missingness-aware imputation and applies protein-centric multiple-testing
  correction.
- **PELSA** (Li et al. 2025 Nat Methods, doi:10.1038/s41592-024-02553-7; PXD034606): a high-trypsin-ratio
  "peptide-centric local stability assay" variant.

**What breaks in a generic engine:**

1. **Semi-specific search is mandatory.** It enlarges the candidate space roughly 10–20× over fully tryptic.
   Fully non-specific search would be needed for PK–PK products, but those are rare after trypsin. Our fragment
   index must support semi-specific digestion efficiently, e.g., by indexing prefix/suffix ladders once.
2. **FDR.** Half-tryptic and fully-tryptic peptides have different prior probabilities. Use group-specific
   (stratified) FDR, or at least report the two classes separately. Pooled FDR lets low-prior half-tryptic
   matches ride on the tryptic score distribution.
3. **Quantification at peptide/precursor level, not protein level.** Peptides must be mapped to cut-site
   coordinates, then normalized against the TrP protein abundance from the matched control sample. The
   LiP/TrP pairing must therefore be explicit in experiment metadata (SDRF), and MBR must not transfer between
   LiP and TrP runs. They are different digests, and half-tryptic peptides do not exist in TrP.
4. **DIA libraries.** Libraries for DIA must include semi-tryptic precursors, which inflates library size and
   lowers DIA identification power. A two-stage approach is common: a TrP-style library plus empirically
   observed half-tryptic peptides.

### A3. Co-fractionation MS (CF-MS) and SEC-MS

**Principle.** Native complexes are separated by SEC, IEX or BN-PAGE into 40–100+ fractions, each measured by
LC-MS. Proteins with correlated elution profiles are inferred to interact.

**Landmark studies:**
- Havugimana 2012 Cell (IEX).
- Kristensen 2012 Nat Methods (PCP-SILAC).
- Wan 2015 Nature (metazoan).
- **Heusel et al. 2019** (Mol Syst Biol 15:e8438, doi:10.15252/msb.20188438). SEC-SWATH-MS with **CCprofiler**,
  using complex-centric detection against prior complexes and **decoy complexes for FDR**. Data: PXD007038
  (HEK293) and PXD010288 (HeLa interphase/mitosis).
- **SECAT** (Rosenberger et al. 2020 Cell Syst): differential complex analysis.
- **EPIC** (Hu et al. 2019 Nat Methods 16:737): Emili/Hu lab, ML on elution profiles.
- **PrInCE** (Stacey 2017 BMC Bioinf; R/Bioconductor version by Skinnider).
- **Skinnider & Foster 2021** (Nat Methods 18:806, doi:10.1038/s41592-021-01194-4). Reanalyzed **206 CF-MS
  experiments** uniformly with MaxQuant 1.6.5, giving 11 M protein-abundance measurements. They benchmarked
  designs and network inference, and recommended co-elution scores and classifier setups. The processed
  resource became **CFdb** (Nat Commun 2023, doi:10.1038/s41467-023-44139-5; PXD042664).
- **Skinnider et al. 2021 Cell 184:4073** (PCP-SILAM, 7 mouse tissues, >125k interactions; PXD012009).

**Why naive pipelines break.**
1. **Match-between-runs across fractions.** MBR transfers an identification into an adjacent fraction because
   the m/z and RT match, even if the peptide was absent there. That widens and smooths elution peaks and
   creates false co-elution, which is exactly the signal used for PPI inference. Adjacent SEC fractions have
   highly similar composition, so MBR "succeeds" in matching chromatographic features that are really another
   peptide or noise.

   *This is a mechanistic argument; we did not find a paper that quantifies it.* The general MBR false-transfer
   problem is documented in the IonQuant MBR-FDR paper (Yu et al. 2021 MCP). There, 44% of yeast proteins had
   transfers into yeast-free samples before filtering.

   What the engine needs:
   - (a) per-fraction provenance flags (MS2-identified vs transferred);
   - (b) **fraction-aware alignment**: only transfer between *replicates of the same fraction index*, or within
     a ±k fraction window when the elution-profile shape supports it, with transfer FDR (IonQuant-style decoy
     transfers) estimated separately for CF-MS;
   - (c) the option to export unimputed matrices. PrInCE/EPIC handle missingness themselves.
2. **Normalization.** Per-fraction total-intensity normalization destroys the elution profile, because
   fractions differ in total protein amount by design. Normalize per replicate series, or not at all, and use
   iRT/spike-ins for RT alignment only.
3. **Protein inference / proteoforms.** Peptides from one gene can elute differently (proteoforms or
   subcomplexes). CCprofiler/SECAT do peptide-level profile clustering. The engine should output
   peptide × fraction matrices, not only protein rollups.
4. **FDR at the PPI/complex level.** This needs decoy complexes (CCprofiler) or reference-set precision
   (PrInCE/EPIC). The engine's role is to supply well-calibrated peptide-level q-values per fraction, plus a
   global q-value across the whole series. Running series-level FDR, rather than per-run FDR, avoids
   accumulating false IDs across 50–100 runs, the same "master FDR" problem as repository-scale meta-analysis.

### A4. Thermal proteome profiling (brief)

TPP was introduced by Savitski 2014 (Science 346:1255784), with the protocol in Franken 2015 (Nat Protoc 10:1567;
data **PXD013774**, verified). Variants are TPP-TR, TPP-CCR, 2D-TPP and PISA. Bacterial TPP: Mateus 2018 Mol
Syst Biol (PXD009495). Analysis uses NPARC (Childs 2019 MCP).

Engine requirements:
- TMT10/16/18-plex reporter-ion quantification with isotope-impurity correction and co-isolation (ratio
  compression) handling: MS3/SPS, or high-resolution MS2 with a purity filter.
- Per-temperature protein summarization without global normalization that would remove the melting signal.
  Normalization should use a stable-protein subset (TPP-TR convention).
- Peptide-level output, so that domain-level stability can be analysed.

The search itself is ordinary DDA, so the main requirement is correct TMT quantification. Sage already supports
TMT, and PXD016766 is its TMT benchmark.

### A5. Glycoproteomics (brief)

Intact glycopeptides are a mass-offset search problem with glycan-specific fragments:
- Oxonium ions trigger detection.
- Y-ions (peptide + partial glycan) dominate HCD.
- c/z ions localize sites under EThcD.

Key tools:
- **MSFragger-Glyco** (Polasky 2020 Nat Methods 17:1125): offset search, with glycan composition assignment
  and glycan FDR in later FragPipe versions.
- **pGlyco3** (Zeng 2021 Nat Methods 18:1515).
- **O-Pair** (Lu 2020 Nat Methods 17:1133).

The HUPO community study (Kawahara et al. 2021 Nat Methods 18:1304, doi:10.1038/s41592-021-01309-x) used
2 serum files (N- and O-glycopeptides, HCD/EThcD/CID). Its consensus is in GlyConnect reference 2943, and the
best-performing tools were Protein Prospector and Byonic. Lessons:
- Joint HCD+EThcD search is best.
- Narrow glycan databases are more accurate and broad ones give more coverage.
- NeuAc/NeuGc and multi-fucose are frequent misassignments.
- Glycan FDR and peptide FDR must be controlled separately.

Glyco needs the same open/offset index as open search, so it comes nearly free once that index exists.
Benchmark data accessions for glyco were **not verified** in this pass. Glyco is a deferred, nice-to-have item.

---

## Part B — Benchmark and gold-standard datasets

### B0. Statistical background for entrapment

**Wen, Freestone, Riffle, MacCoss, Noble & Keich 2025** (Nat Methods, doi:10.1038/s41592-025-02719-x; bioRxiv
10.1101/2024.06.01.596967; code: github.com/Noble-Lab/FDRBench; data: Zenodo 10.5281/zenodo.15073580).
Definitions: r = |entrapment DB| / |target DB|, N_T = target hits, N_E = entrapment hits.

Estimators:
- **Lower bound:** N_E/(N_T+N_E). This is a lower bound only and cannot validate FDR control.
- **Combined (valid upper bound):** N_E(1+1/r)/(N_T+N_E).
- **"Sample"** method, which only counts entrapment hits rescaled against target hits: shown to be flawed and
  able to err in either direction.
- **Paired:** new, tighter upper bound. It needs r=1 and one-to-one paired entrapment peptides, e.g., shuffled
  from each target peptide with the C-terminal residue kept.

Entrapment databases used:
- Shuffled peptides.
- Foreign proteomes: *Arabidopsis*, yeast, primates (*Macaca*, *Callithrix*, *Papio*), mouse. Shared and
  homologous peptides must be removed.

Findings:
- DDA tools (Tide, Sage, MS-GF+, MSFragger) appear to control peptide-level FDR.
- **DIA-NN, Spectronaut and EncyclopeDIA do not consistently control FDR**. Protein-level results were
  "consistently invalid". Single-cell data was worst, with up to 48% inflation at precursor level.

Datasets used include ISB18, HEK293 **PXD001468**, Astral/QE/TripleTOF/timsTOF/Lumos DIA (**PXD042704**,
Panorama), PXD034525, PXD041421, single-cell PXD023325 (jPOST JPST000971), **PXD028735**, **PXD017703**,
**PXD012988** and MSV000084000.

**ProteoBench** (community benchmark platform; bioRxiv 10.64898/2025.12.09.692895;
proteobench.readthedocs.io) adopts the Wen estimators in a **DIA entrapment module (Astral)**. It provides a
pre-digested peptide-level FASTA with paired entrapments, and only Met-ox is allowed. This makes it a turnkey
external check for our DIA FDR.

Design rule for our engine: entrapment is a *test harness* feature. We should be able to append shuffled or
foreign entrapment sequences with a tag, carry the tag through every level (PSM → precursor → peptide →
protein → PTM site → crosslink PPI), and report lower-bound, combined and paired FDP automatically on every
benchmark run.

### B1. Dataset table

Abbreviations: Sz = total project size (all files); Raw = raw-format summary. "✔" = verified by PRIDE API or
ProteoBench server today. "PX" = verified on ProteomeCentral only, sizes unknown. "–" = not verified.
"EBI-ToU" = "EBI terms of use", the default for older submissions with no explicit CC licence. Reuse for
benchmarking is fine, but cite the source.

#### Entrapment / FDR calibration

| Accession | Study | Instrument | Sz / files (raw format) | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD001077 | Vaudel 2012 JPR, *P. furiosus* "standard designed by evolution" | LTQ Orbitrap Velos | 0.5 GB; **1 raw 0.33 GB** + mgf + mzid | EBI-ToU | Tiny DDA entrapment: Pfu run searched vs Pfu+human (or larger); non-Pfu hits are false | ✔ |
| PXD001468 | Chick 2015 Nat Biotechnol, HEK293 (MSFragger/Sage open-search file; Wen 2025 DDA entrapment) | Q Exactive | 116.9 GB; 24 .raw (3.1–3.9 GB) + 24 mzXML | EBI-ToU | DDA speed vs Sage/MSFragger, open search, entrapment (FDRBench) | ✔ |
| ProteoBench DIA-entrapment-Astral | HeLa 50 ng, 15-min DIA, 3 reps (Wen-style paired entrapment FASTA) | Orbitrap Astral | 3 × 3.1 GB .raw; tar.gz 5.65 GB | ProteoBench (open) | DIA FDR validity, direct leaderboard comparison | ✔ (server) |
| PXD042704 | Heil 2023 Astral DIA eval (Panorama Public) | Astral, Lumos | n/a | PX | Used in Wen 2025 for DIA entrapment | PX |
| PXD023325 / JPST000971 | single-cell DIA chip | LTQ Orbitrap (metadata) | n/a | PX | Low-input FDR stress test (Wen 2025: worst inflation) | PX |
| PXD041421 | A549/K562 dda/diaPASEF with batch effects | timsTOF Pro | 532 GB; 96 .d.zip (2.2–9.2 GB) | CC0 | Bruker DIA entrapment, batch effects | ✔ |
| PXD012988 | Mouse myoblast DIA reproducibility | QE HF | 112.7 GB; 15 .raw (2.8–3.2 GB) | CC0 | DIA entrapment (Wen) | ✔ |

#### PTM localization ground truth

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD000138 | Marx 2013 Nat Biotechnol, >100k synthetic peptides + phospho counterparts | LTQ Orbitrap Velos (HCD/CID/ETD) | 69.4 GB; 193 .raw (0.13–0.6 GB) + mgf/xml/csv | EBI-ToU | Gold-standard phosphosite FLR; used by Locard-Paulet 2020 (22 pipelines) | ✔ |
| PXD007058 | Ferries 2017 JPR (doi:10.1021/acs.jproteome.7b00337), synthetic phosphopeptide pools + U2OS TiO2 | Orbitrap Fusion | 32.6 GB; 96 .raw (**pools 0.17–0.22 GB**; U2OS ~1.5 GB) + PD .msf | CC0 | 8 fragmentation schemes (HCD/EThcD/ETciD, OT vs IT): FLR per method | ✔ |
| PXD009449 | Zolg 2018 MCP, ProteomeTools 21 PTMs (phospho, acetyl, methyl, citrulline, ubiquitinyl GG, …) | Fusion Lumos | 150.1 GB; 100 .raw (0.32–0.75 GB) + zips | EBI-ToU | Multi-PTM localization and open-search mass-offset truth | ✔ |
| PXD014525 | Bekker-Jensen 2020 Nat Commun DIA phospho (synthetic phosphopeptides spiked into yeast/HeLa) | Q Exactive HF(-X) | 699 GB; 426 .raw (0.28–2.7 GB) | CC0 | **DIA** phosphosite localization with ground truth | ✔ |
| PXD004732 | ProteomeTools part I (Zolg 2017 Nat Methods) | Orbitrap Fusion | 1.09 TB; 1,460 .raw (0.37–1.07 GB) | EBI-ToU | Unmodified synthetic truth, spectrum-prediction training | ✔ |

#### Multi-species quantitative (HYE / PYE)

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD028735 | Van Puyvelde 2022 Sci Data (doi:10.1038/s41597-022-01216-6) HYE, log2FC H 0 / Y −1 / E +2 | TTOF5600/6600, Synapt, QE HF-X, timsTOF Pro | **3.43 TB**; 845 files: .wiff.scan, .raw(.zip), .d.zip | CC0 | Multi-vendor DDA/DIA quant (ProteoBench DDA-QE uses the Alpha series, 6 × 3.5–3.9 GB) | ✔ |
| ProteoBench DDA-QExactive | Subset of PXD028735 (A/B × 3) | QE HF-X | tar.gz 19.4 GB | CC0 | DDA LFQ ratio accuracy and leaderboard | ✔ |
| ProteoBench DDA-Astral | HYE A/B × 3, 15-min | Orbitrap Astral | 6 × 2.6 GB; tar.gz 15.9 GB | ProteoBench | DDA LFQ on Astral | ✔ |
| ProteoBench DIA-Astral | HYE A/B × 3, 2-Th windows | Orbitrap Astral | 6 × 3.7–3.9 GB; tar.gz 14.1 GB | ProteoBench | DIA LFQ | ✔ |
| ProteoBench DIA-diaPASEF | HYE A/B × 3 (also PXD062685: timsTOF SCP, 44 GB .7z) | timsTOF | tar.gz 43.6 GB | CC0 | Bruker DIA LFQ | ✔ |
| ProteoBench DIA-plasma (JPST003358) | Distler et al. PYE: plasma 1:1, yeast 1:3, E. coli 2:1, 6 reps | timsTOF | **12 .d.zip, 0.96–1.29 GB each**; tar.gz 14.1 GB | ProteoBench | Small Bruker DIA quant with high dynamic range | ✔ (server) |
| PXD070049 | Van Puyvelde 2026 "Generation Beta" HYE (preprint doi:10.64898/2026.01.29.702266), short gradients | Astral, timsTOF Ultra/Ultra 2, ZenoTOF 7600+/8600 | **3.11 TB**; 2,174 files (.raw 433 GB, .d.zip 933 GB, wiff) | CC0 | Modern-instrument, high-throughput quant (tier L) | ✔ |
| PXD002952 | Navarro 2016 Nat Biotechnol LFQbench HYE | TripleTOF 5600/6600 (SWATH) | 470 GB; .wiff/.wiff.scan | EBI-ToU | Historical DIA quant (SCIEX format, low priority) | ✔ |
| PXD017703 | Meier 2020 Nat Methods diaPASEF (incl. TwoProteome HeLa/yeast) | timsTOF Pro | 826 GB in 16 big zips (11–93 GB) | CC0 | Bruker DIA (poor subsetting: monolithic zips) | ✔ |
| PXD014777 | Prianichnikov 2020 MCP MaxQuant-timsTOF (HYE ddaPASEF, plasma, HeLa fractions) | timsTOF Pro | 204 GB; 226 .d.zip (0.06–12 GB) | CC0 | Bruker DDA ID/quant | ✔ |
| PXD003881 | Shen 2018 (IonStar) human proteins spiked in E. coli at 5 levels, N=20 | Orbitrap Fusion | 40.4 GB; 20 .raw (~1.5 GB) + mzXML | EBI-ToU | DDA LFQ spike-in accuracy | ✔ |
| PXD001819 | Ramus 2016, UPS1 in yeast (9 levels × 3) | LTQ Orbitrap Velos | 60.7 GB; 27 .raw (~1.7 GB) | EBI-ToU | Classic low-abundance spike-in | ✔ |
| PXD046283 | Guzman/Olsen 2024 Nat Biotechnol nDIA (Astral) | Orbitrap Astral | 346 GB; 27 .raw (1.3–19 GB) | CC0 | Astral DIA throughput/speed | ✔ |
| PXD049412 | Astral single-cell (up to 5,300 proteins/cell) | Astral, Exploris 480 | 527 GB; 357 .raw | CC0 | Low-input DIA | ✔ |

#### Standard DDA speed / ID benchmarks

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD001468 | (see above) | QE | 24 raw ~3.3 GB | EBI-ToU | Sage blog/MSFragger open-search reference file (first HEK293 run) | ✔ |
| PXD016766 | Eclipse TMT benchmark (Anal Chem 2020, doi:10.1021/acs.analchem.9b05685) | Orbitrap Eclipse | 16.7 GB; 48 .raw (0.2–0.37 GB) | CC0 | **Sage's TMT benchmark** (12 mzML); isobaric quant | ✔ |
| PXD010012 | Meier 2018 MCP PASEF | timsTOF Pro | 630 GB in big zips (≥19.9 GB) | CC0 | ddaPASEF (poor subsetting) | ✔ |

#### XL-MS ground truth

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD014337 | Beveridge 2020, Cas9 95-peptide library (DSS ×3 reps, DSBU, DSSO stHCD/MS3/EThcD/CID-ETD) | QE HF-X, Fusion Lumos | **9.9 GB; 9 .raw (0.85–1.48 GB)** | CC0 | Known-truth CSM/ResPair FDR, cleavable and non-cleavable | ✔ |
| PXD029252 | Matzinger 2022, ribosomal 141-peptide library, 6 linkers, HEK-spiked and SEC variants | QE HF-X, Exploris 480 (FAIMS), Eclipse | 112.7 GB; 84 .raw (0.44–2.1 GB) + xi mzid | CC0 | Known-truth FDR with background; FAIMS; inter vs intra | ✔ |
| PXD042173 | Clasen/Lima 2024 Scout recombinant standards (256 proteins, 32 groups, DSSO) | Fusion Lumos | 233.8 GB; 178 .raw (0.3–2.6 GB) | CC0 | PPI-level ground truth at proteome-like scale | ✔ |
| PXD019120 | Lenz 2021, E. coli lysate SEC → BS3/DSSO → SCX × hSAX | QE HF | 330.9 GB; 314 .raw (0.69–1.18 GB) + mgf | CC0 | Realistic proteome-wide heteromeric PPI-FDR | ✔ |
| PXD024148 | Lenz 2021 sSDA (photo-linker) on affinity-enriched E. coli | Fusion Lumos | 37.8 GB; 50 .raw | CC0 | Heterobifunctional search-space stress | ✔ |
| PXD027161 | MaxLynx 2022, BSA DSSO/DSBU | timsTOF Pro | 3.4 GB; 2 .d.7z | CC0 | Bruker XL-MS reader test | ✔ |

#### LiP-MS

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD031627 | Malinovska/Cappelletti 2023 Nat Protoc, "Correction" set (yeast LiP + TrP), with Spectronaut LiP/TrP reports and MSstatsLiP outputs | Orbitrap (SDRF: Fusion Lumos; filenames "QE1") | 110 GB; 32 DIA .raw (~3 GB) + 8 smaller .raw (~0.8 GB) + reports (2.2 GB) | CC0 | LiP/TrP normalization truth-by-reference; semi-specific DIA | ✔ |
| PXD031616 | same protocol paper, second set | Fusion Lumos | 260 GB; 64 .raw (1.3–7.5 GB) | CC0 | LiP DIA/DDA | ✔ |
| PXD015446 | LiP-Quant 2020 (+ PXD018204, PXD019902) | Q Exactive (+ others) | **1.95 TB**; 430 .raw | CC0 | Drug-target dose–response (known targets as truth) | ✔ |
| PXD034606 | PELSA 2025 | Exploris 480, QE HF | 1.35 TB; 559 .raw | CC0 | High-trypsin LiP variant (optional) | ✔ |

#### CF-MS / SEC-MS

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD007038 | Heusel 2019 SEC-SWATH HEK293 (CCprofiler) | TripleTOF 5600 | 171 GB; 189 .wiff/.wiff.scan | CC0 | Reference SEC-DIA complexome (SCIEX format) | ✔ |
| PXD010288 | Heusel HeLa interphase vs mitosis SEC-SWATH | TripleTOF 6600 | 2.58 TB; 579 **mzXML.gz** (0.55–10.7 GB) | CC0 | Differential complexome (readable mzXML) | ✔ |
| PXD036711 | Collins lab 2023 JPR rapid SEC-SWATH re-org | TripleTOF 5600 | 1.33 TB; mzXML + wiff | CC0 | Perturbation CF-MS | ✔ |
| PXD049340 | Lenz C. 2024 mini-Complexome profiling, FDR-controlled (HEK293) | timsTOF Pro 2 | 478 GB; 288 .d.zip (0.9–7 GB) | CC0 | Bruker diaPASEF CF-MS | ✔ |
| PXD059608 | Guo lab 2025/26 thyroid CF-DIA (Nat Commun doi:10.1038/s41467-026-68686-9) | QE HF | 241 GB; 576 .raw (**0.12–0.58 GB**, 30-min DIA) | CC0 | **Small Thermo CF-DIA fractions; ideal for subsetting** | ✔ |
| PXD044083 | Lee 2024 JPR refined CF-MS human brain | Exploris 480 | 124.5 GB; 192 .raw (0.44–0.83 GB) | CC0 | Thermo CF-MS, 96 fractions × 2 | ✔ |
| PXD012009 | Skinnider 2021 Cell PCP-SILAM mouse tissues | (metadata: Bruker) | 155 GB in **one RawData.zip** | CC0 | Not subsettable; tier L only | ✔ |
| PXD042664 | CFdb (Skinnider/Foster 2023), 2 tar files | various | 311.6 GB | CC0 | Harmonized CF-MS reference | ✔ |

#### TPP

| Accession | Study | Instrument | Sz / files | Licence | Purpose | ✔ |
|---|---|---|---|---|---|---|
| PXD013774 | Franken 2015 Nat Protoc TPP-TR/CCR | Q Exactive | 56 GB; 54 .raw.gz (0.7–1.3 GB) + Mascot .dat | CC0 | TMT10 melting curves (known staurosporine/panobinostat targets) | ✔ |
| PXD009495 | Mateus 2018 MSB E. coli TPP | Q Exactive | 248 GB; 306 .raw | EBI-ToU | Bacterial TPP | ✔ |

Accessions in the brief that turned out wrong or doubtful:
- "PXD000138 Marx": confirmed correct.
- "PXD001468 for Sage": correct only for Sage's *open-search* benchmark. Sage's TMT benchmark was PXD016766, and
  its LFQ was not on these data.
- PXD000647 (tried for TPP) is a Legionella dataset, not TPP. **Do not use.**
- PXD042704, PXD034525 and PXD023325 are not on PRIDE. They are Panorama/jPOST-hosted, and PRIDE API calls
  return 404.

### B2. Raw-file reader availability on Linux

- **Thermo .raw.** ThermoRawFileParser (CompOmics) now ships a **self-contained linux build** (release
  `v.2.0.0-dev`: `ThermoRawFileParser-…-linux.zip`, plus a `net8` build). It no longer needs mono and runs
  without a system .NET. Underneath is Thermo's RawFileReader .NET library, whose licence permits
  redistribution. For a native reader we would need to wrap RawFileReader via .NET hosting or keep TRFP as a
  converter to mzML/Parquet. Plan: **one-time conversion to our columnar format, then amortize**.
- **Bruker .d (TDF/TSF).** **timsrust** (MannLabs; v0.6.7, Oct 2026) is a pure-Rust reader used by Sage. The
  Bruker TDF-SDK `libtimsdata.so` also exists for Linux. Several datasets ship .d as `.d.zip` or `.d.7z`
  (PXD027161), and `p7zip` is needed for the latter.
- **SCIEX .wiff/.wiff2.** There is no native Linux reader. The only route is msconvert under Wine/Docker, and
  neither is available here. We therefore **deprioritize SCIEX-only data**. Where uploaders provided
  mzXML(.gz) (PXD010288, PXD036711), use those.
- **Rust toolchain.** `cargo`/`rustc` are not installed on the node, so a user-level `rustup` install is needed.

---

## Implications for our design

**Engine features forced by the assays:**

1. **A shared mass-offset index serves several needs.** A fragment index that supports (a) open/offset search,
   (b) semi-specific digestion and (c) α-first "unknown-partner-mass" scoring also covers open PTM search,
   glyco (offset + oxonium), LiP (semi-tryptic) and non-cleavable XL. MS Annika 3.0's sparse-matrix
   formulation is the right GPU primitive: batch spectra × peptides as SpMM on cuSPARSE, top-k per spectrum,
   then β lookup by mass on CPU.
2. **Cleavable XL as a mode of linear search.** Doublet detection per spectrum gives two linear peptide-mass
   searches, which reuses the standard path.
3. **FDR is stratified and multi-level by construction.** Use group-specific FDR keyed on peptide class:
   - XL: self vs heteromeric × CSM/ResPair/PPI;
   - LiP: half- vs fully-tryptic;
   - PTM: target site vs decoy-amino-acid FLR, following Ramsbottom 2022 JPR pAla decoys;
   - repository scale: run → experiment → global (picked) protein/site FDR.

   One generic "stratified target–decoy at level L" component should serve all of these.
4. **Experiment-design-aware quantification and alignment.** Use SDRF-driven metadata (fraction index,
   LiP/TrP pairing, temperature/TMT channel). MBR must be **opt-in per design**:
   - same-fraction-index only for CF-MS;
   - never between LiP and TrP;
   - transfer FDR always estimated and provenance-flagged.

   Do not apply global intensity normalization to CF-MS or TPP series.
5. **Entrapment as a first-class harness.** Tagged entrapment sequences (shuffled-paired with r=1, or foreign
   proteome) are propagated to every output level. FDP is computed automatically as lower bound, combined and
   paired. Our continuous-integration regression should *fail* when the combined upper bound exceeds nominal
   FDR by more than a tolerance on any S-tier dataset.

### Proposed tiered benchmark suite

Sizes are the downloaded raw files, and all files are Thermo .raw unless noted.

**Tier S: fast dev loop (≈9.6 GB, minutes per run on 4 cores + 3090). Every commit.**

| ID | Dataset/files | GB | Tests |
|---|---|---|---|
| S1 | PXD001077 `Velos005137.raw` (Pfu) | 0.33 | DDA ID, FDR via foreign-proteome entrapment (Pfu+human DB) |
| S2 | PXD007058 `pPeptideLibrary_pool{1..4}_HCDOT_rep1.raw` | ~0.85 | Phosphosite FLR vs known sites; decoy-amino-acid FLR calibration |
| S3 | PXD014337 `QEx-HFX_DSS_R1.raw` | 0.96 | Non-cleavable XL, real FDR from library design (O(n²) path) |
| S4 | PXD029252 `DSSO_mainlib_1_5HEK_ctrl.raw` | 1.35 | Cleavable XL in human background; inter vs intra real FDR (compare to bundled xi mzid) |
| S5 | ProteoBench DIA-plasma PYE `A9_..._R1,R2` + `B9_..._R1,R2` (.d.zip) | ~4.4 | Bruker reader (timsrust), DIA quant ratios (1, 1/3, 2), DIA entrapment by appended shuffled DB |
| S6 | PXD059608: 5 consecutive fractions of one CF-DIA series (smallest ~0.12–0.3 GB each) | ~1.5 | Fraction-aware MBR guard (transfers must stay within fraction index) |
| | **Total** | **≈9.4–9.6** | |

**Tier M: nightly / pre-merge (≈200 GB, hours).** Covers all assay types with replicates and compares against
published numbers.

- **Quant leaderboard:** ProteoBench DDA-QExactive (19.4 GB), DDA-Astral (15.9), DIA-Astral (14.1), DIA-plasma
  full (14.1). Submit results to ProteoBench for external comparison.
- **FDR:** ProteoBench DIA-entrapment-Astral (5.6), plus 4 × PXD001468 HEK files (~13). The HEK files give a
  direct runtime and ID comparison against Sage, MSFragger and the Wen 2025 DDA entrapment results.
- **PTM:** PXD000138 subset of ~40 pools, HCD and CID/ETD (~12). PXD009449 phospho, acetyl, GG and methyl
  pools (~8). PXD014525 synthetic-spike DIA subset (~10).
- **XL:**
  - PXD014337 full set (9.9).
  - PXD029252 DSSO/DSBU mainlib reps plus HEK-SEC fractions 11–15 (~15).
  - PXD027161 Bruker XL (3.4).
  - PXD042173: ~10 files from one batch (~10). Needs the group map from the Scout supplement.
- **LiP:** PXD031627, 4 LiP + 4 TrP DIA runs plus the provided Spectronaut/MSstatsLiP outputs as reference
  (~27).
- **CF-MS:** PXD059608, one full replicate series (~40 fractions, ~15). PXD044083, one 96-fraction series as an
  optional addition (~60).
- **TPP:** PXD013774, one TMT10 experiment (~10).
- **TMT:** PXD016766, Sage's 12 files (~2.5).

**Tier L: release / paper-grade (multi-TB, days; staged via Aspera/FTP to scratch, one dataset at a time).**
Measures repository-scale throughput, global FDR accumulation and multi-vendor generality.

- **Multi-vendor quant:** PXD028735 (3.4 TB), PXD070049 (3.1 TB, includes ZenoTOF; SCIEX parts skipped until a
  wiff route exists).
- **Large-scale ID:** PXD004732 ProteomeTools (1.1 TB). Gives synthetic truth at scale and tests amortized
  indexing across 1,460 runs.
- **PPI-level XL:** PXD042173 full (234 GB) and PXD019120 (331 GB), checking PPI-FDR at proteome scale against
  the 590-PPI reference.
- **Drug-target LiP:** PXD015446 LiP-Quant (1.95 TB). Known drug targets serve as truth.
- **SEC complexomes:** PXD010288 (2.6 TB mzXML) and PXD049340 (478 GB, Bruker).
- **Astral DIA at scale:** PXD046283 and PXD049412, for runtime per run at Astral data rates and single-cell
  FDR inflation.

Acceptance criteria, applied at every tier and measured against the reference tool on the same tier:
- runtime per GB of raw data and per run;
- IDs at 1% FDR;
- combined-entrapment FDP ≤ nominal + tolerance;
- FLR on synthetic phospho;
- real XL FDR at CSM/ResPair/PPI;
- ProteoBench ratio error (median |log2FC error|, CV);
- LiP adjusted-fold-change agreement with MSstatsLiP;
- CF-MS: fraction of transfers crossing fraction index (must be 0 in strict mode).

Reference tools: Sage/MSFragger/DIA-NN for ID, xiSEARCH+xiFDR/MS Annika/Scout for XL, FragPipe+FLiPPR or
Spectronaut+MSstatsLiP for LiP.

### Open items / caveats

- Glyco benchmark accessions are not verified. The Kawahara consensus is in GlyConnect 2943, but the raw-file
  location was not located.
- ISB18 and MSV000084000 (Wen 2025) are not verified.
- The LiP PXD031627 run-to-condition mapping (which LM34xx are LiP vs TrP, and the 8 small runs) must be read
  from the paper/SDRF before use. The community SDRF lists all runs as DIA, trypsin.
- The library→group truth tables for the XL sets live in paper supplements and the IMP-X-FDR GitHub, not in
  PRIDE. Fetch them alongside the raw files and checksum them in DATA.md.
- The MBR-in-CF-MS harm is argued mechanistically here. A quick experiment would quantify it on S6/M and would
  be a small publishable result: MBR on vs off vs fraction-aware, measuring correlation inflation for
  known-non-interacting pairs.
