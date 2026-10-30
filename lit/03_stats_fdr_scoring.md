# 03 — Statistics, FDR control and neural scoring

Scope: what is proven (and under which assumptions) about target-decoy FDR control; ML post-processors and their failure modes; entrapment-based validation; peptide/protein/repository-scale FDR; PTM localization (FLR); open search; XL-MS; neural "direct ranking" scorers and calibration. Compiled 2026-10-26. DOIs were checked against Crossref/PubMed. Claims marked (P) are proven theorems under the stated assumptions. Claims marked (E) are empirical findings only.

---

## 0. TL;DR

1. **TDC controls FDR only under two assumptions**: (a) *equal chance*: a null (incorrect) match is equally likely to be target or decoy; (b) *independence* of that coin flip across hypotheses and from the scores. Under these, TDC with the "+1" ((D+1)/T) controls FDR in finite samples (P) (Barber & Candès 2015; He et al. 2015). **Spectrum multiplicity breaks (b) at the PSM level.** Do FDR at the **peptide or precursor level** with paired target–decoy competition ("double competition").
2. **FDR is an expectation.** The realised FDP of a single list can be well above α, especially for small lists. Tools exist to bound the FDP (TDC-UB, FDP-SD).
3. **Learned rescoring can leak.** If the same decoys are used to train the model and to estimate FDR, finite-sample control is lost (P/E): Percolator can fail, mostly through multiplicity across CV folds. Two fixes are valid: a **static model trained on disjoint data**, or **RESET**, which splits the decoys into training and estimating sets.
4. **Learned scores add a second leak: sequence priors.** A model that learned what real peptides look like can down-score decoys even on null spectra. That violates equal chance and makes the FDR estimate liberal (Danilova 2019). Sequence-only classifiers separate decoys from targets with AUC 0.64 for reversed decoys and up to 0.81 for DIA-NN-style decoys (Reznikov 2026).
5. **Entrapment is the external audit, and most papers do it wrong.** Use the *combined* (upper bound) or *paired* (tighter upper bound) estimator. The common N_E/(N_T+N_E) is only a *lower bound*, so it can show failure but never success (Wen et al. 2025, Nat Methods). With valid estimators, **no DIA tool consistently controlled FDR**: DIA-NN 1.8.1, Spectronaut 18 and EncyclopeDIA 2.12 all failed somewhere, worst at protein level and on single-cell data.
6. **At scale, FDR must be recomputed globally.** Naively merging per-dataset 1% results gave **28% precursor-level FDR** in MassIVE-KB. The classic protein FDR over-counts decoys at scale; use the picked protein FDR, or the **Picked Protein Group FDR** for protein groups. Occam's-razor grouping becomes anticonservative at scale (The et al. 2022).
7. **Each assay has its own unit of inference.** For PTMs, PSM FDR ≠ FLR ≠ site FDR; use decoy amino acids (Ala/Leu) to estimate FLR. For XL-MS, compute FDR separately for self and heteromeric links, at the level of the claim (residue pair or PPI), after merging. A 5% CSM-level FDR gave 35–60% actual PPI error. Open search and mass offsets need group-aware FDR (Group-walk, CONGA). Naive per-group TDC is *not* rigorous.

---

## 1. Target-decoy competition: what is actually proven

### 1.1 Basic procedure
- **Elias & Gygi 2007**: concatenated target+decoy search. Each spectrum keeps its best match, and FDR is estimated as D/T above a threshold. Nat Methods. https://doi.org/10.1038/nmeth1019
- Estimator with the finite-sample correction: FDR̂(s) = (D(s)+1)/max(T(s),1). q-values are the monotone (cumulative-minimum) version.
- **He et al. 2015** (arXiv:1501.00537) proved three things (P):
  - concatenated TDS gives a *conservative estimate* of FDR at any fixed threshold, but cannot *rigorously control* FDR without modification;
  - with the +1 correction, **peptide-level FDR is rigorously controlled**;
  - **spectrum-level control is difficult**.
- **Levitsky et al. 2017** (JPR, https://doi.org/10.1021/acs.jproteome.6b00144) reached the same correction by a different route. They recommend confidence intervals on the FDP for small lists and very low thresholds.
- **Barber & Candès 2015** knockoff+ (Ann Stat, https://doi.org/10.1214/15-AOS1337) is the general proof: competition-based FDR with "+1" controls FDR when the null target/decoy labels are i.i.d. fair coins independent of the scores. **Rajchert & Keich 2023** (Stat Prob Lett, https://doi.org/10.1016/j.spl.2023.109819; arXiv:2204.13248) studied whether the +1 is needed. Keep it.

### 1.2 Why PSM-level TDC is liberal: multiplicity
**Lin, Short, Noble, Keich 2022** (JPR, https://doi.org/10.1021/acs.jproteome.2c00282) confirmed He et al. empirically.
- If many spectra come from one analyte absent from the database, they all match the *same* wrong peptide. That one coin flip is then counted many times, which violates independence.
- Extreme example: 100 spectra of an absent peptide land on a target with probability ½, and all 100 pass at 1% FDR.
- Three peptide-level schemes were compared:
  - "PSM-only": PSM TDC, then best PSM per peptide, then peptide TDC. This is the most common and the least powerful.
  - "peptide-only": separate target and decoy searches, then peptide-level pairwise competition.
  - **"PSM-and-peptide" (double competition)**: PSM-level competition, then each target peptide competes head-to-head with *its own paired decoy*. This gave **+17% peptides at 1% FDR** over PSM-only with control (E, entrapment). It requires a stored target↔decoy pairing.
- Caveat they state: equal chance is questionable when the target DB contains many "close neighbours" (peptides with highly similar theoretical spectra).

### 1.3 FDR vs FDP
- **Ebadi et al. 2023** "Bridging the False Discovery Gap" (JPR, https://doi.org/10.1021/acs.jproteome.3c00176) shows that TDC controls only the *expected* FDP. A particular list can exceed α. Two remedies:
  - **FDP-SD** (step-down) controls the FDP at a chosen confidence level. Luo et al. 2023 Biometrics, https://doi.org/10.1111/biom.13830; R package `stepdownfdp`.
  - **TDC-SB/TDC-UB** give upper prediction bounds on the FDP of TDC's list (arXiv:2302.11837).
- **Decoy-induced variability.** Keich, Tamura & Noble 2019 (JPR, https://doi.org/10.1021/acs.jproteome.8b00802) measured the spread in accepted PSMs across ten decoy DBs:
  - 4.7% at 1% FDR on a 15k-spectrum set;
  - 10% for 1,000 spectra at 5% FDR.
  - **aTDC** averages over several shuffled decoy DBs. Its improved variant aTDC1+ *empirically* controls FDR with more power. **There is no proof (E only).**
- **Multiple competition.** Emery, Hasam, Noble, Keich (RECOMB 2020, https://doi.org/10.1007/978-3-030-45257-5_4; arXiv:1907.01458) is a framework with d decoys per target that **rigorously controls FDR in finite samples** (P, assuming the null target and its d decoys are exchangeable). It gave **up to 50% more peptides at small FDR thresholds** (E). A related method for multiple knockoffs: Emery & Keich, arXiv:1911.09442.

### 1.4 Separate search, calibration and mix-max
- **Keich, Kertesz-Farkas, Noble 2015** (JPR, https://doi.org/10.1021/acs.jproteome.5b00081) showed that the two "separate target-decoy search" protocols are biased. The unbiased TDC variant discards a random subset of high-scoring targets.
  - **Mix-max** gives unbiased FDR estimates *provided the scores are well calibrated*: a decoy score and an incorrect-target score must have the same distribution for every spectrum. It also estimates π0 and does not discard winners.
- **Calibration matters.** Keich & Noble 2015 (JPR, https://doi.org/10.1021/pr5010983) showed:
  - an uncalibrated score (e.g., raw XCorr) ranks spectra inconsistently across the whole dataset;
  - calibrating it raises IDs at fixed FDR;
  - calibration also makes mix-max less conservative.
- Calibrated score families:
  - exact XCorr p-values: Howbert & Noble 2014, MCP, https://doi.org/10.1074/mcp.O113.036327;
  - res-ev, a high-resolution exact p-value: Lin et al. 2018, bioRxiv https://doi.org/10.1101/290858;
  - MS-GF+ spectral E-values: Kim & Pevzner 2014, Nat Commun, https://doi.org/10.1038/ncomms6277.
- Relevance to us: any **learned score must be calibrated per spectrum** (comparable across spectra, runs and instruments). Otherwise global ranking, and mix-max-style estimators, degrade.

### 1.5 Groups, subsets and side information
- **Applying TDC separately to each group does not rigorously control FDR** (P, counter-example). **Group-walk** (Freestone et al. 2022, Bioinformatics, https://doi.org/10.1093/bioinformatics/btac471) is a version of AdaPT for competition that controls FDR while using group structure. Gains over pooled TDC at 1% FDR: +4% grouped by charge, +3.6% by peptide length, **+26% by modification mass difference** (E). R package `groupwalk`.
- **Subset search** (searching all peptides but reporting only some): Sticker et al. 2017 (Nat Methods, https://doi.org/10.1038/nmeth.4338); response by Noble & Keich 2017 (https://doi.org/10.1038/nmeth.4339). Lin et al. 2021 (JPR 20(8):4153) proposed subset-neighbour search, which needs neighbour handling to avoid liberal FDR.
- **Decoy-free mixture models** exist (Peng et al. 2020, Bioinformatics, https://doi.org/10.1093/bioinformatics/btaa807). They are parametric and offer no finite-sample guarantee. Use them only as diagnostics.

### 1.6 Decoy construction biases
- **Danilova et al. 2019** (JPR, https://doi.org/10.1021/acs.jproteome.8b00991): equal chance "is often violated in practice, even by popular methods". *Biased scoring functions* favour targets, which makes FDR estimates liberal. *Correlated spectra* favour decoys, which makes them conservative.
- **Reznikov et al. 2026** (bioRxiv, https://doi.org/10.64898/2026.03.27.714819) trained sequence-only neural classifiers to tell targets from decoys:
  - reversed decoys: AUC ≈ 0.64;
  - shuffled and DIA-NN-style decoys: ≈ 0.70–0.81;
  - most ESM2-generated decoys were harder to separate than reversed ones.
  - **Reversed decoys remained a strong end-to-end baseline**, and decoy realism was "only weakly informative" for search performance. Estimated FDR tracked entrapment FDP for reversed and ESM2 decoys (yeast sample, maize entrapment).
- **Implication.** Any scorer with a learned *peptide prior* (a decoder language model, a detectability or RT predictor trained on real peptides) can turn this AUC>0.5 into a target preference. That is a direct violation of equal chance.

---

## 2. Semi-supervised post-processors: Percolator, mokapot and their leaks

- **Percolator**: Käll et al. 2007, Nat Methods, https://doi.org/10.1038/nmeth1113. Iterative SVM: decoys are negatives, high-confidence targets are positives. **3-fold cross-validation** (Granholm, Noble, Käll 2012, BMC Bioinf, https://doi.org/10.1186/1471-2105-13-S16-S3) is meant to keep the train and test target/decoy labels apart.
- **Percolator 3.0**: The et al. 2016, JASMS, https://doi.org/10.1007/s13361-016-1460-7. Scales to millions of PSMs (25M spectra of the draft human proteome). Protein inference is by best-scoring peptide over groups of proteins with identical theoretical peptide sets.
- **mokapot**: Fondrie & Noble 2021, JPR, https://doi.org/10.1021/acs.jproteome.0c01010. A flexible Python re-implementation with the same semi-supervised, cross-validated design, so it inherits the same caveats.
- **PEP vs q-value**: Käll et al. 2008, JPR, https://doi.org/10.1021/pr700739d.

**Documented failure (E plus argument):**
- Freestone, Noble, Keich 2024 (JPR 23:1907, https://doi.org/10.1021/acs.jproteome.3c00902) found that Percolator and PeptideProphet can fail to control FDR. Decoys serve two roles: training the classifier and estimating FDR. *Multiplicity* (near-identical PSMs from the same analyte in different CV folds) lets training "peek" at the test labels. This is "essentially guaranteed when multiple runs are combined."
- The same paper rebuts Kong et al.'s claim that narrow-window FDR is generally compromised. After accounting for chimeric spectra and for the different numbers of candidates in narrow vs open searches, standard TDC showed no evidence of failure.

**The fix: RESET / Percolator-RESET** (P under TDC assumptions).
- References: Freestone, Käll, Noble, Keich, RECOMB 2024 (LNCS 14758:448); JPR 2025 24(5):2266, preprint https://doi.org/10.1101/2023.10.26.564068.
- Procedure:
  - Randomly assign each decoy (with p = ½) as either a *training* decoy or an *estimating* decoy.
  - Train only on training decoys (negatives) against "pseudo-targets" (targets plus estimating decoys, unlabeled).
  - Rescore, discard the training decoys, then run competition-based selection (SeqStep+) with FDR̂ = (c/(1−c))·(D_est+1)/T.
  - Single-decoy mode: c/(1−c)=2, because half the decoys were spent on training.
  - Two-decoy mode (two decoys per target; one set trains, one estimates): c/(1−c)=1. Training-phase constants are 3 and 2 respectively.
- Results: single-decoy mode is marginally less powerful than Percolator. Two-decoy mode is marginally *more* powerful and less variable. Slated for Percolator v4.0.
- **Peptide-level training** (one entry per peptide) removes the multiplicity leak.

**Static models.** Fondrie & Noble 2020 (JPR, https://doi.org/10.1021/acs.jproteome.9b00780) learn a Percolator model on a large external dataset and apply it frozen to a small experiment. This increased yield and removed model-induced variability. For FDR it matters because a model fixed before seeing the test data cannot leak the test labels. Equal chance under the learned score still has to hold.

**Rescoring tools built on predicted spectra and RT:**
- Prosit: Gessulat 2019, Nat Methods, https://doi.org/10.1038/s41592-019-0426-7.
- Prosit rescoring for immunopeptides: Wilhelm 2021, Nat Commun, https://doi.org/10.1038/s41467-021-23713-9.
- MS²Rescore: Declercq 2022, MCP, https://doi.org/10.1016/j.mcpro.2022.100266; v3.0 Buur 2024, JPR.
- Oktoberfest: Picciani 2024, Proteomics, https://doi.org/10.1002/pmic.202300112.
- inSPIRE: Cormican 2022, MCP, https://doi.org/10.1016/j.mcpro.2022.100432.
- MSBooster: Yang 2023, Nat Commun, https://doi.org/10.1038/s41467-023-40129-9.
- AlphaPeptDeep: Zeng 2022, Nat Commun, https://doi.org/10.1038/s41467-022-34904-3.
- CHIMERYS: Frejno 2025, Nat Methods, https://doi.org/10.1038/s41592-025-02663-w.

All of these push features into Percolator or mokapot, so they inherit §2's leaks. Evidence:
- **Winkelhardt, Berres, Uszkoreit 2026** (bioRxiv, https://doi.org/10.64898/2026.03.04.709532; FDRBench entrapment over 7 engines):
  - Including all features (notably −log MS-GF+ E-value) in MS²Rescore and Oktoberfest produced "nearly perfect separation between target and decoy PSMs, clearly violating the assumptions." Decoys stopped modelling incorrect targets. Removing that one feature fixed it.
  - Prediction-based rescoring "exhibit[s] a greater tendency to slightly underestimate the true FDR."
  - **Lesson:** feature sets must be audited, and near-perfect target/decoy separation is a red flag, not a success.
- **Castaño & Beaudry 2025** (Proteomics, https://doi.org/10.1002/pmic.202400225): +40–53% peptides and +64–67% PSMs from rescoring MaxQuant output. Up to 75% of the *lost* peptides carried PTMs.

---

## 3. Entrapment: the external audit

### 3.1 Origin
Granholm, Noble, Käll 2011 (JPR, https://doi.org/10.1021/pr200517n) used samples of known content plus foreign sequences to check calibration.

### 3.2 Wen et al. 2025 — the reference framework
Wen, Freestone, Riffle, MacCoss, Noble, Keich, *Nat Methods* 22:1454–1463, https://doi.org/10.1038/s41592-025-02719-x; preprint https://doi.org/10.1101/2024.06.01.596967; software **FDRBench** https://github.com/Noble-Lab/FDRBench.

Setup:
- Search against target ∪ entrapment, with the tool's own decoys added on top.
- Take the tool's list at its reported q ≤ α.
- N_T = original-target discoveries; N_E = entrapment discoveries; r = |entrapment DB| / |target DB|.

| Estimator | Formula | Status |
|---|---|---|
| **Combined** | FDP̂ = N_E·(1 + 1/r) / (N_T + N_E) | **Estimated upper bound**, under an equal-chance analogue: a false discovery is at least as likely to be entrapment as its share r/(1+r). Valid for showing control, but under-powered. |
| **Lower bound** (widely misused) | N_E / (N_T + N_E) | **Lower bound only**. Can show *failure* of control, never success. |
| **Sample** | N_E·(1/r) / N_T | **Invalid**: usually underestimates but can overestimate. Do not use. |
| **Paired** (r = 1, each target paired with a unique shuffled entrapment) | FDP̂* = (N_E + N_{E≥s>T} + 2·N_{E>T≥s}) / (N_T + N_E) | **Tighter upper bound**. N_{E≥s>T}: discovered entrapments whose paired target scored below the cutoff s. N_{E>T≥s}: both discovered, with the entrapment scoring higher. A k-matched generalisation (r = k) is in their Suppl. Note 2.3. |

- Their Table 1 lists published tool papers and benchmarks that "established" FDR control with the lower bound: Lancaster (Spectronaut), Scott (GPS), Zhang (Mzion), Lou et al. 2023 (Nat Commun DIA benchmark, https://doi.org/10.1038/s41467-022-35740-1). Those conclusions are therefore unsupported.
- **Entrapment construction**:
  - Shuffle each target peptide with its C-terminal residue fixed (lengths 7–35).
  - Make each shuffled entrapment unique against all targets and earlier entrapments; FDRBench has `-fix_nc c`, `-I2L` and `-check`.
  - Alternatively use foreign-species peptides, with any match to a target removed. Which species to use is an open question ("what is the right evolutionary distance").
  - Large r moves away from real usage, so prefer **r = 1 (paired) or r = 1–2 (combined)**. With r ≥ 5 the combined and lower-bound estimates differ by ≤ 20%.
- **Protein-level entrapment**:
  - Build a paired shuffled protein per target protein. A peptide shared by several proteins must be *consistently* swapped for the same entrapment peptide in every protein.
  - A protein group containing ≥ 1 original target counts as a target.
- **Findings**:
  - DDA, peptide level, paired estimator: Tide+Percolator-RESET, Sage (LDA), MS-GF+ and MSFragger+MSBooster were all close to the diagonal (control supported). MaxQuant passed at protein level.
  - DIA, versions DIA-NN 1.8.1, Spectronaut 18.7, EncyclopeDIA 2.12.30, ten datasets: **no tool consistently controlled FDR**.
    - Single-cell precursor-level *lower bounds* at 1%: DIA-NN 2.3–4.7%, Spectronaut 3.8–7.6%.
    - Protein level: Spectronaut and EncyclopeDIA consistently invalid; DIA-NN mostly invalid. EncyclopeDIA reached 6.7–9.1% on human-lumos.
    - On single-cell data, DIA-NN's discovery list is inflated by up to 48.3% (precursor) and 44.2% (protein) relative to entrapment-guided cutoffs.
  - Newer DIA-NN versions (2.x) were not tested in that paper.

### 3.3 Other validation designs
- **PyViscount** (Madej & Lam 2025, JPR, https://doi.org/10.1021/acs.jproteome.4c00743) builds a quasi ground truth by randomly partitioning the search space, using unaltered spectra. It complements entrapment.
- **PIP-ECHO** (Solivais et al. 2025, JPR, https://doi.org/10.1021/acs.jproteome.5c00065) is a protocol for controlling and evaluating the FDR of match-between-runs (MBR).
  - MBR accounts for up to 40% of IDs, and up to 75% in single-cell data.
  - Only PIP-ECHO controlled MBR FDR at 1% across benchmarks; FlashLFQ, IonQuant and MaxQuant did not.
  - IonQuant's mixture-model MBR-FDR: Yu et al. 2021, MCP, https://doi.org/10.1016/j.mcpro.2021.100077.

---

## 4. Peptide → protein → repository scale

### 4.1 Protein FDR
- **Picked protein FDR**: Savitski et al. 2015, MCP, https://doi.org/10.1074/mcp.M114.046995; built on about 19,000 ProteomicsDB runs.
  - Each target protein competes with *its own* decoy; only the higher-scoring of the two survives. The protein score is the best peptide q-value (or PEP).
  - Why classic protein TDC fails at scale: "false target PSMs may arise from both correct and false target proteins, whereas decoy PSMs only arise from decoy proteins". Decoy proteins accumulate roughly linearly with data size while true targets saturate, so classic FDR over-predicts false proteins.
- **Picked Protein Group FDR**: The, Samaras, Kuster, Wilhelm 2022, MCP, https://doi.org/10.1016/j.mcpro.2022.100437. Adds three things:
  - (i) **Picked group competition**: walk down groups sorted by score; eliminate any lower-scoring group containing a counterpart (decoy↔target) of a leading protein.
  - (ii) **Rescued subset grouping**: regroup after a 1%-protein-FDR-equivalent PSM filter, which cut the share of discarded peptides from 21% to 15%.
  - (iii) Score = best Percolator PEP. Products of MaxQuant PEPs were anticonservative because those PEPs are poorly calibrated.
  - Validation by entrapment (shuffled 50% or 4% of peptides, 5× database) and by simulations combining up to 500 experiments: **razor/Occam's-razor peptide assignment became increasingly anticonservative**. It reached "no protein groups at 1% entrapment FDR" beyond about 300 experiments on Swiss-Prot+isoforms+TrEMBL.
  - ProteomicsDB reanalysis (≈410M PSMs, 19.8k runs): 18,000 protein groups at 1% group-level FDR.
- **MAYU**: Reiter et al. 2009, MCP, https://doi.org/10.1074/mcp.M900317-MCP200. The first to show that protein FDR inflates with dataset size even at a fixed PSM FDR, and to model it via the protein-size-dependent probability of decoy hits.
- **DIA at scale**: Rosenberger et al. 2017, Nat Methods, https://doi.org/10.1038/nmeth.4398 (PyProphet). Defines run-specific, experiment-wide and **global** contexts, and proposes a *global analyte constraint*: an analyte must pass a global FDR as well as a run-level FDR, to stop false positives accumulating across runs.

### 4.2 Repository-scale aggregation
- **MassIVE-KB** (Wang et al. 2018, Cell Syst, https://doi.org/10.1016/j.cels.2018.08.004): 227 datasets, 658M spectra.
  - Merging per-dataset 1% PSM-FDR results naively would give **28% library-level precursor FDR**.
  - Fix: global 1% precursor FDR, stratified by peptide length (sliding window of 500 precursors), followed by **picked protein FDR** at 1%.
  - Novel proteins follow the HPP guidelines (≥ 2 non-overlapping uniquely mapping peptides ≥ 9 aa; Deutsch et al., HPP MS Guidelines 3.0, JPR 2019).
- **ProteomicsDB**: Picked Group FDR (above).
- **quantms**: Dai et al. 2024, Nat Methods, https://doi.org/10.1038/s41592-024-02343-1. Cloud reanalysis of PRIDE.
- **PeptideAtlas**: MAYU-tuned PSM thresholds.
- **Principle.** Every claim at a given level (PSM, precursor, peptide, protein group, site, PPI) needs its own competition over the full aggregated set, using the best per-entity evidence across runs. Run-level q-values are not additive.

---

## 5. DIA and neural classifiers inside search engines

- DIA-NN (Demichev et al. 2020, Nat Methods, https://doi.org/10.1038/s41592-019-0638-x) trains an ensemble of neural nets per run to separate target and decoy precursors.
- AlphaDIA (Wallmann et al. 2025, Nat Biotechnol, https://doi.org/10.1038/s41587-025-02791-w) uses a semi-supervised target/decoy classifier with options for TDC or "fragment competition".
- Both use decoys for training *and* estimation within the same data, which is the structure RESET was designed to fix. Wen et al. measured the result: DIA tools' failures are concentrated in **low-signal (single-cell) data** and at the **protein level**.
- Pre-trained DIA models aim at zero-shot scoring that does not depend on per-run training:
  - DIA-BERT (Liu et al. 2025, Nat Commun, https://doi.org/10.1038/s41467-025-58866-4): trained on 276M precursors; claims +51% proteins over DIA-NN.
  - DIA-CLIP (Liao et al. 2026, arXiv:2602.01772): dual-encoder contrastive model; claims +45% proteins and 12% fewer entrapment IDs.
  - **Treat their identification gains as unverified** until they are tested with the paired/combined estimators at protein level.

---

## 6. PTM localization: FLR and site-level FDR

- **Separate the quantities.** PSM/peptide FDR says the modified peptide is right. **FLR** is the fraction of reported sites placed on the wrong residue. **Site-level FDR** applies after collapsing PSMs into (protein, residue) sites, possibly across datasets.
- **Localization tools**:
  - Ascore: Beausoleil 2006, Nat Biotechnol, https://doi.org/10.1038/nbt1240.
  - phosphoRS: Taus 2011, JPR, https://doi.org/10.1021/pr200611n.
  - Strategies review: Chalkley & Clauser 2012, MCP, https://doi.org/10.1074/mcp.R111.015305.
  - **LuciPHOr**: Fermin et al. 2013, MCP, https://doi.org/10.1074/mcp.M113.028928. A modified target-decoy approach: decoy phosphopeptides carry the modification on *non-candidate* residues and compete with real site patterns to estimate FLR.
  - **PTMProphet**: Shteynberg et al. 2019, JPR, https://doi.org/10.1021/acs.jproteome.9b00205. A Bayesian mixture model giving per-site probabilities; global FLR is computed from those probabilities, so it is model-based.
- **Decoy amino acid FLR**: Ramsbottom et al. 2022, JPR, https://doi.org/10.1021/acs.jproteome.1c00827.
  - Allow phosphorylation on a residue that cannot carry it, and count those localizations as false. This gives an independent global FLR.
  - **Alanine and leucine were the most reliable decoys.** The choice of decoy residue changes the estimated FLR substantially.
- **Aggregating across datasets**: Camacho et al. 2023, JPR, https://doi.org/10.1021/acs.jproteome.2c00823. Collapses PSM → peptidoform-site → multi-study while tracking decoy-AA FLR. On 8 rice datasets: 6,368 sites vs 4,687 by conventional thresholds with unknown FLR.
- **Why it matters at scale**: Kalyuzhnyy et al. 2022, JPR, https://doi.org/10.1021/acs.jproteome.2c00131. Sites in PhosphoSitePlus/PeptideAtlas supported by a single piece of evidence have estimated FDRs of **84% (pSer), 98% (pThr) and 82% (pTyr)**.
- **Ground-truth sets**:
  - Marx et al. 2013 (Nat Biotechnol, https://doi.org/10.1038/nbt.2585): >100,000 synthetic peptides plus phospho counterparts.
  - Ferries et al. 2017 (JPR, https://doi.org/10.1021/acs.jproteome.7b00337): synthetic phosphopeptide library across HCD/EThcD/OT/IT.
  - ProteomeTools PTM (Zolg et al. 2018, MCP, https://doi.org/10.1074/mcp.TIR118.000783): ~5,000 peptides covering 21 PTMs, 11 fragmentation modes, diagnostic ions.
  - Locard-Paulet et al. 2020 (JPR, https://doi.org/10.1021/acs.jproteome.9b00679): 22 pipelines; localization score distributions differ so much that the same threshold gives different FLRs.
- **Caveat on synthetic sets**: they contain no unmodified-isomer competition from a real proteome and lack biological co-isolation. Use them for calibration curves, and use decoy-AA FLR as the in-situ check.

---

## 7. Open search and mass offsets

- Open search (MSFragger: Kong et al. 2017, Nat Methods, https://doi.org/10.1038/nmeth.4256; Chick 2015, Nat Biotechnol, https://doi.org/10.1038/nbt.3267) widens the candidate set. Unmodified peptides that would be found in a narrow search can be "out-competed by a candidate from the wide window with a randomly high score", which costs power.
- Mass-offset groups also have very different null and alternative distributions, so pooled TDC is badly balanced across groups. Per-offset TDC is not rigorous.
- **CONGA** (Freestone, Noble, Keich 2024, JPR, https://doi.org/10.1021/acs.jproteome.3c00399) controls FDR while matching narrow-search power for unmodified peptides and still finding unexpected PTMs. Steps:
  - (i) **Neighbour filtering**: drop lower-scoring PSMs whose peptides share more than τ = 5% of b/y ions with a higher-scoring match.
  - (ii) **Dynamic-level competition**: peptide level, or "stem" level (unmodified sequence) when variable mods are present.
  - (iii) Groups: one narrow group plus data-driven open groups formed by mass-difference binning and KS-test merging.
  - (iv) **Group-walk**.
- Mass-offset annotation (PTM-Shepherd: Geiszler 2021, MCP, https://doi.org/10.1074/mcp.TIR120.002216) is *post hoc* and does not itself control FDR per offset.

---

## 8. Cross-linking MS

- **Basic estimator**: Walzthoeni et al. 2012, Nat Methods, https://doi.org/10.1038/nmeth.2103. FDR = (TD − DD)/TT, where TD counts target–decoy and DD decoy–decoy pairs.
- **Level matters**: Fischer & Rappsilber 2017, Anal Chem, https://doi.org/10.1021/acs.analchem.6b03745; tool xiFDR.
  - FDR must be set at the level of the claim (residue pair or PPI). Missing this gave **8.4% actual error at a 5% target**.
  - Pre-filtering PSMs and peptide pairs increased residue pairs by 33% with no loss of accuracy.
- **Separate self and heteromeric links**: Lenz et al. 2021, Nat Commun, https://doi.org/10.1038/s41467-021-23666-z.
  - The heteromeric random search space was about 10.6× the self-link space.
  - Controlling FDR on all CSMs and then selecting heteromeric ones enriches false positives: **5% CSM FDR gave 35–36% false PPIs**, and 45–60% under entrapment, wrong-crosslinker-mass and shifted-precursor controls.
  - The fix is to **merge CSMs into PPIs first, then estimate PPI-level FDR on heteromeric links only**. This gave 4.8–4.9% actual error at a 5% target and 590 PPIs at 1% in *E. coli*.
  - Validation used several designs: SEC fractions that cannot crosslink, human/yeast entrapment, a wrong crosslinker mass, and shifted precursors.
- **Ground truth**:
  - Beveridge et al. 2020 (Nat Commun, https://doi.org/10.1038/s41467-020-14608-2): synthetic crosslink library. At an estimated 5% FDR, the **actual false rate was 2.4–32%** depending on the pipeline.
  - Matzinger et al. 2022 (Nat Commun, https://doi.org/10.1038/s41467-022-31701-w): synthetic ribosomal complex.
  - Yugandhar et al. 2020 (Nat Methods, https://doi.org/10.1038/s41592-020-0959-9): structure-based validation can drastically underestimate error rates in proteome-wide XL-MS.

---

## 9. Neural scoring for direct ranking (Percolator-free)

### 9.1 Landscape
- **Casanovo-DB**: Ananth et al. 2024, Bioinformatics, https://doi.org/10.1093/bioinformatics/btae218. The Casanovo de novo transformer (Yilmaz 2024, Nat Commun, https://doi.org/10.1038/s41467-024-49731-x; trained on ~30M MassIVE-KB PSMs) is re-purposed as a PSM score.
  - At 1% FDR (shuffled decoys, TDC): **+31–40% peptides vs XCorr, +52–69% vs Hyperscore, +88–102% vs Andromeda**. Percolator adds further gains.
  - Leakage control: test runs were post-2019, and E. coli and yeast test sets avoid human-peptide overlap.
  - **No entrapment check.** The paper does not test whether the decoder's sequence prior favours targets (equal chance).
- **Tesorai Search**: Burq et al., JMB 2026, https://doi.org/10.1016/j.jmb.2026.169940; preprint https://doi.org/10.1101/2024.08.19.606805.
  - A pretrained spectrum–peptide model with "no decoys in training" and no Percolator, trained on ~16,000 runs from 22 studies.
  - +12%/+9%/+21% peptides vs FragPipe/PEAKS/Proteome Discoverer. The authors report entrapment checks, but the estimator details need verification.
- **Bi-encoders**:
  - yHydra (Altenburg et al. 2026, Bioinformatics, https://doi.org/10.1093/bioinformatics/btag297): joint spectrum/peptide Transformer embeddings for GPU open search. Retrieves the top 50 by embedding, rescores with a conventional score, and runs PSM-level TDC.
  - GLEAMS (Bittremieux 2022, Nat Methods, https://doi.org/10.1038/s41592-022-01496-1): spectrum-only embeddings.
- **De novo FDR**:
  - Sanders, Noble, Keich 2025 (bioRxiv, https://doi.org/10.1101/2025.09.12.675837): a competition-based FDR procedure "at least as powerful as database search".
  - Winnow (Mabona et al., arXiv:2509.24952): a calibrator plus decoy-free FDR estimator. Empirical only, with no finite-sample guarantee.

### 9.2 What makes a learned score FDR-valid (synthesis)
Given TDC theory (§1), a learned score S(spectrum, peptide) followed by TDC at peptide level controls FDR if all of the following hold:
1. **S is fixed before the test labels are seen.** Either it was trained on disjoint data (a static model), or the test-set decoys used for training are excluded from estimation (RESET).
2. **Equal chance under S**: for null spectra, P(best target score > best decoy score) = ½. Learned *peptide priors* violate this. Sources include LM decoders (Casanovo), detectability, RT/IM predictors trained only on real peptides, and features like "is this peptide in the training library". Decoys get out-of-distribution predictions, so their scores are systematically worse even when the spectrum is noise.
3. **Independence**: compete at peptide or precursor level (double competition), not PSM level.
4. **No train/test overlap at the PSM or peptide level.** Reprocessing PRIDE with a model trained on PRIDE identifications is the default *failure* case. The model can memorise previously reported targets, including the ~1% false ones. Decoys were never positives, so their scores are biased down.
5. **Symmetric candidate generation.** Any retrieval stage must treat targets and decoys identically: ANN index, coarse filters, top-k pruning and score caching. If decoys have lower retrieval recall than targets (or the reverse), equal chance fails before scoring even starts. **This is specific to GPU/ANN pipelines and is easy to miss.**

### 9.3 Calibration and conformal views
- TDC/knockoff+ is closely related to conformal inference.
  - Conformal p-values need a set of null calibration scores exchangeable with the test nulls: p_j = (1 + #{null cal ≥ s_j})/(n_cal + 1). Applying BH to them controls FDR (Bates et al. 2023, Ann Stat, https://doi.org/10.1214/22-AOS2244).
  - Decoy PSMs can serve as the calibration nulls. This is essentially a separate target–decoy search, so it relies on the same exchangeability (equal chance) and gains no robustness.
  - **Its practical benefit is calibration**: p-values comparable across runs and instruments enable global ranking, mix-max/π0 gains and repository-level aggregation.
- Per-spectrum calibration (Keich & Noble 2015) remains the key to a single global ranking across thousands of runs. A learned score should be checked for calibration *across* runs, instruments and charge states. Diagnostic: compare decoy-score distributions per stratum.

---

## 10. Implications for our design

### 10.1 FDR engine (defaults)
- **Unit of inference**: precursor level for DIA and peptide level for DDA. Use **paired double competition** (PSM then peptide, each target with its own decoy; Lin 2022). Report PSM-level q only as an auxiliary output.
- Estimator: (D+1)/T with monotone q-values. Optionally use multi-decoy (Emery 2020, rigorous) with d = 2–4 decoys to reduce variance; this matters at strict α and for small runs.
- **Decoys**:
  - pseudo-reverse with protease C-term kept; I/L treated as equal;
  - delete decoys that equal any target (and, for DIA, any target with the same mass/fragment pattern);
  - store an explicit target↔decoy pairing id through the whole pipeline (needed for double competition, picked protein, RESET and paired entrapment);
  - fix the decoy-generation seed list and record it in DATA.md.
- **Groups**: when stratifying (charge, length, open-search mass offset, PTM type, XL type), use **Group-walk**, never independent per-group TDC.
- **FDP reporting**: next to each headline list, report a TDC-UB (or FDP-SD) 95% upper bound on the FDP. Report the variability across ≥ 5 decoy seeds.

### 10.2 Scoring model rules (neural/direct ranking)
- **R1, static model with project-disjoint folds.** Train the scorer on PXD projects partitioned into K folds. Every evaluated or reprocessed project is scored only by a model that never saw it: no spectra, and for benchmarks no peptides, from that project. Keep a frozen held-out benchmark (fixed list of PXDs, instruments and species), including non-human species with no peptide overlap with training.
- **R2, no decoy reuse.** If any on-data fitting uses decoys (fine-tuning, calibration, a Percolator-like layer), use **RESET**:
  - single-decoy: train on a random half of the decoys; estimate with the other half and factor 2;
  - two-decoy mode: decoy set A trains, set B estimates, factor 1.
  - Fold splitting must group by peptide (or stem), never by PSM.
- **R3, equal-chance gates** before any FDR claim:
  - (a) *Sequence-only decoy detectability*: train a peptide-only classifier on our targets vs decoys. Report its AUC. If AUC > 0.6, investigate. (0.6 is our threshold; reversed decoys reach ~0.64 per Reznikov 2026.) Also check that the scorer's *peptide-only* component (e.g., LM log-prior) does not separate them.
  - (b) *Null-spectrum test*: score null spectra (spectra searched with precursor masses shifted by +7–50 Da, or wrong-species spectra) against target+decoy. The target win fraction must be 0.5 ± binomial 99% CI, stratified by length and charge.
  - (c) The target–decoy separation must not be near perfect. Near-perfect separation means a leaking feature (Winkelhardt 2026).
- **R4, symmetric retrieval.** Targets and decoys live in the *same* index with identical ANN parameters. Measure top-k retrieval recall separately for targets (on known-correct PSMs) and decoys (exhaustive vs ANN on decoy-favoured spectra), and require the difference to be within 1%. Re-run every entrapment test after *any* retrieval or pruning change.
- **R5, calibration across runs.** Check calibration per run, instrument, charge and length. Use the decoy score quantiles per stratum (aim for KS-statistic comparability) before pooling for global FDR.

### 10.3 Entrapment protocol (mandatory for every release candidate)
- **E1, peptide/precursor level, paired** (primary):
  - FDRBench shuffled entrapment, r = 1, `-fix_nc c -I2L -check`; one entrapment per target, distinct from all targets *and from our decoys* (different seed and shuffle).
  - Search target ∪ entrapment with our own decoys added.
  - Decoys excluded, take the list at q ≤ α ∈ {0.1%, 0.5%, 1%, 5%}.
  - FDP̂_paired = (N_E + N_{E≥s>T} + 2·N_{E>T≥s}) / (N_T + N_E).
- **E2, combined, foreign-species** (guards against learned sequence priors, because shuffled entrapment looks like decoys to a prior-biased model):
  - For human samples, add a natural proteome absent from the sample (plant, e.g., *Arabidopsis*), then remove every peptide identical to a target after I/L equivalence.
  - Keep r ≈ 1 by peptide count, and also run r = 2 as a check.
  - FDP̂_comb = N_E(1 + 1/r)/(N_T + N_E).
  - **Always report the lower bound N_E/(N_T+N_E) as a failure detector**, never as evidence of control.
- **E3, protein level**:
  - Paired shuffled entrapment *proteins*, swapping each shared peptide with the same entrapment peptide in every protein that contains it.
  - Count any group with ≥ 1 original target as a target.
  - Evaluate after picked-group FDR, at both single-run and aggregated (100s of runs) scale.
- **E4, data**: at least one dataset of each type. DDA HeLa/HEK; DIA bulk; **low-input or single-cell DIA** (the hardest case in Wen 2025); timsTOF; immunopeptidomics (non-specific, large search space); open search; phospho.
- **E5, decision rule**:
  - "Control supported" at α if the paired (or combined) upper-bound curve stays ≤ α at every tested α on every dataset.
  - "Control violated" if the lower bound exceeds α.
  - Anything else is "inconclusive", which is reported as such and never as a pass.
  - Repeat with ≥ 3 entrapment seeds and ≥ 3 decoy seeds.
- **E6, scale test**: concatenate N ∈ {1, 10, 100, 1000} runs or projects. Show that the global peptide and protein-group entrapment FDP does not grow with N. This is the MAYU/Savitski/The scenario.
- **E7, external validators**: PyViscount random partitioning as a second opinion; PIP-ECHO protocol for any MBR/PIP.

### 10.4 Repository-scale FDR ("master" lists)
- Aggregate the best evidence per precursor or peptide across all runs, then do **global** paired peptide-level TDC.
- Then do **Picked Protein Group FDR** with rescued subset grouping. Score protein groups by best calibrated peptide PEP. **No razor/Occam assignment** for FDR.
- Report entities that pass **both run-specific and global q ≤ 1%** (global analyte constraint, Rosenberger 2017).
- Novel or "missing" protein claims additionally need HPP-style evidence: ≥ 2 non-nested unique peptides ≥ 9 aa.
- Publish expected absolute false counts (q × N) for every master list, plus FDP upper bounds.
- MBR/PIP transfers carry their own FDR (PIP-ECHO-style) and are never counted as identifications.

### 10.5 PTMs
- Report three separate numbers: peptidoform FDR, **FLR** and site-level FDR.
- FLR is estimated in situ with **decoy amino acids (Ala and Leu)**, added as allowed phospho acceptors during localization only.
- Calibrate localization scores on Marx 2013, Ferries 2017 and ProteomeTools-PTM synthetic sets. Hold out one set entirely for testing.
- Aggregate across datasets following Camacho 2023: PSM → peptidoform-site → protein-site → multi-study, tracking decoy-AA FLR at each step.
- Flag single-evidence sites, given the 80–98% FDR estimates of Kalyuzhnyy 2022.

### 10.6 Open search, XL-MS and other specialised assays
- **Open/mass-offset search**: CONGA-style pipeline. Narrow plus open search, neighbour filtering, stem-level competition, mass-offset groups, Group-walk. Never per-offset independent TDC.
- **XL-MS**:
  - Use TT/TD/DD with (TD−DD)/TT.
  - Merge to the reported level *before* FDR. Compute FDR separately for self and heteromeric links. Report FDR at the residue-pair and PPI levels.
  - Validate on the Beveridge 2020 and Matzinger 2022 synthetic libraries plus the Lenz-style entrapment, wrong-linker-mass and SEC controls.
- **Other assays (CF-MS, LiP-MS)**: any learned or co-elution score that is thresholded needs either a decoy construct with paired competition or a held-out gold standard. Pooled FDR across heterogeneous strata must go through Group-walk.

### 10.7 Things we must not do
- Claim FDR control from N_E/(N_T+N_E), the "sample" estimator, or from consistency with another tool's list.
- Tune features, hyperparameters or score cut-offs on the entrapment results of the *test* datasets. Entrapment is evaluation only; tuning uses separate datasets.
- Train on decoys and estimate on the same decoys. Use PSM-level CV folds. Pool CV folds whose models are not comparably calibrated.
- Report the PSM-level FDR of merged multi-run data as if it were controlled.
- Use classic (non-picked) protein FDR or razor peptides at repository scale.
- Let ANN/GPU pruning treat decoys differently from targets.
- Change decoy generation, entrapment construction or estimator formulas after seeing results. Fix them in DESIGN.md beforehand, together with the seed list.
