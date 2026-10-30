# DESIGN.md: Vulcan, a repository-scale MS proteomics engine

**Version 0.4 (2026-10-27).** v0.4 moves the hardware of record to the A100/32-core node (§2, D5).
**Previous: 0.3.1 (2026-10-26).** v0.3.1 folds in the stats-implementation issues from results/2026-10-26_stats_sims/design_issues.md.
- v0.2 followed three independent critiques of v0.1: reviews/design_v0.1_{stats,hpc,domain}.md.
- v0.3 follows the round-2 critique: reviews/design_v0.2_round2.md.

Change logs are at the bottom.

Every claim below is provisional until an experiment in §10 validates it or kills it. Literature: LITERATURE.md, lit/*.md.

---

## 1. Goal, deliverable, and the bar

**Goal.** Reprocess PRIDE-scale data with throughput limited by hardware rather than by the algorithm, at no loss of:
- identification sensitivity at *empirically verified* FDR;
- PTM localization accuracy (FLR);
- quantitative precision.

These are measured against Sage, MSFragger and DIA-NN on the same data. MSFragger and DIA-NN are benchmark-only comparators.

**Minimum credible deliverable (MCD).** This is what we commit to. Everything else is a gated research arm. The MCD is an open, GPU-accelerated **DDA** engine with:
1. Ingest that reads all centroids and **infers parameters per run** (calibration, tolerances, analyzer/activation, labels).
2. Closed, chimeric, mass-offset/open and semi-specific search, with GPU matching.
3. **Static-first scoring.** A frozen, project-disjoint, low-capacity model per instrument class (linear or GBM over label-blind and predicted features). Because it is never fit on the evaluated data, one score serves TDC, Group-walk, picked protein and L6. RESET rescoring on the evaluated data is the E4 arm.
4. Exhaustive in-search localization with decoy-amino-acid FLR, validated on held-out synthetic sets.
5. MS1 LFQ and TMT quantification.
6. Picked protein groups plus **global repository FDR shown to stay flat as the number of runs N grows**.
7. Cleavable XL-MS and semi-specific LiP-MS, as demonstrations of stratified FDR and design-aware cross-run logic. CF-MS is a demonstration of transfer provenance and transfer FDR.

**Research arms, each with a pre-stated go/no-go:**
- frozen cross-encoder scorer (E6);
- bi-encoder retrieval (E5);
- spectrum clustering with propagation;
- DIA, starting with Astral nDIA;
- diaPASEF;
- non-cleavable proteome-wide XL.

**The bar.**
| Axis | Reference | Note |
|---|---|---|
| DDA speed | Sage 0.14.7 (MIT, built locally) | Measured on the node of record at 32 threads (and 4 threads for M1 continuity), **and** reported per core-second and GPU-second, plus a cost-normalized comparison. Beating CPU Sage with an A100 is expected and is not by itself the claim. |
| DDA sensitivity | Sage+LDA, Sage+mokapot, Comet+mokapot; MSFragger if licensed | Compared at entrapment-verified FDP and at matched FDP. |
| Literature gains (unverified by entrapment) | Casanovo-DB +31–40% vs XCorr; DDA+ chimeric +15–57% | Treated as hypotheses, not bars. |
| DIA (arm) | DIA-NN 2.7 Academia (local) | Wen 2025 showed lower-bound FDR violations for DIA-NN **1.8.1**, not 2.x. |
| Global FDR | MassIVE-KB: merged per-dataset 1% gave ≈28% global FDR | Our global FDP must stay flat in N. |

**"No loss" (pre-registered non-inferiority).**
- **Test:** on each sealed set, the lower 95% bootstrap CI of the ratio (our IDs / comparator IDs) must be ≥ **0.97**. The bootstrap resamples runs × decoy shuffles. Both methods are averaged over the **same 5 shuffled-decoy FASTAs**; Sage accepts external decoys.
- **Where it is computed:** at each method's own q ≤ 1%, and at matched entrapment FDP.
- **When the comparator's control is not "supported":** if its verdict is "violated" or "inconclusive", only the matched-FDP comparison is used. Analogous tests apply to FLR (site recall at matched true FLR), quant (ε and CV with CIs) and speed.

## 2. Hardware of record and the systems objective

**Node of record (from 2026-10-27, decision D5).** Host gpu-node:
- CPU: 32 physical cores (cgroup affinity, scattered over NUMA nodes 2–7) of a 2-socket AMD EPYC 9354 (Zen 4, **AVX-512 incl. VPOPCNTDQ/BF16**), SMT off.
- GPU: 1× NVIDIA A100 80GB PCIe (sm_80, HBM2e ≈1.9 TB/s, FP16 tensor 312 TFLOPS, 40 MB L2), PCIe Gen4.
- RAM: 1.5 TB. Local NVMe `/tmp`: 862 GB (626 GB free). Project on Lustre (shared-fs).
- `os.cpu_count()` returns 64 but the affinity mask has 32 CPUs, so every thread pool is sized from `os.sched_getaffinity(0)` (`vulcan.hw.n_threads()`), never from `cpu_count`.

**Earlier hardware.** M0–M1 ran on a 4-core EPYC 7742 slice plus an RTX 3090 (24 GB). Results from that node (E1 ingest, stats sims, first baselines) stay valid as recorded, with that hardware named. They are not mixed with A100-node numbers.

**Portability rule.** The allocation can change between sessions. Every result records `bench/hwinfo.py` output (CPU model, affinity count, GPU, driver). A speed comparison is valid only if all of its arms ran on the same node in the same session. CUDA code is built for sm_80 and sm_86 (and PTX), so it runs on both GPUs seen so far.

**Objective.** Maximise **runs/hour end-to-end per node**. We still report **CPU core-seconds per run** and **GPU-seconds per run**, because these are what transfer to other allocations and to a cost model.

What changed with the bigger node:
- **CPU is no longer the binding resource.** With 32 cores, 4-core-era rules (3 decode + 1 orchestrate; every byte of matching on the GPU) relax. Work goes to whichever device finishes the whole pipeline soonest. Ingest and decode can run 8–16 files concurrently on CPU while the GPU searches.
- **GPU memory is no longer a slicing problem for common modes.** With 80 GB, the closed tryptic, phospho and semi-specific human indexes (a few to ~30 GB at 4 B/posting) fit resident. Slicing and streaming stay only for open/non-specific/HLA (30–86 GB) and for cohort-shared superset indexes.
- **Comparator fairness.** Sage and Comet get the **same CPU budget as our CPU side** (all 32 affinity cores), and are also reported at 4 threads, which keeps continuity with M1. Beating CPU Sage with an A100 is expected and is not by itself the claim. Claims are made per node, per core-second and per GPU-second, with a cost-normalized comparison (1 A100 ≈ N cores, stated assumption).
- **CPU kernels.** Zen 4 has AVX-512, so the CPU fallback and oracle paths (Sage, numpy/numba) may use it. The CPU path becomes a real second tier (E2c: CPU-only vs GPU vs hybrid), not just an oracle.

**Reporting.** Every speed number reports:
- wall-clock time, CPU core-seconds, GPU-seconds and H2D bytes;
- spectra/s per device;
- the hwinfo fingerprint;
- whether it was measured from raw or from the store, and whether entrapment was on. Entrapment doubles the search space.

Extrapolation to other hardware is labelled as such.

## 3. Architecture

```
raw ─► L0 ingest (direct, all centroids) ─► Spectral store (Arrow IPC, staged on local disk)
                                                │
                     L0.5 pilot pass per run: recalibration, tolerances, analyzer/activation,
                          reporter check, mass-offset histogram, charge/isotope model
                                                │
FASTA+params ─► L1 proteome index (content-hashed; compact peptidoform records; GPU tiled inverted index)
                                                │
L2 candidate generation (GPU K1): precursor-hypothesis table × index → pair-closed top-K
                                                │
L3 scoring (GPU K2): exact rescoring, label-blind features; predicted-spectrum/RT features (cached);
   localization operator (exhaustive isoforms, decoy-AA)            [arm: cross-encoder reranker]
                                                │
L4 statistics: paired double competition → CONGA/Group-walk strata → RESET (if learned on data)
   → picked protein groups → FLR → XL levels
                                                │
L5 cross-run: MS1 features, alignment, design-aware MBR with transfer FDR, quant (LFQ/TMT)
                                                │
L6 repository: global evidence aggregation → global paired TDC → picked groups → master lists
```

**Principles.**
- Each layer is a batch operator over Arrow tables with fixed schemas.
- No per-spectrum Python.
- **Label blindness:**
  - L0–L3 never receive `is_decoy`, the pair side, or the entrapment flag. These live only in L4 join tables.
  - Enforced by a test: swapping the labels of all pairs, **and** permuting DB/pair storage order with hash-assigned ids, leaves every L2/L3 output bit-identical. This catches positional leaks: targets stored first, index-order tie-breaks, and `rev_` prefix parsing.
  - Pair closure is the one sanctioned exception, and it is label-symmetric by construction.

## 4. Layers

### L0: Ingest and the spectral store
- **Readers:**
  - Thermo: direct RAW→Arrow if E1 validates it (ThermoRawFileParser 2.0 Parquet, or in-process RawFileReader reading the centroid stream). Otherwise RAW→mzML→Arrow.
  - Bruker: timsrust.
  - mzML: a fast reader.
- **Store contents:** **all MS1 and MS2 centroids**, nothing discarded. Peak preprocessing (top-k per 100 Th, deisotoping at high resolution only, precursor removal) is a **search-time view** computed on the GPU (K0). It can be ablated and never requires re-ingest.
- **Spectrum table:**

  | Group | Fields |
  |---|---|
  | Identity | u64 global spectrum id, run id, scan, ms level |
  | Timing | rt, injection time |
  | Precursor | precursor m/z, instrument charge (nullable), isolation window and offsets, master scan, SPS ions |
  | Acquisition | activation (HCD/CID/ETD/EThcD) and supplemental energy, NCE, analyzer (OT/IT/TOF), resolution, FAIMS CV, scan filter string |
  | Bruker | IM / frame–scan map |

- **Peak table:** f32 m/z, f32 intensity (f16 allowed in the GPU wire format), using LargeList or chunked offsets.
- A separate block holds TMT reporter peaks, extracted before any top-k.
- **I/O:**
  - Arrow IPC, uncompressed or lz4.
  - Files are staged to local `/tmp` and read sequentially. **We never mmap on Lustre.**
  - The GPU wire format is fixed-stride blocks of (u32 log-ppm bin, f16 intensity) in pinned, double-buffered host memory. "Zero-copy" is claimed only for the IPC→host path.

### L0.5: Pilot pass (per run)
On a seeded 5–10% subsample, run a wide-tolerance closed search with common variable modifications. It outputs:
- precursor and fragment mass-error curves vs m/z and RT, stored as an offset (never baked into the index);
- empirical tolerances;
- analyzer and activation class, cross-checked against metadata;
- presence of TMT/iTRAQ reporters and SILAC/dimethyl labels;
- a label-blind mass-offset histogram for choosing modifications;
- charge/isotope-error priors and an iRT→RT map.

The pilot's IDs are not reported. The pilot uses its own decoys at q ≤ 1% only to pick calibration anchors, and does so label-symmetrically. Experiment **E0** measures IDs with vs without the pilot pass.

### L1: Proteome index
- **Digest and decoys:**
  - Digest → compact peptidoform records: (protein offset, length, mod mask, mass), ~12 B each.
  - Fragments are generated on the fly from residue-mass prefix sums. **A dense fragment matrix is never materialized.**
  - Decoys: pseudo-reverse with the protease C-term residue fixed; I/L equivalent. A target–decoy collision is replaced by a shuffle seeded by hash(peptide, seed). Pair ids are stored.
- **Entrapment.** Entrapment peptides get their own paired decoys: DB = T ∪ E ∪ dec(T) ∪ dec(E).
- **Semi-specific** decoys come from the decoy proteins, so class counts match.
- **GPU tiled inverted index (K1 layout).** CSR over (log-ppm fragment bin × precursor chunk of 8,192 peptidoforms) with u16 local ids, at 4 B/entry. It is sliced by precursor mass so each slice fits in memory. Semi-specific and phospho slices stream.
- **Fragment ion tables per activation class:**
  - HCD: b/y, a2, immonium, and H₃PO₄ loss in scoring;
  - CID-IT: low-resolution binning, −98 loss;
  - ETD/EThcD: c, z•, z+1, and b/y;
  - fragment charge up to z−1.
- **Content hash** = (FASTA bytes, enzyme, **mod superset**, ion types, decoy recipe and seed, entrapment recipe and seed, code version).
  - Project mod sets are masks.
  - Run-level tolerances and calibration are excluded.
- **Family tags and family-max entrapment:**
  - Index entries carry a family tag: T, E_s, dec_d(T), dec_d(E_s), foreign, dec_d(foreign).
  - K1/K2 keep the top-2 per family per hypothesis, so every (decoy seed, entrapment seed) replicate's competition is reconstructed exactly from **one** search.
  - This is validated once against separate searches.
  - Ids are assigned by hash, not storage order, and ties are broken by hash(sequence, seed).

### L2: Candidate generation (GPU)
- **Query unit: the precursor hypothesis.** One row per (spectrum, charge, isotope error, co-isolated MS1 feature).
  - Charge comes from the MS1 feature or the instrument. When unassigned, z = 2–4 is tried.
  - Isotope error is −1…+3.
  - All hypotheses of a spectrum compete together in L4.
- **MS1 features:** a batched isotope-envelope feature finder in the biosaur2/Dinosaur style. It supplies co-isolated precursors for chimeric search, plus LFQ, MBR and TMT purity.
- **K1 kernel:**
  - chunk-major grid with a warp per (hypothesis, chunk);
  - u8 count accumulators in shared memory;
  - threshold append, then a fused warp top-K;
  - **identical thresholds for every candidate, whatever its label.**
- **Pair closure:** each shortlist (top-K, α-first XL, lazy-prediction survivors, isoform sets) is closed under target↔decoy pairing. If a target is scored, its paired decoy is scored with the same features, and vice versa.
- **Modes:**
  - closed;
  - chimeric, with candidates restricted to MS1-feature precursors in the isolation window;
  - mass-offset (Unimod list plus the pilot histogram);
  - open (shifted and unshifted ions);
  - semi-specific;
  - non-specific (arm, for HLA).
- **Batching:** RAM-sized micro-cohorts of 50–200 runs, with hypotheses sorted by mass within each batch. A global cohort sort happens only if E3 shows >20% end-to-end gain including the barrier.

### L3: Scoring
- **K2 exact rescoring** is warp-per-spectrum. Sorted peaks sit in shared memory, and fragments for the top-K are generated on the fly. Features:
  - hyperscore, matched counts and fractions, ion ladders, ppm statistics;
  - Δscore to the **next candidate regardless of label**;
  - charge, length, missed cleavages, mass shift;
  - label-free chimeric rank, taken from the MS1 feature intensity.
- **Predicted features:** fragment intensities (spectral angle, entropy) and RT/IM residuals from an Apache-licensed predictor (AlphaPeptDeep, or MS2PIP/DeepLC through Koina, or our own).
  - Tryptic libraries are **precomputed once** per (proteome, instrument class).
  - Open and PTM search predict lazily for pair-closed survivors, with a per-run cap and a measured cache hit rate.
  - NCE is chosen empirically from pilot IDs.
  - Strata without predictions (XL, unknown offsets) are handled explicitly as their own strata.
- **Localization operator.** Separate from retrieval and exhaustive. For each top peptidoform:
  - enumerate **all** site configurations, including decoy-amino-acid configurations (Ala primary, Leu secondary);
  - score site-determining ions, aware of neutral losses and of c/z ions under EThcD;
  - add predicted-intensity evidence.

  Scoring is **blind to residue identity given masses**: the same ions and losses for every acceptor. Decoy residues enter only here, never in identification.
- **Default scorer: static-first.** A frozen, low-capacity model (linear or GBM over label-blind and predicted features), one per instrument class (OT-HCD, IT-CID, TOF).
  - Trained on project-disjoint training data, with yeast held out entirely (memorization arm M).
  - It is never fit on the evaluated data, so plain paired TDC, Group-walk, picked protein and L6 all share one score.
- **Arm E4: RESET semi-supervised rescoring.** It replaces the default only if it wins at verified FDP under §6.4.
- **Cross-encoder arm (true split-then-joint):**
  - **Split stage:** shared spectrum and peptide token encoders. Their pooled outputs *are* the E5 bi-encoder embeddings (128-d fp16 ≈ 2.5 GB per 10⁷ peptides, cacheable).
  - **Joint stage:** 1–2 cross-attention layers over the token states of the top-K candidates only. Token states are recomputed rather than cached, which avoids a 50–500 GB cache.
  - A MaxSim/late-interaction head would be a bi-encoder and is used only as an ablation.
  - Sized to the pipeline rate: about 2M parameters at 5×10⁴ spectra/s and K=20, or 20–40M parameters at ingest-bound 3×10³ spectra/s. E2b fixes the size before training.
  - Training negatives are **incorrect targets**: hard, wrong-species and shifted-precursor. Decoy-recipe sequences are never used as labeled negatives.
  - Data are deduplicated by raw sha256 and spectrum hash.
  - Folds are project-disjoint, with the fold→project map committed.
  - It must win E6 **within each instrument class**, and pass G2 and G4, before it can become the default.

### L4: Statistics (see §6 for the pre-registered protocol)
- **Paired double competition:**
  1. All hypotheses of a spectrum compete.
  2. PSM-level target–decoy competition (TDC).
  3. Peptide-level (DDA) or precursor-level (DIA) competition between each target and its own paired decoy.

  The FDR estimate is (D+1)/T, made monotone. Ties are broken by hash(sequence, seed), which is permutation-invariant. Decoys are counted first at equal scores.
- **Strata membership** must be identical within a pair. Group-walk is used; pooled TDC is the fallback. Strata and merge rules are listed below.
- **Open and chimeric search: the CONGA pipeline.** Narrow plus open groups, a label-blind neighbour filter (drop any PSM sharing >5% of b/y ions with a higher-scoring PSM of the same spectrum), and stem-level competition.
- **RESET** whenever the score is learned on the data being evaluated (§6.4).
- **Strata (≤2 factors per mode; all other factors become features):**
  - closed: specificity × label-free chimeric rank;
  - open: CONGA stem/offset bin;
  - PTM: composition.

  The minimum size (200 pair-competitions, counted label-blind) is a power rule, not a validity rule. Pre-registered merges:
  - adjacent offset bins by mass;
  - rare PTM classes into "other";
  - chimeric rank ≥3 pooled.

  The partition is frozen before estimation and may use only label-blind counts or dev data.
- **Proteins:** Picked Protein Group FDR. The group score is the **best peptide score**, not a product of PEPs. Rescued subsets are validated by protein-level entrapment. Razor and Occam assignment are never used for FDR.
- **PTM:** peptidoform FDR (localization-agnostic), FLR (§6.5) and site-level FDR are three separate numbers.
- **XL-MS:**
  - TT/TD/DD with FDR̂ = (TD+1−DD)/TT ("+1" variant), clipped at 0. Always reported alongside the conservative (TD+1)/TT. A level is reported only if TT ≥ 50.
  - Merge to the reported level (CSM → residue pair → PPI) before computing FDR at that level.
  - Self vs heteromeric links are separated by FDR, not by Group-walk. A link involving decoy-A of protein A counts as self.
  - Every CSM must have fragment evidence for **both** peptides (min-of-two features).
  - Linear, mono-link and loop-link explanations compete in the same search.

### L5: Cross-run
- Design-aware, driven by an SDRF/design table.
- **MBR is opt-in.** Every transfer carries a **transfer q-value** (PIP-ECHO-style decoys: decoy peptides plus RT- and m/z-shifted decoy peaks) and a provenance flag. Transfers are never counted as identifications.
- **CF-MS:**
  - a run-specific detection q-value per fraction;
  - a provenance flag for library-borrowed values;
  - **profile-agnostic, fraction-local transfer FDR**. There is no profile-continuity prior: it would inflate exactly the co-elution correlation that CF-MS scores;
  - no global normalization.

  E10 measures false co-elution using spiked foreign-proteome profiles.
- **LiP-MS:**
  - never transfer between LiP and TrP runs;
  - half-tryptic calls exclude protein N-termini, signal/transit peptides and Met removal, and flag in-source fragments;
  - output peptide × run matrices with cut sites for MSstatsLiP. We do not reimplement its correction.
- **Quant:**
  - LFQ: MS1 feature intensities, then a Top-N / weighted-median MaxLFQ / directLFQ roll-up.
  - TMT: impurity correction, reporter tolerance per plex, SPS-MS3 linkage and isolation purity.
  - A QuantUMS-style variance model is an arm.
  - SILAC is out of scope for the MCD.

### L6: Repository scale
- The best evidence per peptide (or precursor) across runs is pooled into a global paired TDC, then picked protein groups, then global FLR.
- Each project is scored by the **fold model that never saw it**. The map is committed.
- **PTM sites at repository scale:**
  1. Best evidence per (protein, site) across runs.
  2. Global decoy-amino-acid site competition, gated by global peptidoform q.
  3. Report expected false sites.
  4. Flag single-evidence sites separately.
- **Amortized index:** one canonical superset index per (proteome, enzyme, mod superset). Each project's mod set is applied as a **mask** over the records' mod field, so projects share an index. Reuse is measured on the T4 sampling frame.
- Calibration and G4 are decided on **training or dev decoys only**. The pooling partition is frozen before estimation, so estimation decoys never choose the groups.
- **Outputs:** a master list carrying a global q, per-run lists carrying a run q, and expected false counts q×N plus FDP bounds.

### DIA (arm, after the MCD)
- Data order: Astral nDIA (ProteoBench DIA-Astral plus the Astral entrapment module) → QE variable-window DIA → diaPASEF.
- Tier 1: GPU spectrum-centric pre-search. Its recall is measured against DIA-NN IDs *before* tier 2 is built.
- Tier 2: peak-group extraction with an interference-aware fragment-selection operator.

## 5. Implementation strategy

- **Python** orchestrates and does statistics and ML (polars/pyarrow, PyTorch).
- **Sage (MIT, Rust)** supplies digestion, mods, the CPU index and the scoring oracle. It is reused through sagepy/PyO3, or as a forked crate extended with top-K plus feature export, pair ids and entrapment. Sage is not reimplemented in numba.
- **CUDA C++** (torch C++ extension or CuPy RawKernel) is used for the irregular kernels K0–K2. Triton or PyTorch is used for dense work: embeddings and the reranker. Numba is allowed only for glue and prototypes.
- Unmodified external baselines: Sage 0.14.7, Comet 2026.01, and DIA-NN 2.7 Academia (version recorded). MSFragger only if licensed.

## 6. Statistical protocol (pre-registered)

Frozen at git tag `prereg-v0.4` (created 2026-10-27; the earlier `prereg-v0.2` tag referenced here was never created. Since §6 text is unchanged since v0.3.1 and no sealed file has been opened — SEALED_LEDGER empty — this is the first binding freeze). Any change needs a dated note in RESEARCH_TREE.md written *before* the affected results are seen.

### 6.0 Amendment A1 (2026-10-27, dated before any re-run of the affected gates; tag `prereg-v0.4.1`)
This amendment follows the independent review in reviews/m3_independent_review.md. Each change makes the protocol stricter or adds a missing specification. None of them relaxes a verdict.

1. **Frozen T1 scorer.** Exactly one scorer goes into T1: `static_gbm_pred.pkl` (sha256 31bfe685…, full hash in results/2026-10-27_M3_scoring/models). Every other variant is exploratory. D0/D1 count as a *selection set*. Before T1, the default's dev numbers must be confirmed on runs that were not used for selection (the D8 stream from E9, its first 10 human runs by manifest order).
2. **Tail-G2 (all learned or predicted-feature scorers, and RESET).** The test uses clean null spectra: those whose best target or decoy passes the method's own 1% threshold. Pass requires all of:
   - n_tail ≥ 200;
   - a one-sided 95% upper bound on the target fraction ≤ 0.60;
   - a one-sided 95% upper bound on (the learned scorer's tail target fraction − the hyperscore pipeline's, on the same nulls) ≤ 0.10.

   With n_tail < 200 the gate is "not evaluable", and that blocks FDR claims for that scorer. The original PSM-level G2 is kept unchanged.
3. **Null construction.**
   - Species nulls are homolog-filtered: any spectrum whose matched peptide (target or decoy) is within Hamming distance ≤ 2 (I/L collapsed) of, or shares a 6-mer with, the host proteome or contaminants is removed. This applies to both labels.
   - Species-null retention times are quantile-mapped onto the host gradient before mixing.
   - At least 3 null projects are used.
   - Offset nulls (n+0.5 Da) are kept as specified.
   - G2, tail-G2 and G5 are reported with and without the RT features.
4. **Gates G1, G3 and G4 are required** (§6.4) before any FDR claim. G3 is run on null spectra binned by score.
5. **Arm M is required for predicted-feature scorers.** Their features come from a peptide-token model, so they have sequence capacity. On dev, this means D7 yeast (PXD031616) searched against yeast + human entrapment at r ≈ 1. T5 is used only at the milestone.
6. **CI convention for entrapment.** Monte-Carlo simulation (§6.7, ≥10⁴ replicates) chooses between (a) bootstrap variance plus between-seed variance and (b) bootstrap only. The rule fixed now: adopt whichever has 95% CI coverage closest to nominal without going below 93%. Until that simulation exists, both are reported and the more conservative one decides. "Consistent with control" is never used to rank scorers.
7. **Arm-C expected lower bound** α·r/(1+r) is reported as a diagnostic only and never changes a verdict.
8. **Re-ranking of the top-K by a static scorer** is allowed only if tail-G2, G1 and G3 pass.
9. **Comparators** include a predicted-feature rescorer on the comparator's own PSMs: Sage with its RT prediction on, plus the same AlphaPeptDeep features fed to mokapot on Sage PSMs. Gains are quoted against the strongest comparator.

### 6.0b Amendment A2 (2026-10-27; adopted before any DN2 file is listed; tag `prereg-v0.4.2`)
The normative text is `reviews/amendment_A2_draft.md` at the tagged commit, adopted **verbatim**. The implementation is commit ce53769 (scripts/a2, vulcan/gates/a2_*, config/a2.yaml), with its stricter-reading choices listed in scripts/a2/IMPLEMENTATION_NOTES.md. The PI resolved three open points before the tag:
1. **G5 evaluability on FS:** FS is evaluable for G5 only with ≥200 target winners **and** ≥200 decoy winners. This mirrors the A1 n_tail ≥ 200 convention and is fixed before any data. Below that, FS G5 is "not evaluable" and, per A2.4.3, does not block. On A and on A∪B, G5 keeps the A1 rule unchanged.
2. **A1 keyword list:** kept verbatim, without word boundaries (it over-excludes, which is the stricter reading).
3. **T4 exclusion of DN2:** enforced in scripts/data/select_t4.py (dn2.tsv). bench/sealed_run.py has no T4 assertion to extend.

### 6.0c Amendment A3 (2026-10-28; adopted before DN2 step 5 and before any T1 look; tag `prereg-v0.4.3` once A3.1 preconditions are met)
The normative text is `reviews/amendment_A3_draft.md` (revision 1, commit f549748), adopted **verbatim**, with its independent review `reviews/amendment_A3_review.md`. It adds v2_sep (gated, claimable) and ce_s (gated, licenses nothing) to the single DN2 step-5 run in two tiers, and adds v2_sep as a pre-registered **secondary** T1 arm with its own error family. The v1 primary (A1) and v1's A2 analysis are unchanged and enforced by the import-closure guard. Step 5 waits for all A3.1 preconditions and the tag; if v1's step 5 started earlier, A3.2 and A3.3 would lapse. Adoption disclosure (A3.1.1) is in RESEARCH_TREE D22.

### 6.0d Amendment A5 (2026-10-28; adopted before any T1/T5 look; it relaxes the claim conditions knowingly)
The normative text is `reviews/amendment_A5_draft.md` (revision 1, commit ae21c30), adopted **verbatim** with its independent review `reviews/amendment_A5_review.md`. **A5 relaxes the conditions for FDR claims**: the T1 look may proceed although v1 FAILED its only completed null-spectrum test (A1; the contamination explanation is unconfirmed) and A2 was not evaluable. A4 was drafted, reviewed twice, not adopted, and not run (D28). Every FDR-filtered number carries the A5 labels and the gate table, and the fixed choices of D28/A5 apply whatever the C9 pilot shows. The PI's choices on the A5 open points are in RESEARCH_TREE D29.

### 6.0e Amendment A6 (2026-10-28; adopted before any HN1 listing and before any T1 look)
The normative text is `reviews/amendment_A6_draft.md` (revision 1, commit dfb36e5), adopted **verbatim** with its independent review `reviews/amendment_A6_review.md`. HN1 is a single human removed-peptide null-spectrum test of tail-G2, G1 and G3 for **v1 only**, at α/2 under A5 S6. Its clock and tag-feasibility rule are part of the text; the status file at the fixed path is the only interface to the A5 look. The PI's choices on its open points are in RESEARCH_TREE D32.

### 6.1 Seeds
| Purpose | Seeds |
|---|---|
| General | `SEEDS = [11, 23, 47, 83, 131]` |
| Decoy shuffles | `DECOY_SEEDS = [1001…1005]` |
| Entrapment, dev | `ENTRAP_SEEDS_DEV = [2001, 2002, 2003]` |
| Entrapment, sealed | `ENTRAP_SEEDS_SEALED = [2004, 2005, 2006]` |
| Calibration decoys | `CALIB_DECOY_SEED = 3001` |
| Pilot subsample, project folds | 11 |

Replacement shuffles are hash(peptide, seed), so they do not depend on DB order. Every entrapment peptide must be unique against all targets and all decoys of all decoy seeds.

### 6.2 Decoy-variability arm
The default pseudo-reverse is deterministic. A real variability arm uses 5 C-term-fixed shuffled-decoy databases (DECOY_SEEDS). What gets replicated is always stated.

### 6.3 Entrapment and the decision rule
**Arms:**
- **(P) Paired shuffled.** FDRBench-style, r = 1, one entrapment per target, with its own decoys.

  FDP̂_P = (N_E + N_{E≥s>T} + 2·N_{E>T≥s}) / (N_T + N_E)

- **(C) Combined, natural foreign proteome.**
  - Human samples use *Arabidopsis*, after removing I/L-identical, single-substitution and isobaric (N↔GG, Q↔AG) neighbours.
  - r is computed on the effective search space (unique, I/L-collapsed peptides after the length/mass filters) and reported per mass decile.

  FDP̂_C = N_E(1+1/r)/(N_T+N_E)

- **(M) Memorization.** A high-resolution sample from an organism absent from training, searched with human peptides as entrapment, seeded-subsampled to r ≈ 1.
- **Lower bound** LB = N_E/(N_T+N_E). It is used only as a failure detector.
- **Required arms:**
  - Static low-capacity scorers (including predicted-feature linear/GBM): **P and C**.
  - Models with sequence capacity (cross-encoder, bi-encoder, any peptide-token model): **P, C and M**.
  - Fixed hand-crafted scores: P or C.

**Statistic.** The mean FDP̂ over replicates (entrapment seeds × decoy seeds), with a 95% CI from bootstrapping target–entrapment pairs combined with between-seed variance.

| Verdict | Condition |
|---|---|
| Supported | mean FDP̂ ≤ α **and** upper 95% CI ≤ 1.5α, for each required arm |
| Supported-noisy | mean ≤ α but the upper CI > 1.5α |
| Violated | the lower CI of the mean LB > α |
| Not evaluable | α·N_T < 20 |
| Inconclusive | anything else |

- **Expected behaviour of the upper-bound estimators.** The paired and combined estimators overstate FDP. In simulation, a correctly calibrated TDC list with true FDP 0.0104 at α = 1% has a paired mean of 0.0124, so it is "inconclusive". A well-calibrated method will therefore often be **inconclusive** rather than supported. In that case we report "consistent with control" if the lower CI of the paired FDP̂ is ≤ α. This is a weaker claim, and it is always labelled as such.
- **Gated levels:** α ∈ {1%, 5%}. At 0.1% and 0.5% (peptide level) results are **reported, not gated**.
- The same rule is used for the development regression gate.
- **Replicates via family-max (§4 L1).** In dev: 3 entrapment seeds with pseudo-reverse, plus *Arabidopsis*, in one pass. At milestones: the full 5 decoy × 3 entrapment.

**Units.** Entrapment is evaluated at peptide/precursor, peptidoform, protein-group and global levels.
- **Protein level:** paired shuffled entrapment *proteins*, with shared peptides swapped consistently. A group with ≥1 original target counts as a target.
- **Reporting:** FDP comes from the augmented search; sensitivity comes from the production search, which has no entrapment.

### 6.4 RESET and learned scores
**RESET (arm E4)** is required whenever a score is fitted on the data being evaluated. The default static-first scorer is not, so it needs no RESET. The exact estimator, including the +1 placement, is written in `vulcan/stats/reset.py` and verified in §6.7.

| Mode | Training/estimation split | FDR estimate | Training constant |
|---|---|---|---|
| Single-decoy | winning decoys, after peptide competition, split by a seeded half | FDR̂ = 2(D_est+1)/T | 3 |
| Two-decoy | decoys from two independent shuffles (exchangeable) | factor 1 | 2 |

- Estimating decoys enter training as unlabeled pseudo-targets.
- There are no CV folds. If folds are ever added, they are grouped by stem across the cohort.
- **Fallback:** if RESET training cannot start (too few pseudo-targets at α_train = 1%), the current score is kept and estimated with factor 2. This is valid. One retry with α_train = 5% is allowed. The fallback rate is always reported.
- **RESET plus Group-walk and RESET plus picked protein are unproven.** Until our simulations (§6.7) and entrapment verify an adapted version (c = 2/3, decoy counts ×2), RESET mode reports **only pooled peptide-level FDR**. Protein FDR comes from a fixed-score pass.

**Gates for learned or predicted components** (all must pass before any FDR claim):
- **G0, label blindness.** The label-swap invariance test in §3.
- **G1, diagnostic only.** The AUC of a sequence-only target-vs-decoy classifier is reported. **The gate:** on null spectra, the scorer's output must not correlate with that classifier's logit in the top score decile (|Spearman ρ| < 0.05, 99% CI).
- **G2, differential null test.**
  - **Pass condition:** the learned pipeline's null target-win rate may exceed the fixed-score pipeline's, on the *same* null competitions, by ≤ 0.01 overall and ≤ 0.02 in the fixed-score top decile. Absolute p is a diagnostic only.
  - **Strata:** instrument class only.
  - **Budget:** one-sided 95% upper bound; ≈ 120k null competitions per class for worst-case independence (`g2_sample_size` gives the exact n). MC must show ≥ 90% power at p = 0.5 and that the gate fires at a +0.02 excess.
  - **Null data:**
    - **Mixture design.** A real human run plus null-provenance spectra go through production, and win rates are computed on null-provenance spectra only.
    - **Null spectra.** Two sources:
      - precursor offsets of n + 0.5 ± 0.1 Da (non-physical mass defects), applied after recalibration;
      - wrong-species spectra against a human DB, after a homolog filter (I/L, single-substitution, isobaric) against the foreign proteome plus contaminants.
- **G3, retrieval symmetry on nulls.** For any pruning, ANN or lazy stage, the recall of the best target and the best decoy relative to exhaustive search, on null spectra, binned by score, must agree within 1%.
- **G4, calibration across strata.** Decoy score quantiles per (instrument, activation, charge) must be comparable, with the KS threshold fixed on dev data. G4 is decided on **training or dev decoys only**. The partition is frozen before any estimation.
- **G5, leak alarm.** On G2 nulls, the AUC for separating target winners from decoy winners must stay ≤ 0.55. Above that, the pipeline halts and the cause is investigated.

**Re-gating.** G2, G3 and dev entrapment are re-run whenever a retrieval, pruning or scoring parameter changes.

### 6.5 FLR
- FLR̂ = n_pA · (N_STY / N_A) / n_pSTY, where N_STY and N_A count candidate acceptor residues in the localized peptides. Ala is primary; Leu (counting I as L) is secondary. Mono- and multi-phospho are reported separately. Sites are also stratified by whether another candidate acceptor lies within 2 residues.
- **E8 (FLR̂ vs true FLR) is a gate** before any FLR claim.
- **Dataset roles:**
  - **Dev and calibration:** PXD007058 (Ferries), plus a dev subset of PXD000138 and PXD009449.
  - **Sealed for FLR claims:** held-out PXD009449 (ProteomeTools PTM) pools, whose per-pool peptide lists are published, plus a held-out PXD000138 subset.
  - PXD007058 lacks open ground truth, because its site list is in closed-access SI. If that SI is obtained it may be added as a second sealed set, but only by a dated note.
- Pools that contain positional isomers of one sequence are excluded.

### 6.6 Dev/sealed split (fixed now; the DATA.md `role` column must match this)
**Dev sets** may gate merges. Their numbers are never reported as claims, only as development curves.

| Id | Data | Purpose |
|---|---|---|
| D0 | PXD001468, 2 HEK files | DDA dev |
| D1 | PXD001077 | Pfu smoke test |
| D2 | PXD000138 dev subset + PXD009449 dev pools + PXD007058 | PTM / FLR |
| D3 | PXD014337 | DSS, non-cleavable XL arm |
| D3b | one PXD029252 DSSO file **other than** T3's | cleavable XL dev |
| D4 | PXD016766 | TMT |
| D5 | ProteoBench DIA-plasma PYE | DIA arm |
| D6 | PXD059608 CF fractions | DIA-CF arm |
| D7 | a DDA LiP/TrP pair, or a DDA CF-MS set | chosen by the data agent; recorded |

**Sealed sets:**

| Id | Data | Purpose |
|---|---|---|
| T1 | ProteoBench DDA-QExactive | IDs, quant, entrapment |
| T2 | held-out PXD009449 pools + held-out PXD000138 subset | FLR |
| T3 | one PXD029252 DSSO + HEK file | cleavable XL |
| T4 | stratified PRIDE sample (below) | headline repository test |
| T5 | a modern high-resolution **yeast** DDA run | memorization arm M |

- **T4 sampling frame:**
  - DDA on Thermo or Bruker; LFQ, TMT or phospho; UniProt reference proteome;
  - projects not in the training sets of the predictors we use;
  - ≤ 8 seeded runs per project, ≤ 300 runs, ≤ 1 TB;
  - replication: 1 decoy × 1 entrapment seed, with full replication on 5 projects;
  - stream-and-delete via `/tmp`;
  - a metadata-only dry run of the selector comes first.
- **T5:** yeast is held out of all training, and the run comes from a lab not in the dev sets. Searched against yeast plus human entrapment.
- **Later:** T6 = ProteoBench DIA-Astral plus the Astral entrapment module; T7 = PXD049412; T8 = LiP truth (PXD015446 known targets); T9 = HLA-I.

**Sealing mechanics.**
1. Sealed files live in `data/sealed/` (`chmod 000`).
2. They are unlocked only by `bench/sealed_run.py`. The script requires a clean tree, a milestone tag and a ledger line written to `results/SEALED_LEDGER.md` before it runs.
3. The sealed entrapment seeds are 2004–2006.

**Predictor leakage check.** Before the prereg tag, every dev and sealed accession and the T4 frame are cross-checked against the training manifests of the third-party predictors (AlphaPeptDeep, MS2PIP, DeepLC, Koina models). Overlaps are excluded or the predictor is retrained. The results are recorded in DATA.md.

### 6.7 Estimator verification
Monte-Carlo simulations (≥10⁴ replicates; pass if mean − 2.326·SE ≤ α) for TDC, paired double competition, Group-walk, RESET and the picked group must give mean FDP ≤ α. A biased simulation (null target-win p = 0.55) must make G2 fire. Both live in `tests/stats/`.

### 6.8 Metrics
- **IDs:** peptides/precursors and protein groups at q ≤ {0.1, 1, 5}%, each with its entrapment verdict. PSM counts are labelled "not FDR-controlled".
- **Localization:** FLR at nominal levels; site recall at matched true FLR.
- **Quant:**
  - median CV with the number quantified and completeness;
  - ProteoBench ε (|log2 FC error|);
  - DE pAUC at empirical FDP.
- **Speed:** as listed under Reporting in §2.

## 7. Benchmarks
Datasets are listed in §6.6 and DATA.md.
- **Dev loop:** D0–D4, minutes per iteration.
- **At milestones:** sealed T1–T5, plus the scale test E9.
  - E9 runs on N ∈ {1, 10, 100, ~1000} **human** runs, so global paired entrapment has one proteome. Other species get per-species master lists.
  - Replication at N=1000 is reduced to 1 entrapment seed × 1 decoy seed (pre-registered).

## 8. Engineering and reproducibility
- **Code layout:** scripts in `scripts/`; the package in `vulcan/`.
- **Benchmarks:** `bench/` writes `results/<date>_<exp>/{metrics.json, config.yaml, env.txt, git_sha}`. Results are committed; raw data are not.
- **DATA.md** lists source, version, sha256, license, purpose and dev/sealed role for every dataset, FASTA and model weight. It includes the training accessions of third-party predictors.
- **Tests** cover masses, digest, decoys and pairs, label-swap invariance, the index against brute force, estimator simulations, and FLR/XL formulas against hand-computed cases.
- **Regression gate on dev:** a paired comparison over decoy shuffles. A change is rejected if IDs drop by more than max(1%, 2× the seed SD), **or** if the §6.3 verdict on dev worsens.
- **Metrics, labels, estimators and decoy recipes are never changed to improve results** (§6 preamble).

## 9. Milestones (re-sequenced)
- **M1, harness.**
  - Toolchain ✓.
  - Dev and sealed data; DATA.md.
  - E1 ingest.
  - Sage and Comet (+mokapot) baselines on D0/D1/D4 with paired entrapment and full speed accounting.
  - Estimator simulations.
- **M2, DDA core.**
  - Full-centroid store.
  - Pilot pass (E0).
  - MS1 features, charge/isotope, ion tables.
  - K0–K2 on GPU (E2).
  - Paired TDC and picked protein.
  - Target: parity with Sage at verified FDP on dev, and higher runs/hour.
- **M3, rescoring and global FDR.**
  - Static-first scorer with predicted features; RESET arm (E4).
  - L6 global FDR and the E9 scale test, brought forward.
  - Gates G0–G5.
  - First sealed look on T1 and T4 (subset).
- **M4, PTM/open.**
  - Mass-offset and open search (CONGA).
  - Localization and FLR (E8).
  - Chimeric search (E7).
  - Sealed T2.
- **M5, assays.**
  - Cleavable XL (T3).
  - DDA LiP.
  - CF provenance and transfer FDR (E10).
- **M6, quant.** TMT (D4) and LFQ on ProteoBench DDA. Can run in parallel with M4–M5.
- **M7, arms.** Cross-encoder (E6), bi-encoder (E5), clustering, Astral DIA.
- **M8:** full T4 run, scale test and REPORT.md.

## 10. Key experiments
| Exp | What it measures | Decision it drives |
|---|---|---|
| **E0** | Pilot pass on/off: IDs vs Sage at verified FDP, on dev | Whether the pilot pass is kept |
| **E1** | Ingest: core-s per 10⁵ MS2 by stage, for TRFP→mzML vs TRFP Parquet vs in-process vs timsrust; 1–4 concurrent files; Lustre vs `/tmp`; QE/Exploris/Astral/Lumos-IT; store bytes per spectrum | The lowest-core-s path. If it is above 30 core-s per 10⁵ MS2, the headline is ingest-bound and GPU work is reprioritized |
| **E2** | Search from the store, identical index and spectra: Sage (32 and 4 threads) vs K1+K2 (GPU) vs our CPU path (E2c), per mode. Wall, core-s, GPU-s, H2D, postings per spectrum. Top-K recall vs brute force for targets and decoys (\|Δ\| ≤ 1%) | A mode moves to the GPU if recall ≥ 99.5% and core-s are freed |
| **E2b** | Reranker throughput vs parameter count (K = 10/20/50, FP16/INT8); predictor peptides/s; cache hit rate | Fixes the model size |
| **E3** | Per-run vs micro-cohort vs global sort, at index sizes 1.6, 30 and 60 GB: runs/hour including barriers, time to first result, peak RAM | Whether to sort globally |
| **E4** | Value of predicted features: ours vs Sage+LDA vs Sage+mokapot at verified FDP | Whether predicted features are kept |
| **E5** | Bi-encoder recall@K vs K1, plus G3 on nulls. **Targeted where indices are huge:** open, non-specific and HLA search (fragment indices of 30–86 GB), not closed tryptic search | Go/no-go for bi-encoder retrieval |
| **E6** | Frozen cross-encoder vs RESET rescoring, per instrument class | Go/no-go for the cross-encoder (needs G2/G4) |
| **E7** | Chimeric gain, rank-stratified, with CONGA filtering | Whether chimeric search stays default |
| **E8** | FLR̂ vs true FLR | Gate for any FLR claim |
| **E9** | Global FDP vs N: naive merge vs our global picked procedure | Validates L6 |
| **E10** | CF-MS transfer modes. Measures: PR-AUC of CORUM intra- vs inter-complex pairs, elution-peak width inflation, entrapment in fractions | The CF-MS transfer policy |

## 11. Risks
- **Ingest-bound throughput.** E1 first; ingest is reported separately.
- **Learned-scorer FDR leakage.** G0–G5, three entrapment arms, RESET.
- **Scope.** The MCD comes first, and arms are gated.
- **A single GPU.** On the A100 (80 GB), closed, phospho and semi-specific indexes fit resident; open/non-specific/HLA still stream. The allocation can change between sessions (§2 portability rule), so kernels must not assume 80 GB.
- **Licensing.** No closed code is embedded.

## Change log v0.1 → v0.2
- **HPC review:**
  - B1: the objective is now CPU core-seconds per run, and all matching moves to the GPU.
  - B2: SpMM and the dense fragment matrix are replaced by a tiled inverted index with on-the-fly fragments.
  - B3: ingest is now direct, with full centroids and views.
  - M1: Sage plus CUDA C++ replaces numba.
  - M2: micro-cohorts replace the global sort.
  - M3: no mmap on Lustre; staged reads.
  - M4: schema completed.
  - M5/M6: reranker and predictor budgets added.
- **Domain review:**
  - B1: MS1 and metadata are now stored.
  - B2: pilot pass and E0 added.
  - B3/stats B4: dev/sealed split.
  - B4: the MCD plus gated arms.
  - M1–M3: hypotheses, MS1 features, ion tables, preprocessing views.
  - M5: RESET is the default.
  - M6: DIA order is Astral first.
  - M8/M9: exhaustive localization and FLR roles.
  - M11–M13: XL, LiP and CF fixes.
  - M14/M15: quant specified, T4 PRIDE sample, TMT and HLA added.
  - M16: fair speed reporting.
- **Stats review:**
  - B1: the Δ-to-decoy feature is removed and label blindness is enforced.
  - B2: three entrapment arms are required for learned scorers.
  - B3: a CI-based decision rule and a real decoy-variability arm.
  - B5: G1 demoted to a diagnostic; G2 is now an equivalence test; G4 and G5 added.
  - M1: pair-invariant strata.
  - M2/M3: CONGA.
  - M4/M5: pair closure and label-blind routing.
  - M6: entrapment decoys and protein-level entrapment.
  - M7: RESET constants and restrictions.
  - M8: negatives and deduplication.
  - M9: residue-blind localization and an FLR formula.
  - M10: the picked-group score.
  - M11: non-inferiority definition.
  - M12–M14: repository calibration, XL counts, hard regimes.

## Change log v0.2 → v0.3 (round-2 review)
- **B1:** G2 is now a differential test per instrument class (≈110k null competitions per class).
- **B2:** the null data are non-physical offsets plus homolog-filtered wrong-species spectra, in a mixture design.
- **B3:** static-first scoring; RESET is the E4 arm.
- **B4:** a role column and sealing mechanics; T5 named (yeast).
- **B5:** FLR truth roles swapped (sealed = held-out PXD009449/PXD000138).
- **M1:** G4 is decided on training or dev decoys only.
- **M2:** ≤2 Group-walk factors per mode, with pre-registered merges.
- **M3:** the verdict table is fixed (upper CI ≤ 1.5α); gating at 1% and 5% only; non-inferiority margin 0.97 on shared decoy FASTAs.
- **M4:** the label-swap test also permutes storage order, with hash ids and hash tie-breaks.
- **M5:** dev DSSO and D7 (DDA LiP or DDA CF) added.
- **M6:** no CF continuity prior.
- **M7:** a predictor leakage check.
- **M8:** a true split-then-joint cross-encoder.
- **Simplifications:** family-max single-pass entrapment; T4 capped; E9 on human runs only.
- **Brief gaps closed:** a canonical superset index with mod masks; repository PTM-site FDR; E5 targeted at huge indices.
