# HPC/GPU review of DESIGN.md v0.1

Reviewer: HPC/GPU systems engineer. Inputs: PROMPT.md, DESIGN.md v0.1, lit/01, with lit/02–03 grepped. I read no other reviews. 2026-10-26.

## 0. Verdict

The statistics are strong. The systems design rests on three wrong premises:
1. It treats the GPU as the scarce resource. On this node the **4 CPU cores** are the scarce resource.
2. It proposes a **dense fragment matrix / SpMM** for candidate generation. That does 10²–10³× more work than an inverted-index join, and it doesn't fit in memory for PTM spaces.
3. It puts **numba re-implementation** effort into closed search, which on this hardware is neither the bottleneck nor unsolved: Sage exists and is MIT.

The real end-to-end bottleneck is **Thermo decoding on the CPU**. Second is **predicted-spectrum generation on the GPU** in open/PTM mode. Fragment matching is third or lower.

## 1. Facts about the node that the design omits (checked today)

- Affinity is cores 112–115 (no SMT, NUMA node 7, likely one Zen 2 CCX: 16 MB L3). **There is no AVX-512**, contrary to lit/01 R3. `os.cpu_count()` returns **128**, so unpinned numba, torch or OpenMP pools oversubscribe 32×.
- The GPU is NUMA-local, **PCIe Gen4 x16**: ~24 GB/s pinned, ~10 GB/s pageable.
- There are 8 NUMA nodes, so "501 GB" is ~62 GB local; larger CPU indices spill to remote memory.
- 1.7 TB of local disk at `/tmp`; stage there.
- The relevant GPU peaks are the **tensor** ones: FP16 (FP32 acc) 71 TFLOPS, INT8 284 TOPS (dense). "35 TFLOPS fp32" is irrelevant here.

## 2. Back-of-envelope numbers

**Assumptions.**
- Human Swiss-Prot, trypsin, ≤2 missed cleavages, 7–50 aa, Met-ox ≤2: ≈3.3×10⁶ target peptidoforms, ≈7×10⁶ with decoys.
- ~28 indexed b/y ions each (z=1, 150–2000 m/z), 8 B per entry (f32 m/z + u32 id).
- A run is 10⁵ MS2 (2 h Orbitrap), 150 peaks kept, 20 ppm (≈0.02 Th) fragment tolerance.

From these, the fragment density is ≈2×10⁸ / 1850 Th ≈ 10⁵ per Th. So **each peak's tolerance window holds about 2,000 postings across the whole DB**, scaled by the fraction of the DB inside the precursor window.

### 2.1 Index sizes
| Space | Peptidoforms (T+D) | Entries | 8 B/entry | 4 B/entry (u16 local id + u16 m/z offset) | Fits 24 GB? |
|---|---|---|---|---|---|
| Tryptic + Mox | 7×10⁶ | 2×10⁸ | 1.6 GB (MSFragger 1.6–2.9) | 0.8 GB | yes |
| + entrapment r=1 | 1.4×10⁷ | 4×10⁸ | 3.2 GB | 1.6 GB | yes. **Entrapment doubles every search cost**; report speed with and without it |
| Semi-tryptic | ~1.5×10⁸ | 4×10⁹ | ~32 GB (MSFragger 55.8 with mods) | 16 GB | only compressed, or sliced |
| Phospho STY≤3 + Mox | 2–3×10⁸ | 6–8×10⁹ | 50–65 GB (MSFragger 86.6) | 25–33 GB | no: slice and stream |
| Open ±500 Da | = tryptic | — | 1.6 GB | — | yes. The cost is on the query side |

The **dense theoretical-fragment matrix** (L1 "GPU form") for phospho would be 2.5×10⁸ × ~100 slots × 4 B ≈ **100 GB**. It must never be materialized. Store each peptidoform as (protein offset, length, mod mask) in ~12 B, which is 3 GB, and **generate fragments on the fly** from residue-mass prefix sums. GPU arithmetic is free; GPU memory is not.

### 2.2 Candidate generation per spectrum
Postings touched ≈ 150 peaks × 2,000 × (fraction of DB in the precursor window).

| Regime | Candidates/spectrum | Postings/spectrum | CPU, ideal (index probe only) | CPU, Sage-class all-in | GPU tiled index (§4 K1) |
|---|---|---|---|---|---|
| Closed narrow | ~4×10³ | ~200 (cost is the 150 probes) | 5–20 µs | **1–3 ms per vCPU** (lit/01: ~10⁴ spectra/s per 32 vCPU) | ~10⁶/s |
| Closed + full-window chimeric (1.6 Th × z 2–4 ≈ 14 Da) | ~7×10⁴ | ~3×10³ | 20–60 µs | ~2–5 ms | ~10⁶/s |
| Semi / phospho (narrow) | ~10⁵ | ~4×10³ | 30–100 µs, plus NUMA/TLB misses on a 30–65 GB index | — | 3×10⁵/s if resident |
| Open ±500 Da | ~3×10⁶ (30% of DB) | ~9×10⁴ scatter-adds into a 3M-entry accumulator (6–12 MB, which spills past L2) | 0.5–2 ms | Sage ≈10²–10³ spectra/s per 32 vCPU [derived] | 10⁵–5×10⁵/s |

Consequences:
- **Closed search is not free on 4 cores.** At Sage's all-in constant a 10⁵-spectrum run costs 100–300 core-s (chimeric up to ~500), the same order as ingest (§2.5). The gap to the ideal 5–20 µs is parsing, preprocessing, charge/isotope hypotheses and full scoring, not the index probe. E2 must decompose it.
- **Open search** takes 4–40 min per run on the CPU and ~1 s on the GPU. This is the clearest GPU win.

### 2.3 Dense/SpMM alternatives for open search (why I reject them)
| Formulation | Work per spectrum | Throughput on 3090 |
|---|---|---|
| Inverted-index join (K1) | ~9×10⁴ postings | 10⁵–5×10⁵/s |
| Dilated bitset (log-ppm bins, 25 KB in smem) × on-the-fly fragments | 3×10⁶ × 28 = 8.4×10⁷ smem lookups (~10¹²/s) | ~1.2×10⁴/s |
| cuSPARSE SpMM (sparse peptide × dense binned spectra) | same 8.4×10⁷ gathers at ~10¹¹/s | ~10³/s |
| INT8 tensor-core GEMM, 1-Th mass-defect bins (D=2,000) | 1.2×10¹⁰ ops at ~150 TOPS effective | ~1.2×10⁴/s, and only a low-resolution prefilter |
| Tensor-core GEMM at high resolution (D≈10⁵) | 6×10¹¹ ops | ~250/s. Dead |

Tensor cores are the wrong tool for high-resolution matching, because specificity needs D≈10⁵. Use them for **embeddings** and the **reranker**.

### 2.4 PCIe
Per run: spectra ~120 MB (~5 ms), top-50 results 40 MB back, the tryptic index 70 ms once, a 60 GB phospho index 2.5 s per pass. PCIe matters in one failure mode only: re-streaming a >24 GB index per small batch (60 GB per 10⁴ spectra caps throughput at ~4×10³ spectra/s). Sorting hypotheses by mass within a run or micro-cohort fixes that; no global sort is needed.

### 2.5 Ingest on 4 cores (estimates; E1 must replace these)
Ingest is TRFP RAW→mzML (zlib/base64) at 1–3 core-min per 120k-scan file, then a Python mzML parse at 2–5k spectra/s/core. That is **≈1.5–4 core-min per run**, so 4 cores sustain ≈1–3 runs/min, or **≈2–5×10³ MS2 spectra/s for the whole node**. A direct RAW→Arrow path (TRFP Parquet, or RawFileReader in-process reading only the centroid stream, no XML) plausibly costs 3–5× less [unverified]. At PRIDE scale the download (50–200 MB/s, so 10–40 s per 2 GB file) is in the same league.

### 2.6 Neural reranker on one 3090
- Budget: ~20 TFLOPS effective FP16 on short-sequence transformers.
- At **50k spectra/s with K=20**, that is 10⁶ pairs/s, or 4×10⁸ FLOP per spectrum (2×10⁷ per pair). A transformer costs ≈2·params·tokens.
- Spectrum encoder: half the budget over 64 peak tokens ≈ **1.5M parameters**.
- Joint head: 20 peptide tokens per pair, with cross-attention to the spectrum ≈ **0.25M parameters** (d=128, 1–2 layers).
- Total **~2M parameters**, about AlphaPeptDeep-MS2 size. A Casanovo-size encoder (~28M parameters, 150 tokens) costs ~8 GFLOP per spectrum, which is ~2.4×10³ spectra/s for the encoder alone.
- A **late-interaction (MaxSim) head** costs 20×64×128 ≈ 1.6×10⁵ MAC per pair, ~60× cheaper than cross-attention. That frees the budget for a ~10M-parameter encoder at K=50.
- Size the scorer to the binding rate, though. From raw files that rate is ingest, ~2–5×10³ spectra/s, which affords a **20–40M-parameter** model. 50k/s matters only for reprocessing from the store.
- **Predicted intensities** (AlphaPeptDeep, ~4×10⁴ peptides/s): open/PTM runs with ~10⁶ unique candidate peptidoforms cost **~25 GPU-s per run**, 10–25× the search kernel. Tryptic: precompute the whole library once (7.9M precursors in ~10 min) instead of predicting lazily.

### 2.7 End-to-end per run (10⁵ MS2, steady state, pipelined)
| Stage | Resource | Estimate |
|---|---|---|
| Download | network | 10–40 s |
| Ingest, TRFP→mzML→Arrow | CPU | 90–240 core-s (**30–80 s wall** on 3 free cores) |
| Closed + chimeric search, Sage-class CPU | CPU | 100–500 core-s |
| Same on GPU (K1+K2) | GPU | <1 s |
| Open search, CPU / GPU | CPU / GPU | 300–3,000 core-s / 1–3 s |
| Predictions: tryptic cached / open lazy | GPU | ~0 / ~25 s |
| Reranker, 2M params | GPU | ~2 s |

**Bottleneck:** CPU decode, then CPU search if it stays on the CPU, then GPU predictions for PTM/open. The GPU sits >90% idle in closed mode unless it takes work off the CPU.

## 3. Issues

### Blockers
**B1. The CPU/GPU split criterion is wrong.** "GPU only where it wins" (L2, E2) compares kernel wall time. On a 4-core node the objective is **CPU core-seconds per run**, because the cores must decode RAW, feed the GPU and run Python. Closed + chimeric search on the CPU costs 100–500 core-s per run (§2.7), as much as ingest or more. Moving it to even an "inefficient" GPU kernel roughly **doubles node throughput**.
**Fix:** maximise runs/hour with the CPU as the binding resource. Give 3 cores to decode and 1 to orchestration. Run all matching, preprocessing and scoring on the GPU. Keep Sage's CPU index as oracle and fallback.

**B2. Dense fragment matrix / SpMM as the GPU search form (L1/L2).** It does 10²–10³× more work than the index join (§2.3), and phospho needs ~100 GB.
**Fix:** a GPU tiled inverted index (K1), fragments generated on the fly for rescoring (K2), and dense/tensor-core work only for embeddings, top-K similarity and the reranker.

**B3. Ingest double-serializes and freezes parameters into the store.** TRFP→mzML→converter encodes and then decodes XML, base64 and zlib on the binding resource. Top-N=150, precursor removal and deisotoping "once at ingest" mean any later change re-ingests the most expensive artifact at repository scale.
**Fix:**
- Store **all centroids** (or top-N ≥ 400 plus a rank column) and apply preprocessing as a GPU view.
- Go direct RAW→Arrow: TRFP Parquet, or in-process RawFileReader (mzPeak's Rust+.NET path; verify it), reading only the FTMS centroid stream.
- Run E1 before writing any engine code.

### Major
**M1. The language strategy wastes the first months.** E2's "numba index vs Sage" answers a question that doesn't matter (§2.2). "Port to Rust later" rarely happens and leaves two diverging code bases.
**Fix:**
- Use Sage (MIT, already built in `tools/sage`) as a Rust crate from day one: digest, decoys, CPU index, oracle. Extend it in Rust with top-K+features export, pair ids, entrapment and full-window candidates.
- Python for orchestration, statistics and ML.
- **CUDA C++** (torch extension or CuPy RawKernel) for irregular GPU kernels. Triton is poor at variable-length postings, smem atomics and warp top-K; it is fine for the reranker.
- Pin numba, OMP, Rayon and torch thread counts to ≤4.

**M2. Cohort "sort once, sweep once" is oversold.**
- It pays only when the index exceeds device memory (semi, phospho), and then it saves ~2.5 s of PCIe per run against 30–80 s of ingest. It buys nothing for a resident tryptic index, and nothing for open search, where every slice meets every spectrum.
- It costs a global ingest barrier: no ingest/search overlap, no incremental runs, no per-run failure isolation.
- The query unit is the **precursor hypothesis**, not the spectrum: z × isotope × co-isolated features gives 3–15 rows per spectrum, each carrying its run-calibrated tolerance.

**Fix:** use RAM-sized **micro-cohorts** (50–200 runs) with a sorted hypothesis table and per-run outputs. Go global only if E3 shows >20% end-to-end gain including the barrier.

**M3. The zero-copy claims are overstated.**
- Zero-copy holds only for **uncompressed** Arrow IPC from the page cache or local disk into CPU memory; LZ4/zstd and Parquet need decoding. The GPU always needs an H2D copy.
- **mmap on Lustre is pathological** (page-fault-sized RPCs, lock contention).

**Fix:**
- Stage to local `/tmp` or read sequentially into RAM.
- Define a GPU wire format: fixed-stride (u32 log-ppm bin, f16 intensity) blocks in pinned, double-buffered memory.
- Use `LargeList` or chunking: 1.5×10¹⁰ peaks is past Arrow's 2³¹ offset limit.

**M4. Schema gaps** for chimeric search, MBR and assays: an MS1 feature table and a spectrum→co-isolated-precursor table (DDA+); analyzer (IT/FT), activation, FAIMS CV, injection time, isolation offsets; a TMT reporter block kept outside top-N; a u64 global spectrum id; the timsTOF frame/scan map.

**M5. Reranker sizing and the "reuse" claim.** No FLOP budget is stated, and §2.6 shows how tight it is: ~2M parameters at 50k/s. Caching peptide encodings (≈5 KB each) for 10⁷–10⁸ peptidoforms takes **50–500 GB**, and fetching them costs about as much as recomputing.
**Fix:** use pooled or late-interaction encodings, and tie the FLOP/spectrum budget to the pipeline rate.

**M6. Predicted-feature cost is unbudgeted for open/PTM** (~25 GPU-s per run). Precompute the tryptic library. For open/PTM, cap unique peptidoforms per run and measure the cross-run cache hit rate.

**M7. GPU residency for semi/phospho.** These indices are 30–65 GB, so they need 4 B/entry encoding plus slicing. On the CPU, use huge pages and NUMA-local allocation, because the index spills past node 7's ~62 GB.

### Minor
- m1. The top-K selection over 3×10⁶ candidates must be fused into the kernel (threshold append, then a warp top-K). Writing scores out costs 12 MB per spectrum.
- m2. Use log-ppm bins (bin = ⌊ln(m/z)/ln(1+tol)⌋) so the tolerance is uniform in ppm, and dilate spectrum peaks by ±1 bin instead of doing two-sided searches.
- m3. f32 m/z gives 0.06 ppm, so the claim stands; intensity can be f16.
- m4. Report from-raw and from-store throughput separately, with thread pinning and entrapment setting.

## 4. Kernels I would build

- **K0: GPU preprocessing.** Top-N, precursor removal, deisotoping and log-ppm binning over the stored centroids.
- **K1: GPU tiled inverted index** for chimeric, semi, phospho, open and DIA tier-1 search.
  - Layout: CSR over (log-ppm fragment bin × precursor chunk of 8,192 peptidoforms), with u16 local ids.
  - Grid: chunk-major, so a chunk's ~0.9 MB of postings stays L2-resident while spectra stream past it.
  - **Warp per (hypothesis, chunk)**, with u8 count accumulators in smem (8 KB per warp) and smem atomics.
  - Candidates above k_min are appended to a buffer, with the same thresholds for targets and decoys (G3).
  - Chunk size trades smem against the peaks × chunks offset lookups.
- **K2: warp-per-spectrum exact rescoring and localization.**
  - Sorted peaks (or a dilated bitset) sit in smem.
  - For top-K candidates, fragments are generated from prefix sums, including mod masks, open-search Δmass (shifted and unshifted ions), positional isoforms and decoy-amino-acid sites.
  - It emits all Stage-A features and site evidence in one pass, at ~10⁴–10⁵ ops per spectrum.
- **K3: tensor-core work.** Embedding GEMMs with fused top-K (clustering, bi-encoder), the reranker, and top-K predicted-vs-observed similarity.
- **K4: Sage CPU index** as the oracle and fallback.

## 5. What E1–E3 must measure to be decisive

- **E1 (ingest).**
  - Measure core-s per 10⁵ MS2, broken down by stage (centroid vs profile reads, XML write/parse, Arrow write), and the store size per spectrum.
  - Compare TRFP→mzML, TRFP Parquet, in-process direct and timsrust, at 1–4 concurrent files, on Lustre vs `/tmp`.
  - Cover QE/Exploris, Astral and Lumos with ion-trap MS2.
  - **Decision:** take the lowest core-s path. If even that is above 30 core-s, the headline is ingest-bound, and GPU work should be re-prioritised to offload the CPU.
- **E2 (search).**
  - Run from the store, with an identical index and spectra, on Sage (4 pinned threads) and on K1+K2: closed, chimeric, semi, phospho, open.
  - Report wall, **CPU core-s**, GPU-s, H2D bytes, postings per spectrum and Sage's per-phase breakdown.
  - Report top-K recall against brute-force K2 on 2k spectra, separately for targets and decoys (|Δ| ≤ 1%).
  - **Decision:** move a regime to the GPU if recall is ≥99.5% and it frees core-s without making the GPU the bottleneck. Drop the numba index.
- **E3 (cohort).** Compare per-run sorted, micro-cohort and global sorting at index sizes of 1.6, 30 and 60 GB. Report runs/hour **including the ingest barrier**, time to first result, H2D bytes and peak RAM.
- **Add E2b:** reranker spectra/s against parameter count (K = 10/20/50, FP16/INT8), plus predictor peptides/s and the cache hit rate. This fixes the model size before any training.

## 6. What's good

- Batch operators over fixed-schema Arrow tables make CPU/GPU swaps testable; the content-hashed index is reused across runs.
- Precursor-mass partitioning under a memory budget is the right primitive (add compression).
- Retrieve-then-rerank, with the bi-encoder as a side channel behind G3, is the right order.
- Chimeric search by default, lazy cached predictions, ingest as its own line item.
- An honest speed protocol: spectra/s per device on the hardware of record, no extrapolation stated as measurement.
- Pre-registered statistics (paired decoys, entrapment decision rule, G1–G3). Apply the same rigour to systems claims: decide on measured core-seconds.
