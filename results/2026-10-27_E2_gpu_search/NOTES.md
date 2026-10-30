# E2: GPU closed search (K0+K1+K2) vs Sage 0.14.7 (dev data only)

These are dev numbers (D0 = PXD001468, 2 HEK files; D1 = PXD001077, 1 Pfu file). They are development curves, not claims (DESIGN §6.6). No sealed data was touched.

## Hardware
- Node gpu-node.
- CPU: AMD EPYC 9354, 32 affinity cores (2 sockets, AVX-512).
- GPU: A100 80GB PCIe (sm_80), driver 580.178.04, nvcc 12.9.
- The full fingerprint is in `metrics.json:hw`; every `timings.jsonl` record also carries it.
- **The node was shared during these runs.** Host loadavg was 14–22 from other agents and users (MS1/biosaur2, TRFP, mokapot, FDRBench). Timed runs were serialized behind `flock /tmp/vulcan_bench.lock`, but unlocked jobs (Java, Lustre I/O) still ran alongside.

## What was built
- **L1 index** (`vulcan/index/`):
  - Digestion and variable-mod enumeration use sage-core 99407db6 through the `vulcan_core` PyO3 crate. The crate now also returns Sage's own f32 masses and hard-checks the permutations.
  - Decoys: pseudo-reverse with the C-term fixed; on a collision, a hash(peptide, 1001, k) shuffle. Entrapment E_s: a hash(target, s, k) C-term-fixed shuffle, with the M1e uniqueness checks.
  - Both recipes are the `scripts/entrapment/peptides.py` ones. `tests/index` checks they are identical.
  - Ids are xxhash64 of the canonical peptidoform. Rows are sorted by (Sage f32 mass, id). Each row stores `partner`, its symmetric pair row.
  - Labels (family T/E/dec(T)/dec(E), tdc_pair, entrap_pair) live only in `labels.arrow`, for L4.
  - Build time: human+contaminants with E2001 is 14.1M peptidoforms in 2 min (CPU); prod (T+dec) is 7.06M in 40 s.
- **K0** (torch, GPU): range filter, global top-150 (Sage) or top-k per 100 Th window, mass sort, deterministic TIC.
- **K1** (CUDA, `vulcan/gpu/kernels.cu`):
  - Layout: CSR over (chunk of 8192 rows × 20-ppm log bin). Each posting is a u32 (u16 local id | u16 fine log-ppm offset), 4 B/posting.
  - For D0 prod: 175,429,652 postings, which equals Sage's fragment count. The index uses 1.5 GB on the GPU and builds in 0.53 s.
  - Kernel: one warp per (hypothesis, chunk tile). Shared-memory u32 counters; a binary search for the local-id range inside each posting list.
  - Then a threshold append (count ≥ 1, the same for every candidate), top-K per spectrum (ties go to the higher row, never to a label), and pair closure.
  - Slicing path: rows are cut at chunk boundaries into slices of < 2³¹ postings, each with its own CSR. All slices are resident today; streaming slices from pinned host memory is left for open search.
- **K2** (CUDA): one block per spectrum, peaks in shared memory. Sage's `score_candidate` is replicated in f32 with `-fmad=false`.
  - Features: hyperscore, matched b/y, summed intensities, intensity-weighted ppm, longest b/y run, delta_next and delta_best regardless of label, Poisson, delta mass, and others.
- **Built for** sm_80 + sm_86 + PTX by `scripts/m2/build_kernels.sh`. `vulcan/gpu/_lib.py` rebuilds automatically when the source hash changes.

## Correctness (`correctness.json`, `topk_sweep.json`)
- **A. G0 invariance.** We shuffled the FASTA entry order, permuted the pair list, and swapped the target/decoy roles of every pair (mods enumerated on the decoy side). All 9 index-side arrays came out bit-identical. The D1 PSM table (top-5, 33,088 rows) is also bit-identical. **PASS.**
  - A nondeterministic float-atomic TIC was found and fixed along the way.
- **B. K1 index vs plain-f32 brute force** (3,000 random spectra each; D1 and D0, E2001 index): the candidate sets, the counts, and the per-hypothesis Sage statistics (scored_candidates, matched peaks) are identical. **PASS.**
  - Bug found and fixed: the quantized lower bound used `ceil`, which dropped exact f32 ties at the tolerance edge (Sage counts them).
  - Quantization is exact by construction, because a quantum (3e-10 relative) is much smaller than the f32 spacing.
- **C. Recall vs exhaustive K2** (every candidate in the precursor windows), K = 500:

  | Set | Target recall | Decoy recall | \|Δ\| |
  |---|---|---|---|
  | D0 prod | 99.90% | 99.78% | 0.12% |
  | D0 E2001 | 99.85% | 99.62% | 0.23% |
  | D1 | 99.99% | 99.91% | <0.1% |

  **PASS** (≥ 99.5% and ≤ 1%).
  - At Sage's K = 50 the test **fails** on D0: 96.8/93.3% (Δ 3.5%), and 95.2/90.0% with entrapment. A K-limited shortlist loses decoys first, so we default to K = 500. The extra cost is +7% of the GPU ms.
  - Recall plateaus at about 99.9% for any K. The residual misses are candidates with K1 count 0 (matched only by b1/b2/y1/y2, which are not indexed). Sage behaves the same way.
- **D. Sage oracle** (same spectra, index and parameters):
  - In K = 50 mode, the top-1 peptidoform is identical for 99.83% of spectra (D1) and 99.23% / 99.38% (D0 files).
  - When the peptide agrees, the hyperscore differs by at most 7e-15. matched_peaks, longest_b/y, scored_candidates and poisson are identical.
  - Where the peptide differs, the cases split into:
    - exact hyperscore ties: 18 / 209 / 177;
    - ours higher: 1 / 56 / 51;
    - Sage higher: 0 / 6 / 6 (tie-breaking at the K = 50 boundary).
  - With K = 500, Sage is never higher, and ours is higher on about 2% of D0 spectra (mostly low-scoring ones).
  - delta_next can differ because pair closure adds candidates that Sage never scores. This is by design.

## E2 speed (`metrics.json:speed`; median of 3 locked reps)
Index: prod (T + dec(T)). Our spectra come from the Arrow store; Sage reads an MGF written from the same store.

| | D0 (88,636 MS2) wall / core-s | D1 (15,374 MS2) wall / core-s |
|---|---|---|
| **Ours, one process** (import + CUDA init + index load + GPU index build + runs) | 7.5 s / 4.9 | 7.4 s / 4.2 |
| of which Python import | 4.2 s | 4.1 s |
| of which runs, warm (store load → PSM parquet) | 0.54 s | 0.065 s |
| of which search only (K0..K2 + features) | 0.24 s | 0.033 s |
| Sage 32 threads | 29.9 s / 237 | 4.5 s / 63 |
| Sage 32 threads, search phase only (its log) | 2.98 s | 0.146 s |
| Sage 4 threads | 46.5 s / 137 | 3.9 s / 10.1 |

- **D0, ours vs Sage-32:** 4.0× faster wall, 48× fewer core-s, 12.5× faster on the search phase alone (2.98 s vs 0.24 s warm).
- **D0, ours vs Sage-4:** 6.2× faster wall, 28× fewer core-s.
- **D1 (one small file): Sage wins on wall time.** Our process cost is dominated by the 4 s Python/torch import from the Lustre-hosted conda env, plus 1 s to load the index. We still use 15× fewer core-s than Sage-32.
- **D0 data volumes:**
  - H2D: 145 MB of peaks per dataset, plus 356 MB for the index upload.
  - GPU peak memory: 9 GB.
  - Per D0 file, GPU-event ms are: K0 16, hypotheses 56 (CPU numpy), K1 12, top-K 11, K2 6, features 17.
- **Pipelined multi-run** (`vulcan/gpu/batch.py`; D0 stores ×5 = 10 runs, 443k spectra, one process):
  - Steady state: 1.8–2.1 s, i.e. 215–250k spectra/s end-to-end from the store, about 0.12 s per 44k-spectrum run, and 0.75 core-s per 10⁵ MS2.
  - Two locked runs with `--io-threads 1` (`D0x5_io1*`) stalled for about 190 s in setup and about 6 s in their first search. During the stall the process was observed in `cl_sync_io_wait` (Lustre). An unlocked rerun right afterwards was normal (6.2 s setup, 2.07 s steady). We treat those two records as environmental and exclude them.
- **Hybrid CPU (E2c): not built.** A K1+K2 run takes about 20 ms of GPU per 44k-spectrum file. The per-run critical path is already on the CPU side: Arrow read and flatten, numpy hypotheses, and feature assembly. A CPU matching tier would not help closed search. What the 32 cores should do instead is I/O and pre/post-processing, overlapped with the GPU (the pipelining above). Sage-32 stands in as the CPU-tier reference.

## First end-to-end IDs (`ids.json`; dev only)
Peptidoform-level paired double competition (`vulcan.stats.tdc.double_competition`, hash tie-break seed 11), q ≤ 1%. The score is the raw hyperscore, with no rescoring.

| | Ours (hyperscore) | Sage PSMs, same estimator | Sage's own peptide_q (LDA) |
|---|---|---|---|
| D0 | **19,601** | 19,605 | 23,637 |
| D1 | **5,758** | 5,758 | 6,037 |

- The same score through the same estimator gives the same IDs. The roughly 17% gap on D0 to Sage's own number is its LDA rescoring; matching it is M3.
- **Entrapment diagnostic** (E2001 index, one entrapment seed, one decoy seed, no CI, peptide-sequence pairing). The paired FDP̂ at 1% is 0.94% on D0 (LB 0.47%) and 1.09% on D1 (LB 0.57%). There are too few replicates for a verdict.

## Caveats
- Sage's `finished in` time (17 s on D0-32) is well below its process wall (29.9 s). We report the process wall, measured by GNU time.
- Our "GPU-s" are CUDA-event-bracketed phase times. They include any CPU work inside a phase (hypotheses are numpy), so they are an upper bound on GPU busy time.
- Speed comparisons are valid only within this node and session (DESIGN §2).

## Next optimization targets
1. **Process startup.** The 4 s import comes from the Lustre-hosted env. Options: a persistent worker or service, or staging the env and index on local NVMe. The L1 index load (1 s) should move to a zero-copy IPC mmap on /tmp.
2. **Move work to the GPU.** Hypotheses (56 ms numpy) and feature assembly (17 ms) to GPU or fused kernels, and K0 to a segmented top-k kernel. Together these could take the 0.12 s per run to about 0.05 s.
3. **Wire format.** u32 bin + f16 intensity (DESIGN L0) would halve the 145 MB H2D. Use pinned, double-buffered transfers.
4. **Family top-2 retention in K1/K2** for multi-seed entrapment from one search, plus the streaming slice path for open, non-specific and HLA search.
5. **M3 rescoring** to close the ID gap to Sage+LDA, and a peak-view ablation (global top-150 vs top-k per 100 Th).
