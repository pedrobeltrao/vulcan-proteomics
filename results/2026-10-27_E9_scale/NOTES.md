# E9 scale test and L6 global FDR (tree 3d; dev-grade; dataset D8)

All numbers are **development curves on dev set D8, not claims** (DESIGN §6.6). No sealed file was touched, no T4-selection project was
used, and §6 / the estimators / the labels were not changed. Node gpu-node: 32 Zen 4 cores, A100 80GB, 1.5 TB RAM
(`metrics.json:hw`). The node was **shared and loaded** throughout: other agents' Sage, Comet and Java jobs, vLLM, loadavg 50-117, and
up to 40 GB of foreign GPU memory. Long E9 runs did not take the bench lock. Speed numbers are from this one session.

Files:
- `metrics.json`: everything below. It is built by `scripts/e9/summarize.py` from the analysis JSONs.
- `e9_n100.json`: N = 1/10/100, three score variants, protein level, A1.1 confirmation.
- `e9_n100_hs.json`: hyperscore reference.
- `e9_famdiag.json`: fair-coin diagnostic, N ≤ 100.
- `e9_n1000.json`: N = 1000, frozen score and hyperscore.
- `c_homology.json`, `download_probe*.json`, `config.yaml`, `git_sha`.

## What was built
- **L6, `vulcan/stats/global_fdr.py`.** The steps:
  1. Spectrum competition within each run.
  2. **Best evidence per peptidoform across runs**, kept per side.
  3. **Global paired competition**, then (D+1)/T, giving the **global q**.
  4. `run_competitions` gives every run's own double competition, so each run's list carries a **run q**.
  5. `naive_union` is the E9 comparator.
  6. `global_outputs` produces the **master list** (pf_id, best run, score, q_global), the **per-run lists** (run, pf_id, q_run, q_global)
     and the **expected false counts** (q×N and D+1 at the threshold). The outputs are in `/tmp/vulcan_e9/l6/N*/`, not committed.
     Example at N=1000 and 1%: 360,392 accepted, 3,601 decoys at the threshold, so about 3,602 expected false. Each run keeps on average
     9,471 globally accepted peptidoforms.
  7. Protein level uses the existing picked protein groups on the global pair winners.
  8. `score_comparability` is the G4-style decoy KS. `runner_up_stats`/`calibrate` is a label-blind per-run calibration.

  Tests are in `tests/stats/test_global_fdr.py` (4 pass). The Monte-Carlo shows global control at N = 1 and 20 while the naive union
  reaches 15% FDP at N = 50. The tests also check equality with the per-run procedure for one run, and permutation invariance under ties.
- **E9 stream (`scripts/e9/`).**
  - `select_e9.py` is the D8 sampler.
  - `dl_bench.py` measures download throughput.
  - `stream.py` downloads, then runs the cohort with arm P, then the GPU search with arm C, then deletes the RAW. It has backpressure and is
    resumable. Downloads are size-checked **and checked against the PRIDE SHA-1**.
  - `score_worker.py` is a persistent L3 feature and frozen-GBM process (2-3 concurrent; claims batches atomically). It deletes stores.
  - `stage_pred_env.sh` stages the pred venv on NVMe.
  - `analyze.py` and `c_homology.py` do the analysis; `summarize.py` builds `metrics.json`.
  - Small additive change to `vulcan/cohort.py`: `--skip-failed`. A run whose ingest or search fails is recorded in `failed` instead of
    aborting the cohort. The default behaviour is unchanged.
- **Data, D8 (DATA.md).**
  - 1,000 human DDA OT-HCD runs from 142 PRIDE projects, drawn from the T4 crawl frame.
  - Excluded: T4 selection, every D*/T*, the 13 M3 training projects and the predictor-training accessions. A re-check found 0 overlaps.
  - Seed 11, nested round-robin order: N = 1, 10, 100, 1000 span 1, 10, 100 and 142 projects. 1.08 TB.
  - Manifest `data/manifests/e9_d8.tsv`, log and sha256 (`e9_d8.sha256`, 1,000 files).
- **Search and score.**
  - The arm-P augmented index is human+contaminants + paired shuffled entrapment E_2001 + decoys (1 entrapment seed × 1 decoy seed, as
    pre-registered for E9).
  - The arm-C index is human + Arabidopsis foreign peptides (r = 1.02) + decoys.
  - Top-5 PSMs per spectrum are re-ranked by the **frozen default `static_gbm_pred`** (D7, A1.1). The hyperscore is kept as a fixed
    reference.
  - Unit is the peptidoform. Ties are broken by hash(pf_id, 11).
  - FDP̂ comes from the augmented searches. No production search was run, so the IDs below are from the augmented search.

## 1. Download throughput (`download_probe*.json`, `metrics.json:stream_sessions`)
All traffic leaves through the cluster HTTP proxy (`proxy.service.consul:3128`); a direct connection is refused. **Aspera is unavailable:**
there is no `ascp` binary, although PRIDE lists `prd_ascp@fasp.ebi.ac.uk` locations.

| | 1 stream | 4 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|
| HTTPS (GB/s), 30-40 s probe | 0.003 | 0.017 | 0.068 | 0.120 | 0.108 | 0.113 |
| FTP via proxy (GB/s) | 0.006 | 0.030 | 0.077 | 0.137 | 0.087 | 0.194 |
| **Sustained in the stream** (HTTPS, 64 streams, 69 min, 811 GB) | | | | | **0.197 GB/s = 668 runs/h** | |

- A single stream gets about 3-6 MB/s. The aggregate saturates at **about 0.1-0.2 GB/s** at 32+ streams, which is the proxy/WAN ceiling.
- That is **about 670 runs/h** at the D8 mean of 1.08 GB per run, and about 1.6 h per 1,000 runs.
- Download **binds**. It is about 7× below the decode capacity measured here on the loaded node (~4.8k runs/h at 24 workers) and about
  40× below M2b's warm-cache figure (~28k runs/h).
- **Integrity problem.** 48 of about 650 verifiable downloads (7%) had the right size but the wrong PRIDE SHA-1. Every one came from a
  transfer that had to be resumed (`curl -C -`) through the proxy.
  - 34 of them were the bulk of the "ingest failed" runs. The other 14 ingested without error, which means they were silently corrupt.
  - Fix: a mismatch now deletes the file and re-downloads it from scratch. The 126 outputs of affected runs were quarantined and redone.
    The final pass caught and fixed 5 more.
  - **350 of the 1,000 D8 files have no PRIDE checksum**, so they cannot be verified. At the 7% rate, about 25 of them may be corrupt.
    For FDP this only costs power, because TDC and entrapment are label-blind to corrupt spectra. Fix for next time: never resume through
    the proxy, or download checksum-less files twice and compare.

## 2. End-to-end runs/hour and per-stage costs (`metrics.json:stage_costs`)
- **End to end** (download → decode → GPU search P + C → L3 features + GBM for both arms → delete):
  - Main session: **773 runs in 79.6 min, about 580 runs/h.** Download bound, with scoring close behind once 3 scorer workers ran.
  - All sessions together, including an OOM incident, the corrupt-download redo and analysis pauses: 1,000 manifest runs → 982 complete in
    about 3 h of streaming wall.
  - **N = 1000 took about 1.6-2 h of streaming**, against the < 12 h budget. The estimate made before starting was 2-2.5 h.

| Stage (per run, mean) | Cost | Capacity on this node |
|---|---|---|
| Download (1.08 GB) | 0.2 GB/s shared link | **~670 runs/h (binding)** |
| SHA-1 + sha256 hashing | 3.7 s on 1 core | ≫ |
| Decode (24 persistent workers; RAW copy to NVMe, store, MS2 prep) | **4.7 core-s**, 16-20 s wall per worker | ~4-5k runs/h on a loaded node |
| Cohort process per batch (setup + decode + search, makespan) | 1.4 s/run, 7.2 core-s/run including worker start | |
| GPU search arm P (T+E+decoys, 14.1M rows) | **0.19 GPU-s** | ~19k runs/h |
| GPU search arm C (+ store reload) | 0.22 GPU-s, 0.55 s wall | |
| L3 features + GBM, per arm | **6.6 s wall** (features 5.2, of which AlphaPeptDeep misses 2.4; GBM apply on CPU 1.15) | 2 arms × 3 workers ≈ 830 runs/h (**co-binding**) |
| L6 analysis at N=1000, per arm and score | load 86 s, global competition 29 s, per-run competitions 65 s, global picked protein 26 s, 982 per-run protein groupings 35 s (24 procs) | |

- **Disk.** The pending RAW budget was capped at 150 GB with 3 unscored batches; spot checks showed up to 190 GB in `/tmp/vulcan_e9`, under
  the 400 GB limit. What remains is 47 GB of PSM, feature and score tables plus L6 outputs, on `/tmp` only.
- **Next levers.** The link is the limit: more nodes or a mirror, or Aspera/Globus if it can be installed. After that comes scoring:
  - a GPU GBM (FIL-style);
  - a pre-warmed prediction library per NCE (cache misses were 37% of feature time);
  - skipping arm C in production.

## 3. E9 result: global FDP̂ vs N, naive merge vs the L6 global procedure
Each cell is the mean FDP̂ with its 95% CI (pair bootstrap, 1 replicate). Columns marked "IDs" count target peptidoforms (protein: target
groups) in the augmented search. "Runs" is the number of runs that completed. Arm P is paired shuffled entrapment (r = 1). Arm C is the
Arabidopsis combined estimator (r = 1.02). Protein level uses arm P (§6.3).

**Default scorer `static_gbm_pred`, α = 1%**

| N (runs / projects) | Level | Arm | **Global FDP̂ %** | Global IDs | **Naive-union FDP̂ %** | Naive IDs |
|---|---|---|---|---|---|---|
| 1 (1/1) | peptide | P / C | 1.03 [0.87, 1.19] / 0.99 [0.84, 1.15] | 30,925 | same (N=1) | 30,925 |
| 10 (10/10) | peptide | P / C | 1.03 [0.92, 1.14] / 0.99 [0.89, 1.10] | 64,969 | 1.50 [1.37, 1.63] / 1.29 [1.17, 1.41] | 67,218 |
| 100 (97/97) | peptide | P / C | 1.21 [1.14, 1.29] / 1.11 [1.05, 1.18] | 203,742 | 3.30 [3.19, 3.41] / 3.02 [2.92, 3.12] | 210,972 |
| 1000 (982/140) | peptide | P / C | **1.16 [1.11, 1.21] / 1.21 [1.16, 1.26]** | 358,164 | **9.33 [9.20, 9.46] / 9.55 [9.43, 9.68]** | 388,460 |
| 1 | protein | P | 0.55 [0.23, 0.87] | 3,616 | same | 3,616 |
| 10 | protein | P | 1.12 [0.77, 1.46] | 5,786 | 2.47 [1.96, 2.98] | 7,024 |
| 100 | protein | P | 1.06 [0.84, 1.29] | 11,179 | 5.91 [5.44, 6.39] | 16,426 |
| 1000 | protein | P | **1.00 [0.80, 1.20]** | 13,402 | **19.13 [18.50, 19.76]** | 28,323 |

**α = 5%, peptide, P / C:**

| N | Global | Naive |
|---|---|---|
| 1 | 5.19 / 4.81 | same |
| 10 | 5.16 / 4.75 | 7.52 / 6.86 |
| 100 | 5.69 / 5.23 | 16.1 / 15.1 |
| 1000 | 5.75 / 5.94 | **36.4 / 37.9** |

At protein level and α = 5%, the global list is at 1.03% (the protein q saturates) and the naive union at 19.4%.

**Fixed hyperscore reference, α = 1%, P / C:**

| N | Global peptide | Naive peptide | Global protein | Naive protein |
|---|---|---|---|---|
| 1 | 1.06 / 0.95 | – | 0.93 | – |
| 10 | 1.03 / 1.01 | 1.56 / 1.58 | 1.34 | 2.63 |
| 100 | 1.04 / 1.05 | 3.25 / 2.97 | 1.17 | 6.25 |
| 1000 | **1.14 / 1.09** | **8.90 / 8.47** | **0.99** | **18.4** |

The full table, including LB, CIs, verdicts and the 0.1% level, is `metrics.json:e9_table`.

**Reading.**
- **The DESIGN claim holds.**
  - The naive union of per-run 1% lists inflates steadily with N: from 1% to **9.3% (peptide) and 19% (protein) at N = 1000**, and to
    36% at α = 5%.
  - §6.3 verdicts for the naive union: every peptide arm is "violated" from N = 100 on (lower CI of LB > α). Protein is violated at 1%
    from N = 100 and at both α at N = 1000.
  - For the global procedure no level is ever violated. Peptide lists are "inconclusive", which is the expected status of the
    upper-bound estimators. Protein lists are "supported" at N = 1000.
  - The L6 global procedure stays at **1.0-1.2% (peptide) and 1.0% (protein)** for N = 1 to 1000, with both scorers. At 5% it stays at
    4.8-5.9%.
  - The global list keeps **92%** of the naive union's peptide IDs at N = 1000 and 47% of its protein groups. Most of the union's extra
    proteins are false.
- **The residual 10-20% excess of global FDP̂ over α at N ≥ 100** is in the range where the upper-bound estimators overstate (§6.3: the
  paired estimator gives 1.24% for a list whose true FDP is 1.04%). Two diagnostics (`metrics.json:fair_coin_family_counts`) split it:
  1. **Arm P: an entrapment-share effect, not a broken coin.** Above the global threshold, entrapment winners and their own decoys are
     balanced: E vs dec(E) = 0.497 at N=1000 (z = −0.4). But dec(E) take a growing share of all decoys: 0.49 at N=1, 0.59 at N=100,
     0.63 at N=1000 (hyperscore: 0.54). Shuffled entrapment peptides are mass- and composition-matched to *present* targets, so across
     many runs they attract more null hits than absent targets do. LB and the paired FDP̂ therefore overstate when N is large. The TDC
     decoy count of the list (3,601 at 1%) is consistent with control.
  2. **Arm C: a small target-favouring tilt of the predicted-feature scorer.** Natural Arabidopsis nulls beat their decoys 54.4% of the
     time at the 1% threshold (z = 5.6, N = 1000) and 52.1% at 5% (z = 6.3). Hyperscore gives 51.6% (z = 1.6) and 50.5%. The tilt is
     **constant in N** (0.536 / 0.543 / 0.538 / 0.544 at N = 1 / 10 / 100 / 1000), so it is a per-competition effect that becomes visible
     only with repository-scale n.
     - **Size.** If absent human targets share it, the FDR at α = 1% is underestimated by about 15-20%, giving roughly 1.2%. That matches
       the C column.
     - It is within the A1 tail-G2 margin (difference to hyperscore ≤ 0.10) but much larger than the M3 G2 overall difference (+0.0001).
     - The 6-mer homology check (`c_homology.json`) does not explain it: accepted foreign peptides share a human 6-mer at the background
       rate (0.905 vs 0.925 random). That criterion has little power, though, because 92% of random foreign peptides share some 6-mer.
- **Calibration and comparability (G4-style, decoys only, dev).** The frozen GBM is **not** comparable across runs.
  - Per-run decoy KS against all other runs has median 0.18 and 90th percentile 0.30. Same-run half-split noise is 0.01-0.03, so the
    median excess is 0.15. Per instrument the KS is 0.04-0.09, and 0.35 for the 8 Exploris 240 runs.
  - Run-level 1% cutoffs span an IQR of 1.08-2.74 score units. The global 1% cutoff is 4.74.
  - Two label-blind per-run calibrations were tried (both valid by construction):
    - runner-up-candidate location/scale: KS 0.31, and only 18k IDs instead of 65k at N=10;
    - low-quantile bulk: KS 0.15, but −13% (N=10) to −68% (N=100) IDs.
  - Both lose power heavily, so they were **not adopted**. The pooled score is the raw frozen score.
  - Validity does not need comparability: best-of-N over exchangeable pairs stays a fair coin, and hyperscore stays flat. Power does. A
    decoy-tail calibration from a separate calibration-decoy set (CALIB_DECOY_SEED 3001; it would need a third search) is the principled
    route, and is listed as next.
- **A1.1 confirmation (first 10 D8 runs, not used for selection; `metrics.json:A1_1_confirmation_first10`).** Per-run lists at 1%,
  static_gbm_pred vs hyperscore on the same augmented searches:
  - IDs: 106,173 vs 62,372 peptidoforms summed over runs (**+70%**).
  - Mean per-run FDP̂: P 0.95% vs 1.39%, C 0.92% vs 1.08%.
  - Verdicts at 1% (10 runs): supported P 3/10 and C 8/10 for GBM, against 2/10 and 2/10 for hyperscore. One run per scorer is
    not-evaluable (180 IDs). There are no "violated" runs.
  - These are dev numbers on unselected runs; the per-run verdicts are inconclusive or supported.

## 4. Failures and issues
- **18 of 1000 runs did not complete; 140 of 142 projects remain.**
  - PXD074511 (8 runs, Exploris 480, 2026): RawFileReader throws in its constructor. This is likely a newer file format than our
    RawFileReader. Needs a reader update.
  - PXD072218 (8 runs, "QEP2 … VQ"): 0 PSMs, so these are probably not DDA MS2 despite the metadata. The scorer error on empty tables is
    recorded.
  - PXD000695 E9_0052: .NET exception in centroid conversion.
  - PXD023001 E9_0735: RawFileReader error.
  - E9_0175 has arm P but not arm C.
- **A GPU OOM on the shared GPU cost a re-download of about 200 runs.** Two score workers grew to 16 and 38 GB of torch cache, so cohort
  index builds failed while the RAWs were already deleted. Fixes: whole-cohort retry with backoff before any RAW is deleted, feature chunk
  of 50k, `empty_cache` per run, `expandable_segments`.
- Corrupt downloads through resumed proxy transfers; see §1.
- **Arm P/C entrapment at repository scale.** Shuffled-entrapment LB/paired estimators lose their "E ~ false target" assumption as N grows
  (the share effect above). Report the fair-coin diagnostic next to the verdict. Possibly propose (not adopted; §6 frozen) an E9-specific
  note in a dated amendment.
- **For 3b / A1, a flag on the default scorer:** the arm-C tilt (54% natural-null target wins vs 51.6% for hyperscore) should be looked at
  with A1 tail-G2 and arm M before T1. It does not break global flatness, but it is the largest scorer-specific signal found so far. It
  agrees with D12, recorded in parallel by M3e: tail-G2, G1 and G3 FAIL for static_gbm_pred, and FDR claims stay blocked. E9 is
  independent evidence of the same direction, on natural foreign nulls at repository scale.
- No production (entrapment-free) search was run, so IDs are from the augmented search. A replicate (seeds 2002/2003, decoy seeds
  1002-1005) was not run at N ≤ 100, though it would be cheap.
- **Future dev/test selections from the PRIDE frame should also exclude `e9_d8.tsv` accessions.** Other agents are now drawing g2null and
  d11/d12 sets.
