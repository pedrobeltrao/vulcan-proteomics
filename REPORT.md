> **Notice.** This work and this report were produced by an AI agent team (Claude Opus 5.5 / Sonnet 5.5). They have not been reviewed or verified by human experts.

# Vulcan: final report (session ending 2026-10-28)

Status: written from committed material (PROMPT, RESEARCH_TREE D0–D39, DESIGN, DATA, LESSONS, results/2026-10-2*/NOTES+metrics). Section 3 reports the single sealed A5 look (T5 + T1), opened at `milestone-M7-2026-10-28` (results/SEALED_LEDGER.md). The node allocation was expected to end about 22:00 on 2026-10-28, which bounded the last day (D31); the node was in fact wiped around 00:30 on 2026-10-29, during the look (D38). Clock times are CET (UTC+1); the sealed ledger is in UTC.

**Conventions.** Every number is **dev** (development sets D0-D12, DN*, D8/D8h; see DATA.md) unless it says *sealed*. "Nominal" = the FDR threshold we asked for. "Verified FDP" = a false-discovery proportion measured against a truth set (library design, spike-in) on dev data; entrapment gives an estimate (FDP̂), not a verified FDP, and entrapment verdicts use the pre-registered words (supported / inconclusive / violated / not evaluable). No FDR claim in this report has passed the null-spectrum gates; see section 4. Hardware of record: gpu-node, 32 affinity cores (2x EPYC 9354), A100 80 GB, shared node (load typically 14-110; timings under `flock` bench lock where stated, comparators often unlocked). Core-s = user+sys CPU seconds; they exclude GPU time, which Vulcan uses and the CPU-only comparators do not. Sage figures are 32-thread unless noted.

## 0 Summary: accuracy and speed vs comparators (all dev data; no held-out test-set comparison)

**The project is missing its held-out test-set comparison.**
- The pre-registered single sealed look (T5 + T1) was started on 2026-10-28.
  - T5 finished.
  - T1 was killed in its final scoring stage when the compute node was wiped (around 00:30 on 10-29).
- It was closed as 'failed (unterminated)' with all outputs deleted unread (ledger commit 742ba48). Nobody has seen any T1 or T5 number.
- T1–T5 therefore remain unseen. A5 allows at most two retries from the same frozen commit (4fe2e2a); one attempt has failed, so up to two retries remain. A retry needs about 13 h of uninterrupted node time, including a check of the rebuilt mokapot env (D38/D39).
- Every number below comes from development sets, several of them used for tuning, as marked in sections 2 and 7. No FDR claim has passed the null-spectrum gates (section 4); error rates are entrapment-estimated on dev.

| Workflow | Accuracy vs comparator (dev) | Speed vs comparator (dev, 32-core EPYC + A100, shared node) |
|---|---|---|
| Standard DDA search (frozen v1) | 0.990× [0.989, 0.992] peptidoforms vs Sage + predicted features + mokapot at matched 1% FDP (D0); +2.8% vs plain Sage | search stage about 23–40× fewer CPU core-s than Sage at 4–32 threads (measured single-process 6.0 core-s) (CPU only; excludes GPU time); RAW → proteins about 5,700 runs/h at about 9 core-s/run vs Sage 70–120 core-s (4–32 threads, from mzML) + about 250 for RAW→mzML |
| Chimeric search (v2_sep, candidate) | 1.150× [1.069, 1.203] peptidoforms vs Sage chimera + predicted features + mokapot on D8h (10 dev runs held out from scorer selection; not a sealed test set); protein groups 0.999× | about 3,900 runs/h at about 15 core-s/run |
| Ion-trap CID / TMT | 1.034× (D4) and 1.093× (D10cf) vs Sage + Prosit features + mokapot | about 0.9 core-s per fraction vs Sage about 54 core-s per file (16 threads, unlocked, indicative) |
| Label-free protein abundance | parity with Sage on the HYE benchmark (matched-set error 0.118 vs 0.1145) | about 1,200–1,280 runs/h at about 15 core-s/run vs Sage search+LFQ 266 core-s/run (+ about 250 for mzML conversion) |
| Phospho + site localization | Ala-decoy FLR point estimates near true; stricter upper-bound rule not met on PXD009449 (n = 183); mixed vs pyAscore; not run-to-run reproducible on about 0.7% of PSMs | about 3,200 runs/h; Sage 3,087 core-s for 2 runs, without localization |
| Open search | −4% peptide stems vs Sage-open | search stage about 65× fewer core-s than Sage-open (32 threads) |
| Cross-linking | DSSO (D3b): about 73% of xiSEARCH heteromeric pairs, at lower true FDP | DSS (D3): about 460× less wall time than xiSEARCH for the whole pipeline, about 8,900× for the search stage |
| Repository-scale FDR (L6) | global peptide FDP̂ 1.16% (arm P) / 1.21% (arm C) and protein 1.00% at N = 1,000 (982 runs completed), vs 9.3% (peptide) and 19.1% (protein) for the naive union | streaming from PRIDE is download-bound, about 580 runs/h |

**Example: 1 TB of RAW files for protein abundance (LFQ), one machine** (estimate, not a measured end-to-end run). Assumes about 1 GB per run, so about 1,000 runs.

| | Vulcan (32 cores + 1 A100) | Sage + mzML conversion (same 32 cores, CPU only) |
|---|---|---|
| CPU cost per run | about 15 core-s | about 266 core-s (search + LFQ) + about 250 core-s (RAW→mzML) |
| Wall time for about 1,000 runs | **about 0.8–1 h** (1,200–1,280 runs/h measured on 96-run cohorts from local disk) | **about 4.5–7 h** (516k core-s over 32 cores at ideal scaling is 4.5 h; 2.3–5 h if mzML already exists) |
| Ratio | about 5–7× less wall time per machine; about 17–34× fewer CPU core-s | — |

Caveats on the 1 TB estimate:
- Vulcan needs an A100-class GPU, so the per-machine comparison favours a GPU node.
- The two runtimes come from different runs: Sage's 266 core-s is on Astral D9, and Vulcan's cohort mixes QE runs with D9 quant.
- Runs/h varied by about ±25% on the shared node.
- Reading from network storage or PRIDE (about 580 runs/h) rather than local NVMe would make I/O the bottleneck for both.
- Accuracy is at parity, and error-rate claims are not gate-supported.

## 1 Approach

**Architecture (DESIGN §3-4).** L0 ingest (direct RawFileReader/Bruker readers to an Arrow store) -> L0.5 per-run pilot (calibration, analyzer/activation, reporter checks) -> L1 content-hashed proteome index (sage-core digestion, paired C-term-fixed decoys, entrapment, GPU tiled inverted index; closed, phospho 184.5M peptidoforms/34 GB, semi-specific 33M, open) -> L2 GPU candidate generation (K1 CUDA, pair-closed top-K=500) -> L3 exact GPU rescoring (K2, label-blind features, cached predicted-spectrum/RT features, localization operator) -> L4 statistics (paired double competition, CONGA/stratified TDC, picked protein, FLR, XL levels) -> L5 cross-run (MS1 features, alignment, MBR with transfer FDR, LFQ/TMT) -> L6 repository (best evidence per peptidoform, global paired TDC, picked protein, master lists).
**Label blindness** is enforced by a test (G0: swapping all target/decoy labels and permuting storage order leaves L2/L3 outputs bit-identical); it passes for closed, chimeric v2_sep, IT-CID (with the pinned Koina store), XL and semi-specific search.

**CPU/GPU co-design.** D0 assumed 4 cores (all matching on GPU); D5 moved to 32 cores + A100, so work goes to whichever device finishes the pipeline first. GPU: index lookup, scoring kernels, feature building, GBM/linear apply, (opt-in) MS1 hill finding. CPU: RAW decode in supervised persistent worker processes (RAW staged to NVMe; mmap of Lustre stalled 5 s -> 77 s), picked protein, greedy MS1 acceptance. One long-lived `vulcan.cohort` process keeps the index resident. Thread pools come from the affinity mask (`vulcan.hw`). Process lessons are in LESSONS.md.

**Statistics protocol and pre-registration history (told as it happened).**
- DESIGN §6 was frozen at `prereg-v0.4` before any sealed access. Entrapment arms P (shuffled) and C (natural foreign), decision rule, FLR decoy (Ala primary), dev/sealed split fixed up front.
- **A1** (D9, after an independent review showed the first gates were weak): froze v1 `static_gbm_pred`, added tail-G2, required G1/G3/G4/arm M, strongest comparator. v1 then **failed** tail-G2, G1 and G3 (M3e). A post-hoc explanation (human-contaminated *S. aureus* null, offset nulls that are shifted human spectra) exists but was found after the failure and is unconfirmed; the failure was recorded as final (D12).
- **A2** (fresh nulls, label-symmetric admission, adopted verbatim before any listing): the walk admitted 6 projects, below the minimum of 10, so A2 is **not evaluable** (D24). One earlier walk had 0 admissions because of an implementation bug (D16), re-run once.
- **IT-CID nulls (DNIT)** were futile (581 usable spectra vs 120k needed).
- **A3** (v2_sep as gated/claimable plus secondary T1 arm): adopted (D22), but its DN2 tiers are "expected not evaluable" by its own power check, and DN2 step 5 was never run.
- **A4** (one new null set for v1 and v2_sep) was drafted, reviewed twice, **not adopted and not run** (D28): simulated P(v1 minimum) 0.08–0.43 (a clean v1 would pass with probability about 0.06–0.4).
- **C9 pilot** (label-free, simulation only): human removed-peptide nulls GO, P(v1 minimum) 0.97-1.00 at <=200 runs (D30).
- **A5** (adopted D29): relaxes the claim conditions knowingly, so the single T1+T5 look proceeds although v1 FAILED and A2 was not evaluable; every FDR-filtered number then carries A5 labels and the gate table; one look is enforced in code.
- **A6/HN1** (adopted D32): one human removed-peptide null test of tail-G2/G1/G3 for **v1 only**, at alpha/2. **Not tagged: infeasible** by its own pre-registered tag rule (D34).
  - Measured admission throughput on the shared node was about 116 runs/h, which gives P(v1 minimum) = 0.63 against the 0.8 required.
  - Its statistics acceptance test also found that the percentile cluster bootstrap behind G2-original undercovers at J ≈ 37 clusters (0.923–0.944 vs 0.95).
  - No candidate was listed, so the pool is unspent, and A6 is ready to run on a quiet node.
- **Final freeze and tags** (D32–D36): G0 for the full v2_sep pipeline passes at the final freeze commit 06ce624, and the D1 dry run is byte-exact for every v1 artefact. The D0 dry run is cited from the rehearsal (disclosed deviation). The closure guard passes. Tags `prereg-v0.4.3`, `prereg-v0.4.4` and `milestone-M7-2026-10-28` were set at 4fe2e2a, and the single A5 look (T5 then T1) opened at 14:23 CET (ledger 13:23:46Z).
- **T2 (FLR), T3 (XL) and T4 (repository sample) stay sealed.** The T2 driver and its dev dry run are committed and ready (results/2026-10-28_T2_driver_dryrun). T2 was not run alongside T1, because of the quiet-node rule and GPU memory (D36).

## 2 Best methods (accuracy and runtime; all dev)

"Runs/h" is RAW to the stated output in `vulcan.cohort`, page-cached unless noted. Comparators were run at their own thread settings; core-s is the transferable number.

### 2.1 Default scorers per instrument class

| Class / method | Accuracy (dev, vs comparator) | Runtime | Source |
|---|---|---|---|
| **OT-HCD DDA, v1 `static_gbm_pred`** (frozen; GBM over label-blind + AlphaPeptDeep features; per-run NCE) | D0 1%: 24,512 peptidoforms vs Sage-LDA 23,703 under the same estimator; fair lead after Sage RT-on is +2.8% at 1%, +8.5% at 0.1%. **Against the strongest comparator (Sage+AlphaPeptDeep features+seeded mokapot) at matched arm-P FDP 1%: D0 0.990 [0.989, 0.992], D1 0.998 [0.997, 1.000]** (all D1 comparator ratios lie in 0.98–1.004), i.e. non-inferior; the gain over Sage is the predicted features. D8h holdout vs hyperscore 1.55 [1.41, 1.71] per run (not vs the strongest comparator). Entrapment at 1%: P inconclusive (1.0-1.2%), C supported on D0/D8h; nothing violated. | Search kernel 0.24 s vs Sage 2.98 s (32 thr). D0 end-to-end search 7.5 s / 4.9 core-s vs Sage 29.9 s / 237 core-s (32 thr), 46.5 s / 137 core-s (4 thr): **about 48x (32 thr) and 28x (4 thr) fewer CPU core-s** (CPU only; excludes GPU time). The M3 erratum notes that 7.5 s / 4.9 core-s were derived; the measured single-process figures are 8.59 s / 6.0 core-s (40x and about 23x). Scoring after speed work 0.54 s / 0.53 core-s / 0.087 GPU-s per run (was 3.89 s / 4.60 core-s), outputs bit-identical. Cohort RAW->protein groups **5,658 runs/h** (6.8 core-s/run, 192 runs, load 60-70); pure search 7,485/h; closed cohort without scoring 20-26k runs/h page-cached. | results/2026-10-27_M3_scoring, M3e_hardening, T1driver_dryrun, E2_gpu_search, scoring_speed, M2b_cohort |
| **IT-CID/TMT, `itcid_linear_pred`** (Prosit CID/TMT features via pinned Koina store) | vs Sage+Prosit features+seeded mokapot at matched arm-P FDP 1%: D4 1.034 [1.030, 1.037], D10cf 1.093 [1.088, 1.097] (sequence units). Entrapment inconclusive, none violated. **FDR claims 'not gate-validated'** (DNIT futile). | 0.92 s / 0.91 core-s per D4 fraction (warm cache); Sage 68 s / 643 core-s for 12 files at 16 threads (unlocked, indicative); first contact ~16 s Koina calls | results/2026-10-27_3c_itcid (+followup) |
| Cross-encoder `ce_s` (1.66M params) | 0.983 [0.968, 0.999] vs v1 on D8h: **NO-GO** as default | 0.56 GPU-s / 0.68 core-s per run | results/2026-10-28_M7_encoders |

Open/PTM/semi search in the same engine (dev, hyperscore pipeline): **open [-150,+500]** 17.2 s / 17.2 core-s vs Sage-open 72 s / 1,153 core-s (4.2x wall, ~65x core-s; GPU time excluded) with 20,488 stems vs 21,400 (-4%); mass-offset mode 2.7 s, +14% stems over closed (M4_ptm). **Semi-specific LiP** search 0.91 core-s vs Sage semi 446-450 core-s per run (M5_lip_cf). Phospho index (184.5M peptidoforms) searches as fast as closed (0.2 s/file); its CPU label build takes 446 s (not optimised).

### 2.1b End-to-end speed for standard abundance and PTM workflows (E12, last day, opt-in/exact, not yet merged)
96-run cohorts, run back to back under the bench lock (node load 33–40); before → after in the same session. Every change is bitwise-identical, or opt-in with default off. (results/2026-10-28_E12_speed_abundance_ptm on branch worktree-agent-a0df6277cef40ec85; merge after the looks, D35.)

| Workflow (RAW → output) | runs/h | core-s/run | Comparator |
|---|---|---|---|
| closed search + scoring + proteins | 4,570 → 5,753 | 9.9 → 9.1 | Sage 70–120 core-s/run + 250 for RAW→mzML |
| LFQ, MBR off / on | ~560 → 1,280 / 1,194 | 19.7 → 15.3 / 15.5 | Sage search+LFQ 266 core-s/run (D9 Astral; not like-for-like data) |
| chimeric v2_sep | 2,567 → 3,867 | 16.2 → 14.9 | not timed |
| phospho + localization + Ala-decoy FLR | 2,362 → 3,172 | 10.7 → 9.8 | Sage 3,087 core-s for 2 runs incl. index build, no localization |

- **Stage speed-ups:** directLFQ 22–28 → about 5 s per experiment (identical files); phospho index build 444 → 81 s (byte-identical, 184.5M peptidoforms); scorer ready 18–23 s instead of 37–52 s.
- **Next limits:** decode (6–7 core-s/run, vendor reader) and GPU-stage Python.
- **Noise:** absolute runs/h varied about ±25% between rounds on the shared node, so compare within a row.

### 2.2 Chimeric search v2_sep (candidate, not claimed)
- Method: each co-isolated MS1 feature becomes a query unit, CONGA neighbour filter, per-stratum TDC (`sep`; strata = chimeric rank 1, 2, >=3); linear scorer on the v1 feature set retrained on 40 projects. Group-walk was rejected (stratum >=3 violated on D0, 3.1%).
- **Accuracy vs the fair comparator (Sage chimera R=2 + predicted features + seeded mokapot), matched arm-P FDP 1%, peptidoforms: D0 1.149 [1.147, 1.152], D1 1.137 [1.120, 1.140], D8h holdout 1.150 [1.069, 1.203].** Against the non-chimeric comparator it was 1.17-1.22. **Protein groups on D8h: 0.999 [0.951, 1.039], i.e. no protein gain.** At 0.1% the gain shrinks (D8h 1.096; D1 own-q 0.964 [0.959, 0.971]). Sage chimera's overall arm-P FDP is near nominal (0.91-1.09% at 1%), but its rank>=2 IDs run at 4-13% FDP at nominal 1%, so ratios favour the comparator. No stratum is violated for v2_sep on D0, D1, D8h.
- **Runtime per D0 run (exact ports):** 16.9 s / 32.5 core-s -> 9.3 s / 15.5 core-s; cohort chimeric **2,987 runs/h (steady 4,643), 19 core-s/run** vs closed 5,433 in the same session; MS1 is 79% of core-s. (results/2026-10-28_E7_chimeric_v2, E7b_fair_speed)
- **Opt-in GPU MS1 (E11), bitwise-exact on 13 dev runs:** MS1 12.4 -> 1.8 core-s per D0 run; whole chimeric pipeline 14.3 -> 4.4 core-s; 96-run cohort 3,307 -> 3,647 runs/h incl. setup (1.10x), steady 5,718 -> 11,892 runs/h. Merged opt-in (default cpu) so the frozen path is unchanged. (results/2026-10-28_E11_gpu_ms1)
- Cost vs v1: ~2.5 -> 9.3 s wall per D0 run. v2_sep never passed a null-spectrum gate, so it is labelled 'candidate'.

### 2.3 Global FDR (L6)
Best evidence per peptidoform -> global paired TDC -> picked protein (`vulcan/stats/global_fdr.py`). On D8 (982 dev runs, frozen v1, 1 entrapment x 1 decoy seed, nominal 1%): entrapment-estimated global peptide FDP arms P/C 1.03/0.99% (N=1), 1.21/1.11% (N=100), **1.16/1.21% (N=1000)**; protein 1.00% at N=1000. Naive union of per-run 1% lists: **9.3% peptide, 19.1% protein at N=1000**. So global control stays near nominal where the union inflates; the 0.2% excess has two identified causes (entrapment-share effect on P; target-favouring tilt on C, see 5). Streaming from PRIDE is download-bound: ~580 runs/h (0.1-0.2 GB/s proxy ceiling); GPU search 0.19 GPU-s/run, decode 4.7 core-s/run. (results/2026-10-27_E9_scale)

### 2.4 PTM / FLR
Ala-decoy FLR (PSM level, D8): estimated vs true FLR of accepted sites, mono-phospho, dev truth. At nominal 1%/5%: PXD000138 0.29%/3.7%, PXD009449 0.87%/4.4%, PXD007058 proxy 0.0%/2.3% (calibrated or conservative). The Leu decoy under-estimates slightly (1.31% at nominal 1% on PXD000138); secondary only. Site recall at matched true FLR 1% vs pyAscore: ahead in 6/6 cells (+0.6 to +17.5 points; PXD007058 0.758 vs 0.613, untuned); at 5% pyAscore is ahead in 4/6. PXD000138/PXD009449 numbers are tuned on dev; PXD007058 is a single-engine proxy truth. Localization runs at 5.8–23k PSMs/s (CPU prep dominates), against pyAscore's 7–16k/s on one thread. (results/2026-10-27_M4_ptm)

**Last-day revisions** (results/2026-10-28_T2_driver_dryrun, E12):
- **Stricter T2 rule.** Under the T2 driver's rule (a cell passes only if the one-sided 95% upper bound of true FLR is ≤ α), the dev E8 gate is **not passed**. PXD000138 Ala passes at 1% and 5%. PXD009449 is inconclusive at 1% (0.86%, upper bound 1.34%) and at 5% (4.44%, upper bound 6.93%), because only 183 single-Y peptides carry that set.
- **pyAscore, paired bootstrap.** We are higher only on PXD009449 at 1%, show no difference on PXD000138 at 1%, and are lower on PXD000138 at 5%. The '6/6' above was based on point estimates.
- **[!] Localization is not run-to-run reproducible.** GPU atomic sums break exact score ties at random, and 313 of 45,845 PSMs on PXD009449 changed their best site configuration between two identical runs. Over 11 repeats, no dev verdict changed, and M4's values lie inside the observed range in 30 of 36 fields.
  - An opt-in deterministic tie rule is reproducible but accepts 1.7% fewer sites at 1% FLR on PXD009449. It is not adopted, because the T2 path stays as pre-registered (D35).

### 2.5 Quantification
- **LFQ (D9 Astral HYE, MBR off, after M6d fix):** precursor median |eps| 0.122, CV 4.6%; matched-set eps 0.118 vs Sage 0.1145; dLFQ pAUC 0.656. Honest summary: at parity with Sage on accuracy (not clearly better), better CV earlier (5.0 vs 6.1%). D11 same-precursor eps 0.100 vs 0.108 before. Cost +5.8 core-s per Astral run over RAW->PSM (MS1 4.5 + quant 1.3); Sage search+LFQ 266 core-s/run, FlashLFQ 153 core-s/run (locked, 32 thr). (M6_quant, M6d_quant_accuracy)
- **MBR** is opt-in. False-transfer FDP on D9 H-runs 0.63-1.23% at nominal 1%; fresh D12 0.18-0.79%, but the transfer decoy under-counts ~2.4x there (see 7).
- **TMT (D4):** reporter values equal Sage's (r = 0.9999996), 0.5 s/run; Sage 48.6 s / 792 core-s for 12 files.

### 2.6 XL / LiP / CF
| Assay | Result (dev) | Runtime | Source |
|---|---|---|---|
| XL DSSO D3b | hetero residue pairs 334 at true FDP 1.2% (nominal 1%), 456 at 4.8% (nominal 5%); PPI 169 at 1.2%. Submitter xiSEARCH+xiFDR: hetero 459 at 2.0%. **Sensitivity ~73% of xiSEARCH.** | search 4.7 s, +6.2 s CPU site localization; xiSEARCH extrapolated ~27 h / 1.48M core-s | M5_xl |
| XL DSS D3 | 0 false heteromeric links; self residue pairs at nominal 5%: 119 at 9.2% true FDP (conservative estimator 4.7-4.9%); xiSEARCH 129 at 0.8% is better on self pairs | ours 0.43 s search; xiSEARCH 3,812 s / 53,671 core-s | M5_xl |
| LiP D7 | 26,976 LiP peptidoforms at nominal 1% (10,928 half-tryptic); MSstatsLiP: 20/20 significant TrP-adjusted SNCA peptides lie in the fibril core 37-99. No entrapment on the semi class yet. | 0.91 core-s vs Sage semi ~450 core-s | M5_lip_cf |
| CF D10cf | fraction-local transfer (default): 0% off-peak values; naive global MBR 9% off-peak, but highest CORUM PR-AUC on shared proteins (0.088 vs 0.082), so the design hypothesis is only partly supported. Null set underpowered (13-17 profiles). | - | M5_lip_cf |

Homo-family XL pairs are reported as a separate class without FDR estimation (D10), because TT/TD/DD cannot estimate them; this costs sensitivity.

### 2.7 Ingest / MS1
Direct Thermo ingest 17.8 core-s/1e5 MS2 (E1, old node; 19.2 on the new node) vs 531 for TRFP->mzML (E1; cold-Lustre inflated, 235 on the node of record, hw_adapt); 6.5 in a persistent worker. MS1 finder vs biosaur2 (D0): recall 0.85, precision 0.65, 28-30 core-s/run vs 1,350-1,520.

## 3 Test-set results

**No held-out test-set result exists.** The single A5 look opened on 2026-10-28 at 14:23 CET (ledger 13:23:46Z; `milestone-M7-2026-10-28` @ 4fe2e2a).
- **T5** (arm M) completed with exit 0.
- **T1** reached its final scoring stage (00:20), but the compute node was wiped before it finished. It was closed with `sealed_run.py --recover` as 'failed (unterminated)': outputs deleted unread, no number seen (results/SEALED_LEDGER.md, 742ba48).
- **Retry:** under A5 a retry re-runs both T5 and T1 from the same frozen commit; at most two retries are allowed, and none has been used, so up to two retries remain.
  - The retry needs the node-local mokapot env rebuilt. The core versions match, but the frozen digest does not (D38). The rebuilt env's full package list is now recorded (results/2026-10-30_mokapot_env_rebuild).
  - It then needs a dev check that the comparators stay within the D33 band, plus the look itself; D39 puts the total at about 13 h of uninterrupted node time.
  - This session ended before that. T1–T5 remain sealed and unseen.

**HN1/A6:** not tagged: infeasible (D34). No null-spectrum gate result exists for this session beyond A1 (FAILED) and A2 (not evaluable).

T2, T3 and T4 remain sealed (D36), and no DIA result exists; see section 7.

## 4 Statistical status

| Gate / attempt | Applies to | Outcome (literal) | Data |
|---|---|---|---|
| G0 label-swap/permutation bit-identity | v1, v2_sep (at final freeze 06ce624), ce_s, IT-CID (pinned store), XL, semi | pass | dev |
| G2-original (PSM-level differential) | v1, v2, ce_s, IT-CID | pass; nearly vacuous (winner changes on 0.08% of nulls). **Caveat (D34):** the cluster percentile bootstrap undercovers with few clusters (0.923–0.944 at J=37), so these passes are weaker than nominal | dev nulls |
| **A1 tail-G2** | v1 | **FAILED** (n_tail 751; target fraction UCB 0.742 > 0.60; diff UCB 0.124 > 0.10); contamination explanation unconfirmed | DN, 6 null projects, dev |
| **A1 G1** | v1 | **FAILED** (rho 0.048 [0.019, 0.076]) | DN, dev |
| **A1 G3** | v1 | **FAILED** on clean-null subset (one bin 1.03 pp, ~100 per bin); passes on all nulls | DN, dev |
| G4 (decoy KS across strata) | v1; v2_sep partition | pass (weak; tau 0.54 / 0.298) | training decoys, dev |
| G5 (leak AUC) | v1, v2, ce_s | pass (0.504-0.506) | dev nulls |
| Arm M (yeast + human entrapment) | v1; v2_sep | inconclusive, not violated (1.05% / 1.04% at nominal 1%) | D7, dev |
| **A2 / DN2** (fresh nulls, tail-G2/G1/G3) | v1 (+v2_sep, ce_s via A3) | **A2 not evaluable**: 6 projects/18 runs vs >=10; \|A\| 24,581 vs 120k; tail proxy 22 vs 400; step 4 not run, no gate quantity computed | DN2 |
| **DNIT** (IT-CID nulls) | IT-CID scorer | **DNIT futile**: 581 usable spectra vs 120k; tail-G2/G1/G3 not evaluable | DNIT |
| A3 DN2 tiers for v2_sep | v2_sep | **expected not evaluable** by power check; step 5 never run; v2_sep has never passed a null-spectrum gate | synthetic/dev |
| A3 MC (sec 6.7) | v2_sep lists, union, picked protein | pass in all merged-stratum lists (1e4 reps). This verifies the estimator on simulation, not the scorer. | simulation |
| **A4 / DN3** | v1, v2_sep | **A4 not run** (drafted, reviewed twice, not adopted; P(min) 0.08-0.43) | - |
| C9 pilot | design | GO by simulation; not a gate | dev, label-free |
| A5 | T1/T5 look | adopted; relaxes claim conditions; single look opened 2026-10-28 14:23 CET (see section 3) | sealed T5, T1 |
| **A6/HN1** | v1 only | **not tagged: infeasible** (tag rule: P(v1 minimum) 0.63 < 0.8 at measured throughput); pool unspent | human removed-peptide nulls (dev dry run only) |
| Entrapment verdicts, per-run/global | all | no arm "violated" on any dev set for v1/v2_sep; arm P mostly 'inconclusive' at 1% (estimator overstates by ~20%); natural arm C tilts 54% | D0, D1, D8, D8h |
| Sage chimera rank>=2 (comparator) | comparator | 'violated' on D8h (5.3%) | D8h |
| E8 FLR (Ala) | PTM | point estimates calibrated or conservative on dev; under the stricter T2 upper-bound rule the dev gate is not passed (PXD009449 inconclusive at 1%); sealed T2 not run | D2, dev |
| MBR transfer decoy | quant | under-counts ~2.4x on D12 Lumos; MBR opt-in | D12 |
| XL | XL | heteromeric true FDP ~1.2% at nominal 1% on D3b; self pairs 9% at nominal 5% (D3) | D3, D3b |

Net: **FDR claims for v1, v2_sep and IT-CID are not gate-supported.** v1 FAILED the only completed null-spectrum test; no null-test of v2_sep or IT-CID exists. Claims rest on entrapment (dev) and the A5 labels.

## 5 Negative results
- **E5 bi-encoder retrieval: NO-GO.** Open-search recall@1000 79.5% (shifted 57.2%), no faster than the inverted index (0.28 vs 0.24 ms/spectrum), target/decoy asymmetry at K=10.
- **E6 cross-encoder: NO-GO as default** (0.983x v1 on D8h; -24% at 0.1%). Hybrid ce+static not run.
- **Group-walk over chimeric strata:** stratum >=3 violated on D0 (3.1%); recorded final; v2_sep per-stratum TDC replaced it.
- **IT-CID nulls (DNIT) futile**; **A2 not evaluable**; **A4 not run** (section 4).
- **Spectrum-level neural scoring** adds nothing over static predicted features at our operating point (D19).
- **19 richer label-blind v2 features:** no held-out gain (-0.13% linear).
- **Per-run score calibration (G4-driven):** costs 13-72% of IDs, not adopted.
- **Tight tolerances (E0):** -0.8 to -5.8% IDs with Sage; recalibration alone gives 0%.
- **Decoy realism (3g):** natural foreign nulls beat reversed decoys 54.4% (frozen score) vs 51.6% (hyperscore); formal verdict 'mixed'; predictor naturalness (H2) not supported; real homology/plants carry part. No remedy adopted.
- **Chimeric protein depth:** vs the fair comparator no protein-group gain on D8h.
- **Global MBR for CF-MS** adds artefacts (9% off-peak) though raises PR-AUC; fraction-local stays default.
- **Open search** loses ~1,400 narrow stems vs closed and 4% vs Sage-open.
- Mokapot 0.10 ignores `--seed`; Koina predictions are not bit-deterministic (2.8% entries differ after rebuild): both handled by wrappers/pinned store.

## 6 Follow-ups
1. **Retry the T1+T5 look** from 4fe2e2a after the D1 mokapot-env check (D38/D39; about 13 h of uninterrupted node time).
2. **T2 look** (driver ready, about 20 min on an idle node); then T3 and T4 on a quiet node.
3. **A6/HN1 on a quiet node.** It needs ≥ 192 admitted runs/h or an earlier start, and the G2 bootstrap must first be replaced, because the percentile type undercovers at J ≈ 37.
4. **Deterministic localization tie rule:** re-validate it on dev under a new pre-registration before it becomes the default.
5. **Merge E12** (abundance/PTM speed-ups) after the looks; make GPU MS1 the default after a repeated bitwise-equality check.
6. v2_sep null-spectrum test with multi-precursor nulls (needs its own amendment); IT-CID nulls designed for unit-resolution MS2.
7. Speed: decode (5 core-s/run) and protein stages are next after GPU MS1; GPU isotope grouping; MS1 shared once with LFQ; phospho label build (446 s); CONGA vectorisation; multi-core directLFQ.
8. XL static scorer and top-K choice (self pairs); semi-class entrapment for LiP; spiked-foreign CF set.
9. Open-mode narrow-stem displacement; predicted-intensity localization evidence.

## 7 Limitations
- **Dev vs sealed.** All numbers above are dev; several (D0/D1 scorer selection, PXD000138/PXD009449 localization settings, M6d quant choices on D9) are tuned on those sets. D0/D1 are a selection set.
- **FDR is not gate-validated** (section 4). Entrapment estimators overstate by ~20%, so many verdicts are 'inconclusive'; natural arm C sits at 1.0-1.2% at alpha 1%.
- **DIA arm not done.** Only the D5 data are present; no DIA result.
- **XL sensitivity below xiSEARCH** (~73% hetero on D3b; self pairs weaker on D3); homo-family pairs not FDR-controlled.
- **MBR decoy under-count** (~2.4x on D12 Lumos; exact on D9 Astral); MBR stays opt-in.
- **Shared-node timings.** Load 14-110 on all measurements; some comparator timings were unlocked (xiSEARCH 2 h hold, Sage in 3c) and are indicative. Cohort throughputs are page-cached unless streamed; PRIDE streaming is ~580 runs/h, download-bound. Per-session numbers (E7b 5,433 closed vs M2b 20k+ vs scoring-speed 5,658) are not comparable across rows.
- **Frozen path vs reported speeds.** Throughputs quoted for chimeric v2_sep and E11 GPU MS1 are for opt-in/candidate paths, not the frozen v1 T1 path.
- **Comparators:** MSFragger and DIA-NN were not run (MSFragger unlicensed); XL comparator on D3b is a subset extrapolation; E7b's Sage chimera comparator is narrow-precursor only.
- **Scope:** ion-trap and Astral are out-of-class for the OT-HCD model (Astral scored out of class on D9). Transfer to Bruker, Exploris, Astral DDA scoring is untested; human paralog 'neighbours' unmeasured.
- PXD007058 has only a single-engine proxy truth; T2, T3 and T4 remain sealed.
- **Localization is not run-to-run reproducible** (GPU atomic ties, about 0.7% of PSMs on PXD009449). Dev verdicts are stable across 11 repeats, but individual site calls are not.
- **Weak G2-original.** It uses a percentile cluster bootstrap that undercovers at about 37 clusters.
- **Concurrent-agent costs.** Several agents shared one node and one GPU. /tmp filled to 100% three times, killing runs, and OOM events forced retries. Timing numbers carry this load.

## 8 Remaining ideas
Hybrid cross-encoder + static features (new pre-registration); abundance-aware quant window; sibling values at peptide level only; clustering arm; DIA arm; per-rank-controlled Sage chimera comparator; Astral-class scorer and paired entrapment on D9; GPU GBM and a seed-3001 calibration decoy for global tails; spiked-proteome CF null; bi-encoder only if a faster ANN index than the inverted index is shown.
