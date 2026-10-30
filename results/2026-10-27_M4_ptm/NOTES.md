# M4: PTM / mass-offset / open search, localization operator, E8 on dev (D2)

Dev data only (DESIGN §6.6). Sealed T2 is not touched: no file under `data/sealed/` and no file under
`data/truth/sealed/` is read by any M4 script. The PXD000138 library FASTA lives in `data/truth/sealed/`, so it is
not used. Searches use the dev per-pool FASTAs (`data/truth/PXD000138/pool_*.fasta`) plus human background.

## 0. E8 truth definitions (written and committed BEFORE any localization result was computed)

These definitions are fixed now. If one changes, a dated note goes here and in RESEARCH_TREE.md before results are
looked at again. Estimators, labels and metrics follow DESIGN §6.5 and `vulcan/stats/flr.py` unchanged.

Common rules
* Unit: a PSM that passes peptidoform-level paired double competition at q ≤ 1% (localization-agnostic, DESIGN L4),
  is the rank-1 PSM of its spectrum, and carries ≥ 1 phospho. Localization is then done by the L3 operator.
* A reported site is *true* iff its residue position equals the truth position for that peptide. Sites of PSMs whose
  peptide has no truth entry are excluded from true FLR (they still count in FLR̂ "all"); FLR̂ is reported both on all
  localized PSMs and on the truth-eligible subset, so FLR̂ vs true FLR is compared on identical site sets.
* Mono- and multi-phospho are reported separately (§6.5). The E8 calibration table uses mono-phospho.
* §6.5: pools containing positional isomers of one sequence are excluded; checked programmatically below (a sequence
  with two different truth positions is dropped).

PXD000138 dev pools 58 (E10), 60 (E12), 66 (F6), 72 (F12): **design truth**
* Each pool entry is a concatenation of equal-length synthetic units (backbone + variable positions). The unit length
  is the shortest period P such that the concatenation splits into units with identical residues at every position
  except a few variable positions. The *acceptor slot* is the variable position whose residue set is exactly {S,T,Y};
  the other variable positions take ~20 residues (including S/T/Y confounders, by design).
* Truth: a PSM is truth-eligible iff its stripped sequence equals one unit exactly; the true phospho site is that
  unit's acceptor slot. (Junction peptides spanning two units are artifacts of the concatenation, not synthesized.)
* If the design inference does not give exactly one {S,T,Y} slot per pool, the pool is dropped from E8 (recorded).

PXD009449 dev pool Ymod_Phospho (BA4, 3xHCD file): **kit-design truth**
* The ProteomeTools PTM kit carries a phospho-Tyr on every modified peptide. Truth-eligible peptides: sequences in
  the submitter MaxQuant `evidence.txt` of this pool (any raw file of the pool) that contain **exactly one Y**; the
  true site is that Y. Peptides with 0 or ≥ 2 Y are excluded from true FLR.

PXD007058 (Ferries 2017, pools 1–4 HCDOT rep1): **proxy truth (declared here before use)**
* No design table is available (closed-access SI) and the submitter files are a single engine (MaxQuant 1.5.8.0,
  `txt/` per pool), so the "≥ 2 engines" idea in DATA.md is not possible. The proxy is:
  1. From the pool's `txt/msms.txt`, take MaxQuant PSMs with exactly one phospho (STY), not reverse, not contaminant.
  2. Per (pool, stripped sequence): the MaxQuant site = argmax of `Phospho (STY) Probabilities`. The sequence is
     truth-eligible iff (a) ≥ 1 PSM has max site probability ≥ 0.99, (b) ≥ 90% of its mono-phospho PSMs with max
     probability ≥ 0.75 agree on the same argmax site, and (c) no second site of that sequence reaches max
     probability ≥ 0.99 in any PSM (positional-isomer exclusion, §6.5).
  3. Exclude sequences whose MaxQuant `Proteins` map only to wells of the sealed PXD000138 T2 pools
     (C1, G4, B12, F1): the Ferries library IDs carry Marx-library well suffixes, and we keep T2-related sequences out
     of dev entirely.
  4. Search DB for this set: human canonical + contaminants + one FASTA entry per stripped sequence in the four
     submitter `peptides.txt` files (the library FASTA `SF_080317_Phospholibrary_v2.fasta` is not public).
* Known bias: the proxy keeps only sites MaxQuant localizes confidently, so it favours easy sites and true FLR
  measured against it is optimistic. PXD007058 numbers are a development curve and never gate on their own; the
  E8 gate reading uses PXD000138 and PXD009449 (design truth) first.

---------------------------------------------------------------------------------------------------------------------
(Everything below was written after §0 was committed.)

## Hardware and conditions
gpu-node, EPYC 9354 (32 affinity cores), A100 80GB PCIe, driver 580.178.04, CUDA 12.9 (fingerprint in
`metrics.json:hw` and in every `timings.jsonl` / `stats_*.json`). Shared node (loadavg 14–30, M3 and baselines agents).
All timed runs (ours and Sage) were taken under `flock /tmp/vulcan_bench.lock`; the phospho index CPU build (446 s)
was not (it ran next to other jobs). Ours: from the Arrow store, warm (2nd pass in one process); Sage: process wall
from GNU time on an MGF of the same store. No entrapment in any M4 search.

## 1. What was built
* **Variable-mod index (L1).** `vulcan_core` (Rust) now supports a variable peptide-N-term mod (Sage `^`) carried in
  mask bit 63 (`NTERM_BIT`); mass, canonical id and b/y follow Sage (b starts at the N-term delta, y at mass − delta).
  STY phospho + M-ox work through the existing one-mod-per-residue table. Engine row packing went from 27 to 28 bits
  (up to 268M peptidoforms with `chunk_size` 32768). Closed search on D0 01A is **bit-identical** to the E2 PSMs.
* **K2 shift mode** (`vk_k2_shift`): per-candidate delta = observed − candidate mass; every b/y ion is matched
  unshifted and shifted, and the delta is placed on the single residue that maximises the hyperscore (prefix sums,
  O(L) per candidate; MSFragger-style localization-aware open scoring). The first version (union of shifted and
  unshifted, no localization) lost narrow IDs to inflated open candidates: open 19,562 → 20,488 stems after the fix.
* **Modes** (`vulcan/ptm/search.py`): closed hypotheses + one hypothesis per mass offset (33 offsets: Unimod list
  + D0 pilot-histogram peaks, `scripts/m4/offsets.py`) + one open hypothesis per spectrum×charge (delta in
  [−150, +500] Da). All compete in one top-K / one ranking, label-blind. Open K1 runs in spectrum batches with
  per-batch top-K (exactly the one-shot top-K) and a bounded append buffer; the K1 kernel is unchanged.
* **CONGA** (`vulcan/ptm/conga.py`): top-5 PSMs/spectrum → label-blind neighbour filter (>5% shared b/y) → stem
  competition (paired double competition, unit = stem) → Group-walk over narrow + 1-Da shift bins merged by the
  pre-registered `merge_adjacent` (≥200 label-blind pair competitions). Pooled TDC is reported alongside.
* **L3 localization operator** (`vulcan/ptm/localize.py`): exhaustive C(n,k) configurations over S/T/Y + decoy
  residue (A primary, L/I secondary), b/y at z ≤ min(z−1, 2), −H3PO4 on every fragment holding a relocalized mod for
  every acceptor alike; the scorer receives only masses/positions (no letters). Site score = S(best) − best S of a
  configuration without that site. GPU (torch) batched.
* **Tests:** `tests/ptm/test_localize.py` (exhaustive enumeration, true-site recovery, **decoy-AA G0 cases**: A vs S
  with equal masses → bit-identical config scores; letter permutation among acceptors → identical; neutral loss
  applies to A/Y exactly as to S), `tests/gpu/test_ptm_modes.py` (var-mod masses/masks = Sage incl. `^`; K1 = f32
  brute force on the var-mod index; **G0 role swap + permutation bit-identical for closed, mass-offset and open**).
  `tests/index`, `tests/gpu`, `tests/ptm`: 19 passed after the last kernel change.

## 2. Phospho / variable-mod closed search (task 1)
Index d2_phospho_prod: human SwissProt + contaminants + D2 dev libraries; STY phospho, M-ox, pep-N-term acetyl,
≤3 var mods, trypsin 2 missed. **184.5M peptidoforms** (92.3M targets; 26× the closed index), 5.17G postings,
**34.3 GB resident on the A100** (10 slices), GPU index build 6.1 s, L1 load 2.6 s, CPU build 446 s (unlocked;
321 s of it is Python label-string assembly, an easy target). Peak GPU memory 57.8 GB.

| Search (warm) | D0 01A (41,820 MS2) | D0 01B (46,816) | per-run core-s |
|---|---|---|---|
| closed, M-ox index (7.06M) | 0.21 s | 0.14 s | ≈ wall |
| phospho index (184.5M) | 0.20 s | 0.21 s | ≈ wall |
| D2 9 runs, phospho index | 0.65 s total (Ymod 55k MS2: 0.17 s) | | |

Narrow windows keep the phospho search as fast as the closed one: index size costs memory, not time.

## 3. Mass-offset and open search on D0 (task 2)
| | wall (2 files) | core-s | stems q≤1% Group-walk | pooled TDC | narrow / shifted | 1-Da shift bins (≥5 stems) |
|---|---|---|---|---|---|---|
| ours closed (CONGA pipeline) | 0.35 s | 0.35 | 17,714 | 17,714 | 17,714 / 0 | 0 |
| ours mass-offset (33 offsets) | 2.7 s | 2.7 | 20,188 (+14.0%) | 13,015 | 18,500 / 1,688 | 10 |
| ours open [−150, +500] | 17.2 s | 17.2 | 20,488 (+15.7%) | 7,832 | 16,270 / 4,218 | 17 |
| Sage open, same pipeline, 32 thr | 72.0 s | 1,153 | 21,400 | 13,579 | 17,538 / 3,862 | 10 |
| Sage open, 4 thr | 155.9 s | 574 | (same PSMs) | | | |
| Sage's own LDA q (rank 1) | | | 24,834 peptides | | | |

* Speed: open search is **4.2× faster in wall and ~65× fewer core-s than Sage-32** on the same node; mass-offset
  mode is 27× faster in wall than Sage open. Literature (not run; MSFragger unlicensed): MSFragger open 500 Da on
  the same 41,820-spectrum HEK file took 5.4 min on 4 cores (Kong 2017); on PXD001468 Sage open 44 s vs MSFragger
  225 s on 32 vCPU, 24,682 vs 24,911 PSMs (Sage paper/blog). Hardware differs, so these are context only.
* IDs: Group-walk is essential (pooled TDC loses 35–62% of the gain). Ours open finds more distinct offsets
  (16, 43 carbamyl, 128/−128 K, −17, 156/183, −89, 32, 80…) and fewer ±1/±2 Da (we handle isotope errors in narrow
  hypotheses), but **4% fewer stems than Sage-open** through the identical pipeline: ours loses ~1,400 narrow stems
  when open candidates join (17,714 closed → 16,270 narrow in open). Mass-offset mode keeps narrow IDs (18,500).
  Sage's own LDA number is not comparable (no CONGA, no stem competition, learned score).

## 4. E8 on dev (tasks 3–4), `e8.json`
IDs: localization-agnostic peptidoform paired double competition at q ≤ 1% (hyperscore). Localized phospho PSMs:
PXD000138 3,858 (3,827 mono; 2,277 truth-eligible), PXD009449 Ymod 45,845 (26,483 eligible), PXD007058 506 (444
mono; 194 proxy-eligible). Truth reports (design inference: 1 {S,T,Y} slot per pool, 0 positional isomers;
183 single-Y Ymod peptides; 72 proxy sequences, 6 sealed-well sequences excluded) are in `e8.json:truth_reports`.

**Calibration, mono-phospho: FLR̂ threshold from all localized PSMs → true FLR of the STY sites accepted among
truth-eligible PSMs** (n = accepted eligible STY sites)

| Set | decoy | ours @FLR̂ 1% | ours @FLR̂ 5% | pyAscore @1% | pyAscore @5% |
|---|---|---|---|---|---|
| PXD000138 | A | 0.29% (1,739) | 3.7% (2,155) | 0.0% (1,046) | 2.0% (2,161) |
| PXD000138 | L | **1.31%** (2,131) | 3.6% (2,228) | 0.89% (1,919) | 2.1% (2,226) |
| PXD009449 | A | 0.87% (23,678) | 4.4% (25,924) | 0.38% (22,460) | 4.9% (25,826) |
| PXD009449 | L | 0.87% (23,429) | 4.5% (26,039) | 0.54% (22,881) | 0.54% (22,881) |
| PXD007058 (proxy) | A | 0.0% (138) | 2.3% (172) | 2.4% (165) | 2.3% (172) |
| PXD007058 (proxy) | L | 0.0% (87) | 1.7% (173) | 0.0% (84) | 2.4% (168) |

Whole-list (no threshold) FLR̂ vs true on eligible: A: 5.3 vs 3.7% (138), 5.6 vs 4.4% (9449); L: 2.6 vs **3.6%**
(138), 3.9 vs **4.5%** (9449). **The Ala estimator is conservative-to-calibrated everywhere; the Leu estimator
under-estimates on both design-truth sets** (1.31% at nominal 1% on PXD000138). Reading: A passes E8 on dev
at 1%/5%; L does not (by a small margin). Not a claim; the gate itself is on sealed T2.

**Site recall at matched true FLR** (correct sites / eligible PSMs, same PSMs for both tools)

| Set | decoy | ours @1% | pyAscore @1% | ours @5% | pyAscore @5% |
|---|---|---|---|---|---|
| PXD000138 | A | 0.880 | 0.874 | 0.911 | 0.930 |
| PXD000138 | L | 0.890 | 0.855 | 0.943 | 0.957 |
| PXD009449 | A | 0.886 | 0.866 | 0.936 | 0.927 |
| PXD009449 | L | 0.877 | 0.859 | 0.939 | 0.859 |
| PXD007058 | A | 0.758 | 0.613 | 0.866 | 0.943 |
| PXD007058 | L | 0.737 | 0.562 | 0.876 | 0.918 |

Ours is ahead at 1% on all six cells (+0.6 to +17.5 points), behind pyAscore at 5% on 4 of 6 (pyAscore's best-site
accuracy is higher on PXD000138/7058; its score has many ties at 0, which hurts it at strict thresholds). The plain
identification search's own rank-1 isoform (STY only) is correct for 97.1% / 93.8% / 90.7% of eligible PSMs, so
the A/L candidates cost accuracy by design (adjacent A/S with only tiny b-ions between them, e.g. E10 `A[STY]X…`).

Dev tuning disclosure: `peaks_per_window` {4, 8, 12, 20, global-150} × weight {rank, intensity}, neutral loss on/off,
max fragment charge {1, 2} were tried on PXD000138 + PXD009449 (truth used); chosen: 12 per 100 Th, intensity
weight (mean recall at true FLR 1% over the 4 cells 0.883 vs 0.874 for the initial rank/8). PXD007058 was not used
for tuning. pyAscore was run at its documented defaults (bin 100, top 10, 0.02 Da) with the same STYA/STYLI groups
and the same residue-blind neutral loss.

## 5. Runtime per stage (task 5)
* Ingest (store): D2 RAW → store with vulcan.io.convert, 4 threads, unlocked (`/tmp/vulcan_m4/logs/convert_*.json`).
* Index: phospho CPU build 446 s (digest 4.7, expand 21.4, arrays/labels 321); GPU index build 6.1 s; L1 load 2.6 s.
* Search: table §2/§3; per-phase GPU ms in `stats_*.json` (open: K1+top-K 8.3 of 8.6 s per file; K2 shift ≈0.1 s).
* CONGA (CPU, Python): neighbour filter 5–8 s per D0 dataset (Python set intersections; vectorise next).
* Localization: ours 5.8k–23k PSMs/s (GPU scoring + per-PSM numpy prep, which dominates), pyAscore 7k–16k PSMs/s
  (1 CPU thread). Peptidoform FDR + spectra load: < 15 s per accession.

## 6. Open issues
1. Open-mode narrow loss (−1,400 stems vs closed) and −4% vs Sage-open: candidates' shifted scoring still beats some
   correct narrow PSMs. Options: score open candidates with a delta-mass prior feature in the M3 static scorer, or
   keep narrow and open as separate CONGA groups before the spectrum competition (CONGA's original two-pass form).
2. Leu decoy estimator is mildly anti-conservative on dev (1.31% at 1%, 3.6% vs 2.6% whole list on PXD000138; I/L
   sit at the variable X positions of the Marx design). Ala stays the primary per §6.5; keep L as secondary only.
3. N-term acetyl is peptide-N-term (Sage `^`): protein-N-term needs digest positions; it inflates the index
   (184.5M rows). Engine limit is now 2^28 rows.
4. Predicted-intensity evidence for localization (DESIGN L3) is not added: D4 bars AlphaPeptDeep/DeepLC on D2;
   Prosit/MS2PIP via Koina is the allowed route.
5. Localization host prep is per-PSM numpy; Python label assembly dominates the phospho index build.
6. PXD007058 proxy truth is single-engine and small (72 sequences, 194 PSMs). The declaration is in §0
   (commit 97daf98) and still needs to be recorded in RESEARCH_TREE.md by the lead.
7. Proposal (NOTES only, §6 unchanged): E8 reads FLR̂ at the PSM level with PSM-level acceptor counts; §6.5 does not
   say whether N_STY/N_A count per PSM or per distinct peptide. Recommend fixing PSM-level in a dated note before T2.
