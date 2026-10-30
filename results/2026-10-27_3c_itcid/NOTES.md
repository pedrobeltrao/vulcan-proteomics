# 3c: IT-CID instrument-class static-first scorer (dev only)

Status: done (dev only; no sealed file touched; §6, metrics, labels and estimators unchanged; OT-HCD frozen scorer untouched).
The two "pre-stated" sections below were committed before any IT-CID model was trained or applied to dev data (commits 92d7109, 6615a53).
Files: `metrics.json` (everything), `eval.json` (D10cf), `eval_D4.json`, `sage_eval.json`, `gates.json`, `check_da.json`, `koina_repeat.json`, `bench.json`, `timings.jsonl` (Sage/mokapot), `models/`, `configs/`, `config.yaml`, `git_sha`. Code: `scripts/m3c/`, `vulcan/predict/koina.py`, `vulcan/rescore/itcid.py`, Da mode in `vulcan/gpu/{kernels.cu,host.py,engine.py,_lib.py}`.

## Pre-stated model-choice rule (fixed 2026-10-27, before training finished and before any dev evaluation)
- Candidates: the four frozen variants trained by `scripts/m3c/train_static_itcid.py` with the M3 recipe:
  `itcid_linear_base`, `itcid_gbm_base`, `itcid_linear_pred`, `itcid_gbm_pred`.
- Criterion: **training CV only**. For each variant, the leave-one-project-fold-out held-out count of rank-1 target PSMs at
  PSM q <= 1% (TDC on the held-out training runs), summed over the 3 folds, at that variant's CV-selected hyper-parameter
  (`cv_heldout_psms_at_1pct` in `models/training_log.json`).
- Eligibility: the variant's held-out count must exceed hyperscore's on the same held-out runs.
- Choice: the eligible variant with the largest held-out count. Tie-break: any variant within 0.5% of the best counts as
  tied; among tied variants prefer linear over GBM (lower capacity), then base over pred (fewer inputs).
- The chosen model is frozen (sha256 in `metrics.json` and `models/*.json`) and is the only IT-CID scorer whose dev numbers
  are headline. The other three are exploratory and never replace it on the basis of dev results.

## Pre-stated evaluation (dev only: D4 = PXD016766 TMT11 SPS-MS3, all 12 fractions; D10cf = PXD015406 SEC, 29 fractions)
- D1 (PXD001077) is not IT-CID: its store has 15,374/15,374 MS2 as FTMS + HCD, so it is excluded.
- IDs at q <= 0.1 / 1 / 5 %: peptidoforms (paired double competition, pooled over the fractions of each set as one
  experiment) and picked protein groups, for hyperscore and the chosen scorer (others exploratory); Sage 0.14.7 LDA and
  Sage+mokapot 0.10.0 by their own peptide q (I/L-collapsed sequences, M1e units) on identical spectra.
- Entrapment §6.3: arm P (seeds 2001-2003, pseudo-reverse decoys) and arm C (Arabidopsis, r from M1e), CI convention (a).
- Gates not depending on the pending A2 null set: G0, G2-original, G4, G5. Arm M only if an IT-CID yeast run exists among
  dev (none does: D7 is OT-HCD).

## Rule applied (training CV only; recorded and committed before any dev evaluation)
Held-out rank-1 target PSMs at q <= 1%, 3 project folds, best hyper-parameter (hyperscore on the same held-out runs: 84,091):

| variant | hp | CV held-out | refit train |
|---|---|---|---|
| itcid_linear_base | C 0.01 | 110,727 | 113,041 |
| itcid_gbm_base | n 150 | 106,480 | 106,393 |
| **itcid_linear_pred** | C 0.01 | **127,081** | 127,760 |
| itcid_gbm_pred | n 300 | 125,542 (-1.2%, not tied) | 126,972 |

**Chosen and frozen: `itcid_linear_pred`** (sha256 667ce39cb653026baed40bb18b25f5cbe196b192e2d5a7f2f8e7aa3618c23b89).
The other three are exploratory.

## What was built
- **Engine (IT-CID support).** The K1/K2 kernels took ppm tolerances only (M6 used +-600 ppm as a stand-in). New
  `SearchParams.fragment_unit = "da"` -> Sage `Tolerance::Da` bounds (m + tol) in K1, the K1 brute-force reference and K2
  (`tol_bound`, a `__constant__` flag; ppm arithmetic unchanged op for op). Low-res binning = coarser K1 log bins
  (`bin_ppm` 200; a pure layout knob: PSMs identical to 20 ppm bins, K1 3.2x faster).
  - OT-HCD unchanged: D0 (2 runs) and D1 top-5 PSM tables re-searched **bit-identical** to M3.
  - Da mode on D4 A1 (`check_da.json`): K1 == brute force on 2M candidates (0 mismatches); top-K(500) recall of the
    exhaustive best target/decoy 100% / 99.8%; Sage oracle (M1e D4 run, TRFP mzML, +-0.5 Da): same top-1 peptide on
    **99.7%** of Sage's q <= 1% target scans (54% over all scans: low-score noise matches, centroids differ).
- **Predictor.** AlphaPeptDeep has no CID/ion-trap MS2 model (HCD instruments only), so Prosit 2020 CID (label-free) and
  Prosit 2020 TMT with fragmentation "CID", CE 35 (TMT), plus Prosit iRT / iRT-TMT, via Koina (`vulcan/predict/koina.py`,
  cached by the existing content-keyed PredCache). D4-compliant; training projects of these models are not in any 3c set.
  Peptides > 30 aa / non-standard residues: `pred_missing` = 1 (1-4% of PSMs; up to 18% in the G2 mixtures).
- **Features** (`vulcan/rescore/itcid.py`, a subclass so the OT-HCD/scoring-speed path is not edited): fragment matching
  +-0.5 Th, TMT N-term on b ions, fixed CE (no NCE grid), missing-prediction handling, the same 9 predicted features + `pred_missing`.
- **Training data** (DATA.md, role train; `data/manifests/train_itcid*.tsv/.sha256`, fold map committed): 13 runs / 13 PRIDE
  projects (8 human, 5 mouse): 6 TMT SPS-MS3 (Fusion/Lumos), 7 label-free (LTQ Orbitrap / Velos / Elite), selected by a
  seeded rule committed before download; IT-CID verified on the spectra (20 walked candidates rejected: 17 had Orbitrap
  MS2, 3 TMTpro). Excluded: every manifest accession (dev, sealed, T4 selection, M3 training, predictor training, DN, D8,
  ...), non-human/mouse species, and (PI instruction) a mechanical DN2 guard: 0/58 overlap with every accession the A2 DN2
  walk has seen, no ALLOW-reachable organism. scripts/a2 was read, never run.
- **Training** (`scripts/m3c/train_static_itcid.py`): the M3 recipe unchanged (Percolator-style labels on training runs,
  3 iterations, leave-one-project-fold-out hp, refit, freeze).

## Dev results (all 12 D4 fractions pooled; all 29 D10cf fractions pooled; peptidoforms by paired double competition;
"seq" = I/L-collapsed sequences, the unit of the Sage rows, which use Sage's / mokapot's own peptide q on identical spectra)

| D4 (TMT11 SPS-MS3) | 0.1% | 1% | 5% | seq @1% | protein groups @1% | verdict 1% (P / C -> overall) | verdict 5% |
|---|---|---|---|---|---|---|---|
| hyperscore | 29,725 | 43,437 | 54,665 | 42,410 | 6,734 | P 0.94% sup / C 1.06% inc -> **supported** | inconclusive |
| **itcid_linear_pred (chosen)** | 56,134 | **60,873** | 64,827 | **59,082** | **7,858** | P 1.11% inc (cwc) / C 0.97% sup -> inconclusive | inconclusive (P cwc) |
| itcid_gbm_pred (expl.) | 56,308 | 59,343 | 63,071 | 57,521 | 7,832 | P 1.07% inc (cwc) / C 1.09% inc | inconclusive |
| itcid_linear_base (expl.) | 45,454 | 55,589 | 62,043 | 54,154 | 7,586 | inconclusive | inconclusive |
| itcid_gbm_base (expl.) | 45,287 | 54,127 | 59,733 | 52,590 | 7,547 | inconclusive | inconclusive |
| Sage 0.14.7 LDA (seq, own q) | 47,154 | 55,151 | 60,029 | 55,151 | 7,608 | P sup / C inc -> inconclusive | inconclusive |
| Sage + mokapot (seq, own q) | 47,906 | 55,606 | 60,051 | 55,606 | 7,662 | P sup / C inc -> inconclusive | inconclusive |

| D10cf (SEC CF-MS, label-free) | 0.1% | 1% | 5% | seq @1% | protein groups @1% | verdict 1% | verdict 5% |
|---|---|---|---|---|---|---|---|
| hyperscore | 25,864 | 40,547 | 53,666 | 37,435 | 3,933 | P 1.18% / C 1.30% -> inconclusive | inconclusive |
| **itcid_linear_pred (chosen)** | 63,407 | **72,160** | 78,451 | **65,786** | **5,394** | P 1.11% / C 1.19% -> inconclusive (no cwc: P lower CI 1.02%) | inconclusive |
| itcid_gbm_pred (expl.) | 64,374 | 69,071 | 73,890 | 62,604 | 5,417 | P 1.16% / C 1.18% -> inconclusive | inconclusive |
| itcid_linear_base (expl.) | 44,576 | 55,884 | 66,065 | 51,584 | 4,879 | P 1.17% / C 0.89% sup | inconclusive |
| itcid_gbm_base (expl.) | 45,742 | 53,534 | 62,417 | 48,884 | 4,850 | P 1.19% / C 0.94% sup | inconclusive |
| Sage 0.14.7 LDA | 45,054 | 55,640 | 63,795 | 55,640 | 5,009 | P inc / C inc -> inconclusive | inconclusive |
| Sage + mokapot | 44,273 | 55,447 | 63,647 | 55,447 | 5,019 | P inc / C inc -> inconclusive | inconclusive |

Reading (dev curves, not claims; the model was chosen on training CV only, so D4/D10cf were not used for selection):
- **Depth.** The chosen scorer gives +40% (D4) / +78% (D10cf) peptidoforms over hyperscore at 1%, and in sequence units
  **+6.2% (D4) / +18.6% (D10cf) over Sage+mokapot** (+7.1% / +18.2% over Sage LDA); +16% / +31% at 0.1%; proteins +2.6% /
  +7.5% over Sage+mokapot. M6's D4 gap (-24% vs Sage, hyperscore only) is closed. Caveat: our count uses the paired
  double competition, Sage's its own q; A1.9's strongest comparator (Sage + the same predicted features in mokapot) was not run.
- **Entrapment.** Nothing is violated anywhere (lower CI of LB <= alpha for every scorer, both sets, both alphas). At 1% the
  learned scorers are "inconclusive" with the usual ~1.1-1.2% paired estimate; D4 keeps "consistent with control" on P;
  D10cf does not (P 1.11%, lower CI 1.02%), and its arm C is 1.19% [1.07, 1.31]% (r = 1.02), i.e. above nominal by the
  upper-bound estimator; hyperscore is also inconclusive there (P 1.18%, C 1.30%), so the set itself (29 sparse SEC
  fractions pooled, few IDs per fraction) is hard for the estimators. Protein-level arm P: supported at 1% on D10cf,
  inconclusive on D4. This is the least comfortable number and is reported as is.
- **Model choice.** The rule picked linear_pred (CV 127,081 vs gbm_pred 125,542). On dev the two are within 2-4%; no switch.

## Gates (frozen IT-CID scorers; `gates.json`)
| gate | itcid_linear_pred | notes |
|---|---|---|
| G0 label blindness, fresh prediction cache per side | LF (D10cf): **pass**; TMT (D4): **FAIL** (strict) | PSM tables bit-identical on both; on TMT 5 fragment features differ by <= 3e-5 because Koina returned different floats for 38 of 30,748 identical TMT peptidoforms (<= 2.4e-4, server-side; same inputs, label-independent). 2 x 20k repeated queries in identical batches gave 0 differences (`koina_repeat.json`), so the nondeterminism is load/batch dependent |
| G0, shared content-keyed prediction cache (predictor queried once per peptidoform, as production) | **pass** LF and TMT | everything bit-identical incl. all 4 scores |
| G2 original (class IT-CID) | **pass** | n = 504,807 null competitions (449,105 offset, 55,702 yeast after the homolog filter removed 10,731) >= budget 107,570; overall diff -0.0018 (UCB -0.0007), top decile UCB -0.0014; passes by source; winner changed on 21% of nulls. All 4 variants pass |
| G4 (training decoys, charge strata) | **pass** | tau 0.324, KS vs rest 0.18-0.25; all variants and hyperscore pass |
| G5 leak AUC on G2 nulls | **pass** | 0.5013 [0.4997, 0.5029], top decile 0.502 (hyperscore 0.5026) |
| arm M | **not run** | no IT-CID yeast run exists among dev (D7 is OT-HCD); the yeast IT-CID runs fetched here (DN3c) are G2 nulls only |
| tail-G2, G1, G3 (A1) | **pending** | depend on the A2 fresh null set; note the A2 DN2 walk selects OT-HCD projects, so an IT-CID null set for the IT-CID class is still needed |
| A1.9 strongest comparator | **not run** | Sage + Prosit-CID features in mokapot |

So, per A1.2/A1.4/A1.8, FDR claims for the IT-CID scorer stay blocked (tail-G2, G1, G3 and arm M pending; G0 strict fail on TMT with fresh caches).

## G2 null set DN3c (deviations)
- Yeast IT-CID runs, seeded rule committed before download; PI check: no DN2 candidate used (guard, 0/58 overlap;
  S. cerevisiae is excluded from DN2 in config/a2.yaml). Pre-specified 3 projects; **only 2 accepted** (PXD008785,
  PXD010849; 3 walked projects had Orbitrap MS2) because the next one (PXD002607) stalled at ~13 kB/s on PRIDE; decided
  before any null was searched. Hosts: the 29 D10cf fractions (label-free; TMT-class nulls not covered).
- Offset nulls dominate (89%): unlike OT-HCD (595 usable in M3), n + 0.5 Da offsets of IT-CID host spectra mostly find
  candidates within 10 ppm here. Not investigated; G2 passes on each source separately.

## Runtime (A100 + 32 Zen 4, shared node, load 25-35)
- Production config, persistent process, warm prediction cache, under `flock` (bench.json; median of 4 runs, setup excluded):
  **D4 0.92 s / 0.91 core-s per fraction** (20k MS2: search 0.11 s, features 0.62 s, linear apply 0.03 s);
  D10cf 0.43 s / 0.20 core-s per fraction (3.4k MS2). First run after process start pays the prediction-cache load (10-12 s).
- First contact (cold cache, unlocked): Koina prediction dominates, 16 s (D4) / 17 s (D10cf) per fraction, network-bound
  (~4-6k predictions/s over 8 HTTP streams); a local CID predictor or a precomputed library would remove this.
- Sage 0.14.7 at 16 threads, unlocked, from MGF: D4 68 s / 643 core-s (12 files), D10cf 438 s / 4,138 core-s (29 files);
  mokapot adds 34 / 79 and 299 / 528 (s / core-s). Not a locked comparison: indicative only.

## Caveats / open
- Koina is a remote service: throughput and float determinism are outside our control (G0 above). Production should
  precompute a content-keyed library (one prediction per peptidoform) or host the Prosit CID ONNX locally.
- RT for TMT comes from Prosit iRT-TMT; RT MAD 2-11 min on training runs.
- D10cf arm C at 1% is above nominal by the upper-bound estimator for every learned scorer (and hyperscore): needs a look
  at per-fraction pooling before any CF claim.
- Training is dominated by one fold (fold 0 = 56.6k of 84.1k hyperscore PSMs); linear_pred beats gbm_pred in folds 0 and 1 (83,561 vs 82,340; 24,070 vs 23,284), gbm_pred wins fold 2 (19,918 vs 19,450).
