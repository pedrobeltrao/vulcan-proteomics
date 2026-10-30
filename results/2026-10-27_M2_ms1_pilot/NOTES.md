# M2: MS1 features, precursor hypotheses, L0.5 pilot pass, E0 first read (2026-10-27)

- **Data:** dev only. D0 is PXD001468 01A and 01B (QE HCD-OT); D1 is PXD001077 Velos005137 (HCD-OT); one D4 file is used only as the reporter positive control. Sealed data was never touched.
- **Node:** gpu-node, 32 Zen 4 cores plus an A100 (`hwinfo.json`).
- **Timing:** all timed runs went through `flock /tmp/vulcan_bench.lock` (`timings.jsonl`).
- **Files:**
  - `metrics.json` holds the curated numbers. `pilot_e0.json` is the full per-run dump (pilot report, hypothesis stats, E0 arms).
  - `validation_*.json` compares against biosaur2.
  - `<run>.calib.json` is the per-run calibration (offsets plus tolerances plus priors).
  - `<run>.ms1.json` holds feature-finder timings and hypothesis stats.
  - `config.yaml` records every parameter.

## 1. Feature finder (vulcan/ms1): review and fixes

The recovered WIP worked but cost **~650 core-s per D0 run** and reported 2.5× more features than biosaur2. The changes:

| Change | Effect (01A) |
|---|---|
| Hill linking rewritten as an O(n) two-pointer merge per scan pair (numba, parallel over scans). It is bit-identical to the old global searchsorted+lexsort (tested). Root labelling is one pass, because pred < i always. | link 21.5 s → 0.5–1.3 s wall |
| Valley splitting is parallel over hills, using suffix maxima | 4 s → 2 s |
| Isotope search visits hills in m/z order and filters scan overlap with contiguous m/z-sorted start/end arrays. The count pass only tests isotope 1. | isotopes 35 s × 32 thr → 3.4 s; total 650 → **~30 core-s** |
| The trace cosine uses full-hill norms (biosaur2) | precision vs biosaur2 0.43 → 0.46 |
| Averagine "explained ≥ 0.5" rule and the charge-harmonic ban (both from biosaur2) | small gains |
| **Removed the "overlap < 3 scans ⇒ accept without trace test" shortcut.** It was the main noise source. | features 1.02 M → 0.67 M; precision 0.46 → 0.65; recall 0.75 → 0.85; mono errors (ID-shifted) 10.4% → 7.6% |

Also fixed in `hypotheses.py`: the window-scan output buffer was capped at `n_ms2 × 64`, which allowed out-of-bounds writes in numba. It is now an exact two-pass allocation.

**GPU:** the torch hill linker is 0.26 s versus 0.55 s on CPU (8 threads, 1.7 core-s) for 50 M peaks. But torch init takes about 10 s and linking is under 10% of the runtime, so the GPU does not pay. The CPU path is the default; `device="cuda"` is kept.

### Validation against biosaur2 0.3.4 (D0, m/z 350–1500)

Matching is any-match: same charge, ≤ 10 ppm, overlapping RT and apex within 0.5 min.

| | 01A | 01B |
|---|---|---|
| features ours / biosaur2 | 625,764 / 401,659 | 587,833 / 370,943 |
| recall (biosaur2 features found) | **0.847** (Q1 0.82 … Q4 0.88 by intensity) | **0.849** |
| precision (ours found by biosaur2) | **0.646** | **0.635** |
| ≥ 3-isotope features only: recall / precision | 0.825 / 0.546 | 0.833 / 0.564 |
| mass diff, matched (median / robust SD / p95 \|·\|) | 0.00 / 0.09 / 0.92 ppm | 0.00 / 0.09 / 0.89 ppm |
| apex RT diff p95 | 0.10 min | 0.11 min |
| log10 mono-apex intensity Pearson r | 0.990 | 0.992 |
| unmatched biosaur2 features that match ours one isotope off (ours lower / higher) | 12.9 k / 6.7 k | 12.2 k / 6.1 k |
| **ID support:** Sage q ≤ 1% target PSMs on a feature with the exact mono — ours vs biosaur2 | **0.872 vs 0.864** | **0.868 vs 0.859** |
| ID support: on a feature only after an isotope shift — ours vs biosaur2 | 0.076 vs 0.077 | 0.080 vs 0.078 |
| runtime, core-s / run (wall) | **29.6–30.6** (11 s at 32 thr; 16–17 s at 4 thr, 27 core-s) | **28.5** (10.8 s) |
| biosaur2 runtime, core-s / run (wall) | 1523 (458 s, plus mzML parse) | 1350 (422 s) |

Runtime includes the Python import, the store read, features and the hypothesis table (0.26 s). Velos takes 4–5 core-s.

Precision against biosaur2 is a definitional comparison, not a truth measure. The ID-based support is the better
check: we recover slightly more identified precursors with the correct mono and have the same isotope-shift rate, at about 50× fewer core-s.
The remaining "our mono is one lower" set (about 13 k) consists mostly of 2–3-isotope features with a weak M−1 hill.

Synthetic truth (`tests/ms1/test_features.py`) gives ≥ 97% recall and precision with uniform noise, the correct charge for z=2 envelopes (no z=1 harmonic), valley splitting, missing-mono handling, and rejection of heavy 2-peak patterns.

## 2. Precursor-hypothesis table (vulcan/ms1/hypotheses.py, SCHEMA.md)

The table has one row per (spectrum, precursor group, isotope error −1…+3). Groups come from three sources:
- co-isolated MS1 features (charge 2–6, rank ≤ 5, fraction ≥ 1%);
- the instrument (charge known, no agreeing feature);
- a z = 2–4 fallback.

The schema is pinned by a test. The table carries raw and calibrated m/z, `neutral_mass`, a per-source isotope prior from the pilot, the feature id and intensity, window intensity and isolation fraction, the label-free chimeric rank, and instrument agreement with its isotope offset.

| | 01A | 01B | Velos |
|---|---|---|---|
| MS2 | 41,820 | 46,816 | 15,374 |
| rows / spectrum (groups × 5 isotope errors) | 10.2 | 10.1 | 7.8 |
| precursor groups / spectrum | 2.04 (1: 41%, 2: 30%, 3: 18%, ≥4: 12%) | 2.02 | 1.57 |
| groups by source: feature / instrument / fallback | 78.1 k / 7.3 k / 0 | 85.8 k / 8.7 k / 0 | 23.6 k / 0.5 k / 0 |
| MS2 with ≥ 1 feature group | 91.5% | 90.4% | 97.9% |
| MS2 with ≥ 2 feature groups (chimeric candidates) | 55.0% | 53.1% | 41.2% |
| instrument charge+mono agrees with a feature | 82.5% | 81.3% | 96.8% |
| of those, instrument mono off by ±1..3 isotopes (the feature corrects it) | 20% (+1: 3.6 k, −1: 1.9 k) | 20% | 2.4% |
| top-ranked feature charge = instrument charge | 92.6% | 93.1% | 99.2% |
| **coverage of baseline IDs:** Sage q ≤ 1% PSMs whose (z, mass) is a row | **99.2%** (feature 91.2%, instrument 8.0%) | 98.7% | 99.4% |
| isotope error of the explaining feature row (0 / +1 / −1) | 92.3% / 6.2% / 0.8% | 92.2% / 6.2% / 0.8% | 97.4% / 1.8% / 0.5% |

Every MS2 in these runs has an instrument charge, so the fallback path is exercised only in tests.
Pilot-derived isotope priors for feature rows on 01A: P(0) = 0.89, P(+1) = 0.087, P(−1) = 0.016.

## 3. L0.5 pilot pass (vulcan/calib) and E0 first read

**Pilot setup:**
- A 10% seeded subsample (seed 11) is searched with Sage 0.14.7 at ±50/±50 ppm, isotope −1..3, M-ox + protein N-term acetyl, through the `vulcan.calib.searcher.Searcher` interface (`SageSearcher` now; the GPU search implements `search(batch, cfg, workdir)`).
- Anchors are all rank-1 PSMs with Sage spectrum q ≤ 1%, decoys included (a score threshold, label-symmetric). There were 2142 / 2468 / 1123 anchors, of which 20 / 23 / 10 were decoys.
- Calibration is ppm = f(RT) + g(m/z), with piecewise-linear robust knots. A component is kept only if it lowers the MAD.
- Fragment errors come from Sage `--annotate-matches`.

**Tolerances (μ ± 4 SD of the sigma-clipped anchor-error core)** before → after recalibration:

| run | precursor (ppm) | fragment (ppm) | model used |
|---|---|---|---|
| 01A | [−4.2, 4.9] → [−4.5, 4.4] (width −2%) | ±6.9 → ±6.7 (−3%) | m/z only (no RT drift); precursor MAD-SD 1.36 → 1.32 |
| 01B | [−1.8, 6.8] → [−4.1, 4.1]: a +2.5 ppm offset removed. The symmetric half-width needed drops from 6.8 to 4.1 (−39%). | [−5.1, 9.1] → ±6.8 (−2 ppm offset removed) | m/z |
| Velos | [−12.5, 11.5] → [−7.0, 7.0] (**−41%**) | [−13.3, 14.3] → [−10.6, 11.0] (−21%) | RT drift −3…+4 ppm; precursor MAD-SD 3.1 → 1.9 |

Relative to the prod settings (±10 / ±20 ppm), the empirical tolerances are 2.2–2.9× tighter on the QE runs.

Only 83–93% of anchors fall inside μ ± 4 SD. The rest are a flat wide-window tail plus a negative-ppm shoulder, which is deamidation (0.984 Da) matched as isotope error +1 (0.019 Da off, i.e. −19 ppm/kDa).

**Run diagnostics:**
- **Analyzer/activation:** metadata (OT/HCD) agrees with the data on all three runs: ppm-scale fragment errors, b/y matched-intensity share 0.28 against c/z 0.004–0.008, and no 1/3 low-mass cutoff.
- **Reporters:** none on D0/D1. Positive control D4 dq_00095: TMT detected in the SPS-MS3 spectra (90% of spectra against 1.2% at control offsets) and correctly absent in its IT-CID MS2.
- **Label-blind offset histogram** (open search −150…+500 Da, rank-1 PSMs, no labels):
  - D0: unmodified, isotope ±1/±2/+3, oxidation, **carbamylation +43.006 (urea)**, formyl +27.995, unannotated +53.919, +301.99 and +16.998, and ammonia loss.
  - D1: unmodified, deamidation, oxidation.
  - No TMT, iTRAQ, SILAC or dimethyl offsets in any run.

**E0 first read.** Sage on the full MS2 of each run; Δ peptides at peptide q ≤ 1% versus arm A:

| arm | 01A | 01B | Velos |
|---|---|---|---|
| A raw spectra, prod tol (±10/±20) | 11,005 pep / 20,071 PSM | 13,255 / 23,312 | 6,038 / 10,870 |
| C calibrated spectra, prod tol | +0.01% | +0.05% | +0.02% |
| B calibrated spectra, pilot tol (core 4 SD) | **−5.0%** (PSM −5.7%) | **−5.8%** | −0.8% |
| D raw spectra, pilot tol | −4.8% | −5.1% | +0.2% |
| F calibrated, pilot precursor, prod fragment | −3.7% | −3.8% | −0.6% |
| G calibrated, prod precursor, pilot fragment | −1.8% | −1.9% | −0.1% |
| H calibrated, pilot precursor, fragment 6 SD | −3.9% | −4.6% | −0.6% |
| E calibrated, cov99 tol (prec ±12–13, frag ±17–20) | **+0.7%** (PSM +1.1%) | **+0.4%** (+0.8%) | +0.4% (+0.7%) |

**Read:**
- Recalibration alone is neutral for Sage on these instruments: its ±10 ppm window already covers the 2.5 ppm offset and the Velos drift.
- **Tightening tolerances to the calibrated error core costs Sage 4–6% of peptides on D0.**
  - Most of the loss comes from the precursor (F): the lost PSMs are the deamidation/isotope-confusion shoulder and the heavy tail, which Sage's ±10 ppm plus isotope errors 0..2 currently absorb.
  - A tight fragment window (G) costs about 2%. Sage's hyperscore and `min_matched_peaks` benefit from low-intensity fragments with larger errors.
- Slightly *wider*, coverage-based windows (E) gain 0.4–0.7%.
- So the pilot's value for a Sage-like scorer is in calibration and diagnostics, not in narrow windows. A narrower window only pays if the search explicitly models the deamidation/isotope offsets (mass-offset mode) and the decoy competition benefits from the smaller search space. That is not observed here.
- These are dev numbers at Sage q-values without entrapment-verified FDP. They are not a claim.

**Pilot cost:** 7–46 s wall per run. Sage at 32 threads on a 4k-spectrum subsample takes about 167 core-s for the wide search and 150 core-s for the open search on D0. That is about 10× the feature finder and wasteful (LESSONS: 32-thread Sage on small inputs spins). Run it at 4–8 threads, or move it to the GPU search once that lands.

## 4. Bugs found in passing

- **Sage tolerance sign convention.** Sage 0.14.7 accepts peptide masses in [exp + lo, exp + hi], i.e. error (exp − calc) in [−hi, −lo].
  - The first E0 run passed asymmetric error windows unconverted. On 01B (+2.5 ppm offset) arm D collapsed to 5,030 peptides (−62%). Fixed in `SageSearcher` (windows are now negated and swapped), and the run was repeated; the bad run is not reported.
  - Sage's fragment windows are centred on the experimental peak in the fragment index but on the theoretical ion in scoring, so asymmetric fragment windows are ill-defined. They are now symmetrised.
  - **Relevant for the GPU/Sage-parity work:** asymmetric tolerances should not be fed to Sage.
- Sage writes `isotope_error` in Da (× 1.00335), not as an integer.

## 5. Open issues

1. E0 needs entrapment-verified FDP (DESIGN §6) and a scorer that can use narrow windows. Re-read E0 with the GPU search.
2. About 13 k features per run whose mono is one isotope lower than biosaur2's. Consider an M−1 plausibility test (averagine ratio of M−1/M0).
3. Isotope tolerance is fixed at 8 ppm. biosaur2's per-run isotope-error calibration (σ ≈ 1.5 ppm here) could tighten it.
4. Not implemented: the pilot's iRT→RT map, SILAC/dimethyl detection beyond the offset histogram, ion-trap (low-res) calibration in Da, and Bruker/IM features (`im` is ignored).
5. Hypothesis rows: 10 per spectrum on D0 (5 isotope errors × 2 groups). Feature rows almost never need −1/+3 (0.8% / 0.07%). Pruning +2/+3 for feature rows (pilot prior < 0.01) would cut about 35% of rows. That decision belongs to the search (E2/E7).
6. The biosaur2 timings were taken under the lock, but this agent's own untimed exploration ran concurrently. Its core-s is robust; its wall time is an upper bound.
