# M5a: cross-linking MS (dev only), 2026-10-27

Dev data only: D3 is PXD014337 `XLpeplib_Beveridge_QEx-HFX_DSS_R1` (DSS, Cas9 95/96-peptide library in 12 groups) and D3b is PXD029252 `20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2` (DSSO, ribosomal main library). **Sealed T3 was not touched.** Nothing under `data/sealed/` or `data/truth/sealed/` was read. The D3b library table and FASTAs were re-fetched from their public URLs into `/tmp/vulcan_m5/truth`, and their sha256 values are identical to the manifest. These are dev numbers, not claims (§6.6).

Code:
- `vulcan/xl/`: `linkers.py`, `doublets.py`, `search.py`, `sites.py`. These are new modules; no PSM column or kernel was changed.
- `scripts/m5/`: `truth.py`, `evaluate.py`, `run_m5.py`, `exhaustive_d3.py`, `comparators.py`, `summarize.py`.
- `tests/gpu/test_xl.py`: G0 for XL plus synthetic recovery.

Settings are in config.yaml and full numbers in metrics.json; `eval_D3*.json` holds every level/class/α cell.

## What was built
- **Cleavable (DSSO, D3b).**
  - Doublets are detected on the GPU over all spectra: pairs 31.972/c apart, an isotope check for c ≥ 2, top 300 peaks.
  - Doublets that have a *complementary* doublet (P_i + P_j + 158.004 = M_k) are ranked first, then by intensity, keeping 10 per spectrum. Ranking by intensity alone lost true precursor doublets to stub-carrying fragment ions, which also show the 31.97 spacing.
  - Candidate peptide masses are the doublet masses (plus a −1 isotope variant), their precursor complements (isotope error 0–3), or both.
  - Two linear-like searches then run on the L1 index: every row in each mass window is scored by K2 localization-aware shifted scoring, with the alkene stub as the shift and again with the thiol stub. The evidence is the better of the two.
  - The top 5 rows per window are kept, closed under target↔decoy pairing.
  - Pairs are formed under the precursor constraint (10 ppm, isotope error 0–3).
- **Non-cleavable (DSS, D3).**
  - α-first open search: one open hypothesis per spectrum. K1 counts unshifted fragments and keeps the top-20 α per spectrum, pair-closed via `engine.select_topk`.
  - For each isotope error, the β mass is looked up by mass. A decoy has the same mass as its target, so the window is pair-closed by construction.
  - α and β are both scored by K2 with shift = M_k − m(peptide).
- **CSM score** = min(evidence α, evidence β). This puts fragment evidence for both peptides into the score (DESIGN L4). The best pair per spectrum is the CSM, with ties broken by pf_id hash.
- **Competition:** each spectrum is also searched closed (linear) and with mass offsets for mono-links (+L+H2O, +L+NH3) and loop-links (+L), using M4 `search_modes`. A CSM is kept only if its XL score beats the best competitor. On D3b, of 31,340 rank-1 XL explanations, 10,493 win; the losers go to linear (5,702), mono-link (10,481) or loop-link (4,664). On D3: 1,781 of 4,222 win.
- **Link sites:** re-localized on the CPU over K/S/T/Y/N-term only (K2 places the shift on any residue).
- **Label blindness:** nothing in `vulcan/xl` reads labels. `tests/gpu/test_xl.py` extends G0 to XL. Under role swap plus pair-list permutation, the CSM tables (both linkers, including competitor columns, the homo flag and sites) are bit-identical. Synthetic cross-link recovery is 100% for DSS and 89% for DSSO (test thresholds 0.8 / 0.9 doublet hit rate); 3/3 pass. `tests/stats/test_flr_xl.py` passes (6/6). No kernel changed, so tests/index and the existing tests/gpu were not re-run for this milestone.

## Finding: homo-peptide pairs cannot be estimated by TT/TD/DD
A pair whose two sides are the same peptide, up to target↔decoy pairing (α = β, or α with its own decoy), can only win as TT or DD. The TD candidate scores min(e_t, e_d), which never beats the TT or DD candidate, so TD never wins. In the merged self class, these DD then *subtract* in (TD + 1 − DD)/TT. Before the fix, D3b self peptide pairs at nominal 1% had a true FDP of 16%.

Fix:
- Such CSMs are flagged label-free in L3 (`homo_family`: sequence(α) ∈ {sequence(β), sequence(partner(β))}).
- They are excluded from the TT/TD/DD estimation and reported separately as **not FDR-controlled**.
- An exploratory DD/TT curve is shown for them; it is not preregistered and not used for any claim.

The estimator and §6 are unchanged. This is a reporting-class decision made on dev data, proposed as tree decision D10 and to be confirmed by the reviewer under D9. Library homodimers are numerous in these designs (D3b: 500 true of 525 homo TT; D3: 132 of 133), so this costs sensitivity until a sound homo-pair estimator exists.

## Results: true FDP from the library designs (primary estimator (TD+1−DD)/TT; R = reportable, TT ≥ 50)
Peptides that are not in the library (human background, contaminants, non-library Cas9 or ribosomal peptides) count as false. So do links between groups. A peptide is in the library if its sequence equals a library peptide or one contains the other (I = L).

**D3b (DSSO).** Nominal 1% / 5%: accepted TT entities, with true FDP in parentheses.

| level | self 1% | self 5% | hetero 1% | hetero 5% |
|---|---|---|---|---|
| CSM | 116 (0.86%) R | 159 (2.5%) R | 775 (0.90%) R | 1457 (5.15%) R |
| peptide pair | 0 | 47 (2.1%) n.r. | 312 (1.28%) R | 420 (5.0%) R |
| residue pair | 0 | 46 (2.2%) n.r. | 334 (1.20%) R | 456 (4.8%) R |
| PPI | 0 | 0 | 169 (1.18%) R | 192 (2.6%) R |

- Homo-family CSMs (not controlled): 525 TT, 32 DD, 0 TD; 500 true.
- The conservative (TD+1)/TT is at or below nominal everywhere (e.g. hetero residue pairs at 1%: 308 accepted, true FDP 0.65%).

**D3 (DSS, α-first).** All library links are Cas9–Cas9, so they fall in the self class. The hetero class is Cas9×contaminant or contaminant×contaminant, i.e. all false.
- **Hetero:** 0 TT accepted at 1% and 5% at every level (no false heteromeric link was reported).
- **Self CSM:** 118 at 1% (0.0%) R; 276 at 5% (4.35%) R.
- **Self peptide pair:** 0 at 1%; 110 at 5% (**10.0%**) R.
- **Self residue pair:** 0 at 1%; 119 at 5% (**9.2%**) R.
- **Conservative estimator:** 82 peptide pairs (4.9%) and 86 residue pairs (4.65%) at 5%.
- **Homo-family:** 133 TT, 2 DD; 132 true.

So the primary estimator **exceeds nominal 5% for D3 self peptide and residue pairs** (about 11 false of 110). The false entries are high-scoring between-group Cas9 links and a few non-library Cas9 or contaminant self-links. This is not a mechanical artefact we could identify: the self search space here is one protein, and the number of decoys is small.

## Search-space reduction, non-cleavable (D3)
Index: 108,650 rows (targets plus decoys). Search space over the 5,502 XL-eligible spectra:
- All unordered pairs: 5.9e9.
- Mass-matching pairs (exhaustive; 10 ppm, isotope error 0–3): 1.50e9, median 226k per spectrum, involving 1.34e8 (spectrum, row) scorings.
- The exhaustive reference scores all of those rows with K2 (9.1e8 K2 evaluations, including per-isotope rescoring): 358 s wall. That run was GPU-shared and *not* under the lock (the lock was held for about 1 h by the xiSEARCH comparator), so its timing is informational only.

| search | pairs considered | reduction vs exhaustive (pairs / K2 evals) | wall | CSM self 1% / 5% | peptide pair self 5% (FDP) | residue pair self 5% (FDP) |
|---|---|---|---|---|---|---|
| α top-5 | 4.9e5 | 3,082× / 976× | 0.57 s | 125 / 227 | 98 (9.2%) | 104 (8.7%) |
| α top-10 | 9.3e5 | 1,613× / 503× | 0.41 s | 114 / 254 | 105 (8.6%) | 112 (8.0%) |
| **α top-20 (default)** | 1.8e6 | **829× / 256×** | 0.43 s (locked run) | 118 / 276 | 110 (10.0%) | 119 (9.2%) |
| α top-50 | 4.6e6 | 326× / 102× | 1.5 s | 128 / 285 | 90 (5.6%) | 94 (4.3%) |
| exhaustive | 1.5e9 | 1× | 358 s | 0 / 262 | 77 (2.6%) | 86 (2.3%) |

Reading:
- α-first accepts **as many or more** CSMs and links than exhaustive at the same estimator. The exhaustive space admits about 1000× more random pairs, which raises the decoy-pair competition: rank-1 agreement with exhaustive is only 26–39% among its top-1000 CSMs, and the remainder are mostly random high-scoring pairs.
- **Open issue:** pair-level true FDP falls as the space grows (9% at top-5–20, 4.3% at top-50, 2.3% exhaustive). Pruning is label-symmetric (K1 thresholds are blind, the shortlist is pair-closed, and β windows contain both members of a pair), so this is unlikely to be a leak. A plausible mechanism is that the decoy-pair *density* in the pruned space under-represents false same-protein pairs. There are only 5–11 false entities per cell, so it is not resolved.
- Candidate dev default: top-50 (326× fewer pairs, 1.5 s). Not switched here, so the timed run stays at top-20; to be decided with the reviewer.

## Comparators
The same truth code (`scripts/m5/comparators.py`) was applied to every comparator. TT only; "non-homo" excludes same-peptide links, as our controlled set does.

**Installed:**
- xiSEARCH 1.8.12 and xiFDR 2.4.2 (Apache-2.0), built from GitHub source with Maven.
- Scout, MS Annika and pLink 3 were checked by the comparator agent; none was a clean Linux install (see the agent report).

**D3 (DSS), our xiSEARCH 1.8.12 + xiFDR 2.4.2 run** (same DB and tolerances: 10/20 ppm, K/S/T/Y/N-term, hydrolysed/amidated mono-links; the lab's templates have no loop-links and add a non-covalent pseudo-linker; boosting off; 32 threads, under the bench lock):
- xiFDR's **default** self/between grouping is degenerate here. The self group has 0 TD and 3 DD, so 1% and 5% are identical: 265 residue pairs at 8.3% true FDP.
- The comparator agent recommends xiFDR's **single-pool** run (no grouping), which we use as the comparator:
  - residue pairs: **131 at 1% (true FDP 2.3%)** and **173 at 5% (9.8%)**;
  - of these, self residue pairs are 129 at 0.8% and 157 at 0.6%;
  - the excess comes from Cas9×contaminant heteromeric links (2 at 1%, 16 at 5%), all false. Pooling self and hetero lets them in, which is the failure our separate classes avoid.
- **So on D3, xiSEARCH is clearly better on self-link pairs than we are.** We have 0 residue pairs at 1% and 119 at 9.2% at 5%; it has 129 at 0.8%. Our score (min of two hyperscores) is the likely weak point.
- **Runtime: 3,812 s wall / 53,671 core-s** (0.69 s and 9.8 core-s per spectrum; 96 GB heap). Ours: 0.43 s search plus 0.64 s sites, i.e. **about 8,900× less wall for the search stage and about 50,000× fewer core-s for search plus sites**. The whole pipeline (ingest 6.8 s + index 0.35 s + search + sites ≈ 8.3 s) is about 460× lower in wall time.

**D3b (DSSO), submitter xiSEARCH 1.7.6.5 + xiFDR 2.1.5.5 result** (PRIDE mzid for this exact file, "1% link + CSM"):
- 1008 TT CSMs (1.3% true FDP). Residue pairs: 546 (1.8%); hetero 459 (2.0%); self 87 (1.1%). PPI hetero: 237 (1.7%).
- Non-homo: 504 residue pairs at 2.0%.
- Ours at nominal 1%: hetero residue pairs 334 (1.2%); PPI 169 (1.2%).
- So we have about 73% of xiSEARCH's heteromeric residue pairs at lower true FDP. Their DB differs slightly (it is the submitter's DB; truth is sequence-based).

**D3b runtime comparison:** a full xiSEARCH run was projected at many hours. It was replaced by a timed seeded 2,000-spectrum subset (seed 11), which took **4,490 s wall / 67,947 core-s** (2.25 s and 34 core-s per spectrum). Linear extrapolation to the full 43,596-spectrum file is **about 27 h wall and 1.48M core-s** (labelled "subset extrapolation"; the node was contended). Ours: 4.7 s search. Ours: 4.7 s search (+ 6.2 s sites, + 9.0 s ingest).

**Published** (not re-run): Beveridge 2020 reports actual FDR at nominal 5% of pLink 11.3%, StavroX 5.2%, xi 2.5%, Kojak 2.4–32%, MeroX 4.9% and XlinkX 0% (pooled designs, not exactly this file). Matzinger 2022: all six engines exceeded nominal 1%.

## Runtime
Measured on gpu-node (EPYC 9354, 32 affinity cores; A100 80 GB, shared with other agents' GPU jobs) under `flock /tmp/vulcan_bench.lock`. The hwinfo fingerprint is in metrics.json. The search is in-process, 3 passes, with bit-identical output across passes.

| stage | D3b (43,596 MS2) | D3 (5,502 MS2) |
|---|---|---|
| ingest RAW→store (copy to NVMe + decode, child process) | 9.0 s wall / 4.8 core-s | 6.8 s / 4.5 core-s |
| L1 index build | 1.5 s (265k peptidoforms) | 0.35 s (109k) |
| GPU index build | 0.16 s | 0.008 s |
| XL search, cold / warm | 4.9 / 4.7 s wall (≈ core-s; single thread plus GPU) | 0.45 / 0.43 s |
| ↳ warm stages | doublets 0.04, K2 linear-like 0.50, evidence/top-N 0.23, pairing 0.92, competitors 0.18, CSM table 0.17; about 2.3 s of host glue (window construction and chunking) is untimed | K1 α 0.04, β lookup 0.02, K2 α+β 0.03, pairing 0.26, competitors 0.03 |
| link-site localization (CPU, Python) | 6.2 s for 10,493 CSMs, now the largest stage | 0.64 s |
| L4 FDR + truth | about 20 s, Python and not optimized (evaluation only) | about 5 s |

## Issues / next
1. **Sensitivity is below xiSEARCH on both files.** On D3b: hetero residue pairs at 1% are 334 at 1.2% (ours) vs 459 at 2.0%. On D3 self links: 129 at 0.8% (xiSEARCH single-pool) vs 0 at 1% for us. A better XL score (static scorer) is the main next step. The score is a plain min of hyperscores; a frozen static XL scorer (per D7) is the obvious next step.
2. The D3 self peptide/residue-pair FDP exceeds 5% under the primary estimator (conservative is OK). Too few entities to tell bias from noise; needs the other D3 replicates (R2/R3), which are not downloaded.
3. Homo-family pairs need a sound estimator (decision D10 pending review).
4. Site localization and the 2.3 s of host glue should move to the GPU.
5. Truth is group-based (IMP-X-FDR style), so residue positions are not checked. Overlapping library peptides from different groups (e.g. Cas9 KFDNLTK vs FDNLTKAER) make residue-level truth lenient.
6. D3b `mainlib_rep2` has no HEK in the sample. The human proteins act only as a false-target DB background, not as a real background proteome; T3 is the real test.
