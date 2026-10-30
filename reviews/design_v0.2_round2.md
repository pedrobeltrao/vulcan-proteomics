# Round-2 review of DESIGN.md v0.2

Independent reviewer, 2026-10-26. Inputs: PROMPT, LITERATURE, DESIGN v0.2, DATA.md, manifests; v0.1 reviews skimmed. **[B]** blocker, **[M]** major, **[m]** minor. "est." counts must be replaced by M1 measurements.

**Verdict.** v0.2 resolves most round-1 blockers; the architecture is sound. Remaining problems: §6 arithmetic/validity errors (some inherited from round 1); three inconsistent scores (RESET, fixed-score pass, fold models); and a dev/sealed split that DATA.md (declared authoritative) contradicts. Do not tag `prereg-v0.2` until §1–2 are fixed.

---

## 1. Were the round-1 fixes adequate?

**Adequately resolved:** HPC B1–B3; domain B1/B2; stats B1, M2/M3 (CONGA), M6 (entrapment decoys); XL merge-then-FDR; the "4-core Sage vs a 3090" clause.

**Fixes that are wrong, over-engineered or inconsistent:**

**[B1] The G2 sample size is wrong and the gate cannot pass.** At n = 16.6k the 99% CI half-width is *exactly* 0.01. Any p̂ ≠ 0.5 fails, so a perfectly fair scorer passes with probability ≈ 0. For 90% power, n = (4.22·0.5/0.01)² ≈ **44.5k per cell**. The top-decile clause multiplies this by 10, giving **≈ 450k null competitions per stratum**. G2 also tests the decoy recipe and target-side "neighbour" structure, not the scorer.

**Fix:** make G2 **differential**: the learned pipeline's null target-win rate may exceed the fixed-score pipeline's on the *same* competitions by ≤ 0.01 overall and ≤ 0.02 in the fixed-score top decile; report absolute p as a diagnostic. Stratum = **instrument class** only (OT-HCD, IT-CID, TOF). Cost ≤ ≈ 110k null competitions per class (11k in the top decile at margin 0.02 dominates).

**[B2] G2's null data are not null.**
- *Shifted arm:* [+7, +50] minus integers still admits carbamylation (+43.006), Na (+21.98), K (+37.96), +32/+48 oxidation and paralog substitutions (A→T +30.01, L→F +33.98, D→Y +48.02); targets have paralogs, decoys don't, so a *fair* scorer fails. **Fix:** non-physical mass-defect offsets, n + 0.5 ± 0.1 Da, applied after recalibration.
- *Wrong-species arm:* conserved peptides and contaminants are true matches. **Fix:** apply the entrapment homolog filter (I/L, single-substitution, isobaric) against foreign proteome + contaminants to the human DB.
- *Pipeline:* semi-supervised training on a pure-null run trains on noise. **Fix:** **mixture design** — real human run + null spectra through production; evaluate win rates on null-provenance spectra only.

**[B3] The RESET restriction contradicts the MCD and L6.**
- MCD items 3 and 6 promise RESET rescoring *and* picked protein plus global FDR.
- §6.4 says RESET supports only pooled peptide FDR, so protein-level and global lists come from a "fixed-score pass".
- L6 uses fold models, which are static.

Three scores; protein-level non-inferiority would pit our fixed score against MSFragger+Percolator. CALIB_DECOY_SEED also implies a second search of every run.

**Fix ("static-first"):** default = a **frozen, project-disjoint, low-capacity model** (linear/GBM, label-blind + predicted features) per instrument class. TDC holds without RESET, so Group-walk, picked protein and L6 share one score. RESET becomes an E4 arm; if kept, fit calibration and G4 on **RESET's training-decoy half** and drop CALIB_DECOY_SEED.

**[M1] G4 drives grouping from estimation decoys.** "Strata that fail [decoy-quantile KS] stay separate under Group-walk" makes the group partition depend on decoy labels, which violates Group-walk's fixed-partition assumption. **Fix:** decide G4 only on training or calibration decoys, or on the dev data, and freeze the partition before estimation.

**[M2] The Group-walk strata explode and the merge rule is unspecified.**
- Six crossed factors (charge × offset × length × PTM × specificity × chimeric rank) leave most cells below 200.
- "(pep, z)" charge strata conflict with the peptide competition unit: charge requires precursor-level units.

**Fix:** ≤ 2 factors per mode (closed: specificity × chimeric rank; open: CONGA stem/offset-bin; PTM: composition); other factors become features. Pre-register merges: adjacent offset bins by mass, rare PTM classes → "other", rank ≥ 3 pooled, counted on label-blind pair counts. Min size is a power rule, not a validity rule.

**[M3] The verdict table overlaps and the non-inferiority margin is unattainable.**
- "Supported" (mean ≤ α) contains "Supported-noisy". Define Supported = mean ≤ α **and** upper CI ≤ 1.5α.
- With ~4.7% decoy-induced variability at 1% FDR, a lower CI ≥ 0.99× on six runs needs a +3–5% *win*. **Fix:** average both methods over the same 5 shuffled decoy FASTAs (Sage accepts external decoys); margin 0.97.
- If the comparator is "violated"/"inconclusive", compare at matched entrapment FDP only.

**[M4] The label-swap test is too weak.** Swapping flags cannot detect *positional* leaks:
- targets stored first, so ids below N_T are targets;
- K1 top-K ties broken by index;
- `rev_` prefix parsing.

**Fix:** the test also permutes pair/DB order and assigns ids by hash; K1/K2 break ties by hash(sequence, seed).

**[m1]** RESET "training constant 3/2" is undefined: write the estimator explicitly (+1 placement matters at 0.1%) and verify in §6.7. G5's AUC > 0.99 depends on ID rate; compute it on G2 nulls (alarm > 0.55).

---

## 2. New blockers and majors

**[B4] DATA.md contradicts the pre-registration.**
- DESIGN §6.6 says "DATA.md is authoritative", but DATA.md has no dev/sealed column.
- DATA.md lists the *sealed* T2 (PXD007058) for "FLR calibration" (a dev use) and fetches it with the dev data.
- Its seeds differ: DECOY 1001–1003 vs 1001–1005, no sealed entrapment seeds, no calibration seed.
- D2 (PXD000138), D4 (PXD016766), T1 and T5 are not in the manifest. T5 is not even named.
- The FASTA table is the placeholder `FASTA_TABLE`.

**Fix before the tag:** `role ∈ {dev, sealed, train}` in the manifest; sealed files in `data/sealed/` (`chmod 000`), unlocked only by `bench/sealed_run.py`, which requires a clean tree, a milestone tag and a ledger line; name T5.

**[B5] The sealed FLR set has no ground truth.** DATA.md states that the PXD007058 library site list is only in the closed-access JPR SI. The plan's proxy ("sites supported by ≥ 2 engines") is circular for an FLR claim, and "exclude pools with positional isomers" cannot be applied without the list. **Fix:** swap roles. Make PXD007058 dev and seal **held-out PXD009449 (ProteomeTools PTM) pools**, whose per-pool peptide lists are published, plus a held-out PXD000138 subset. If the PXD007058 SI is obtained later, it can return as a second sealed set.

**[M5] Assay demos lack dev data in the MCD's modality.**
- XL: dev D3 is DSS (an arm); the only cleavable data is sealed T3 (one file). Add a dev DSSO file (another PXD029252 file, or PXD019120).
- LiP: no dev set; add a DDA LiP/TrP pair.
- CF-MS: D6 is **DIA** but the MCD is DDA-only; pick a DDA CF set (CFdb/PXD042664) or move CF to the DIA arm.

**[M6] CF "profile-continuity prior" manufactures the artefact it guards against.** Favouring transfers that make an elution profile continuous inflates exactly the correlation that CF-MS analysis scores. **Fix:** profile-agnostic, fraction-local transfer FDR; E10 measures false co-elution with spiked foreign-proteome profiles.

**[M7] Leakage through third-party predictors.** AlphaPeptDeep, MS2PIP and DeepLC training sets include PRIDE projects, likely PXD001468 among them. If an evaluated project was in predictor training, its original false PSMs may have memorized RTs and intensities. **Fix:** before the prereg tag, cross-check every dev and sealed accession, and the T4 sampling frame, against the predictors' training manifests. Exclude overlaps or retrain.

**[M8] The cross-encoder arm is not a cross-encoder.** A late-interaction MaxSim head is ColBERT-style, effectively a bi-encoder: it drops the joint attention that makes a cross-encoder worth its cost. **Fix:** split stage = shared spectrum/peptide token encoders whose pooled output *is* the E5 embedding; joint stage = 1–2 cross-attention layers over top-K token states only. Cache pooled peptide embeddings (128-d fp16 ≈ 2.5 GB per 10⁷), recompute token states for top-K, avoiding the 50–500 GB cache.

**[m2]** MS1 feature finding must go on the GPU or be budgeted in E1. Isotope −1…+3 × z 2–4 = up to 15 hypotheses per spectrum; report inflation in E2.

---

## 3. Practicality: can each gate be computed?

| Gate / level | Needed | Likely available (est.) | Verdict |
|---|---|---|---|
| G2 as written | ≈ 450k null competitions per stratum × many strata | D1: 1 Velos file (~10⁴ MS2); shifted D0: ~10⁵ peptide-level competitions per offset draw | **Infeasible** |
| G2 differential, per instrument class | ≈ 110k | Shifted copies of dev runs (OT-HCD) plus 5–10 foreign-species runs per class | Feasible |
| Peptide entrapment, α = 0.1% | N_T ≥ 20k | D0 (2 QE files) ≈ 25–35k; D1 Pfu ≪ 20k | D0 is evaluable but ~15 entrapment wins (±55% CI), so expect **Supported-noisy at best**; D1 not evaluable |
| Protein entrapment | α·N_groups ≥ 20 | 5–8k groups | Evaluable at 1%, marginal at 0.5%, **not at 0.1%** |
| M arm | N_T ≥ 2k (1%), 20k (0.1%) | Pfu: few k | Only 1%. **Hold yeast out of all training** and use a modern yeast DDA run for T5 to reach 0.1–0.5% |
| XL T3 | Library-truth false counts | 1 file, ~1–3k CSMs | Report at 1% and 5%; PPI level only at 5% |
| Replicates | 5 decoy × 3 entrapment × 3 arms = 45 augmented searches per dev set | — | Breaks "minutes per iteration" |

**Simplifications that keep the rigor:**
1. **Family-max single-pass entrapment.** Tag index entries by family (T, E_s, dec_d(T), dec_d(E_s), Arabidopsis); K1/K2 keep top-2 per family per hypothesis, so every replicate's competition is exact post hoc. Dev: pseudo-reverse × 3 entrapment seeds + Arabidopsis in one pass; full 5 × 3 at milestones; validate once against separate searches.
2. **α levels by unit.** Pre-register 0.1% as peptide-only and "report, don't gate". Gate on 1% and 5%.
3. **Require all three arms only for sequence-capacity models.** Static low-capacity scorers need P + C. The cross-encoder and bi-encoder need P + C + M.
4. **T4 within MCD scope.** T4 samples SILAC/DIA/XL strata the MCD cannot process; 30–50 whole projects ≈ 1–5 TB, and replicated comparator baselines cost hundreds of Sage-hours. **Fix:** frame = DDA, Thermo/Bruker, LFQ/TMT/phospho, UniProt proteome, not in predictor training; ≤ 8 seeded runs per project, ≤ 300 runs, ≤ 1 TB; 1 decoy × 1 entrapment seed, full replication on 5 projects; stream-and-delete via `/tmp`; dry-run the selector on metadata now.
5. **E9 on human-only runs.** Global paired entrapment needs one proteome. Run heterogeneous species as per-species master lists.

---

## 4. Brief (PROMPT.md) coverage gaps

- **Split-then-joint cross-encoder and encoder reuse:** mis-specified (M8).
- **Bi-encoder coarse-to-fine:** E5 does not say *where* it should win. Target it at open, non-specific and HLA search, where fragment indices reach 30–86 GB. Closed tryptic search is not a target, because Sage-style indexing already wins there.
- **Amortized indexing:** the content hash includes pilot-chosen mods, so projects rarely share an index. **Fix:** a canonical superset index per (proteome, enzyme, mod superset), with project mod sets as masks over the records' mod field; measure reuse on the T4 frame.
- **Repository-scale PTM-site FDR:** L6's "global FLR" has no site aggregation. Specify best evidence per (protein, site) across runs, global decoy-AA site competition gated by global peptidoform q, and expected false sites; otherwise single-evidence sites accumulate (82–98% FDR, LITERATURE §2.6).
- **DIA:** a "standard" workflow in the brief, but arm-only. Consider moving Astral nDIA tier 1 (spectrum-centric, reuses the chimeric DDA path) into the MCD for a cheap DIA-NN comparison.
- **"Synthetic PTM localization sets":** blocked by B5.

---

## 5. The five things to do first

1. **Reconcile DATA.md with §6 and implement sealing (B4, B5, M5, M7):** roles, seeds, named T5, swapped FLR truth, dev DSSO/LiP/DDA-CF sets, predictor-overlap check, sealed dir + runner + ledger. Then tag.
2. **Rewrite §6 gates before coding them (B1, B2, M1–M4).** Add Monte-Carlo checks that G2 passes at p = 0.5 with ≥ 90% power and fires at p = 0.52.
3. **Adopt static-first scoring (B3);** RESET becomes the E4 arm.
4. **M1 harness with family-max entrapment** in the Sage fork; Sage/Comet baselines on D0/D1/D4; publish the **measured feasibility table** (N_T per α, entrapment wins, null competitions per class, core-s) and fix gated α levels from it before engine work.
5. **E1 ingest + T4 metadata dry run:** core-s per 10⁵ MS2 per reader, and total T4 bytes/runs — showing whether the headline is ingest-bound and whether T4 fits ~1 TB and a few node-days.
