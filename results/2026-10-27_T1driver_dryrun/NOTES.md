# 3f: pre-registered T1 evaluation driver, dry run on dev D0/D1 (2026-10-27)

Only dev data were used: D0 = PXD001468, 2 HEK runs; D1 = PXD001077, 1 Pfu run. No file under data/sealed was read. §6, A1, the metrics, the labels and the estimators are unchanged; the driver only calls `vulcan.stats`.

Files:
- per dataset: `D0/`, `D1/` with `metrics.json`, `tables.md`, `driver.log`, `run_meta.json`, `timing.jsonl`
- `config.yaml`, `git_sha`
- `t4_selector/` (M7)

## What the driver does (`bench/t1_driver.py`, `bench/t1/`)
It is one command per dataset. Every stage resumes from its own outputs.

1. **db.** One multi-replicate L1 index holds the families T, F, E_2001..3 and the shuffled decoys dec_1001..1005:
   - F is the arm-C Arabidopsis set for human. It comes from the M1e neighbour filter.
   - The decoys are C-term-fixed hash shuffles.

   From it the driver writes **25 replicate DBs**, each as an index directory plus a peptide FASTA:
   - prod_j: 5 production DBs;
   - P_j_s: 15 arm-P DBs (5 decoy × 3 entrapment seeds);
   - C_j: 5 arm-C DBs.

   **Every method searches exactly these peptide sets**, so the external decoys and entrapments are identical by construction. §6.1 is verified on the final sequences: no sequence falls in two families or decoy sets, there are no I residues, and the partners are C-term-fixed shuffles. 0 violations on D0 and D1.
2. **ours.** `vulcan.cohort` does the ingest and the prod search. The engine then searches all 25 replicates with report_psms 5. The engine's prod output is **bit-identical to the cohort's** (checked on D0 and D1). Then come the L3 features and the frozen `static_gbm_pred.pkl` (sha256 31bfe685… asserted). The estimator is paired double competition.
3. **Comparators** run on TRFP mzML:
   - Sage 0.14.7 (LDA, `predict_rt: true`);
   - Sage+mokapot;
   - **Sage+pred+mokapot**: Sage's own top-1 PSMs, plus Vulcan's 9 AlphaPeptDeep/RT L3 features appended to Sage's PIN, then mokapot. This is the MS2Rescore-like strongest comparator (A1.9).
   - Comet 2026.01 + mokapot.
4. **evaluate.**
   - **Units, identical across methods:**
     - peptidoform = (I/L sequence, Met-ox positions);
     - sequence = q is the minimum over the sequence's peptidoforms;
     - picked protein groups = `picked_group_fdr` over each method's sequence-level scores, with proteins mapped from the same FASTA for every method.
   - Each method uses its own q (Sage peptide_q, mokapot peptide q, our double competition).
   - **Verdicts:** §6.3 for P (15 replicates) and C (5 replicates), under both A1.6 conventions. (a) decides, because 3e's coverage simulation chose (a); the worse of (a) and (b) is reported as well.
   - **§1 ratio:** run-level IDs summed over runs and averaged over shuffles. The 95% CI comes from a percentile bootstrap (2000 draws, seed 11) that resamples runs × decoy shuffles.
   - The ratio is computed (i) at own q (this counts only if the comparator's verdict is "supported") and (ii) at matched arm-P FDP of 1% and 5% (first crossing on a 160-point grid of the method's own q, recomputed in every draw).
   - Vectorized counts are asserted equal to `vulcan.stats.entrapment` on every replicate.

## Dry-run results (dev; D0/D1 are a *selection set*, so these are development curves, not claims)
Ratios are ours / comparator. Brackets are the 95% CI over runs × decoy shuffles. "**counts**" = the comparator is "supported" at 1%, so its own-q comparison applies.

**D0 (2 runs).** At 1%, every method is "inconclusive" at peptide level except Comet+mokapot peptidoform (supported). The arm-P means are 1.01–1.20%; arm C is 0.83–1.02%.

| vs | level | own q 1% | matched FDP 1% | matched FDP 5% |
|---|---|---|---|---|
| Sage (RT on) | peptidoform | 1.025 [1.022, 1.028] | 1.028 [1.026, 1.029] | 1.006 [1.002, 1.010] |
| Sage (RT on) | protein | 1.013 [1.011, 1.016] | 1.014 [1.012, 1.016] | 1.014 |
| Sage+mokapot | peptidoform | 1.021 [1.016, 1.025] | 1.026 [1.023, 1.027] | 1.001 [0.998, 1.006] |
| **Sage+pred+mokapot** | peptidoform | **0.990 [0.988, 0.994]** | **0.990 [0.989, 0.992]** | 0.992 [0.990, 0.994] |
| **Sage+pred+mokapot** | sequence | 0.992 [0.989, 0.995] | 0.991 [0.990, 0.993] | 0.993 |
| **Sage+pred+mokapot** | protein | 0.995 [0.993, 0.997] | 0.995 [0.993, 0.997] | 0.995 |
| Comet+mokapot | peptidoform | 1.037 [1.025, 1.048] **counts** | 1.034 [1.026, 1.039] | 1.008 [1.002, 1.013] |
| Comet+mokapot | protein | 0.997 [0.989, 1.004] **counts** | 0.997 [0.988, 1.004] | 0.997 |

D0 IDs at 1% (mean of 5 shuffles):

| Method | Peptidoforms | Sequences | Protein groups |
|---|---|---|---|
| ours | 24,328 | 22,981 | 6,429 |
| Sage | 23,727 | 22,434 | 6,320 |
| Sage+mokapot | 23,822 | 22,512 | 6,337 |
| Sage+pred+mokapot | 24,576 | 23,180 | 6,465 |
| Comet+mokapot | 23,443 | 22,163 | 6,267 |

At 0.1%, at the peptidoform level, ours vs Sage is 1.093 and vs Sage+pred+mokapot is 0.989 [0.982, 0.997].

**D1 (1 run, so the CI is over decoy shuffles only).** Every ratio lies in 0.98–1.004, with lower CIs ≥ 0.979, and all pass 0.97. Most methods are "supported" at 1%; Sage+mokapot is inconclusive. Protein level is not-evaluable at 1% (α·N_T < 20). Full table: `D1/tables.md`.

**Reading.**
- **Non-inferiority passes** on every level vs every comparator, on both sets.
- **Ours vs Sage (RT on) is +2.8% at matched FDP 1% on D0.** This reproduces the review's corrected figure (M4), not M3's +3.4%.
- **Against the strongest comparator, the same AlphaPeptDeep features fed to mokapot on Sage's PSMs, we are ~1% lower** at both own q and matched FDP. So on dev, the gain over Sage comes from the predicted features, not from the scorer. A T1 "gain" claim must be quoted against this comparator (A1.9); only "no loss" is supported.
- **The CIs on D0 are narrow and probably too optimistic.** With 2 runs, resampling runs has only 3 distinct outcomes. T1 has 6 runs.

## Sage predict_rt: reconciling the tree 1e conflict
Both statements were true, about different runs:
- **M1e's baselines had RT on.** All `scripts/baselines/run_searches.py` configs set `predict_rt: true`.
- **M3's "Sage-LDA (our TDC)" reference had RT off.** It used the E2 oracle runs, built by `vulcan/gpu/sage_bridge.py:88`, which sets `"predict_rt": False`. That is the run the review (M4) flagged.

The driver uses RT on and checks that the setting takes effect: with predict_rt true, Sage fits its RT model, so `predicted_rt` is non-zero on 98–99% of PSMs. Spearman(aligned RT, predicted RT) on confident targets is 0.92 on D0 and 0.94 on D1. With predict_rt false, predicted_rt is 0 everywhere. The control run with RT off (`--rt-off-control`, prod_1001) shows the cost:
- D0: Sage 23,432 → 23,732 peptidoforms at 1% (+1.3%), and 20,621 → 20,865 at 0.1%;
- D1: +0.1%.

The PIN carries `sqrt(delta_rt_model)` either way; it is constant when RT is off.

## Things found while building it (fixed in the driver, recorded)
1. **mokapot 0.10.0 ignores `--seed`.** The CLI seeds `np.random`, but `brew()` uses `default_rng(None)`. Identical inputs gave 37–68 different peptides at 1% on D1. Sage also writes PIN rows, and numbers their SpecIds, in thread order.
   - **Fix:** `bench/t1/mokapot_seeded.py` seeds `default_rng` from the CLI-seeded global state, and PIN rows are sorted by content.
   - mokapot is now bit-reproducible. All mokapot results here come from this version.
2. **mokapot sometimes writes "mokapot score" with the opposite sign.** Its q-values are right; this was D1 seed 1003. It is the M1e "collapsed protein count" (660 ± 369).
   - The driver orients every method's score so that q is non-increasing in score. This is label-blind and uses the method's own q.
   - With the seeded mokapot it flips 1 of 25 Sage+mokapot replicates on D1 and none on D0 (`n_replicates_score_sign_flipped` in metrics.json).
3. **Comet with `spectrum_batch_size = 0`** writes a PIN identical to the M1e batching (41,068 rows, 0 differences) in about half the wall time. The `-j` peptide index is not used, because its output differed in 4 rows.
4. **Sage 0.14.7 is not bit-reproducible on augmented searches.** Discriminant values differ slightly between identical runs on 2 of 4 checked replicates.
   - A full D1 re-run from a **clean worktree at the committed HEAD**, with a fresh work dir, gave identical production IDs for all 5 methods at every α and level, and identical verdicts.
   - The only differences were ≤1e-4 in arm-P means for the Sage-based methods, and ratios within ±0.0004.
   - Vulcan and Comet are bit-identical.
5. Index row strings exceed 2 GB on human (`take()` goes through large_string). Peptidoforms do not share residue offsets, so sequences are compared as fixed-width bytes.

## Driver timing (unlocked, shared node, loadavg 30–77; operational, not a speed claim)
**D1 end to end** from the clean worktree: 43 min wall. This covers the db stage including arm-C/foreign filtering, ours, comparators and evaluation.

**D0**, per stage:

| Stage | Wall time |
|---|---|
| arm-C filter | 5 min |
| multi-replicate index build (107.7M peptidoforms) | 24 min |
| materialize 25 DBs + verify | 15 + 3 min |
| ours: cohort ingest + prod search | 22 s |
| ours: 25 replicate searches | 66 s |
| ours: features + scoring | 5 min |
| Sage, 26 searches, ~65 s each (Comet's 50 file searches ran alongside, 6 jobs × 4 threads) | 29 min |
| Sage+pred features | 10 min |
| mokapot, 75 runs (6 parallel) | 4 min |
| evaluate | 36 min |

D0 total ≈ 2.2 h. Driver core-s are in `timing.jsonl`.

**T1 projection:** HYE DB is about 1.3× human, and there are 6 QE HF-X runs. TRFP takes ~6 × 4 min in parallel. Expect ~4–6 h, with Sage and Comet scaling in the number of spectra.

## T4 selector (review M7): `scripts/data/select_t4.py`
The selector now mechanically excludes, each with its own reason column:
- `train` (`train.tsv` + DATA.md role=train);
- `dev_manifest:*` (dev_sealed, dev_sealed_extra, tier_s, d9, d4_ref, g2null);
- `e9_d8_stream` (`e9_d8.tsv`);
- `predictor_training` (as before).

`assert_disjoint()` is importable and is also available as `--check`.

These exclusions remove 159 projects from the 14,190-project frame: 142 D8/E9, 13 train, 4 G2-null.

- **Frame mode** (exclude before allocation and shuffling) changes 12 of the 40 selected projects, all in `orbitrap|LFQ|none` cells. The cause is that the shuffled cell lists get shorter; none of the dropped projects was itself excluded. Output: `t4_selector/frame/`.
- **Draw mode** (default) shuffles the pre-registered frame and skips excluded candidates while drawing. It **reproduces the committed `data/manifests/t4_selection.tsv` byte-identically** (40 projects, 300 runs; no skip was needed). Output: `t4_selector/draw/`.
- `assert_disjoint` passes on the committed selection.

**The sealed T4 selection is unchanged. Whether frame mode should replace it is the user's decision.**

The `bench/sealed_run.py` T4∩train assertion (2 lines) was not applied. Also, `--dataset T4` has no `dev_sealed.tsv` rows yet.

## Sealed command (NOT executed)
```
python3 bench/sealed_run.py --dataset T1 --verify --purpose "3f T1 pre-registered ID evaluation (DESIGN §1, §6.3, A1): frozen static_gbm_pred vs Sage (RT on), Sage+mokapot, Sage+AlphaPeptDeep+mokapot, Comet+mokapot on identical external decoys; 5 decoy x 3 sealed entrapment seeds (2004-2006)" -- python3 bench/t1_driver.py --dataset T1 --rt-off-control
```
In T1 mode the driver:
- refuses to run without `SEALED_DATASET=T1` / `SEALED_FILES`, which sealed_run sets;
- forces seeds 2004–2006;
- asserts that the HYE FASTA sha256 equals T1's (`cae44ed2…`);
- refuses any mzML cache;
- deletes the staged stores and mzML of the sealed runs at the end.

Outputs go to `results/<date>_T1_sealed/T1/`.

## Readiness
**The pipeline is ready.**
- Committed: 8656413, 4cadf63, 2c658b1.
- It re-runs from a clean worktree at HEAD.
- Its dependencies (`scripts/e9/stage_pred_env.sh`, A1.6 `CI_CONVENTIONS`, the `e9_d8`/`g2null` manifests) are now committed by their owners.

**The T1 look itself is still blocked by D12** (A1 tail-G2/G1/G3 fail; A2 pending). Running it also needs a clean tree (other agents have uncommitted `vulcan/gpu/*` edits) and a `milestone-*` tag.
