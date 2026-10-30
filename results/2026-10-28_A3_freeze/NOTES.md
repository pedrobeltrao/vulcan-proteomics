# A3 freeze REHEARSAL: A3.1.3, A3.1.4, A3.1.6, A3.1.7 (dev D0/D1 only; nothing sealed)

- **Date:** 2026-10-28. **Agent:** the A3 freeze agent. I did not build E7, E7b, A2 or the A3 gate code.
- **Status: candidate, not final.** D29 (A5) requires A5's code to land before the A3 freeze. D31 adds the cec4f37 GPU MS1 merge. So this is a rehearsal at candidate commit **e0b94aa**. `config/a3.yaml` has the values under `freeze_candidate` (status: candidate). **`freeze.commit` stays TODO**, so `require_freeze_commit()` still refuses for every A3 script and for the sealed T1 driver. **No tag was created.**
- **Not touched:** `bench/sealed_run.py`, anything under T1–T5, DN2 data, `/tmp/vulcan_a2`, A4/A5/A6 drafts and checks. Models, metrics, controls, labels, estimators and existing tests are unchanged. Two existing test lists were extended only to classify new scripts, as the test requires.

## Bottom line

| Precondition | Result at e0b94aa |
|---|---|
| A3.1.6/7 driver | `bench/t1_driver.py --a3` (forced for T1) adds vulcan_v2_sep, C2 = `sage_chim2+pred+mokapot` and `sage_chim2`. **R = 5 is excluded** by E7b's own rule: R = 5 had to win on both D0 and D1, and it lost on D1 (6,238.3 vs 6,238.9). There is no wide-window variant, so the wording is "Sage chimera at narrow precursor tolerance". `a3_secondary` writes the LOO backstop, the merged-stratum verdicts, void flags and the Sage-rank split |
| A3.1.6 dry run, v1 bitwise | **PASS on D1 and D0.** Driver at prereg-v0.4.2 (A1 set) vs e0b94aa (`--a3`): every v1 artefact is byte-identical (D1 51/51, D0 102/102 parquets). Every v1 metrics value is identical on all common keys (D1 and D0; the tag lacks only the later `deciding_conservative` keys), and every A1 comparator artefact recomputed from the same Sage/Comet searches is byte-identical (D1 400/400, D0 425/425 files). v1 is also identical to the independent 3f run on both sets |
| A3.1.6 Sage within tolerance | **Met for the with/without comparison:** delta 0, because the downstream artefacts are identical. **Against an independent earlier run** (3f, other Sage invocations): Sage and Sage+mokapot are within 1e-4 / ±0.0004. **Sage+pred+mokapot exceeds both on D1**: arm-P mean up to 1.6e-4, v1/C1 ratio up to 4.3e-4 (protein). **On D0 the independent-run deltas are larger**: Sage+pred+mokapot arm P up to 4.3e-4, v1/C1 ratio up to 9.1e-4 (protein); Sage+mokapot v1 ratio up to 2.3e-3 (peptidoform, `dryrun_D0.json` by_level). Sage itself stays ≤ 2e-5. The cause is Sage's own non-determinism: all 25 replicate TSVs differ. The documented tolerance came from one D1 rerun pair and is too tight for the mokapot-based comparators on D0 |
| A3.1.4 G0, full v2_sep pipeline | **PASS at e0b94aa** (`g0.json`). On D1 + E2001 vs the role-swapped, row-shuffled index, each side used fresh MS1, its own hypotheses and a fresh prediction cache. MS1 hypothesis and feature tables, chimeric PSM tables (42,833 PSMs), NCE, every v2 feature, and the v1 and v2 scores are all bit-identical. The pilot at c6116ec also passed |
| Closure guard (Q1) | **PASS at e0b94aa** with the new allow-list entries (90 files). **FAIL at HEAD 279945d**: `vulcan/cohort.py`, `vulcan/ms1/features.py` and `vulcan/ms1/gpu.py` were changed by the cec4f37 merge after the dry run, so their diffs are not covered (`closure_guard.json`) |
| A3.1.3 freeze fields | Recorded as candidate: tree hashes, Sage, AlphaPeptDeep, mokapot wrapper, prediction-cache manifest, E7b commit 7a2d6ee, speed-path decision. Final values must be recomputed at the final freeze |

## What was done

### 1. Code (commits c6116ec, e0b94aa)
- `bench/t1/a3.py` is new, and `bench/t1_driver.py --a3` was added. The `evaluate.py` change is additive: it exposes `LAST_CURVES`.
- **Prediction-cache isolation (e0b94aa).** The A3 arm predicts through its own `work/predcache_a3`.
  - The first dry run (c6116ec) shared the cache. The chimeric stage then predicted some Sage+pred peptides first, in different batches. AlphaPeptDeep is not batch-composition invariant: 426 of 296,873 RTs and 8 MS2 vectors differed on D1.
  - As a result, the A1 comparator Sage+pred+mokapot moved by up to 63 peptidoforms at 0.1%. **The A3 additions reached an A1 comparator; this is now fixed.** v1 itself was unaffected, because it is scored first.
- **`vulcan.gates.a3_config.freeze_status` fix.** The config cannot name its own commit, so "at the freeze" now means: HEAD is the freeze commit or a descendant **with identical frozen tree hashes**, on a clean tree. Before this fix the check could never pass.
- **Closure guard made stricter.** It now follows `python -m` subprocess modules (`vulcan.cohort`, which ingests and searches for v1, and `vulcan.ms1`) and the `scripts/e7`, `scripts/e7b` roots. The previous guard did not see `vulcan/cohort.py`, which is on v1's T1 path and was changed by E7b.
- `scripts/a3/g0_freeze.py`, `scripts/a3/dryrun_compare.py` and `tests/a3/test_a3_freeze.py` are new. `tests/a3` passes 63/63.
- IMPLEMENTATION_NOTES notes 33–40 document every choice.

### 2. Dry-run design (`scripts/a3/dryrun_compare.py`)
- **Worktrees.** Clean git worktrees at prereg-v0.4.2 ("without", the A1-tagged driver) and at e0b94aa ("with", `--a3`). Both use the same replicate DBs (the db stage, unchanged since the tag, is symlinked).
- **Fresh in both runs:** cohort ingest, all replicate searches, the prediction cache and v1 scoring.
- **"With" run:** everything is fresh (TRFP mzML, Sage, Comet, Sage chimera, mokapot).
- **"Without" run:** reuses only the "with" run's raw Sage and Comet search outputs. It recomputes the Sage+pred features and every mokapot run.
- **Same commit without `--a3`** (D1): the A1-set metrics are identical to the `--a3` run on every key.
- **Against independent earlier runs** (3f D1 at 2c658b1): v1 `methods.vulcan` and Comet+mokapot are identical. The Sage-based deltas are given above.
- **v2_sep vs E7b** (D1): the v2_sep record is identical. v2_sep / C2 at peptidoform matched FDP 1% is 1.1379 (E7b: 1.1374). C2 now has its own cache and its own Sage run.
- **MS1 tables** of D0 (2 runs) and D1 at e0b94aa are byte-identical to E7's pre-E7b tables.

## Open, blocking or differing for the final freeze
1. **The final freeze must be redone after the A5 code and cec4f37.**
   - **vulcan/ tree:** HEAD's `vulcan` tree differs from the candidate's (cohort.py and ms1 from the GPU MS1 merge), so the guard fails at HEAD.
   - **Re-runs needed:** the D0+D1 dry run (`dryrun_compare.py`; about 1 h for D1 and 3 h for D0 on a loaded node) and G0 (2 min), from a clean worktree at the final code commit. Then append allow-list entries for the new diffs.
   - **Expected:** GPU MS1 is opt-in (default CPU), so the default path should stay bitwise. The dry run must show this.
   - **A5 code:** if it edits `bench/t1_driver.py` or `bench/t1/*.py`, every allow-list entry pinned to those diffs is void, which is intended.
2. **D0 metrics** finished after the first commit and are added (`dryrun_D0.json`). The D0 same-commit-without-`--a3` evaluate was not run, for time; D1 has it. D0 `a3_secondary`: LOO full-ratio checks 105/105, no void flag, 0 runs without MS1 hills.
3. **Sage+pred+mokapot tolerance** (above). Suggestion for the PI: record its independent-run variability, about 1.6e-4 arm P and 4.3e-4 ratio at protein level, as the documented non-reproducibility. Otherwise the independent-run reading of A3.1.6 fails marginally for C1. The with/without reading passes exactly.
4. **Sequencing reminder (A3.3.1, D24).** A2 is not evaluable, so as written there is no T1 look under A3 unless A5 (D29) replaces that clause.
5. **Data provenance.** The dry runs reused the 3f/E7 replicate indexes in `$VULCAN_SCRATCH/vulcan_t1/{D0,D1}`. These are node-local and rebuildable by the unchanged `dbs.py`.

## Files
- `dryrun_D1.json`, `dryrun_D0.json`: the full comparisons.
- `D1/`: metrics and tables of the with, without and same-commit-no-A3 runs.
- `D0/`: v1 artefacts, A1 comparator artefacts, and the with/without metrics and tables.
- `g0.json` (e0b94aa) and `g0_pilot_c6116ec.json`.
- `closure_guard.json`.
- `predcache_manifest.json`. Shard file hashes are not content hashes: the D0 tag and freeze caches differ by file hash, but the PINs built from them are byte-identical.
- `metrics.json`.
