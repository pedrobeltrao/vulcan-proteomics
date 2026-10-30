# A3/A5 FINAL freeze (D29, D31, D32, D33, D34): A5 implementation landed; G0, D1 dry run, closure guard at the final commit

- **Date:** 2026-10-28. **Agent:** the A5 implementer and final-freeze agent. T1/T5 exposure: none (no sealed file, no ledger entry beyond its header, no T1/T5 output read). Dev data only.
- **Freeze code commit F = `06ce624`** (`06ce624e945940738bda757e0f96061e85ab5f4a`). `config/a3.yaml freeze` and `config/a5.yaml hashes/env_manifests` are recorded in the commit after it (configs cannot name their own commit; `a3_config.freeze_status` accepts descendants whose frozen trees are byte-identical).
- **No tag created.** The PI tags `prereg-v0.4.3` (and `prereg-v0.4.4`, after reading the A5 claim-script review).

## Bottom line

| Item | Result |
|---|---|
| A5 implementation list (A5.11.1-A5.11.8) | Implemented and tested (below); claim script reviewed T1-blind: **adopt** (reviews/claims_a5_review.md, e757c30 + addendum c88b5e9; its R1 condition is met: the D7 arm-M end-to-end test now runs on the committed dev output) |
| G0, full v2_sep pipeline, at F (clean worktree) | **PASS** (`g0.json`): MS1 hypotheses and features, chimeric PSM tables (42,833), NCE, v2 features, v1/v2 scores bit-identical; identical to the rehearsal's record at e0b94aa in every field but wall time |
| A3.1.6 dry run, D1, prereg-v0.4.2 vs F (clean worktrees) | **PASS, exact**: v1 artefacts 51/51 byte-identical; v1 metrics columns identical; A1 comparator artefacts recomputed from the same Sage/Comet searches 400/400 byte-identical (`dryrun_D1.json`) |
| A3.1.6 dry run, D0 | **Not re-run (D32 deviation).** Cited: rehearsal at e0b94aa, v1 artefacts 102/102, A1 comparator artefacts 425/425, v1 metrics identical (results/2026-10-28_A3_freeze/D0/). Per-diff argument below |
| Closure guard (v1 closure, both entries) | **PASS at the evidence commit** (`closure_guard.json`); 9 new allow-list entries pinned to the exact diffs at F |
| Dev `--look` dry run (A5.11.8) | **PASS** (`dev_look/`): one joint T5(D7)+T1(D1) look through the real sealed_run/t5_armM/driver/claim script in a clone; exit 0 both, leak test pass, both manifests match the ledger result lines, scratch deleted, claim script exit 0; 49.4 min |
| Sage+pred+mokapot tolerance (D32/D33) | Reported, not a gate: A3.1.6 rests on the exact with/without comparison. Independent-run variability (this D1 run vs the 3f run): below |

## 1. A5 implementation (what landed before F)

| A5 item | Files | Tests |
|---|---|---|
| A5.11.1 one look | `bench/sealed_run.py` (`--look a5`, `--recover`, `--retry`, `--dev`), `bench/env_manifest.py` | `tests/bench/test_sealed_run_one_look.py` (13: second look refused incl. after ledger rewrite; plain T1/T5 refused; ledger-only commit skipping; missing prereg/milestone tag; hash mismatch; C9 mismatch (verdict and sha); GO branch A6 rules; caller `VULCAN_SAGE` dropped; both ledger lines before unlock; console whitelist; unterminated look blocks until recover; SIGTERM failure lines + deletion; retry rules (deviation, new tag, latest only, <=2, not after success); leak test; frames-only summary) |
| A5.11.2 driver hardening, completeness, additions | `bench/t1_driver.py` (`look_guard`, deletion, manifest, exit 3), `bench/t1/completeness.py`, `bench/t1/a5.py` | `tests/bench/test_t1_driver_look.py`, `tests/bench/test_t1_a5.py` |
| A5.11.3 T5 arm M | `bench/t5_armM.py` | `tests/bench/test_t5_armM.py` (17); dev equivalence on D7 exact vs the committed v2_sep (A3 preconditions) and v1 3e tables (`results/2026-10-28_A5_impl/`) |
| A5.11.4 claim script | `bench/t1/claims_a5.py` (labels, gate tables, voids, family S, P1, headline, A5.9.3 state from the A6 status file at the look's base commit, D33 non-reproducibility column) | `tests/bench/test_claims_a5.py` (incl. D0 and D7 end-to-end), `tests/bench/test_a5_config.py` |
| A5.11.5 config | `config/a5.yaml` (look, hashes of 18 committed files incl. the review, 6 env manifests, gates, labels, c9, a6, deviations []) | label texts pinned to the draft |
| A5.11.6 freeze extension | `vulcan/gates/a3_config.FROZEN_TREES` (+ scripts/e7b, m3, m3e, entrapment, e9, baselines) | `tests/a3/test_a3_freeze.py`, `test_a3_config.py` |
| A5.11.7 lint | `bench/t1/wording_lint.py` | `tests/report/test_report_labels.py` |
| A6 hook (D32) | `results/2026-10-28_A6_hn1/a6_status.json` read at the look's base commit (closed enum; absent/unverifiable -> "stopped: node time limit"); `sealed_run` requires it (and `prereg-v0.4.5` only for gate-run statuses) in the GO branch before 18:00 CET | in the tests above |

Implementation choices: `scripts/a5/IMPLEMENTATION_NOTES.md`; claim-script stricter readings: `config/a5.yaml claim_script_stricter_readings`.

Issues found and fixed by the dev look dry run: the env allow-list lacked `LD_LIBRARY_PATH` (no GPU driver in the child); `t5_armM --mode dev` refused sealed_run's stand-in files; the look root must be writable shared storage (the shared filesystem outside the repo is read-only; `$VULCAN_LOOK_ROOT` is NFS-shared); NFS-safe deletion after SIGTERM.

## 2. Final freeze (D32)

### 2.1 What was run, from clean worktrees
- `F` worktree `$VULCAN_SCRATCH/vulcan_a3freeze/wt_final` (clean tracked tree; resource symlinks only).
- **G0**: `scripts/a3/g0_freeze.py` (pred env) -> `g0.json`, PASS.
- **D1 with**: `bench/t1_driver.py --dataset D1 --a3` at F, fresh work (replicate indexes shared, unchanged dbs.py), 2,235 s.
- **D1 without**: prereg-v0.4.2 worktree, fresh `ours` stage (v1 artefacts byte-identical to the rehearsal's tag run as well, 51/51), then comparators+evaluate on the "with" run's raw Sage/Comet searches (`--mzml-cache`), 490 s.
- **Compare**: `scripts/a3/dryrun_compare.py` -> `dryrun_D1.json`.

### 2.2 D0 (D32 deviation), argued per diff
The rehearsal showed byte identity on D0 at e0b94aa. The diffs between e0b94aa and F in the v1 closure are:
1. `vulcan/cohort.py`, `vulcan/ms1/features.py`, `vulcan/ms1/gpu.py` (cec4f37, GPU MS1): the GPU finder is opt-in (`ms1_backend=gpu`), the default is the CPU path, which neither the T1 driver nor step 4b changes; tests/ms1 34/34 bitwise. The D1 run at F reproduces the v1 artefacts and the MS1-dependent v2_sep/G0 tables exactly.
2. `bench/t1_driver.py`: look-mode code runs only with `SEALED_LOOK_ID` set (absent in both dry runs) and the A5 additions run after evaluate, writing a new key; v1 and A1 numbers are computed before them (D1 identical; the dev look shows the same v1 columns with the guard active).
3. `bench/t1/a5.py`, `bench/t1/completeness.py`: new, additive, read-only on outputs.
4. `bench/sealed_run.py`, `bench/env_manifest.py`: imported only by the look guard.
5. `vulcan/gates/a3_config.py`: the frozen-tree list check (no computation path).
None of these touches a D0-specific code path (no dataset-conditional code), so D1 identity plus the rehearsal's D0 identity is the evidence. This is disclosed as the D32 deviation.

### 2.3 Sage independent-run variability (D33: reported, not a tolerance)
This D1 "with" run vs the independent 3f run (2c658b1): v1 (`methods.vulcan`) and Comet+mokapot identical.

| Comparator | max abs delta arm-P mean FDP | max abs delta v1 ratio at gated 1%/5% | max abs delta v1 ratio, any threshold |
|---|---|---|---|
| Sage | 2.2e-5 | (see json) | 1.9e-5 |
| Sage+mokapot | 4.2e-5 | (see json) | 1.8e-5 |
| Sage+pred+mokapot | 2.3e-5 | 4.7e-4 (protein, matched 1%) | 1.4e-2 (sequence, own q 0.1%: 74 IDs) |

Together with D33's D0 values (Sage+pred+mokapot <= 4.3e-4 arm P / 9.1e-4 ratio; Sage+mokapot <= 2.3e-3 ratio), these are the comparator non-reproducibility reported next to T1 comparator ratios by the claim script. The 0.1% own-q ratio is reported, never gated.

### 2.4 Dev look dry run (A5.11.8)
- Clone at ed63785 (look code identical to F except two later refinements in `sealed_run.py`, both unit-tested: A6 tag-less statuses (D34) and the /tmp free-space precondition).
- One `sealed_run.py --look a5 --dev`: T5 = D7 (3 runs, entrap 2001-2003, decoy 1001) through `t5_armM.py --mode dev`; T1 = D1 through the driver with the look guard active. Exit 0/0; leak test pass; manifests == ledger; work and scratch deleted; claim script (`--look-id`) exit 0 (`dev_look/claims.md`).
- **v1 bit-identical with and without the preceding arm-M step**: `methods.vulcan`, `matched_fdp_curves.vulcan` and every v1-side ratio field identical to the plain F run; only Sage-based comparator fields differ (separate Sage invocations).
- Recover/retry/SIGTERM: exercised by the unit tests and, for SIGTERM, by an aborted dev look (failure lines written, outputs deleted).
- Wall 49.4 min.

## 3. For the look (operator)
- Command (one invocation, clean tree, `milestone-*` tag on the base commit, `prereg-v0.4.3` and `prereg-v0.4.4` ancestors):
  `rm -rf /tmp/vulcan_look_env && conda_envs/pe/bin/python bench/sealed_run.py --look a5`
  (no other argument; commands and env from `config/a5.yaml look`). Then, once both result lines exist: `bench/t1/claims_a5.py --t1 <root>/T1/out/metrics.json --t5 <root>/T5/out/armM.json --look-id <id> --out <dir>`.
- **Duration**: dev look 49 min (1+3 small runs). The D0 A3 driver (2 runs) took 3.0 h on a loaded node; T1 has 6 runs and a larger DB, so expect roughly 4-9 h for T1 plus ~0.5 h for T5. `look.time_limit: 9h`. This does not fit 18:00-19:30.
- **Disk**: the look refuses unless node-local `/tmp` has 150 GB free (`look.min_scratch_free_gb`); at 14:10 it had ~55 GB.
- **Clean tree**: `a3_config.freeze_status` needs no untracked files under vulcan/bench/scripts/config/tests and the frozen trees byte-identical to F. `bench/` contains the T2 driver; any T2 commit to `bench/` after F breaks the freeze (a re-freeze would be needed).
- **Status at 14:20: the freeze does NOT hold at HEAD**, because the T2 agent committed `bench/t2*` after F (5696a9a, 7bebe32, f6836b4); `freeze_status` gives `head_tree_hashes_match: False` and the look would refuse. Nothing else in a frozen tree changed (G0, D1 dry run and guard stay valid: bench/t2* is outside every T1/T5/v2_sep import closure). Once T2 stops committing to `bench/`, run `python scripts/a5/refreeze_t2_only.py` (refuses unless only bench/t2* changed, the guard passes and the A5 hashes match), then commit `config/a3.yaml`; then tag.
- Mokapot env is node-local (`/tmp/vulcan_m1e/envs/mokapot`); env manifests are checked at look start.

## Files
`dryrun_D1.json`, `g0.json`, `closure_guard.json`, `predcache_manifest.json` (shard file hashes; T1 builds its caches fresh), `metrics.json`, `D1/` (with/without metrics and timings), `dev_look/` (record, claims output, D7 armM, D1 metrics and manifest).
