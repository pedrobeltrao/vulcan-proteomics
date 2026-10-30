# A5 implementation notes (DESIGN §6.0d; normative text reviews/amendment_A5_draft.md rev 1, ae21c30)

Implementation choices, stricter readings and deviations of the A5.11 implementation list. Every deviation that
loosens nothing is listed here; the claim script's own stricter readings are in `config/a5.yaml
claim_script_stricter_readings`.

| A5 item | Where | Notes |
|---|---|---|
| A5.11.1 one look | `bench/sealed_run.py`, `tests/bench/test_sealed_run_one_look.py` | see 1-9 |
| A5.11.2 driver hardening, completeness | `bench/t1_driver.py` (`look_guard`, `write_manifest`), `bench/t1/completeness.py`, `tests/bench/test_t1_driver_look.py` | see 10-13 |
| A5.11.2.5 additions | `bench/t1/a5.py`, `tests/bench/test_t1_a5.py` | fork notes in the module docstring |
| A5.11.3 T5 arm M | `bench/t5_armM.py`, `tests/bench/test_t5_armM.py`, `results/2026-10-28_A5_impl/` | see 14 |
| A5.11.4 claim script | `bench/t1/claims_a5.py`, `tests/bench/test_claims_a5.py`, `reviews/claims_a5_review.md` | see 15 |
| A5.11.5 config | `config/a5.yaml`, `scripts/a5/fill_hashes.py`, `bench/env_manifest.py`, `tests/bench/test_a5_config.py` | see 16 |
| A5.11.6 freeze extension | `vulcan/gates/a3_config.FROZEN_TREES`, `tests/a3/test_a3_freeze.py` | |
| A5.11.7 lint | `bench/t1/wording_lint.py`, `tests/report/test_report_labels.py` | |
| A5.11.8 dev look dry run | `scripts/a5/dev_look_dryrun.py`, `results/2026-10-28_A3A5_final_freeze/` | see 17 |
| A6 hook (D32) | `bench/t1/claims_a5.py` (`a6_state`), `bench/sealed_run.py` (`a6_ready`) | see 8 |

1. **Ledger format.** Look rows carry `look=<id> attempt=<n> c9=<GO|NO-GO> a6=<state>` in the purpose column; result
   rows carry `result: look=<id> exit=<rc>, wall=<s>s, manifest_sha256=<sha256 of <ds>/out/manifest.sha256>`; the leak
   test writes a `look` row. A look **succeeded** iff both datasets have `exit=0` result rows and the leak row says
   `leak_test=pass`. A look is **unterminated** while a dataset row or the leak row is missing.
2. **History scan.** Every `+|` line of `git log -p -- results/SEALED_LEDGER.md` plus the current file, de-duplicated.
   Any opening row ever written for T1 or T5 (also a pre-A5 plain row, which would get a `legacy-n` id) refuses a new
   look; only `--retry` remains.
3. **Ledger-only commit skipping** (unchanged rule, extended): a `sealed-ledger:` commit that touches only the ledger
   and a failed look's `results/<date>_look_<id>_failure.md` is skipped when looking for the milestone tag.
4. **Environment.** The child environment is exactly `look.env` (explicit values) plus `SEALED_FILES`,
   `SEALED_DATASET`, `SEALED_LOOK_ID`, `SEALED_LOOK_ROOT` (= `<root>/<id>/<ds>`) and `SEALED_LOOK_SCRATCH`
   (= `/tmp/vulcan_look_<id>`). Besides the variables A5.11.5 names, the allow-list carries what `env/activate.sh` sets
   for the GPU path (`CUDA_HOME`, `TORCH_CUDA_ARCH_LIST`, `TORCH_EXTENSIONS_DIR`, `VULCAN_NVCC_GENCODE`), `HOME`, `LANG`,
   `PROJECT_ROOT`, `CONDA_PREFIX`. The staged-tool variables are set to the values the stager itself sets
   (`$VULCAN_LOCAL_ENV/{trfp,vulcanraw/VulcanRaw.dll,dotnet}`), so the driver's post-re-exec assertion sees equality.
5. **Work in node-local scratch.** `--work` may lie under `$SEALED_LOOK_SCRATCH/<ds>` (A5.5.3 allows node-local scratch
   for speed); `--out` must lie under `$SEALED_LOOK_ROOT`. Both must be absent or empty.
6. **Console.** Child stdout+stderr are line-buffered into `<ds>/logs/console.log`; a line reaches the terminal only if
   it matches a `look.console_whitelist` regex. sealed_run's own lines are fixed templates.
7. **Failure.** On a non-zero exit or a failed leak test both datasets have still run (A5.5.2). The summary keeps the
   last whitelisted line, the exception type names and the `file:line:function` frames (paths relative to the repo or
   basenames); no message text. The look root, the scratch and any leaked file are deleted unread. SIGTERM and SIGINT
   terminate the child, relock the files and write `failed (terminated|interrupted)` rows for every dataset without a
   result plus a `leak_test=not run` row. SIGKILL / node loss leaves an unterminated look for `--recover`.
8. **C9 and A6.** `c9_branch` recomputes GO from the committed pilot file at `c9.commit` (sha256 must match): GO iff the
   configured quantity (`summary['f0.20/own/R3/kappa1.0/budget200'].P_v1_minimum`) **and** (stricter) the minimum over
   every walk cell at 200 runs are >= 0.8. In the GO branch, before the A6 time box (D31.3: 18:00 CET 2026-10-28), the
   look additionally needs `prereg-v0.4.5` as an ancestor and a committed `results/2026-10-28_A6_hn1/a6_status.json`
   with a status in the closed enum; after the time box it proceeds and the claim script maps an absent file to
   "stopped: node time limit" (A6.2).
9. **Retry.** Only the latest look, only if it failed, never after any successful look, at most two, only with a dated
   `deviations` entry `{date, retry_of: <id>, ...}` and a milestone tag not used by any earlier look. The other A5.5.6
   retry conditions (T1-blind review, re-passed dry runs, re-freeze) are procedural and are not checked in code.
10. **Driver guard** (`look_guard`, active iff `SEALED_LOOK_ID` is set; after the re-exec): refuses non-default
    `--stages/--methods/--ours/--mzml-cache/--keep-sealed-derivatives`, `--out/--work` outside the look root (or
    scratch, 5) or non-empty, a look id that is not the open ledger line of `$SEALED_DATASET`, missing prereg tags, and
    any `VULCAN_*` variable not equal to `look.env` (`VULCAN_T1_FMATCH_0001`, `VULCAN_T1_NO_LABEL_CACHE` are therefore
    refused). `--rt-off-control` is forced.
11. **Completeness** (`bench/t1/completeness.py`, A5.5.4) checks presence only: all `A3_METHODS` with verdicts and no
    `missing_replicates`, the run list, the 25-replicate design, every A3.3.3 comparator in `meta.a3`, both numerators
    against every comparator, `fdp_reached` and `n_boot_nan` on every matched ratio (own-q bootstraps cannot produce
    NaN: the denominator is floored at 1e-9), one LOO entry per run, and `a5.complete`. The driver also counts every
    "missing output" event of the evaluate loader. On any violation the driver writes `out/INCOMPLETE`, the manifest,
    and exits 3.
12. **Deletion** in the look: the whole `work/` except `index/` (A5.5.3), then `out/manifest.sha256` (sha256 of every
    file under `out/`). The pre-A5 sealed deletion (stage, mzml, ms1) is subsumed.
13. **Dev look mode.** The guard is keyed on `SEALED_LOOK_ID`, not on the dataset role, so the dev dry run exercises it
    on D1.
14. **T5 dev scripts not refactored.** `scripts/m3e/armM.py` and `scripts/a3/dev_armM.py` are left unchanged so their
    committed outputs keep their provenance; `bench/t5_armM.py` carries their design and its dev mode is checked
    against their committed tables (`results/2026-10-28_A5_impl/t5_dev_equivalence.json`). Deviation from A5.11.3.1
    (refactor), recorded; it loosens nothing.
15. **Claim script inputs.** With `--look-id` the claim script reads both manifest sha256s and the base commit from
    the look's ledger rows (refusing unless the look succeeded) and reads the A6 status file with `git show
    <base>:<path>`; sealed inputs without `--look-id` are refused.
16. **Hashes and environment manifests.** `hashes` are the sha256 of the committed bytes of every look code file and
    of the claim-script review; `env_manifests` are `bench/env_manifest.py` digests (Sage binary bytes; Comet, TRFP,
    mokapot, pe and pred envs by their conda-meta / dist-info RECORD / pyvenv / interpreter-link files, which pin every
    installed file's hash; full byte hashing of GBs of CUDA libraries at look time was judged too slow). The staged env
    under `VULCAN_LOCAL_ENV` is a copy of the pe/pred envs (`cp -a`), so it is covered by their manifests; it is rebuilt
    from scratch for the look by deleting `/tmp/vulcan_look_env` before the look.
17. **Dev dry run** in a `git clone --shared` (fake tags, stand-in sealed files, fake A6 status "not completed"; never the
    project ledger). The commands read the real dev raws (D7 for T5, D1 for T1).
