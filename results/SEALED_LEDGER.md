# SEALED_LEDGER

Every use of a sealed dataset (`data/sealed/`, DESIGN.md section 6.6) is appended here by `bench/sealed_run.py`
*before* the files are unlocked. Do not edit by hand. A line without a following `result` line means the run was
interrupted or is in progress.

| date (UTC) | sha | tag | dataset | command | purpose |
|---|---|---|---|---|---|
| 2026-10-28 13:23:46Z | 4fe2e2a6ef35 | milestone-M7-2026-10-28 | T5 | `$VULCAN_ROOT/conda_envs/pe/bin/python bench/t5_armM.py --mode sealed --runs '{files}' --out '{root}/out' --work '{scratch}/T5/work'` | look=20261028T132346Z-4fe2e2a6 attempt=1 c9=GO a6=not tagged; A5 look (DESIGN §6.0d A5.6): arm M on T5 (PXD052597, 2 runs) for v1 (closed, bench/t1/ours.py Scorer) and v2_sep (chimeric, q_sep, merged strata deciding), 15 replicates, sealed entrapment seeds 2004-2006 |
| 2026-10-28 13:23:46Z | 4fe2e2a6ef35 | milestone-M7-2026-10-28 | T1 | `$VULCAN_ROOT/conda_envs/pe/bin/python bench/t1_driver.py --dataset T1 --out '{root}/out' --work '{scratch}/T1/work' --rt-off-control` | look=20261028T132346Z-4fe2e2a6 attempt=1 c9=GO a6=not tagged; A5 look (DESIGN §6.0d A5.5): T1 v1 primary (A1 comparators, §1 no-loss P1) + v2_sep pre-registered secondary arm (family S, 'Secondary arm; not gate-validated') + A5 additions; --a3 forced |
|  |  |  | T5 | result: look=20261028T132346Z-4fe2e2a6 exit=0, wall=1975s, manifest_sha256=6ec740302311bffaaae4e6ba23cd1090828405595128cba25397c7d71170544d |  |
|  |  |  | T1 | result: look=20261028T132346Z-4fe2e2a6 failed (unterminated) |  |
|  |  |  | look | result: look=20261028T132346Z-4fe2e2a6 leak_test=not run (unterminated); outputs deleted unread |  |
