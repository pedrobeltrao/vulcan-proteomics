# T2 driver: dev dry runs and M4 E8 reproduction (2026-10-28)

Dev data only. **T2 exposure: none; T1 exposure: none.** No file under `data/sealed/` was opened. Sealed
truth/library files were only sha256-hashed by `bench/t2_driver.py --check`, and their contents were never printed.
Driver: `bench/t2_driver.py`, `bench/t2/` (choices in `bench/t2/IMPLEMENTATION_NOTES.md`), config `config/t2.yaml`.
The runs of record are `D2/` (look id `dev-d2-4`, git f6836b4) and `D2T_guard_tripped/` (`dev-d2t-3`, git 7bebe32). Both are fresh,
full runs: RAW → store → FASTA → index → search → E8 → claims, and they follow PI decision D35 (frozen localizer
unchanged, plus 10 repeat localizations for the run-to-run band). The working tree also held other agents'
uncommitted A5 edits (sealed_run.py, t1_driver.py, A3 freeze files). The T2 driver imports none of them.

## 1. D2 = M4's E8 set: reproduction (`D2/compare_m4.json`, `bench/t2/compare_m4.py`)
| Stage | vs M4 (results/2026-10-27_M4_ptm and its /tmp/vulcan_m4 work dir) |
|---|---|
| search FASTA | byte-identical (sha256) |
| PSM tables, 9 runs | `pyarrow.Table.equals` true for all 9 |
| rank-1 PSMs, peptidoforms at q≤1%, phospho PSMs, mono/multi, truth-eligible, search rank-1 site accuracy | identical (3 accessions) |
| truth tables (design / single-Y / proxy) | identical dicts (4,800 / 183 / 72 peptides) |
| pyAscore evaluation, A and L, every field | identical |
| our L3 localizer (frozen, non-deterministic) | M4's committed value lies inside the driver's 11-run range in 30 of 36 checked fields. The other 6 sit 1–2 sites or ≤0.0003 recall outside it, and inside the earlier 5-run M4-mode range (`m4mode_spread.json`) in 4 of those 6 |

**Run-to-run variability (D35).** `vulcan.ptm.localize.score_configs` sums peak weights with a CUDA atomic
`index_add_`. Exact ties between configurations (by construction, e.g. S vs A with identical ion sets) are then broken
by last-bit noise. On D2 the 11-run ranges are:
- PXD000138 Ala: 1% n 1739, true FLR 0.29% (constant); 5% n 2155–2158, true FLR 3.38–3.62%.
- PXD000138 Leu: 1% n 2131–2136, true FLR 1.31–1.54%.
- PXD009449 Ala: 1% n 23646–23676, true FLR 0.81–0.91%; 5% true FLR 4.39–4.49%.
- Recall at matched true FLR: up to ±0.002.

In D2 and D2T, no E8 or recall verdict changes across the 11 runs. M4's values can fall just outside an 11-run
range, so the printed band is a range, not a bound; the claim table says so.

A deterministic variant (`torch.use_deterministic_algorithms`) was tried first (`superseded_*_deterministic/`). It
moved three cells outside the noise range (exact ties form one block in the FLR̂ curve: PXD009449 1% 0.32%, PXD000138
Leu 5% jumps). Under D35 it is not used by the driver.

## 2. D2 claim table under the pre-registered stricter rules (`D2/claims.md`)
- **E8 gate (dev): "not passed".**
  - PXD000138 Ala: PASS at 1% (0.29%, UCB 0.63%) and at 5% (3.43% [3.38–3.62], UCB 4.58% [4.51–4.82]).
  - PXD009449 Ala: **inconclusive at 1%** (0.86% [0.81–0.91], UCB 1.34%) and at 5% (4.44%, UCB 6.93%). Only 183
    single-Y peptides carry this component, so the sequence-cluster bootstrap is wide.
  - Leu PXD000138 1%: FAILED (1.54% [1.31–1.54]; secondary, not a gate input).
  - PXD007058: n < 200 in every cell (not evaluable; never a gate input).
  - So M4's dev reading "Ala passes E8" holds on point estimates only. Under the one-sided-UCB rule it holds for
    PXD000138 only.
- **Recall vs pyAscore at matched true FLR (Ala), paired sequence-cluster bootstrap:**
  - PXD009449 1%: **higher** (0.886 vs 0.866, CI [0.011, 0.028]).
  - PXD000138 1%: no difference shown (0.880 vs 0.874, CI [−0.020, 0.058]).
  - PXD000138 5%: **lower** (0.914 vs 0.930, CI [−0.032, −0.003]).
  - PXD009449 5%: no difference shown.
  - So M4's "beats pyAscore at 1% on 6/6" (point estimates) becomes "higher" in 2 of 6 cells (PXD009449 A and L).
- Recall rows at estimated FLR̂ print "void: E8 gate not passed" (the computed values stay in the statistic).
- Completeness test: deleting pyAscore output from a copy of metrics.json gives "look failed: incomplete", exit 1, and
  the comparison rows read "not evaluable (comparator missing) -> not beaten". Re-running `claims.py` on the same
  metrics.json gives byte-identical claims.md.

## 3. D2T = T2-shaped dev analogue (`D2T/`)
The dev PXD000138 pools plus the four dev PXD009449 pools: Ymod_Phospho and the non-phospho Kmod_Dimethyl, Kmod_Formyl
and Rmod_Citrullin. T2's PXD009449 component consists **only** of non-phospho K-acyl pools (DATA.md: Ymod_Phospho was
pinned to dev). The non-phospho pools give 0 phospho PSMs at q ≤ 1% (row `D-nonphospho-PXD009449`). Their truth status
is "none (not a phospho kit)". Results: §5.

## 4. Incidents (shared node), all fixed in the driver
1. GPU OOM at the search while other agents held ~30 GB. The driver now waits for ≥ 64 GB of free GPU memory (logged,
   max 90 min).
2. /tmp filled up (6.9 GB free at one point). Scratch moved to /dev/shm. A D2T attempt was then OOM-killed in the index
   build: tmpfs is charged to the shared job memory cgroup (256 GiB, with all agents at 243–274 GB), and the build peaks
   at **115 GB RSS**. The driver now waits for 170 GB of cgroup headroom before the build. A watchdog stops the driver
   below 6 GB headroom or 2 GB of free scratch. `--check` reports GPU memory, scratch free space and headroom.
3. The leak test walks `data/cache`, which takes about 2 min. It flags any file written under results/, bench/,
   scripts/, config/, data/cache or data/truth during the look, **including files written by other agents**, so the
   T2 look needs the same exclusive node as the A5 look.

## 5. D2T results and runtime
- **`D2T_guard_tripped/` (dev-d2t-3) completed with exit 0**, but its memory watchdog tripped during the index build
  (cgroup headroom 5.0 GB, caused by other agents plus this run's own 115 GB). `interrupt_main` was lost inside the
  native call, so the run carried on. Fix 5696a9a: a tripped guard now fails the run at the next stage boundary. Results
  of that run:
  - E8 gate "not passed". PXD000138 Ala PASS at 1% (0.29%, UCB 0.63%) and at 5% (3.43% [3.43–3.71], UCB 4.58%).
    PXD009449 Ala inconclusive at 1% (0.83% [0.83–0.91], UCB 1.30%) and at 5% (4.22%, UCB 6.64%).
  - Recall at matched true FLR vs pyAscore: PXD009449 1% higher (0.895 vs 0.875); PXD000138 1% no difference shown;
    PXD000138 5% lower.
  - Non-phospho pools: 0 phospho PSMs.
  - No verdict changed across the 11 runs.
- **A clean re-run with the guard fix (dev-d2t-4) was stopped by me at the PI's request** (the T1/T5 look needs the
  node) during its E8 stage, after convert, index and search had finished. Its scratch was deleted and it has no
  outputs. The guard-fix code path is therefore untested end to end.
- **Measured dev duration of the T2-shaped pipeline (8 runs; D2T dev-d2t-3, shared node at load ~40):**
  - convert 6.6 s (8 parallel, 4 threads each);
  - index build 674 s (peak RSS 115 GB);
  - search 10 s;
  - E8 with 10 repeat localizations and the bootstraps 351 s;
  - total 1,057 s wall (17.6 min), 1,334 core-s.
  - D2 (9 runs): total 920 s.
  - Add about 2 min for the leak test in sealed mode, plus any waits for GPU or memory headroom.
  - On an idle node, expect about 20 min. Recommended time limit: ≥ 45 min (2× the dev run plus margin).
