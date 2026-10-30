# T2 driver: implementation notes and stricter-reading choices

Author: T2-driver agent (fresh context), 2026-10-28. **T2 exposure: none. T1 exposure: none.** No file under
`data/sealed/` was opened. Files under `data/truth/sealed/` were only sha256-hashed by `t2_driver.py --check` against
the committed manifests; no content was printed or read by me.

Normative basis: DESIGN §6.5 (FLR), §6.6 (T2), RESEARCH_TREE D4 (predictor leakage), D8 (PSM-level FLR̂ primary,
site-level secondary; Ala primary, Leu secondary), M4 E8 dev pipeline and truth declarations
(`results/2026-10-27_M4_ptm/NOTES.md` §0, commit 97daf98; code `scripts/m4/e8.py`, `scripts/m4/truth.py`), A5 look
mechanics used as a pattern (A5.5.3 root/console/leak test, A5.5.4 completeness).

## 1. What the driver runs (frozen pipeline = M4 E8)
- Search DB: human SwissProt canonical + contaminants + the dataset's libraries (the M4 `d2_fasta` rule; on D2 it is
  byte-identical to M4's FASTA). T2: the four sealed PXD000138 pool FASTAs and the union of the four sealed PXD009449
  pools' `evidence.txt` sequences.
- Index and search: STY phospho, M-ox, peptide-N-term acetyl, ≤3 var mods; GPU closed search, report 1, params of
  `scripts/m4/bench.sh d2_phospho` (config/t2.yaml `search`).
- PSM list: localization-agnostic paired double competition on hyperscore, unit (sequence, composition), q ≤ 1%,
  seed 11 (imported unchanged from `scripts/m4/e8.py`).
- Localization: `vulcan.ptm.localize` with the M4-chosen LocParams (12 peaks / 100 Th, intensity weight, max frag
  charge 2, H3PO4 on every acceptor); candidates STY+A (primary) and STY+L/I (secondary).
- Comparator: pyAscore 1.0.1 at M4's settings, on the same PSMs, pinned by tree hash (persistent copy
  `data/cache/t2_pylib/pylib`; M4's copy lived in node-local /tmp).
- `scripts/m4/e8.evaluate` is called unchanged, so every M4 E8 number is recomputed with the same code; the T2
  additions are in `bench/t2/analysis.py`.

## 2. Ambiguities and the reading taken (stricter unless stated)
1. **T2's PXD009449 pools are K-acyl kits (Kmod_Acetyl/Succinyl/Ubiquitinyl/Crotonyl; DATA.md D2 note: Ymod_Phospho
   was pinned to dev).** They contain no phosphopeptides by design, and the M4 PXD009449 truth rule (single-Y,
   phospho-Tyr kit) applies only to Ymod_Phospho. §6.5's estimator is defined for phospho (N_STY/N_A). Reading:
   the PXD009449 T2 component has **no phospho truth** and is "not evaluable" for E8 and recall. It is still searched
   and reported descriptively (phospho PSMs at q ≤ 1% in pools without phospho, row `D-nonphospho-*`). No K-acyl
   localization/FLR is invented here: that would be a new, unregistered estimator. **Consequence: the T2 E8 gate and
   all T2 FLR claims rest on the four PXD000138 pools only** (76/G4, 25/C1, 24/B12, 61/F1; pool 61 is small). This
   needs the PI's attention before the look (open issue 1).
2. **E8 pass criterion.** §6.5 says "E8 is a gate" without a numeric rule; M4 read "true FLR ≤ nominal" on point
   estimates. Stricter, fixed now (config `stats.e8_rule`): per component, Ala, mono-phospho, PSM-level, threshold
   FLR̂ ≤ α from all localized sites, true FLR on truth-eligible accepted STY sites; **PASS iff the one-sided 95% UCB
   ≤ α**, UCB = max(Clopper-Pearson, sequence-cluster bootstrap 95th percentile; 2000 draws, seed 11); inconclusive
   iff point ≤ α < UCB; FAILED iff point > α; not evaluable iff no truth or n < 200 accepted eligible sites (A1 n ≥ 200
   convention). Gate PASS iff every evaluable gate component PASSes at both 1% and 5% (any FAILED → FAILED; otherwise
   "not passed"). Leu is reported with the same cell rule but is not a gate input (D8).
3. **"E8 is a gate before any FLR claim."** Rows that use FLR̂ thresholds (recall at estimated FLR) print
   "void: E8 gate not passed" unless the gate PASSes. Recall at matched **true** FLR does not use FLR̂ and is not voided.
4. **"Beats the comparator."** M4 compared point estimates. Stricter: paired sequence-cluster bootstrap of
   ours − pyAscore recall (thresholds re-chosen in every replicate); "higher" iff the 95% CI lower bound > 0. A missing
   or failed comparator is "not evaluable → not beaten", and makes the look incomplete (A5.5.4 pattern).
5. **§6.5 "pools that contain positional isomers of one sequence are excluded".** M4 dropped isomeric sequences;
   the literal text excludes the whole pool. The driver excludes the pool (no dev pool is affected; the PXD000138
   design has one slot per pool, so it cannot trigger there).
6. **§6.5 close-acceptor strata** (another candidate acceptor within 2 residues; M4 did not compute them) and **D8
   site-level FLR̂** are computed and reported descriptively (`D-strata-*`, `D-site-*`).
7. **Predicted features / Koina (D4, LESSONS 10).** The frozen E8 pipeline uses no predictor. Adding Prosit/MS2PIP
   features now would be a new pipeline chosen after M4 tuning, so none is used; there is nothing to pin. The driver
   fails if any of koina/koinapy/peptdeep/alphabase/deeplc/ms2pip is imported (`predictor_guard` in metrics.json).
8. **PSM-FDR labels.** The pipeline is neither v1 nor v2_sep, so A5's v1/v2_sep labels (L-S, L-Q, L-R, L-G) do not
   apply literally; A5.1.5 leaves T2 under its own preconditions. Stricter: every row that rests on the q ≤ 1% PSM list
   carries L-PSM (hyperscore TDC, not entrapment-checked on T2, no null-spectrum test), the A5 designation label L-D
   verbatim, the truth label and the dev-tuning disclosure (config `labels`).
9. **One look.** `sealed_run.py --dataset T2` does not enforce one look for T2 (only T1/T5 via `--look a5`). The driver
   refuses unless the ledger file holds exactly one T2 line and git history exactly one `sealed-ledger: T2 at` commit,
   the tree is clean, `$SEALED_DATASET` is T2, and the out dir does not exist.
10. **PXD007058** (dev proxy truth) is part of the D2 reproduction only; it never enters a gate (M4 NOTES §0).

## 3. Run-to-run variability of the frozen localizer (PI decision D35)
`vulcan.ptm.localize.score_configs` sums matched-peak weights with a CUDA `index_add_` (atomic, order-dependent
float64). Exactly tied configurations (e.g. S vs A with identical ion sets) then get scores that differ in the last bits,
so best site and site delta (≈1e-16 instead of 0) vary from run to run (dev evidence:
`results/2026-10-28_T2_driver_dryrun/m4mode_spread.json`; M4's committed values lie inside the 5-run M4-mode range).
I first ran the localizer under `torch.use_deterministic_algorithms(True)`; that moved three dev cells outside the
noise range (exact ties form one block in the FLR̂ curve). **PI decision D35: the look uses the frozen path unchanged.**
So the driver calls `scripts/m4/e8.run_ours` as is. Run 0 is the result of record; 10 further localizations of the
same PSMs (config `stats.localization_repeats`) give a [min, max] band that the claim table prints next to every
number of ours, plus a verdict-stability note ("same verdict in all 11 runs" or the verdict counts). Verdicts are
run 0's; nothing is averaged or selected across runs. pyAscore is deterministic (checked: identical evaluations to M4).
The reproduction check is therefore: deterministic stages bit-identical, and M4's committed numbers inside the band.
`run_ours_deterministic` stays in analysis.py only for that diagnosis; the driver does not call it.

## 4. Look mechanics
- Root outside the repo: `config/t2.yaml look_root` (`<root>/<look-id>/{out,logs}`, mode 0700). The sandboxed agent
  cannot create it (read-only for agents); the operator's `--check` verifies it is writable before the ledger line.
- Node-local work in `/dev/shm/vulcan_look_t2_<look-id>` (staged RAW, stores, FASTA, index ~30 GB, PSMs): off /tmp
  (full twice on 2026-10-28; D35). tmpfs is charged to the shared job memory cgroup (256 GiB) and the index build peaks
  at 115 GB RSS, so the driver waits for 170 GB cgroup headroom before the build (logged, max 90 min), waits for 64 GB
  free GPU memory before the search, and a watchdog stops the driver below 6 GB cgroup headroom or 2 GB scratch free
  space (a failed run, rather than the kernel OOM-killing another agent; a dev D2T attempt was OOM-killed this way).
  Deleted at the end in sealed mode.
- Console: fds 1/2 go to `<root>/<look-id>/logs/driver.log`; the terminal shows only whitelisted templates.
- Leak test (sealed): after deletion, no file newer than the start under `results/` (ledger excepted), `data/cache`,
  `data/truth`, `bench`, `vulcan`, `scripts`, `config`, `/tmp/vulcan_m4`; the scratch dir must be gone.
- Completeness (claims.py, checked by the driver): every component, both decoys, both tools, all runs, every verdict
  cell and bootstrap present; else "look failed: incomplete" and exit non-zero.
- `--check` (operator, before `sealed_run.py`): look root writable, pyAscore hash, sealed truth/library hashes,
  sealed RAW present (existence only), GPU idle, /tmp free space.
