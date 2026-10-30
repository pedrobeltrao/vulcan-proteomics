# A6 (HN1): implementation, dev dry run, a6_power re-run and tag criterion (A6.8.2 steps 2–3)

- **Role:** A6 implementer (fresh context, T1-blind). I did not write A6, I ran no selection, admission or gate step, and I listed no HN1 candidate on PRIDE. Commissioning prompt: `commissioning_prompt_implementer.md` (verbatim).
- **Data:** DEV only: D0 (PXD001468, 2 runs) and the 10 D8h runs (E9_0000–0009). Sealed data: none. Per-label quantities on null candidates: none. Full-index labels were read only by the certification module, only at rank-1 rows of each run's host full-index search (`admit_hn1.LabelView` guards this).
- **Date:** 2026-10-28, 12:24–13:40 CET.

## Result: "not tagged: infeasible" (a label-free stop). No project was listed, so no pool is burnt.

The result has two independent causes. Each one alone means "not tagged".

1. **Throughput rule (A6.8.2 step 3) fails.**
   - N is r × (16:30 − 15:15). The measured best case is N ≈ 145 runs (`throughput.json`).
   - With the dry-run yields, P(v1 minimum within 145 runs) is **0.63 at ×1** (needs 0.8) and **0.003 at ×0.5** (needs 0.5) (`power/tag_criterion.json`).
   - Both conditions need **N ≥ 240 runs**: P is 1.00 at ×1 and 0.57 at ×0.5. For a 15:15 tag that is r ≥ 192 admitted runs/h. The best case measured here is about 116/h effective.
2. **The `a6_stats` acceptance test (A6.8.3) fails item (iii).** The G2-original percentile bootstrap undercovers (`a6_stats_acceptance.json`):

   | Unit, heterogeneity | Coverage (nominal ≥ 0.95, 1,000 replicates, J = 37) |
   |---|---|
   | cluster, heterogeneous | 0.923 |
   | cluster, homogeneous | 0.932 |
   | project, homogeneous | 0.944 |

   - Even the full envelope covers only 0.938 in the heterogeneous cluster case.
   - Items (i) and (ii) pass: exact reuse of `a2_stats`, and reproduction of `a4_power.cap_factors` and `c9_power.tail_iut` to 1e-12.
   - The G1 jackknife-t (0.989–1.000) and the TOST bootstraps (0.954–0.977) pass.
   - Replacing the bootstrap type (for example bootstrap-t) would be a rule-level choice. I did not make it (A6.9).

## What was built (outside every frozen tree; `scripts/a6/`, `tests/a6/`, `config/a6.yaml`)

| File | Content |
|---|---|
| `admit_hn1.py` | Two-engine certification:<br>• hyperscore, per-run paired double competition, PSM q ≤ 0.001, peptidoform q ≤ 0.01, tie-break seed 11;<br>• Sage `tools/bin/sage` with the T1 plain config and DB = `prod_1001.fasta` (prod T + its decoys); rank-1 target, spectrum_q ≤ 0.001, same I/L sequence and charge.<br>**Symmetric I:** pair winners ≥ the run threshold, either label.<br>**\|I\| count assertion:** label-free maxima per pair key; must equal the I set.<br>**Project-wide R_P:** (i) near I, both members and isobaric variants; (ii) ladder f ≥ 0.05 over every row in any certified window; (iii) contaminant FASTA ∪ Panel B. Pair-closed.<br>**Freeze** files plus a sha256 manifest, written before the reduced index exists.<br>**Null stage:** reduced-index production search with a single hypothesis.<br>**Swap test** from the raw PSM tables on an index rebuilt with `swap_labels=True`.<br>**Label guard:** `LabelView`. |
| `reduced_index.py` | Row-subset `L1Index` with `row_map`; `GpuIndex` and `search` unchanged. Exactness test (below). |
| `select_hn1.py` | Exclusions, giving `config/a6_excluded_accessions.tsv`: pool **1,415**, as at drafting.<br>Walk seed from the tag sha; refuses before the tag. PRIDE listing also refuses before the tag.<br>File guard, cluster slots ≤ 5, stopping rule on capped counts, `selection_log.tsv`, `hn1.tsv`. |
| `mixtures_hn1.py` | Self-host mixture: host = full-index certification search; null copies (prov 20) and FS copies (prov 22) on the reduced index; concatenated before features. Label-free winners only.<br>**R12:** no v1 feature counts host PSMs, proteins or charges of a candidate peptide. The run-level features that see the twins (`hs_runz`, the RT residuals, the NCE choice) are listed, direction pro-target. |
| `a6_stats.py` | Envelope rows, cap, IUT precedence, stricter difference clause, subset check (acceptance test above). |
| `gates_hn1.py` | Step-5 one-session runner. The only reader of null-row labels. Refuses unless the tag, a clean tree, the committed step-4 record and its hashes are in place, and it is before 17:55. Writes `a6_status.json` through `status.py`. |
| `status.py` | Closed-enum writer and reader for `results/2026-10-28_A6_hn1/a6_status.json`. Absent or unverifiable reads as "stopped: node time limit". |
| `dry_run.py`, `throughput.py`, `a6_power_dryrun.py` | Dry run, throughput model, and the a6_power re-run with the tag rule. |
| tests | 67 fast tests pass: status, select, label guard, gates, admit with synthetic I. `test_a6_stats.py`: 23 pass, 3 fail, item (iii). |

**GPU and lock discipline.** A6 processes take `/tmp/vulcan_a6.lock` and never take the bench lock. The node has one A100. During the dry run it was shared with a T5 bench and a cohort run that used up to 60 GB, and one A6 run hit CUDA OOM and was retried. So the GPU **cannot be separated**, and under A6.8.2 A6 yields to the other users.

## Checks (dev)

**Reduced-index exactness** (`reduced_index_exactness.json`, D0, 3,000 single-hypothesis spectra, R = 9,576 rows) passes:
- R = ∅ is bit-identical to the full index;
- every returned row maps outside R;
- every hyperscore equals the full index's for the same (spectrum, row);
- the exhaustive table equals the full table minus R.

0.03% of production (spectrum, row) keys are not reported by either full search: new shortlist entries after removal. They have no full-index value to compare against.

**Swap test** passes on 5 distinct dev projects: PXD018572 (twice), PXD014823, PXD001468, PXD069854 and PXD057360.
- Label-free index files are bit-identical.
- I, R_P, A, A_cand and the tail proxy are identical.
- Part 2 of the yield run used `--no-swap` to save time.

**Count assertion** holds on every run. A unit test checks that a target-only I (the pilot's construction) fails it.

**An implementation bug found and fixed during the dry run** (before any tag; label-free):
- The first D8h2 run took the certified charge and isotope error from the unrestricted PSM table, which gave wrong z for some certified spectra. Sage agreement was 59% instead of 99%.
- It is fixed: both are taken from the rank-1 rows. Its log is kept as `dry_run/D8h2_events_buggy_z.log`, and its numbers are not used.

## Dry-run yields (final A6 rules, 12 DEV runs; `power/a6_cert_yield.json`)

| | Certified (2 engines) A | A_cand | Tail proxy |
|---|---|---|---|
| A6 drafting (`a6_cert_yield.json`) | 63,692 | 23,049 | 118 |
| **This dry run** | **63,666** | **16,651** (×0.72) | **63** (×0.53) |

- Sage agreement is unchanged.
- The losses come from the stricter-than-drafting removal that A6 adopted:
  - symmetric I with both pair members;
  - the project-wide ladder over every row in any certified window;
  - (iii) as A1.3-near to non-human contaminant and Panel B proteins. This flags 1.27 M of the 7.06 M index rows, against exact occurrence at drafting.
- The tail proxy, which binds the stop, roughly halves. The walk therefore needs more runs: P(v1 minimum) is 0.80 at 160 runs (×1) and 0.57 at 240 runs (×0.5).

## Full a6_power re-run on the dry-run yields (`power/a6_power_dryrun.{json,out}`; reviews/checks/a6_power.py unchanged; seed 131; 400 walks)

| Cell | Runs at the v1 minimum, p10 / 50 / 90 | J p50 | Drafting p50 |
|---|---|---|---|
| final rules, dry-run yields ×1 | 107 / **136** / 174 | 48 | 102 |
| ×0.5 | 183 / **234** / 293 | 80 | 163 |
| one-engine (sensitivity; drafting counts) ×1 | 74 / 92 / 116 | 33 | 92 |

- The binding criterion is still the capped tail proxy (1.02× its minimum at the stop).
- **Gate behaviour at the stop:**
  - clean tail-G2 pass 0.905 / 0.87 at κ = 1 / 0.5, with invalid pool 0.095 / 0.065;
  - G1 clean pass 0.99; G3 clean pass 1.00.
  - These are within MC error of the drafting figures. The obstacle is the run count, not the power at the stop.

## Throughput (`throughput.json`; pipelined D8h run, 10 parallel PRIDE downloads)

**Per run:**
- download: 3.5 MB/s per stream (0.6–1.5 GB in 137–465 s);
- fill latency to the first admission: 313 s;
- CPU (core-s): decode 6.7, TRFP 69, Sage 83 (1-run projects), near-test pool 90 per pass (×2 with the swap test). Total ≈ 352.

**Bounds on r (admitted runs/h):**
- serial GPU admission (26.7 s wall/run): 135;
- swap-test process (28.9 s/run): 124;
- download (0.197 GB/s aggregate): 593;
- CPU: 327 / 225 / 123 / 20 / 0 if other node users load 0 / 10 / 20 / 30 / 40 cores of 32.

**Load from other users during the dry run:** 23–57 (T5 bench, cohort/speed track, A5 work). The A5 freeze and G0 run overlap the 15:15–16:30 window.

**N for the 75-min window**, after fill latency: about 145 at other-load ≤ 20, 24 at other-load 30, 0 at 40.

**Not measured:**
- 3-run projects: the dry run used 1–2 runs per project. Sage amortises a little; R_P near-I grows with the window rows; D0's 2 large runs took 379 s in an older, slower build.
- A second admission process. It would lift the serial bounds but not the CPU bound.
- The killed part of the D8h run: an external process killed 4 TRFP and 3 curl processes at t ≈ 508 s, so 3 of 10 runs completed.

## What would make A6 taggable (for the PI; nothing here was changed)

1. **Statistics.** A G2-original bound with ≥ 0.95 coverage at J ≈ 37, for example a bootstrap-t or jackknife-t component like G1's, reviewed T1-blind. It would then be re-tested by the acceptance test.
2. **Clock.** N ≥ 240 runs needs about 2.0 h of walk at the best measured r (≈ 124/h), plus a node with other-load ≤ 20 during the walk. The window is only 75 min, so the walk would have to start by about 14:30 with a quiet node. That cannot happen under D31.

## Disclosure, gaps and deviations

- **n_min re-derivation (A6.6.4) not run.** The A4 provisional table is used (D32 default 8); this is a disclosed gap.
- **R_P scope.** R_P is evaluated over rows in the union of certified single-hypothesis windows widened by 0.5 ppm, plus their partners. Rows outside these windows cannot be candidates of any null search, so outputs are identical (exactness test). This is recorded as an implementation note.
- **(iii) reading.** "Occurs in" a human contaminant protein is exact I/L substring. A1.3-near is applied to non-human contaminant and Panel B proteins, as written.
- **Sage DB.** Sage used `$VULCAN_SCRATCH/vulcan_t1/D8h/fasta/prod_1001.fasta` (T1 prod DB of D8h, rev_ decoys).
- **Contaminant row cache.** The contaminant-near flag for all rows is cached at `/tmp/vulcan_a6/contam_rows_b540a3391edd0a12.npy`, with provenance (FASTA sha256s and index content hash) checked on load.
- **Label reads.** Labels were read only at host rank-1 rows. No step 2–4 output carries a per-label key; `test_no_label_reads.py` checks this.
