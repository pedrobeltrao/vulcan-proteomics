# M1e comparator baselines (dev data only), 2026-10-27, node gpu-node

Sage 0.14.7 (built-in LDA), Sage+mokapot 0.10.0, Comet 2026.01+mokapot on D0 (PXD001468, 2 HEK files), D1 (PXD001077, Pfu) and D4 (PXD016766, TMT; 2 of 12 files, see Deviations).
Nothing under data/sealed was touched. All numbers: metrics.json; tables/*.tsv. Fingerprint: timings.jsonl (`hw`), env.txt. Node: 32 affinity cores of 2x EPYC 9354 (Zen 4), 1.5 TB RAM, A100 not used by these tools.
Code state: git_sha (scripts committed before the result files).

## Main tables (alpha = 1%)
Peptides = unique I/L-collapsed sequences at the tool's own peptide-level q <= 1% in the production search (protein FASTA, engine-native decoys, no entrapment). Proteins@1% = vulcan.stats.protein.picked_group_fdr (subset grouping, rescue on) on the tool's peptide scores, target groups at q <= 1% (contaminant entries are counted like any target). Verdicts are the §6.3 rule (static low-capacity scorer, required arms P and C) from the augmented searches: P = paired shuffled (3 seeds 2001-2003), C = combined natural foreign proteome (Arabidopsis for human; human for Pfu). `*` = verdict is "inconclusive" but the lower 95% CI of the FDP estimate is <= alpha (the weaker "consistent with control" claim of §6.3, labelled as such). Wall / core-s = search + rescoring of the production workload (conversion excluded, reported below), at 32 threads and at 4 threads, benchmark lock held, numbers from one session on one node.

#### D0 (PXD001468, Q Exactive, HCD, Orbitrap MS2)

| tool | peptides@1% | proteins@1% | verdict P (pep) | verdict C (pep) | verdict P (prot) | verdict C (prot) | wall s / core-s @32 thr | wall s / core-s @4 thr |
|---|---|---|---|---|---|---|---|---|
| sage | 22507 | 6320 | inconclusive* | supported | inconclusive* | inconclusive* | 34 / 239 | 46 / 140 |
| sage+mokapot | 22556 | 6335 | inconclusive* | supported | inconclusive* | inconclusive | 52 / 387 | 63 / 202 |
| comet+mokapot | 22377 | 6747 | inconclusive* | supported | inconclusive* | inconclusive* | 93 / 1734 | 315 / 900 |

Raw -> mzML (TRFP, counted separately; both Sage and Comet need it): wall 459 s, core-s 501, max RSS 3667 MB.

#### D1 (PXD001077, LTQ Orbitrap Velos, HCD, Orbitrap MS2 (R=7500))

| tool | peptides@1% | proteins@1% | verdict P (pep) | verdict C (pep) | verdict P (prot) | verdict C (prot) | wall s / core-s @32 thr | wall s / core-s @4 thr |
|---|---|---|---|---|---|---|---|---|
| sage | 5904 | 834 | inconclusive* | supported | not-evaluable | not-evaluable | 5 / 55 | 5 / 11 |
| sage+mokapot | 5899 | 837 | inconclusive* | supported | not-evaluable | not-evaluable | 24 / 115 | 19 / 56 |
| comet+mokapot | 5912 | 852 | supported | supported | not-evaluable | not-evaluable | 23 / 217 | 38 / 112 |

Raw -> mzML (TRFP, counted separately; both Sage and Comet need it): wall 38 s, core-s 37, max RSS 594 MB.

#### D4 (PXD016766, Orbitrap Eclipse, TMT11 SPS-MS3, MS2 CID ion trap)

| tool | peptides@1% | proteins@1% | verdict P (pep) | verdict C (pep) | verdict P (prot) | verdict C (prot) | wall s / core-s @32 thr | wall s / core-s @4 thr |
|---|---|---|---|---|---|---|---|---|
| sage | 9969 | 4288 | inconclusive* | inconclusive* | supported | inconclusive* | 31 / 233 | 38 / 108 |
| sage+mokapot | 10109 | 4358 | inconclusive* | inconclusive* | inconclusive* | inconclusive* | 42 / 329 | 48 / 151 |
| comet+mokapot | 10166 | 4579 | inconclusive* | supported | supported | supported | 43 / 401 | 114 / 315 |

Raw -> mzML (TRFP, counted separately; both Sage and Comet need it): wall 136 s, core-s 127, max RSS 604 MB.


Peptide-level counts at 0.1/0.5/5% and PSM counts (labelled not FDR-controlled) are in tables/ids.tsv; FDP-hat, lower bound, CIs and verdicts at every alpha in tables/entrapment_peptide.tsv and tables/entrapment_protein.tsv.

## Verdicts and what they mean
- No tool is "violated" at any gated alpha on any dataset, at peptide or protein level (lower CI of the lower bound never exceeds alpha).
- Arm P (paired, upper-bound estimator) is "inconclusive" for most tool x dataset cells at peptide level, with mean FDP-hat about 1.0-1.1% at alpha = 1% and a lower CI below alpha. That is the behaviour §6.3 predicts for a calibrated list (the paired estimator overstates by ~20%), so these read "consistent with control", not "supported". Only Comet+mokapot on D1 is "supported" at peptide level for P.
- Arm C is "supported" at peptide level for all three tools on D0 and D1, and for Comet+mokapot on D4; Sage and Sage+mokapot on D4 are inconclusive* (mean FDP-hat 1.0-1.1%).
- Protein level: D1 is "not-evaluable" (alpha * N_T = 8 < 20). D0/D4: mixed supported / inconclusive*. Arm C protein-level uses r from the effective peptide space (1.022 for Arabidopsis; 12.94 for human-vs-Pfu), as §6.3 specifies. Sensitivity: the ratio of proteins (1.89 human/Arabidopsis; 8.39 Pfu) would turn several D0/D4 protein-C cells to "supported" (it lowers FDP-hat), so the peptide-space r is the conservative choice and is what the table uses (metrics.json: `sensitivity_r_protein_count`).
- Sensitivity (production, 1%): peptides Sage / Sage+mokapot / Comet+mokapot = D0 22507 / 22556 / 22377, D1 5904 / 5899 / 5912, D4 9969 / 10109 / 10166. Protein groups = D0 6320 / 6335 / 6747, D1 834 / 837 / 852, D4 4288 / 4358 / 4579. Comet+mokapot is highest at protein level on all three (+2 to +7% over Sage-LDA) and on D4 peptides (+2%), lowest on D0 peptides (-0.6%). No superiority claim: one dataset each, a single run or two, entrapment CIs overlap.
- Decoy-variability arm (external shuffled decoys, DECOY_SEEDS 1001-1005, tables in metrics.json `decoy_variability_S_external_shuffled_1001_1005`): peptides@1% mean (sd) D0 Sage 22585 (41), Sage+mokapot 22681 (43); D1 Sage 5906 (4), Sage+mokapot 5916 (8), Comet+mokapot 5928 (11); D4 Sage 9941 (32), Sage+mokapot 10112 (40). These are within 0.6% of the production values. One D1 Sage+mokapot seed gave a collapsed protein count (mean 660, sd 369 over 5 seeds) because the mokapot q-value curve of that seed is non-monotone at protein level; reported as observed.

## Speed (DESIGN §2): per stage, 32 vs 4 threads
Input path: **raw Thermo -> indexed mzML via ThermoRawFileParser 2.0.0.0 (vendor peak picking, zlib), staged on node-local NVMe /tmp, then the tool reads mzML**. Both Sage and Comet need the mzML, so conversion cost is counted separately and not in the table columns above. Conversion, 1 thread, wall / core-s / max RSS: D0 (2 files) 459 s / 501 / 3.7 GB; D1 (1 file) 38 s / 37 / 0.6 GB; D4 (the 2 files used) 136 s / 127 / 0.6 GB. The raw copy to /tmp (cp) is untimed.
MS2 counts, per-tool spectra/s and core-s per 1e5 MS2 are in metrics.json (`speed`). "Rescore" = mokapot (single worker; its BLAS threads follow the thread budget, which is why its core-s at 32 threads is inflated).

| dataset | tool | DB | stage | wall s @32 | core-s @32 | RSS MB @32 | wall s @4 | core-s @4 | RSS MB @4 |
|---|---|---|---|---|---|---|---|---|---|
| D0 | sage | prod | search | 34 | 239 | 5717 | 46 | 140 | 4816 |
| D0 | sage+mokapot | prod | rescore | 18 | 148 | 583 | 17 | 62 | 526 |
| D0 | comet | prod | search | 78 | 1592 | 5409 | 300 | 844 | 5139 |
| D0 | comet+mokapot | prod | rescore | 15 | 142 | 507 | 15 | 56 | 453 |
| D0 | sage | P2001 | search | 71 | 440 | 10114 | 95 | 258 | 8916 |
| D0 | sage+mokapot | P2001 | rescore | 21 | 149 | 587 | 19 | 64 | 534 |
| D0 | comet | P2001 | search | 1384 | 6979 | 5371 | n/a | n/a | n/a |
| D0 | comet+mokapot | P2001 | rescore | 15 | 140 | 502 | n/a | n/a | n/a |
| D1 | sage | prod | search | 5 | 55 | 473 | 5 | 11 | 394 |
| D1 | sage+mokapot | prod | rescore | 19 | 60 | 308 | 14 | 44 | 310 |
| D1 | comet | prod | search | 14 | 168 | 4890 | 28 | 73 | 4626 |
| D1 | comet+mokapot | prod | rescore | 9 | 50 | 306 | 10 | 40 | 305 |
| D1 | sage | P2001 | search | 7 | 71 | 879 | 7 | 16 | 780 |
| D1 | sage+mokapot | P2001 | rescore | 11 | 54 | 310 | 8 | 39 | 310 |
| D1 | comet | P2001 | search | 41 | 301 | 4896 | 41 | 111 | 4635 |
| D1 | comet+mokapot | P2001 | rescore | 11 | 52 | 307 | 9 | 40 | 307 |
| D4 | sage | prod | search | 31 | 233 | 4419 | 38 | 108 | 3743 |
| D4 | sage+mokapot | prod | rescore | 11 | 96 | 454 | 9 | 43 | 394 |
| D4 | comet | prod | search | 32 | 309 | 2749 | 103 | 270 | 2732 |
| D4 | comet+mokapot | prod | rescore | 11 | 92 | 402 | 11 | 45 | 347 |
| D4 | sage | P2001 | search | 64 | 410 | 8839 | 80 | 237 | 7504 |
| D4 | sage+mokapot | P2001 | rescore | 11 | 95 | 452 | 10 | 45 | 396 |
| D4 | comet | P2001 | search | 784 | 3339 | 2723 | n/a | n/a | n/a |
| D4 | comet+mokapot | P2001 | rescore | 10 | 93 | 405 | n/a | n/a | n/a |

Findings:
- Extra threads buy little on these small inputs and cost a lot of cores. Sage on D0 production: 34 s / 239 core-s at 32 threads vs 46 s / 140 core-s at 4 threads (wall 1.4x faster, 1.7x more core-s). On D1 (15k MS2): 5 s either way, 55 vs 11 core-s (5x). On D4 (2 files): 31 vs 38 s, 233 vs 108 core-s. Reason: mzML parsing and index building are serial; the fair per-core comparator figure for small runs is the 4-thread one.
- Comet scales better: D0 production 78 s at 32 threads vs 300 s at 4 (3.8x faster for 1.9x the core-s); D4 32 s vs 103 s. Comet is ~6-7x more core-s than Sage on D0/D4 production searches (4-thread: 844 vs 140 on D0).
- Entrapment doubles+ the search space: Sage D0 P2001 (peptide-level DB, 4.8M targets+entrapment, 9.7M with decoys) costs 440 core-s at 32 threads (1.8x production). Comet on the same peptide-level DB (no-cut enzyme, no protein-level pre-index) costs 6979 core-s / 1384 s on D0 and 3339 core-s on D4, 5-20x production. Treat Comet augmented times as set-up specific (peptide-level DB, searched without cleavage), not a Comet property.
- Timing quality: every timed prod and P2001 run held the benchmark lock, but other agents' unlocked jobs shared the node throughout (load average at start of the locked runs: median 12, range 6-26 at 32 threads; median 16, range 12-20 for the 4-thread pass). Treat wall times as +-20-30%; core-s (user+sys) is less sensitive. The 3 correctness DBs per tool (target, P2002, P2003, C, S*) are timed but unlocked.
- Comet 4-thread runs were done on production workloads for all datasets and on P2001 only for D1 (Comet on a 2.4M-peptide DB at 4 threads is ~80 min/file). Sage 4-thread runs cover prod and P2001 on all datasets.
- All Sage runs used `predict_rt: true` (RT prediction on; LDA uses it), chimera off, `report_psms: 1`, `--batch-size = min(n_files, threads)`.

## Entrapment databases (scripts/entrapment; spec check)
build_dbs.py was reviewed against §6.1/§6.3. Implemented as specified: unique-against-all-targets, C-term-fixed pseudo-reverse decoys (palindromes and collisions replaced by hash(peptide, 1001) shuffles), paired shuffles from FDRBench 1.1.1 (-check, seeds 2001-2003, r = 1), hash-seeded order-independent repairs, arm C with I/L, single-substitution (Hamming-1) and isobaric filters, r on the effective space with per-mass-decile r.
Deviations found and fixed (outputs of the P and C DBs changed; not an output-identical refactor):
1. §6.1 requires every entrapment peptide to be unique against all targets **and all decoys of all decoy seeds** (1001-1005); the old code only checked the seed-1001 pseudo-reverse decoys. Now E is also checked against the hash(peptide, seed, 0) shuffled decoys of every target for the five decoy seeds, and no shuffled decoy of an entrapment may equal a target or another entrapment (peptides.py `shuffle_decoy_set`, `shuffle_hits`). Short peptides collide often (E is a composition twin of its target): human P2001 needed 25,650 repairs (1.06%; 15.5k shuffle-decoy collisions, 10.0k E equal to a shuffled target decoy) and 1,668 pairs could not be repaired within 200 draws and stay as unpaired targets (1,839 in total; 0.08%). Pfu P2001: 1,451 repairs, 40 unpaired.
2. The isobaric filter (N/GG, Q/AG) used an order-dependent greedy replacement that missed e.g. AGG vs QG. Now both sides are expanded N->GG, Q->AG and each G/A run is sorted (slightly broader than the spec, conservative: only removes more foreign peptides). Arm C human: 441 isobaric removals, 67,175 single-substitution removals, 4,533 I/L-identical.
3. A second pass of the repair loop originally trusted pre-computed collisions for replacement candidates; the independent verifier caught candidate-vs-candidate collisions. Fixed (SEadd tracking) and verified.
4. New: scripts/entrapment/verify_dbs.py re-derives every property from peptides.tsv with exact-string keys (20/20 DB files pass: results .../entrapment/verify.log); tests/entrapment (7 tests) cover shuffles, order independence, isobaric key, paired and foreign arms, shuffled-decoy arm. tests/index passes (4 passed): `hash_shuffle`, `pseudo_reverse`, `build_decoys` are unchanged, so decoy outputs are identical to what vulcan/index mirrors; only the entrapment/foreign-peptide selection changed.
5. New arm S (DESIGN §6.2 decoy variability): `build_dbs.py --arms S` writes T + 5 external shuffled-decoy DBs (hash(peptide, seed, k), C-term fixed, re-draw on collision).
Remaining differences from the spec: P entrapment starts from FDRBench's own (Java-seeded, order-dependent) shuffle, with hash-seeded repairs for violations; fully hash-seeded generation would be order-independent throughout (vulcan/index already does that). FDRBench-unpaired targets (175 for P2001 before repair, 1,839 after the stricter rule) have no entrapment partner. The C arm for D1 uses human as the "foreign" proteome (r = 12.9), as built in the previous session.

## Other set-up notes and deviations
- **D4 uses 2 of 12 files** (dq_00095 A1 and dq_00101 B1). Comet on the peptide-level DBs is 5-10+ min per file per DB at 32 threads; 12 files x 5 DBs x 2 thread budgets was not feasible. D4 numbers are not the full dataset; all three tools use the same two files.
- Production searches use each engine's own decoys (Sage reversed peptides via generate_decoys, Comet decoy_search 1) and trypsin on the protein FASTA. The P/C/target/S databases are peptide-level FASTAs with external decoys searched with no cleavage (Sage cleave_at "$", Comet enzyme 11). Sage writes only PSMs with q <= 1, so the decoy tail beyond that is absent from protein-level inference.
- External shuffled-decoy comparison (reviewer B1): the harness supports it now via `--dbs S1001..S1005` (Sage generate_decoys=false with the external decoy FASTA). Done: Sage and Sage+mokapot on D0, D1, D4; Comet+mokapot on D1 only; IDs only (no entrapment verdict on S yet). Not done: Comet on D0/D4 with S; S combined with P/C entrapment (entrapment peptides are built against the five shuffled decoy seeds already, but the searches use pseudo-reverse decoys); the matched-FDP comparison. A T1 driver can reuse run_searches.py (`--dbs`) and build_dbs.py.
- mokapot 0.10.0 writes per-file q-values when given several PINs: Comet multi-file runs use `--aggregate` (joint q-values like Sage's single PIN). The first D0 Comet+mokapot attempt used the default and was discarded.
- mokapot env, entrapment DBs, mzML and search outputs live in /tmp/vulcan_m1e (rebuilt by scripts/baselines/run_all.sh: setup_mokapot_env.sh, scripts/entrapment/run_all.sh, 01_convert.sh); DB manifests/checksums are in entrapment/.
- Sage binary prints 0.14.6 (crate string not bumped); it is the v0.14.7 tag build. Stock binary used (sage-znver4 gave identical results, no speed gain per hw_adapt).
- Timings are not comparable with 2026-10-26_M1e_baselines (4-core node) and are not mixed.
- Not done: Sage 4-thread runs of the non-timed DBs, Comet 4-thread P2001 on D0/D4, per-decile r tables beyond metrics.json, MSFragger (unlicensed), run-level bootstrap (single or two runs per dataset), RESET/Group-walk variants.
