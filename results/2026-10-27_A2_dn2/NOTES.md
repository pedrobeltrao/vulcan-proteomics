# A2 DN2, steps 2–4 (2026-10-27/28): 6 projects admitted, below the minimum; step 4 not run

**Result (label-free).** The single D16 re-walk (seed 23, `select_dn2.py` with the D16 fix d6d199c; operator wrapper
`scripts/a2_ops/run_dn2_walk.py` from a clean worktree at 8fb2330) exhausted all 58 frame candidates and **admitted 6
projects (18 runs)**. The A2.9.3 stopping criteria were **not met**:

| criterion | required | reached |
|---|---|---|
| admitted projects | ≥ 10 | **6** |
| \|A\| | ≥ 120,000 | **24,581** |
| A spectra with an exhaustive candidate | ≥ 72,000 | **9,608** |
| hyperscore tail proxy (≥ 33.875) | ≥ 400 | **22** |

Per D16 point 4 ("fewer projects than the minimum → record and stop; that is the result") **the mixtures were not
built (step 4 not run)**. Steps 4b and 5 are not run by this agent (D22). No per-label quantity, tail statistic or gate
quantity was computed, printed or stored at any point.

## Admitted projects (walk order)
| accession | organism | own proteome | runs | spectra | certified | \|A\| | A with cand. | tail proxy | markers (trigger) |
|---|---|---|---|---|---|---|---|---|---|
| PXD006818 | Streptomyces lividans tk24 | UP000028682 | 3 | 40,150 | 1,983 | 352 | 143 | 3 | B|CSN1S1|9913 |
| PXD018953 | Caballeronia mineralivorans | UP000035963 | 3 | 141,124 | 62,234 | 5,561 | 1,634 | 2 | none |
| PXD044078 | Bacillus pumilus | UP000001355 | 3 | 160,805 | 49,889 | 6,494 | 2,975 | 4 | none |
| PXD011591 | Rhodobacter capsulatus | UP000002361 | 3 | 194,576 | 32,954 | 2,750 | 2,022 | 2 | none |
| PXD062500 | Deinococcus radiodurans | UP000002524 | 3 | 382,856 | 47,923 | 9,037 | 2,771 | 10 | H|ALB|9606 |
| PXD038243 | Emiliania huxleyi | UP000013827 | 3 | 86,929 | 8,516 | 387 | 63 | 1 | B|CSN1S1|9913 |
| **total** | | | 18 | 1,006,440 | 203,499 | **24,581** | 9,608 | 22 | |

Certified spectra removed before A (label-free screens): (b) f ≥ 0.20 111,787; (b) Hamming 813; (d) family 21; (c-pre) A1.3 66,297. B candidates (pending c-final): 981,859. FS donors would be the 24,581 A spectra.
One admitted run (PXD006818 / SL322_set1_rep5) has 0 certified spectra. All of its spectra are B candidates
(DEVIATIONS item 17).

## Walk outcome (58 candidates)
29 protocol text/SDRF (verbatim keyword list, §6.0b point 2); 10 genus slot taken; **6 marker-screen exclusions**
(Fervidobacterium ×4: PXD054269, PXD054285, PXD054267, PXD054272; Chlorella PXD015681; Frankia PXD018045; 2–17
non-keratin H/B members each); 6 admitted; 3 no reference proteome (genus-only or `sp.` organisms, D16 / item 13);
4 too few usable RAWs. Full log: `selection_log.tsv`, `data/manifests/dn2_select.log.json`, `walk_events.jsonl`.

## Checks
- **P1 swap test (PASS):** `admit.swap_test` on PXD006818 (the first admitted project). The gate index was rebuilt with
  T and D swapped. The label-free index files are identical, and the membership sha256 is identical
  (a05b78f9…37a2) and equal to the committed one.
- **Reproducibility across restarts:** all 12 project decisions and all 36 file sha256 values are identical in every
  attempt that reached them (2, 2r, 2r2, 2r3, 2r4).
- **Downloads:** 18 admitted RAWs. 15 match the PRIDE SHA-1; 3 (PXD062500) have no PRIDE checksum (size + sha256,
  identical in 2 independent downloads).
- **Host threshold:** 33.875 (max of the D0 hosts' hyperscore 1% thresholds), identical in the main checkout and the worktree.
- **MS1 availability (A3.1.5):** every one of the 18 admitted runs' stores contains MS1 scans (4,272–9,062 per run;
  `ms1_availability.json`). Stores, RAWs and admission files are retained persistently in `data/cache/a2/`.

## Paths for any later step
- Membership and counts: `results/2026-10-27_A2_dn2/{membership.sha256, admission_counts.json}` (commit in `git_sha`).
- Persistent data (gitignored): `data/cache/a2/{admission/<acc>/*.npz|counts.json|marker.json, store_dn2/*.vspec.arrow(.json),
  raw/dn2/<acc>/*.raw, host_threshold.json, swap_test.json}`. Manifests: `cache_admission_stores.sha256`, `cache_raw.sha256`.
  Working copy: `/tmp/vulcan_a2` (node-local; the `.npz` `store` field points to `/tmp/vulcan_a2/store_dn2/`).

## Deviations
`DEVIATIONS.md` items 1–21. The most consequential: item 11 (run 1 proteome-query bug → D16), item 13 (`sp.` =
genus-only), items 16–20 (clean worktree; zero-certified-run guard; transient-download and ENOSPC resumptions).
