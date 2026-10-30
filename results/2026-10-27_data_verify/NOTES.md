# Data verification, 2026-10-27 (32-core node)

Command: `python3 scripts/data/verify_data.py --out results/2026-10-27_data_verify/verify.tsv` (full table in `verify.tsv`).

| kind | OK | MISSING | PARTIAL | MISMATCH | total |
|---|---|---|---|---|---|
| fasta | 8 | 0 | 0 | 0 | 8 |
| raw (dev) | 44 | 0 | 0 | 0 | 44 |
| sealed raw (stat only) | 17 | 0 | 0 | 0 | 17 |
| truth | 43 | 0 | 0 | 0 | 43 |
| metadata (SDRF) | 2 | 0 | 0 | 0 | 2 |
| all | 114 | 0 | 0 | 0 | 114 |

Before the run (first check this session): 104 OK, 1 MISSING (D7 SDRF), 1 PARTIAL (D7 LM2487_DDA.raw .part, 713 MB of 1356 MB), 3 D7 raw files without a recorded sha256, 8 per-pool FASTAs not in any manifest.

What was done
- Only D7 LM2487_DDA.raw and the D7 SDRF were actually missing. The other interrupted-log items (D2 58.raw, T1 Alpha_02) had in fact been completed by later per-set runs (the old combined logs are stale); size and sha256 are correct. LM2487 was resumed by HTTP Range, and its SHA-1 matches PRIDE.
- Sealed files: 17 of 17 present at the exact expected size, all mode 000, never opened or chmod'ed by the verifier. Their sha256 (sealed.sha256, recorded at download) is not re-verified here.
- Per-pool truth FASTAs regenerated with make_pool_truth.py (byte-identical to before); sha256 now in `data/manifests/derived_truth.sha256`.
- PRIDE SHA-1 cross-check: `pride_sha1.tsv`, 37 OK, 0 mismatch, 28 without a PRIDE checksum.
- Script fixes in fetch_tier_s.py: `--jobs`, PRIDE SHA-1 check after download, locked/atomic cache and manifest writes (safe for parallel processes), manifest rewritten only when changed, partial sealed files made writable before resuming.
- Total on disk (OK rows): 61.3 GB (data/raw 29 GB, data/sealed 29 GB, FASTA 136 MB, truth 172 MB).

Unresolved: see the "Open issues" section of DATA.md (PXD056598, PXD080610, InfinDIA PXD005573 and MSV000093613 mapping). PXD036062 is a MassIVE dataset (explains the PRIDE 404).
