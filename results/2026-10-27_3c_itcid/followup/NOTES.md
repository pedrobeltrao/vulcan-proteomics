# 3c follow-up: pinned predictions (G0), strongest comparator, IT-CID fresh nulls

Status: done (dev only). No sealed file touched. The frozen `itcid_linear_pred` and its training/selection, §6, metrics,
labels, controls and tests are unchanged. Seeds: mokapot 11, bootstrap 11, DNIT walk 47 (§6.1 SEEDS[2]).
Files: `metrics.json` (everything), `g0_pinned.json`, `pinned_{verify,rebuild,manifest}.json`, `pinned_vs_3c_features.json`,
`eval_pinned.json`, `ratios.json`, `comparator_timings.json`, `NULL_RULES.md` (+ dated notes 1–4),
`dnit_walk_events.jsonl`, `dnit_host_threshold.json`.
Code: `vulcan/predict/pinned.py`, `scripts/m3c/{pinned_cache,g0_pinned,sage_pred_itcid,followup_eval}.py`,
`scripts/a2_itcid/{walk_dnit,gates_dnit}.py`, `config/a2_itcid.yaml`.

## 1. Pinned prediction store and G0
- **Store.** `vulcan.predict.pinned` holds the raw float32 Koina outputs, keyed by (model, fragmentation, CE, charge,
  peptide string as sent to Koina). Shards are immutable and named by their sha256. `MANIFEST.json` pins the snapshot,
  and every shard is sha256-verified on load. A key is fetched from Koina at most once (flock'd appends). Mode "frozen"
  cannot call Koina.
  - Location: `data/cache/koina_pinned` (persistent, gitignored; manifest copy committed here).
  - Size: 20.9 M entries in 12 tables, 3.0 GB.
  - Snapshot `416ea7bf…`; verify: all shard hashes and table digests OK.
  - `PinnedPredCache` replaces the node-local pf_id cache for every consumer.
- **Rebuild** (`pinned_cache.py rebuild`): it re-queried Koina for all 20.9 M keys into a new store (34 min).
  **590,146 entries (2.8%) came back with different floats** (max |Δ| 3.1e-4). For charge 5–6, 74–100% of entries differ.
  So Koina is not reproducible across sessions for the label-free CID model either, not only for TMT. The rebuild script
  therefore reports differences and never overwrites the pinned snapshot: the pinned snapshot is the reference.
- **Every dev replicate re-featurized from the store** (D4 and D10cf × prod, E2001–3, C; 205 runs).
  - 0.4–1.7% of feature values differ from 3c's (fresh-Koina) features, by ≤ 4e-4.
  - IDs at 0.1/1/5%, sequences, proteins, arm-P/arm-C FDP and every verdict are **identical to 3c** for hyperscore and
    `itcid_linear_pred` on both sets (D4 60,873 / 59,082 seq / 7,858 proteins; D10cf 72,160 / 65,786 / 5,394).
- **G0 (strict, both sides read the pinned store).** Search both sides, fill the store once with the union of keys,
  then let each side compute every feature and all 4 frozen scores with its own new frozen-mode store instance.
  - LF D10cf: **PASS**. 162,420 PSMs; PSM tables, all 11 features and all 4 scores bit-identical.
  - TMT D4: **PASS**. 77,247 PSMs, all bit-identical. 3c's strict fail is resolved by pinning.
  - Total 184 s.

## 2. Strongest comparator (A1.9): Sage + the same Prosit features + seeded mokapot
- **Setup.**
  - Sage: 3c's Sage 0.14.7 searches, unchanged (M1e DBs prod/P2001–3/C, ±0.5 Da, TMT statics on D4, MGFs from our stores).
  - Features: our 10 IT-CID features (incl. `pred_missing`, label-blind, RT calibrated per run like ours), computed on
    every Sage PSM from the pinned store and appended to the PIN.
  - Rescoring: `bench/t1/mokapot_seeded.py`, seed 11, canonical row order.
  - Also re-run: seeded mokapot on the plain PIN (3c's run was effectively unseeded).
- **Units and CIs.** Sequence units (I/L), as in 3c. Ratio = run-level IDs summed over runs, as `bench/t1/evaluate.py`
  computes it. CI = percentile bootstrap over runs (2000 draws, seed 11). 3c has a single decoy set, so J = 1.
- **Matched FDP.** Each method's own-q threshold is moved until its arm-P paired FDP (sequence level, pair key = target
  sequence, mean of 3 seeds) reaches f.

| ours / comparator | D4 own q 1% | D4 matched FDP 1% | D10cf own q 1% | D10cf matched FDP 1% |
|---|---|---|---|---|
| **Sage + Prosit features + mokapot (strongest)** | **1.036 [1.033, 1.039]** | **1.034 [1.030, 1.037]** | **1.090 [1.085, 1.094]** | **1.093 [1.088, 1.097]** |
| Sage + mokapot (seeded) | 1.074 [1.069, 1.078] | 1.077 [1.072, 1.081] | 1.140 [1.135, 1.146] | 1.128 [1.123, 1.133] |
| Sage (LDA) | 1.081 [1.077, 1.085] | 1.085 [1.082, 1.090] | 1.137 [1.133, 1.142] | 1.126 [1.122, 1.131] |

- **Other thresholds** (ours / strongest):
  - own q 0.1%: D4 1.031, D10cf 1.081.
  - matched FDP 5%: D4 1.047 [1.042, 1.051], D10cf 1.098 [1.093, 1.103].
- **Dataset-level sequences at 1%:**
  - D4: ours 59,082; Sage+pred+mokapot 57,764; Sage+mokapot 55,622; Sage 55,151.
  - D10cf: ours 65,786; Sage+pred+mokapot 63,347; Sage+mokapot 55,393; Sage 55,640.
- **Reading.**
  - The strongest comparator closes most of the gap: 3c's +6.2% / +18.6% over Sage+mokapot shrinks to +3.4% / +9.3% at
    matched arm-P FDP 1%. Most of 3c's apparent gain came from the Prosit features, as D17 found for OT-HCD.
  - Unlike OT-HCD (0.990× on D0), ours stays above the strongest comparator on both IT-CID sets, with every lower bound > 1.
- **Caveats.**
  - The CIs resample runs only (one decoy set), so they are narrower than T1's runs × shuffles CIs.
  - Our arm P uses our E2001–3 indexes, while Sage's uses the M1e P2001–3 DBs (same seeds, as in 3c).
  - At own q 1%, the sequence-level arm-P FDP on D4 is 1.29% for ours vs 1.15% for the strongest comparator; this is why
    the matched-FDP ratio is the fairer number.
  - Sage-side `pred_missing` is 9% of PSMs (peptides > 30 aa).
- **Runtime** (shared node, unlocked, indicative only):
  - Feature PINs: D10cf 161–846 s per DB (one outlier at 2,344 s under contention); D4 56–102 s per DB.
  - Seeded mokapot on the pred PIN: D10cf about 220 s / 450 core-s; D4 about 44 s / 92 core-s.
  - Our dev re-featurization from the warm store: 0.5–2 s per run after first contact. While filling, runs were
    Koina-bound (up to about 100 s per run, about 10k predictions/s).

## 3. IT-CID fresh null set DNIT (tail-G2 / G1 / G3)
- **Rules.** Written and committed (2814f06) before any listing or download: A2 verbatim except the IT-CID adaptations
  in `NULL_RULES.md`.
  - Gate DB: human prod. Hosts: D10cf. Search: ±0.5 Da. Identity screen f(c) at ±0.5 Da. No FS nulls (they cannot work
    at unit tolerance).
  - Organism class: free-living, non-host-associated, excluding every ALLOW genus and every DN2-seen or DN2-reachable
    project. No 3c training or null project.
  - Downloads verified by PRIDE SHA-1 and size.
- **Dated deviations, all before the affected step, none touching gate data:**
  1. Tail-proxy threshold taken over the 24 hosts with a defined 1% threshold (33.14).
  2. PRIDE SHA-1 also taken from `checksum.txt`; the walk was re-run once from the start (only metadata seen).
  3. Futility rule added after the first admitted project.
- **Walk (seed 47).**
  - 28 decisions: 15 protocol/SDRF, 7 with < 3 usable RAW, 3 with no PRIDE SHA-1, 3 admitted. All admitted runs were
    100% IT-CID; marker screen passed.
  - Admitted: PXD014351 *Dinoroseobacter* (n_A 461 of 26,431 certified), PXD008499 *Nannochloropsis* (24 of 7,475),
    PXD011268 *Trichococcus* (96 of 2,991).
  - The identity screen at ±0.5 Da removes 86–89% of certified spectra: at unit tolerance, ≥ 20% b/y coincidence with
    some human candidate in the window is common.
  - **|A| = 581, tail proxy = 0.** Projection at the 16-project cap: about 3,100, against the 120,000 required.
- **Outcome: power insufficient. Stopped by the futility rule (note 3). No mixture built, no gate run.**
  - tail-G2, G1 and G3 for `itcid_linear_pred` remain "not evaluable", so IT-CID FDR claims stay blocked.
  - A tail proxy of 0 also means that no A spectrum reached the host's 1% hyperscore threshold: even a much larger
    certified pool would give n_tail ≈ 0 under these rules.
  - DNIT is spent under P4 (`dnit.tsv`, `dnit.sha256`, `dnit_select.log.json`, DATA.md). A further attempt needs a new
    amendment designed for unit-resolution MS2 (for example an identity screen that does not remove almost every certified
    spectrum) and must report this attempt.
  - The P1 swap test was not run (stratum A is not used).

## Open / caveats
- The pinned store pins *one* Koina draw. The 3c dev numbers reproduce exactly from it, but a rebuild differs in 2.8% of
  entries. Production should ship the pinned library, or host Prosit locally (ONNX).
- Item 2's ratios are dev, sequence-level and Sage-DB-based. The D10cf arm-C excess noted in 3c is unchanged.
