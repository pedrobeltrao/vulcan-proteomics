# Second independent review of Amendment A4 (draft revision 2, commit 138f102)

- **Reviewer:** independent statistician (agent), second reviewer (D26 item 4). I had no role in A1–A4, the first A4 review or the T1-sequencing opinion, the DN2/DNIT walks, E7/E7b, M7, the T1 driver or any gate run.
- **Date:** 2026-10-28.
- **Read:** `reviews/amendment_A4_draft.md` (rev 2); `reviews/checks/a4_r1_walk.{py,json,out}`, `a4_power.py` (functions `dn2_yields`, `dn2_tier0_projects`, `cap_factors`, S0 block of `a4_power.json`), `a4_pool_count.py`; `reviews/amendment_A4_review.md`; `reviews/a4_t1_sequencing_review.md`; RESEARCH_TREE D9–D27; DESIGN §6 (§6.0–§6.0c, §6.6); `reviews/amendment_A2_draft.md` (A2.3, A2.9.6), `reviews/amendment_A3_draft.md` (A3.1, A3.3, A3.6), `scripts/a2/IMPLEMENTATION_NOTES.md`, `scripts/a3/IMPLEMENTATION_NOTES.md` (freeze notes 33–37); `vulcan/gates/a2_stats.py` (`tail_g2_a2`, `g1_gate_a2`), `vulcan/gates/a3_config.py` (`freeze_status`); `bench/sealed_run.py`; `bench/t1_driver.py` (CLI, sealed-mode branch); `bench/t1/a3.py` (header); `config/a3.yaml` (freeze section, all TODO); `results/SEALED_LEDGER.md` (header only, no entries).
- **Label-free DN2 records read:** `results/2026-10-27_A2_dn2/{NOTES.md, admission_counts.json, walk_events.jsonl ("listed" events only), selection_log.tsv}`; `data/manifests/{dn2,dnit,g2null}_select.log.json` (decision reasons); `data/manifests/{g2null,dnit}.tsv` (accession column). Local PRIDE crawl cache `data/cache/pride_meta/` (crawl table, project records, cached file listings: names and sizes only).
- **Not read:** anything under `data/sealed/` or `data/escrow/`; any sealed T1–T5 file; any DN2/DNIT admission array (`data/cache/a2/admission/*.npz` not opened; I took DN2 yields from the committed S0 block of `a4_power.json`); `/tmp/vulcan_a2`; anything under `results/2026-10-28_E7b_fair_speed/` (my copy of the exclusion scan skips that directory). No PRIDE call. No per-label quantity of any null set.
- **Checks I ran** (scratch `…/scratchpad/reviewer2/`; nothing in the repo was run via `main()`, written or overwritten; `git status` of `reviews/checks/` unchanged):
  - `samelab.py`: A4.3.3 item 5 (same-lab exclusion) evaluated on the tier-0 projects and on the pool, with A4.3.4 lab normalisation with and without the `labPIs`/`submitters` person strings.
  - `walk2.py`: re-run of the `a4_r1_walk` walk (its pool, protocol flags and yield model imported unchanged) adding, cumulatively: the |A| cap factor; A4.3.4 person links; item 5 on the pool; item 5 on tier 0; no tier 0. Tier-0 order bug fixed. 300 walks per cell, seed 131+70 as the original.
  - `walk3.py`: the same variants with per-project RAW counts drawn from the **empirical** distribution of DN2's own 21 listed projects (`walk_events.jsonl` `n_raw_records`, usable ~ Binomial(n, 0.8)) instead of the assumed log-normal medians 6/12/24.
  - `peptone.py`: A2 contamination terms dropped by the revision-1 screen, on rev-1 passers and on tier 0.
  - Distribution of RAW counts in the 56 cached PRIDE file listings and DN2's listing events.

---

## Verdict: **do not adopt revision 2 as written. Adopt after a revision 3 that fixes B1–B2 and M1–M6; no change of architecture is needed.**

Revision 1 fixed the first review's three blocking items correctly in substance (word boundaries; a realistic walk; capped/uncapped intersection–union). The statistics A4 adds (outer envelopes, α split, IUT) are sound and stricter. The problems are in the two revision-2 additions:

1. **Tier 0 contradicts A4's own same-lab rule.** Read as written, A4.3.3 item 5 excludes 3 (arguably 4) of the 5 eligible tier-0 projects, and the reason it does so — a facility that produced five marker-contaminated DN2 candidates — is exactly the hazard item 5 exists for. The walk check implements neither item 5 nor the A4.3.4 person links, so the stated evaluability (and the "+0.02 to +0.06" from tier 0) is optimistic.
2. **The escrow's unseal trigger (a) opens T1 when DN3 is "not evaluable"**, which is broader than D26 decision 3 ("until A4 step 5") and reopens the leakage route the sequencing opinion called the strongest: the next null-set amendment designed by people who have seen T1. Under my re-derivation, "not evaluable" is the more likely outcome in realistic cells.

The escrow mechanics also leave real read paths open (driver work directory, output default, manifest sizes), and A4's own code placement collides with the A3 freeze.

---

## BLOCKING

### B1. Tier 0 and A4.3.3 item 5 (same-lab exclusion) are mutually inconsistent, and the strict reading removes most of tier 0 (A4.3.3 items 3 and 5; A4.3.10)
- **Text.** Item 3 carves the six DN2 projects out "of items 1–3 and of the file-level guard". It does **not** carve them out of item 5. Item 5 excludes "every project whose normalised lab (A4.3.4) equals that of … any DN2 candidate that is **not** in tier 0". A4.3.4's normalised lab includes every `labPIs` and `submitters` name.
- **Facts (crawl table + cached records; label-free):**

  | Tier-0 project | Normalised lab / persons | Non-tier-0 DN2 candidate(s) with the same string | Item 5 as written |
  |---|---|---|---|
  | PXD018953 *Caballeronia* | jean armengaud | PXD054269, PXD054285, PXD054267, PXD054272 (all **marker-screen excluded**, 2–17 non-keratin H/B members), PXD018045 (**marker excluded**), PXD054278, PXD054280, PXD052144, PXD062210 | **excluded** |
  | PXD044078 *Bacillus* | jean armengaud | same | **excluded** |
  | PXD006818 *Streptomyces* | tassos economou; konstantinos tsolis | PXD006819 (genus slot taken) | **excluded** |
  | PXD062500 *Deinococcus* | "yongqian zhang"; siyu zhu | PXD027969 *Deinococcus radiodurans*, lab "**zhang yongqian**" (same PI, reversed name order) | passes only because the normalisation is order-sensitive |
  | PXD038243 *Emiliania* | rosalind rickaby; craig dedman | — | passes |
  | PXD011591 | silke leimkühler; beata siemiatkowska (person hit) | — | already out (protocol) |

  A4.3.10 instead says the Armengaud pair "form one cluster" and are expected admitted, and `a4_r1_walk.py` (`pool_r1`, lines 74–101; tier-0 loop lines 211–232) applies item 5 to nobody.
- **Why this matters statistically, not just textually.** Item 5's stated purpose is the "shared-instrument carry-over route", applied to *all* same-lab projects "not only the contaminated ones". The Armengaud facility (the same MSxx-yyy run naming across PXD018953, PXD044078 and PXD054269) produced five DN2 candidates excluded for human/bovine serum markers (ALB, C3, IGHG1, LTF, TF, SERPINA1 …). Carving the two Armengaud projects out of item 5 would admit exactly the configuration item 5 was written to exclude. A4.1.6 ("stricter reading wins") resolves the conflict towards exclusion, which makes A4.3.10's table, its "expected contribution", A4.0 and A4.7.2 wrong.
- **Pool effect not modelled.** Item 5 also removes **32** pool projects (17 of them rev-1 protocol passers, mostly tier I) that `a4_r1_walk` keeps (`samelab.py`).
- **Fix.**
  1. Say explicitly that item 5 applies to tier 0 (my recommendation, and the reading A4.1.6 forces). Then tier 0 = {PXD062500, PXD038243}, and with order-insensitive name normalisation (below) = {PXD038243} only. Record this before the tag, as the draft does for PXD011591. If the PI prefers to keep the Armengaud/Economou projects, that is an explicit carve-out of item 5, which has to be argued against the five marker exclusions and decided in a dated D-note before the tag, not left to the implementer.
  2. Make person-name normalisation order-insensitive (sorted alphanumeric tokens, single-letter middle initials dropped: "zhang yongqian" = "yongqian zhang"; "brook l nunn" = "brook nunn"). Apply it to item 5 and to cluster links.
  3. Re-run `a4_r1_walk` with item 5 on the pool and on tier 0, the person links of A4.3.4, the |A| cap factor (M4 below), and replace every tier-0 figure (A4.0, A4.3.10, A4.7.2, A4.11.1, A4.13).

### B2. Unseal trigger (a) opens T1 on a "not evaluable" DN3, broader than D26 and against the sequencing safeguards (A4.8.4(a))
- **Text.** A4.8.4(a) unseals on "the commit of DN3 step-5 results … **or the commit recording DN3 as not evaluable under A4.3.8**".
- **D26 decision 3** says outputs stay unread "until A4 step 5 or until the project nears its time limit". A not-evaluable DN3 has no step 5. The draft widens the PI's decision.
- **Sequencing opinion, safeguards 1 and 3:** "If the v1 minimum is not met, T1 stays blocked until a further null-set amendment is tagged and its walk reaches the same point", and "the author and reviewer of any later null-set amendment" must be T1-blind. Rev 2 adopts neither (A4.8.3 only covers "A4 agents"). A4.12 says the opinion is "made normative in A4.8". That is not true of these two items.
- **Size of the problem.** Under my re-derivation (table R1 below), P(v1 minimum not met) is about 0.45–0.80 in the cells that DN2's own file counts support. So the most likely path is: DN3 not evaluable → T1 unsealed → A5 (the null-set amendment A4.6 already budgets at α/2) is designed by T1-exposed people. Arm C (*Arabidopsis* natural sequences) on T1 is direct information on whether the scorers prefer natural sequences, which is what the tier-II (land-plant) rows test.
- **Fix.** Delete "or the commit recording DN3 as not evaluable" from (a). Then add:
  - "If DN3 is not evaluable, the escrow stays sealed until the step-5 commit of the next null-set attempt for v1, or trigger (b)."
  - "The author, reviewer and adopter of any later null-set amendment, and its walk and step-4b/5 agents, must be T1-blind, with exposure statements."
  
  If the PI wants unsealing on not-evaluable, that is a new decision that knowingly accepts the T1→A5 route. It must be recorded as such, and A5 must then carry a T1-exposure disclosure.

---

## MAJOR

### M1. The escrow does not contain the T1 driver's outputs (A4.8.2 items 1–2)
- **`T1_OUT` does not exist.** `bench/t1_driver.py` takes `--out` (default `results/<date>_T1_sealed/T1`) and `--work` (default `/tmp/vulcan_t1_<uid>/T1`, lines 134, 187–189). The work directory holds the replicate indexes, every method's PSM tables, mokapot outputs, `predcache`/`predcache_a3` and the chimeric outputs. Sealed mode deletes only "staged copies / stores / mzML" (docstring; line 436). The other files, which are T1 results, stay in a shared `/tmp` that every same-uid agent can read.
- **Manifest sizes leak.** "File names and sizes" (item 2) of ID lists, PSM tables or per-stratum tables are monotone in ID counts. The committed manifest would be a coarse T1 readout for everyone.
- **`chmod 000` is not access control** for agents running as the same uid; it only stops accidental reads. During the run (hours), the 0700 escrow directory is fully readable, and it sits under the repository, where repo-wide `grep -r`/`find` by other agents reaches it.
- **Fix (normative text):**
  - escrow mode passes `--out` **and** `--work` inside one escrow root;
  - that root is on node-local disk outside the repository;
  - at exit the whole root (out + work + caches + logs) is tarred and then deleted;
  - the tar is **encrypted to a public key whose private key is not on the cluster** (held by the PI). `escrow.py open` decrypts only after the trigger commit;
  - the manifest lists file names only (or the tar's total size and sha256 only);
  - a test asserts that, after an escrow-mode D0/D1 dry run, no file newer than the run start exists outside the escrow root, except the ledger and manifest, under `/tmp/vulcan_t1_*`, `results/` and the cache dirs;
  - the console check is a **whitelist of line templates** (stage names, job counters), not a blacklist regex. Sage, Comet and mokapot print counts (for example "N PSMs at q ≤ 0.01").

### M2. "Time limit" trigger (b) is discretionary, and A4 completion is conditional ("if time allows") (A4.8.4(b))
- As written, the PI can declare the final-report phase the day after the escrow run. That turns escrow into an ordinary early look and restores the sequencing opinion's persistence risk: a T1-exposed PI deciding whether to continue a long or failing DN3. "A4 runs to completion if time allows" is an exit that the sequencing opinion's safeguard 4 explicitly closed ("stopping or abandoning A4 for any reason other than its own label-free rules is a recorded protocol violation").
- **Fix.**
  - (b) may not fire before DN3 steps 2–3 are committed (the earliest point the sequencing opinion allowed a look).
  - In the adoption note, the PI fixes a calendar date, or a rule computed from a committed quantity, for "time limit". The sequencing opinion's "pre-commit a deadline and resources for A4" recommendation becomes normative.
  - Replace "if time allows" with: "A4 runs to completion; stopping it other than by its own rules is a protocol violation reported with T1."

### M3. A4's own code placement collides with the A3 freeze; the "where possible / A3.6 deviation route" fallback is not valid (A4.4.3, A4.4.5, A4.8.1.5)
- `a3_config.freeze_status` requires HEAD's tree hashes of `vulcan/`, `bench/`, `scripts/a3/`, `scripts/e7/`, `scripts/e7b/` to equal the freeze commit's (lines 53–89; A3 notes 33). The sealed T1 driver calls `require_freeze_commit()`.
- A4 puts `vulcan/gates/a4_stats.py` under `vulcan/` (A4.4.3), and the escrow changes under `bench/` (A4.8.1.5). Anything landing after the freeze makes HEAD's hashes differ. Then the escrow run (which needs the A4 tag, hence the A4 implementation) and A4 step 4b/5 for v2_sep refuse to run. The same holds for D27's GPU MS1 merge "after the freeze".
- The A3.6 deviation route is for engineering fixes to `gates_a3.py` after 4b. It is not a route for adding new T1 outputs (stratum-1 family-S inputs) or new modules to frozen trees.
- The freeze is not yet taken (`config/a3.yaml` freeze: TODO), so this is cheap to fix now.
- **Fix.** Make the order normative, with no "where possible":
  1. A4 rev 3 adopted (at least A4.8 and A4.5 frozen);
  2. A4 code in `vulcan/`, `bench/` and the A3-frozen script trees committed and tested, including the escrow mode, `escrow.py`, `claims_a4.py` and the stratum-1 outputs of `bench/t1/a3.py` with the A3.1.6 bit-identity dry run re-done;
  3. the A3 freeze commit;
  4. `prereg-v0.4.3`;
  5. `prereg-v0.4.4`.
  
  Alternatively, keep all A4 code outside the frozen paths (`scripts/a4/`). Then say so, and drop `vulcan/gates/a4_stats.py`. D27's merge must wait until after the T1 escrow run and DN3 step 5, or be done on a branch.

### M4. The evaluability numbers are optimistic, and the input that drives them is partly knowable now (A4.0, A4.7.2, A4.11.1)
Corrections to the walk model (`a4_r1_walk.py`):
- (i) **No item 5** on the pool or tier 0 (B1).
- (ii) **Cluster links use the crawl `lab` only** (`walk()` keys `(g, l)`), not the A4.3.4 `labPIs`/`submitters` strings.
- (iii) **The cap uses min(proxy, A_cand)** (`capped()`, line 135), not the three factors of A4.3.6.
- (iv) **Tier 0 is never shuffled.** `sorted(tier0, key=lambda d: random.Random(order_seed).random())` gives every element the same key, so the order is always that of `admission_counts.json` (I verified: identical output for seeds 1000, 1001, 5000). The effect is small.
- (v) **Tier-0 marker screen.** It passes stage 1 with probability 0.9 and stage 2 with the generic e. But DN2's NOTES show 3 of the 5 eligible tier-0 projects already had one non-keratin panel protein at 3 runs (PXD006818 and PXD038243: bovine CSN1S1; PXD062500: human ALB). One more pooled protein triggers exclusion, and stage 2 pools up to 48 runs.

The draft calls the per-project usable-RAW count unknowable. DN2's own 21 listing events (label-free, committed) give n_raw_records quantiles 4 / 6 / **12** / 16 / 36 (p10/p25/p50/p75/p90). Similarly, 3 of DN2's 6 passers carrying one marker at 3 runs argues for e ≥ 0.35 rather than 0.15. Re-derived (300 walks per cell; table R1):

| Model | P(v1 min), e = 0.15 | P(v1 min), e = 0.35 | P(full), e = 0.15 / 0.35 |
|---|---|---|---|
| rev 2 as simulated, empirical DN2 RAW counts | 0.75 | 0.37 | 0.50 / 0.11 |
| + \|A\| factor + A4.3.4 person links | 0.72 | 0.34 | 0.41 / 0.08 |
| + item 5 on pool (tier 0 carved out) | 0.58 | 0.23 | 0.29 / 0.05 |
| + item 5 on tier 0 (literal; PXD062500 kept) | 0.53 | 0.22 | 0.28 / 0.06 |
| full spec, no tier 0 | 0.54 | 0.17 | 0.22 / 0.05 |

- **Fix.**
  - Re-run with (i)–(v).
  - Add the empirical-RAW cell, and state it as the **reference** cell in A4.0 and A4.11.1 ("about 0.2–0.55 for v1, 0.05–0.3 for v2_sep tier 1"), with the 6/12/24 grid kept as sensitivity.
  - Say plainly that tier 0 contributes ≈ 0 under the strict reading.
  - The PI should decide with these figures. This is the "evaluability honesty" item: it does not change any gate, but it changes the expected value of spending 3–4 days and most of the pool.

### M5. The response to the first review's M2/B3 is weaker than A4.12 presents it; stale size figures (A4.3.6–A4.3.7, A4.6, A4.10.6, A4.12)
- **J_eff ≥ 10 is (almost) implied by the 10% cap.** If every capped share is ≤ s, then Σw² ≤ s and J_eff ≥ 1/s = 10. Shares can exceed 10% slightly only because the factor is the minimum across quantities. In every v1-minimum state of my 300-walk runs, the capped J_eff was ≥ 10.0 (minimum 10.01), while the maximum capped share reached 0.145. The criterion never bound. A4.12 lists only "J_eff ≥ 15" as not adopted and presents J_eff ≥ 10 as the M2 fix. In fact the first review's M2 request is not adopted, and the protection is the cap plus the IUT.
- **Concentrated-leak size, rev-2 run (`a4_r1_walk.json` tail_r1).** 0.040 (κ = 1, both leak patterns) against nominal 0.025; MC SE ≈ 0.011, so a 95% interval of about 0.02–0.06. The text still says "≤ 0.037" (A4.6), "0.017–0.037" (A4.10.6) and "0.22–0.33 capped-only" (A4.3.6; the JSON says 0.27–0.29). A4.9's "≤ 0.04" is the only current figure.
- **Fix.** Either:
  - (a) state that J_eff ≥ 10 is redundant with the cap; that the residual is a concentrated-leak size of ≈ 0.04 (CI 0.02–0.06) vs 0.025 nominal at the v1-minimum states; and that the PI accepts it in the adoption note; or
  - (b) adopt J_eff ≥ 12 on capped quantities (binding, unlike 10) and simulate its cost.
  
  Update every size figure to one source.

### M6. A4 amends A3/D12 by "reading as" and leaves the stratum-1 family-S estimator unspecified (A4.8 preamble; A4.8.1.4)
- "A3.3.1's 'no T1 look if v1 fails or is not evaluable' is read as 'no T1 claim'" re-interprets normative A3 text. The sequencing opinion asked for explicit replacement text for D12, A3.3.1, the first row of A3.3.7, and A3.3.2/A3.3.8. A3.3.1's rationale (not spending T1 on a scorer without validated FDR) is reversed: after the unseal, T1 is spent for every future scorer (v3) whatever DN3 shows. The draft never says so.
- The stratum-1-only family-S inputs (A3.3.5) are to be "added to `bench/t1/a3.py` if absent". Their estimator is not defined: matched arm-P FDP on a stratum-1-only TDC list, or counts restricted to stratum 1 with acceptance from the full list as in A3.3.6, plus the LOO minimum. That is an open analysis option in normative text (A3 Q6).
- **Fix.**
  - Replace by quotation: D12's "the T1 look stays blocked" → "T1 may be computed under escrow (§6.0d A4.8); FDR claims stay blocked until v1 passes". A3.3.1 first bullet and A3.3.7 row 1 → "no T1 claim; T1 reported in full under A4.8.5".
  - Add: "Unsealing spends T1 as a confirmatory set for every scorer."
  - Define the stratum-1 list estimator in A4.8.1.4 (my suggestion: the stratum-1 list = v2_sep accepted units of merged stratum 1 at common q_sep, with matched arm-P FDP computed on that list's own paired counts, LOO as A3.3.4). It must be reviewed T1-blind before the escrow run.

---

## MINOR

1. **IUT precedence** (A4.3.6). "Not evaluable on either pool means not evaluable. A fail on either means fail" conflicts when one pool fails and the other is NE or invalid. `tail_iut` (lines 314–319) gives *fail* precedence even over an *invalid* pool. A2 (`tail_g2_a2` lines 136–141) gives invalid precedence over fail within a pool, and a fail is final. **Fix:** any invalid → invalid; else any fail → fail; else any NE → NE; else pass. Apply this to every IUT row, and re-run tail_r1.
2. **Clean invalid-pool rate is not reported.** Under the IUT, a clean scorer at the v1-minimum states is "invalid pool" 0.07 (κ = 1) and 0.10 (κ = 0.5) of the time (`a4_r1_walk.json`). That is an extra not-evaluable route. Report it in A4.7.3 and A4.10.
3. **A2 contamination term dropped.** `peptone from (meat|casein)` (A2.3 step 2) is absent from the rev-1 contamination group (A4.3.5(a); `R1_CONT`, line 109), contradicting "nothing relaxed". It hits 0 current rev-1 passers, so restoring it costs nothing. Add `tryptone` as well: two tier-0 projects show bovine casein at 3 runs.
4. **Final n_min vs fixed stopping thresholds.** A_cand ≥ 120,000 and 450,000 (A4.3.7) are fixed, but final n_min "can only be equal or larger" (A4.5.3). If G1's final n_top exceeds ≈ 12,000, "v1 minimum met" no longer implies G1 is evaluable. Define the thresholds as functions of the final n_min (label-free, before the tag).
5. **File-guard carve-out too wide** (A4.3.3 item 3). Carve out only a tier-0 project's own DN2 files (matches in `dn2.sha256`). Other sha256 matches (for example T-set, A1 or training files) stay unusable.
6. **DN2 npz bit-identity check** (A4.3.10). Compare the recomputed npz sha256 with `results/2026-10-27_A2_dn2/cache_admission_stores.sha256`. Never open the DN2 arrays, because the stores hold T∪D candidate identities.
7. **A1-class check** (A4.5.5). It is gating but sits under "reported, never used to pass". It does not say which pool. Move it, and apply it on both pools (IUT).
8. **Escrow retries** (A4.8.2.4). "Same commit, code" contradicts "may only follow an engineering fix". A fix changes the code, the driver hash in `config/a4.yaml` and possibly frozen trees (M3). Say what is re-hashed, and who reviews (T1-blind). Name the final successful run as the result of record. A failure summary of "stage + exception type" is too thin to debug. Allow the traceback's file:line:function frames without messages or locals, and require debugging by reproduction on D0/D1.
9. **`sealed_run.py` conflicts.** It requires a `milestone-*` tag on HEAD (rule 2) and skips only ledger-only commits. The manifest commit (A4.8.2.2) would break a retry's tag check. Also, the ledger "result" line publishes the T1 wall time (a speed number). Specify the tag rule for escrow, and move the wall time into the escrow.
10. **Speed results under load.** The escrow run "overlaps the walk" (A4.4.8, A4.7.4). T1 speed (a §2 claim) would then be measured under contention from DN3 admission compute. Run the escrow on an otherwise idle node, or label T1 speed "not claims-grade".
11. **Stale or wrong cross-references.** A4.13 says "resolutions in A4.14" (there is no A4.14). "D26 item 4" is the second review, and "thresholds unchanged" is part of D26 item 1. A4.10.5 says "about 8%" while A4.7.3 gives a clean fail rate of 0.05. A4.3.10's "39–67 runs" should be 41–67 (JSON medians). The `a4_r1_walk.py` docstring still describes rev 1 (≤ 3 per cluster, J_eff ≥ 15).
12. **Claims script inputs.** `claims_a4.py` also needs the A3.1.8 dev-precondition verdicts (D23) as inputs. Name them.

---

## Did revision 1/2 resolve the first review and the sequencing opinion?

| Item | Status |
|---|---|
| B1 (`IgE`/digest) | Resolved. Leading `\b`, case-sensitive `\bIgE\b`, regression tests. One A2 term was dropped (minor 3) |
| B2 (optimistic walk) | Largely resolved in method. Inputs still omit item 5, person links and the \|A\| factor, and the RAW distribution is not anchored (M4) |
| B3 (concentrated leaks, cap re-weighting, 1.5α) | Resolved by the IUT (size-correlated leak 0.29 → 0.00) and by withdrawing 1.5α. Residual 0.040 vs 0.025 under concentrated leaks, with stale figures (M5) |
| M1 (envelope) | Resolved (A4.5.2) |
| M2 (cap on three quantities, J_eff, no silent no-cap) | Text resolved. Simulation uses two quantities. J_eff ≥ 10 is non-binding (M5) |
| M3 (multiplicity) | Resolved |
| M4 (organism class) | Resolved |
| M5 (exclusions) | Resolved for accessions. Item 5 added, but it conflicts with tier 0 (B1) |
| M6 (T1 safeguards normative) | Partly. Escrow is normative, but safeguards 1 (blocked after a not-evaluable DN3), 3 (later null-set authors T1-blind) and 4 (no abandonment) are missing or weakened (B2, M2) |
| M7 (open options) | Resolved, except the new open items: the stratum-1 estimator (M6) and the code-placement fallback (M3) |
| Minor 1–12 | Resolved, except stale figures (minor 11 here) |

---

## What should not change
- **The architecture.** Label-free walk, admission, cap and stopping; the membership hashes of both pools committed before any gate statistic; separate walk, 4b and step-5 agents; one step-5 session, v1 first.
- **Margins unchanged.** The α split (α/2 for every re-tested row, shared by v1, v2_sep and ce_s); the outer envelope over the analytic, project and cluster intervals; the capped/uncapped intersection–union for every A- and B-based row; the 10% cap on the minimum of the factors. Validity controls are not hulled.
- **Completed test = step 5 ran.** A2 counts zero. The next attempt is at α/4 after a completed test, α/2 after a label-free stop.
- **The word-bounded protocol screen with its field split**, and the PXD011591 outcome recorded before the tag.
- **Reusing DN2's admitted projects is statistically sound in principle.** Only label-free DN2 information exists, it is walked first whatever its yield, and it spends no α. What must change is the rule conflict (B1), not the principle.
- **Escrow as the T1 mechanism.** One look enforced in code, every branch computed in one run, a pre-hashed claim script, reporting regardless of outcome, labels until DN3 step 5, no external release before unseal.
- **Two-stage marker screen** pooled over all selected runs; the file-level SHA-1 and (name, size) guard; permanent exclusion of every reached candidate; the byte-identical host-threshold check; the P1 swap tests.

## Numbers I re-derived
- **Pool:** 576 exclusion accessions / 351 rows removed; tiers 171 / 258 / 131 (identical). A4.3.3 item 5 removes a further **32** pool projects (17 rev-1 passers; 25 / 5 / 2 by tier) with person strings, or 23 with the crawl lab only.
- **Tier 0 under item 5:** PXD006818, PXD018953 and PXD044078 hit (crawl lab alone suffices); PXD011591 hit via person string (already protocol-excluded); PXD062500 hit under order-insensitive names (PXD027969, "zhang yongqian"); PXD038243 clear.
- **Walk, draft model reproduced** (lab-only, no item 5; 300 walks): P(v1 min) 0.25 / 0.07 / 0.87 / 0.45 / 0.99 / 0.94 for (6, 0.15) / (6, 0.35) / (12, 0.15) / (12, 0.35) / (24, 0.15) / (24, 0.35), against the draft's 0.24 / 0.08 / 0.87 / 0.49 / 1.00 / 0.91.
- **Walk, full spec, same grid (6/0.15, 6/0.35, 12/0.15, 12/0.35, 24/0.15, 24/0.35):**

  | Variant | P(v1 min) | P(full stop) |
  |---|---|---|
  | Tier 0 carved out | 0.11 / 0.04 / 0.71 / 0.29 / 0.99 / 0.84 | 0.00 / 0.00 / 0.44 / 0.10 / 0.95 / 0.62 |
  | Item 5 on tier 0 too | 0.09 / 0.03 / 0.61 / 0.25 / 0.98 / 0.78 | 0.01 / 0.01 / 0.40 / 0.08 / 0.93 / 0.62 |
  | No tier 0 | 0.08 / 0.02 / 0.67 / 0.21 / 0.97 / 0.77 | 0.01 / 0.00 / 0.40 / 0.08 / 0.93 / 0.61 |

  The draft states 0.08–0.24 (6 runs per project), 0.49–0.87 (12) and 0.91–1.00 (24). Under the full spec these are 0.02–0.11, 0.21–0.71 and 0.77–0.99. MC SE is ≤ 0.03.
- **R1 (empirical DN2 RAW counts):** see the table in M4. Scratch outputs are `walk2.json` (grid) and `walk3.json` (empirical).
- **Empirical RAW counts:** DN2 listed projects (n = 21) have n_raw_records p10/p25/p50/p75/p90 of 4 / 6 / 12 / 16 / 36. The 56 cached prior-use listings have usable RAWs (0.3–1.5 GB) of 0 / 0 / 12 / 39 / 100.
- **Capped J_eff at v1-minimum states:** minimum ≥ 10.0 in every cell (non-binding); maximum capped proxy share 0.12–0.145.
- **Tail-G2 (from `a4_r1_walk.json`, not re-run):** concentrated-leak IUT pass 0.040 / 0.040 (κ = 1) and 0.023 / 0.030 (κ = 0.5); capped-only under the largest-cluster leak 0.29 / 0.27; clean pass 0.93 / 0.85, clean invalid 0.07 / 0.10, clean fail 0.00 / 0.05.
- **Tier-0 order in `a4_r1_walk`:** never shuffled (sort key constant).

## What the PI needs to decide (before revision 3)
1. **Tier 0 under item 5 (B1).** Exclude the Armengaud/Economou (and Zhang) projects, as recommended, or carve them out explicitly with a reason.
2. **Unsealing after a not-evaluable DN3 (B2).** Keep sealed, as recommended and per D26's wording, or unseal and accept a T1-exposed A5.
3. **Time limit (M2).** A calendar date or a committed rule for trigger (b).
4. **Evaluability (M4).** With realistic inputs, P(v1 minimum) is roughly 0.2–0.55 and P(full stop) 0.05–0.3. Proceed knowing that a not-evaluable DN3 consumes most of the remaining pool.
