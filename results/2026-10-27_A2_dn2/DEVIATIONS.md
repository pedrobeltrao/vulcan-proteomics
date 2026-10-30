# A2 DN2, steps 2–4: dated deviations and ambiguity resolutions

Operator: the data/admission agent (steps 2–4 only; it does not run step 5). Binding text: DESIGN §6.0b =
`reviews/amendment_A2_draft.md` at tag `prereg-v0.4.2` (40d9fb9), with `scripts/a2/IMPLEMENTATION_NOTES.md`.
Rule (A2.9.5): every note below is written **before** the step it affects; where the text is open, the stricter reading wins.
`scripts/a2/`, `vulcan/gates/`, `config/a2*.{yaml,tsv}` are not modified (the selector's `require_tag` guard checks this).

## 2026-10-27 19:57 CET — before step 2 (no DN2 candidate listed yet)

1. **Download location (operational, no effect on selection).** `select_dn2.PrideIO` downloads to `data/raw/dn2` (Lustre).
   The operator instruction and LESSONS 4 require local NVMe. The walk is run through the wrapper
   `scripts/a2_ops/run_dn2_walk.py`, which calls the unmodified `select_dn2.main()` and replaces only the I/O class
   (`PrideIO` subclass): RAWs go to `/tmp/vulcan_a2/raw/dn2/<acc>/`. Ingest then reads from `/tmp` directly (no copy).
2. **PRIDE SHA-1 verification (stricter; LESSONS 6).** The tagged `PrideIO.download` checks size only. The wrapper also
   computes SHA-1 and sha256 in one pass and compares the SHA-1 with the PRIDE archive API v3 `checksum` field; if that
   field is empty, with the project's `checksum.txt` in the PRIDE FTP directory, if present. A file with neither is
   recorded as `no_pride_checksum` (size + sha256 only). On a mismatch the file is deleted and re-fetched, at most 3
   times; if it still mismatches **the walk aborts** (a project is *not* skipped for a download failure, because the
   rule has no such skip). The abort would be recorded here and reported; it is not resolved by the operator.
   `list_files` additionally carries the PRIDE checksum in each file dict; `usable_raws`/`pick_files` ignore the extra key,
   so the file choice is unchanged.
3. **Restarts after an infrastructure failure are resumptions, not re-draws (P4).** The walk is deterministic (fixed
   order, seeded picks, deterministic searches). If it crashes (network, ingest hang, OOM), it is re-run from the start
   with the same code and inputs: verified RAWs are reused, everything else is recomputed. Every walk decision is
   appended to `/tmp/vulcan_a2/walk_events.jsonl` with a run number; after a restart the decisions of the earlier run
   are compared with the new run's, and **any difference stops the procedure** and is recorded here. No candidate is
   ever added, removed or re-ordered by hand.
4. **Steps 2 and 3 interleave (IMPLEMENTATION_NOTES 19).** The stopping rule needs `admit.py`'s label-free counts after
   each project, so the selector runs the own search, marker screen and admission inside the walk. "Each step is
   committed before the next" is therefore satisfied at commit granularity: commit 1 = step 2 artifacts (manifest
   `dn2.tsv`, `dn2.sha256`, selection log, DATA.md); commit 2 = step 3 artifacts (membership hashes, label-free counts,
   swap test, host threshold); commit 3 = step 4 (mixture manifests). Commit 1 precedes commit 2, which precedes any
   mixture build. No gate statistic exists at any of these points.
5. **Host threshold first.** `scripts/a2/host_threshold.py` runs once before the walk (IMPLEMENTATION_NOTES "first
   step-3 run must start with host_threshold.py"). It reads labels of the two D0 *host* runs only and stores two scalars.
6. **P1 swap test.** `admit.py --swap-test <first admitted project>` is run as soon as that project is admitted
   (concurrently with the rest of the walk; it cannot influence the walk). **If it fails, step 4 is not run** and the
   failure is reported (P1 violated = admission invalid).
7. **Disk (< 300 GB on /tmp).** After the marker screen *excludes* a project, its RAWs and stores are deleted from /tmp
   (their sha256/SHA-1 are kept in the log). Admitted RAWs stay on /tmp until step 4 is done. Persistent copies under
   `data/cache/a2/` (gitignored) are: admitted-run stores, `admission/`, `host_threshold.json`, mixtures and their PSMs,
   each with a committed sha256 manifest. RAWs are not copied (re-fetchable; checksums committed in `dn2.sha256`).
8. **End of candidates (ambiguity).** A2.9.3 says that at the cap of 16 the gates run anyway. It does not say what
   happens if candidates run out first with the criteria unmet. Reading: the walk ends there (nothing is extended or
   re-drawn), the shortfall is recorded as the result, and **step 4 is not run without a PI decision** (operator
   instruction "record it and stop"). If the cap is reached with criteria unmet, step 4 *is* run, as A2.9.3 states.
9. **Label discipline in the wrapper.** The wrapper prints and stores only what the selector already logs (project
   decisions, marker counts, admission label-free counts, stopping state) plus download checksums. It never opens the
   gate index labels. The own-search target/decoy split stays inside `own_search.py` (not stored, as in the tagged code).
10. **Run-time exclusion set.** The selector reads the union of `data/manifests/*.tsv` at start. Other agents may add
   manifests during the walk; the set used is fixed at start, and its size and sha256 are recorded in NOTES.md.

## 2026-10-27 20:20 CET — after the walk (run 1), before any step 3/4 action

11. **The walk admitted 0 projects. The reference-proteome lookup in the tagged selector is broken (implementation
    defect, not a data outcome).** `select_dn2.PrideIO.proteome` queries UniProt with
    `(organism_name:"<genus species>") AND (proteome_type:1)`. The current UniProt REST API returns **no rows** for
    `proteome_type:1`, even for a known reference proteome (`(organism_id:243230) AND (proteome_type:1)` returns nothing,
    while `UP000002524` has `proteomeType = "Reference proteome"` and `reference:true` finds it). So every candidate that
    reached step 5 of the walk (17 of 58) was logged "no reference proteome" and skipped. The code was never run
    before the tag (IMPLEMENTATION_NOTES: "Code only, never run"); `scripts/data/select_t4.py` uses the working filter
    `reference:true`.
    - A metadata-only diagnostic (no file listed beyond what the walk listed, nothing downloaded or searched) with
      `reference:true` finds a reference proteome for **all 17** skipped candidates. They cover exactly 10 genera
      (Fervidobacterium ×5, Streptomyces ×2, Caballeronia ×2, Bacillus, Rhodobacter, Chlorella, Deinococcus, Emiliania,
      Synechocystis, Frankia).
    - Second latent defect: for a genus-only frame organism ("Fervidobacterium", PXD054278/PXD054280) the
      `genus + species` fallback is the genus alone, and with the working filter it would match a *different species'*
      proteome (F. ngatamarikiense). The stricter reading would be "no reference proteome" for genus-only organisms.
    - The other 41 decisions (37 protocol-text/SDRF, 4 too few usable RAW) are independent of the defect and valid.
      Most protocol exclusions come from the verbatim A1 keyword list without word boundaries (`standard`, `dia` in
      "media", `acetyl` and `N-termin` from search settings), as fixed by the PI before the tag (§6.0b point 2).
    - **Decision (stricter reading; operator instruction "record it and stop", "never re-draw, extend or prune DN2
      outside the rule"):** I do **not** patch the query or re-run the walk. Re-running would need a change to the
      tagged A2 code (or an override in the wrapper) after the first walk's outcome is known, which is a protocol
      decision for the PI, not the operator. Steps 3 and 4 are not run (nothing admitted). No own-species search,
      admission, swap test or mixture exists; no gate quantity of any kind was computed.
    - **For the PI:** no spectrum was downloaded or searched, so a corrected re-run cannot have been informed by any
      admission or gate outcome; only metadata (protocol text, file lists, proteome existence) has been seen. Even
      with the fix, at most 10 genera remain, i.e. exactly the `min_projects = 10` minimum, so the stopping criteria
      can be met only if all 10 projects pass the marker screen and jointly reach |A| ≥ 120,000, ≥ 72,000 A spectra
      with a candidate and tail proxy ≥ 400. If any one is excluded, the walk ends below 10 projects (DEVIATIONS item 8).

## 2026-10-27 20:40 CET — PI decision D16 (842c961): fix two selector bugs, re-walk once. Written before the fix and run 2

12. **Scope of the fix** (only the two bugs D16 names). In `select_dn2.PrideIO.proteome`: (i) the query filter is
    `reference:true` (as `select_t4.py`); (ii) the match must be at species level. The query string and the matching are
    moved into two pure functions (`proteome_query`, `match_reference_proteome`) so that each can be unit-tested; this
    changes no other behaviour. Keyword lists, order, seeds and every other rule are unchanged.
13. **Ambiguity: what counts as "only a genus".** A frame organism whose second word is the placeholder `sp.`/`sp`/`spp.`
    (e.g. "Synechocystis sp. pcc 6803") does not name a species, and the old `genus + species` fallback key
    "synechocystis sp." would match whichever *Synechocystis sp.* reference proteome UniProt lists first (possibly
    another strain). **Stricter reading: a placeholder epithet counts as genus-only, so the row has no reference
    proteome** (logged "no reference proteome"), as does a one-word organism (D16). A row with a real species epithet
    matches a hit whose organism name starts with the full frame organism, else with "genus species".
14. **Tag guard.** `select_dn2.require_tag` refuses to run when `scripts/a2` differs from the tag, and D16 does not
    move the tag. The operator wrapper therefore replaces the guard (wrapper only; the tagged guard is not edited) with
    a stricter equivalent: every A2 path (`scripts/a2`, `vulcan/gates`, `config/a2*`) must be identical to the tag
    **except** `scripts/a2/select_dn2.py`, which must be byte-identical to the blob in the committed D16 fix commit, and
    the working tree for those paths must be clean. Anything else aborts.
15. **Run 2 bookkeeping.** Run 1's outputs are preserved in git (2650cb8). Run 2 overwrites `data/manifests/dn2*` and
    appends to `walk_events.jsonl` with `run = 2`. This is the single re-walk D16 allows: if run 2 crashes for an
    infrastructure reason, item 3 (deterministic resumption, decisions compared) applies, and it is not a third draw.
16. **Committed code only (2026-10-27 20:55, before run 2).** The shared checkout has another agent's uncommitted 3c
    edits to `vulcan/gpu/{engine,host,_lib}.py` and `kernels.cu` (a Da fragment-tolerance mode). To keep steps 2–4 on
    committed code, run 2 and everything after it run from a clean git worktree `/tmp/vulcan_a2/wt` at commit df2e380.
    The CUDA library is rebuilt there from the committed `kernels.cu`. Gitignored resources (`conda_envs`, `data/fasta`,
    `data/cache`, `tools/*` except `vulcan_kernels`) are symlinked. Check: `host_threshold.py` re-run in the worktree
    is byte-identical to the run-1 value. Outputs written into the worktree (`data/manifests/dn2*`) are copied back to the
    main checkout for commit.
17. **(2026-10-27 21:35, before resuming run 2) Crash on a run with zero certified spectra.** Run 2 stopped (EXIT 1)
    at PXD006818 (Streptomyces lividans TK24, UP000028682). The project had passed the marker screen. A label-free
    diagnostic gives certified counts of 1,067 / 916 / **0** for its three runs. For the empty run,
    `admit.gpu_exhaustive` → `engine.search_exhaustive` has an empty candidate set and `psm_table` fails on `k2 = None`.
    This is an implementation crash (the GPU path was never run before the tag), not a rule question. The meaning is
    unambiguous: a run with no certified spectrum contributes no A spectrum, and all of its spectra get code 0 (B
    candidates). `admit_project` already computes exactly that for n = 0 once the exhaustive call returns empty arrays.
    **Fix (wrapper only; tagged A2 code and `vulcan/gpu` untouched):** the wrapper wraps `admit.gpu_exhaustive` so that
    a run with zero certified spectra returns empty arrays without a search. Every other call is unchanged. The swap
    test runs through the same wrapper (`run_dn2_walk.py --swap-test`). Run 2 is then **resumed** under item 3 / 15
    (same code otherwise, same order; label `A2_WALK_RUN=2r`). Verified RAWs are reused, and the decisions of the
    resumed run for PXD054269, PXD054285 and the 58-candidate prefix must equal those of the crashed attempt, or the
    procedure stops. This is not a third draw.
18. **(2026-10-28 00:55, before resuming run 2 again) Transient HTTP 403.** Run 2r stopped (EXIT 1) at PXD015681
    (Chlorella vulgaris). The first RAW returned HTTP 403 from ftp.pride.ebi.ac.uk, and the same URL returned 200 a
    minute later. This is an infrastructure failure. Wrapper changes (I/O only): each curl call is retried up to 10 times
    with linear backoff; RAWs of marker-excluded projects are now kept on /tmp until the walk ends, so a resume does not
    re-download them (deleted afterwards; item 7 amended). The walk is resumed under items 3 / 15 as run `2r2`.
    Decisions before PXD015681 must reproduce those of 2r: PXD054269, PXD054285, PXD054267 excluded; PXD006818,
    PXD018953, PXD044078 and PXD011591 admitted with identical counts; PXD054272 as in 2r.
    Clock note: the times given in items 12–17 are approximate. The commits are authoritative: items 12–15 in 0265820,
    before the fix d6d199c; item 16 in 5c90498 and item 17 in ad617f5, each before the run it governs.
19. **(2026-10-28 01:30, before resuming as run `2r3`) Unverifiable re-download.** Run 2r2 aborted (item 2 rule) on the
    first file of PXD054269. One curl HTTP error was followed by repeated size mismatches (the wrapper did not log the
    sizes). A manual `curl` of the same URL gave the full size and sha256 `cfbdc952…bc763`, identical to run 2r's
    record. The file was moved into place. Wrapper changes (I/O only): size mismatches are logged; 6 attempts instead of
    4; **and, stricter, a file without a PRIDE checksum must also match the sha256 recorded for it in any earlier run of
    this walk** (`walk_events.jsonl`). Resumed under items 3 / 15 as run `2r3`; the same decision-equality condition as
    item 18 applies.
20. **(2026-10-28 ~05:35, before resuming as run `2r4`) /tmp full (ENOSPC).** Run 2r3 stopped while ingesting
    PXD018045 (Frankia alni). The node-local NVMe was 100% full: other agents hold about 820 GB, and this walk about
    36 GB. All decisions of 2r3 up to that point equal those of 2r (items 18–19). It also decided PXD062500
    (Deinococcus) admitted and PXD038243 (Emiliania) admitted, both new. Operational changes (I/O only; no rule
    affected): every DN2 RAW except PXD018045's is moved to the persistent `data/cache/a2/raw/dn2/<acc>/` on Lustre and
    symlinked back to `/tmp/vulcan_a2/raw/dn2/<acc>`. On resume, the wrapper re-hashes the file through the symlink and
    ingest finds the existing NVMe store (`conv` returns "exists"), so no RAW is decoded from Lustre (LESSONS 4). The
    partial store `Q20060_Frankia-11.vspec.arrow.partial` is deleted. This also covers A3.1.5 (RAWs of admitted runs
    retained persistently).
21. **(2026-10-28 ~08:00) End of the walk; step 4 not run.** Run 2 (attempts 2, 2r, 2r2, 2r3, 2r4) exhausted all 58
    candidates with **6 projects admitted**. The stopping criteria are unmet (projects 6 < 10, |A| 24,581 < 120,000,
    A with candidate 9,608 < 72,000, tail proxy 22 < 400). Every decision and every file sha256 is identical in all
    attempts that reached it (checked from `walk_events.jsonl`). Per D16 point 4 and item 8, this is the result:
    **mixtures are not built (step 4 not run)**, and per D22 / A3 steps 4b and 5 are not run by this agent.
    Housekeeping: during the ENOSPC incident (item 20) the RAWs of PXD054269 and PXD054285 (excluded) were deleted
    from /tmp before item 20 was written. They were re-downloaded in 2r4 and matched their recorded sha256. RAWs and
    stores of all excluded projects are now deleted (checksums in `selection_log.tsv`). Admitted RAWs, stores
    (with MS1) and admission files are persistent in `data/cache/a2/` (A3.1.5).
