> **Notice.** This fact-check was produced by an AI agent team (Claude Opus 5.5 / Sonnet 5.5). It has not been reviewed or verified by human experts.

# Fact-check of the Vulcan reports (2026-10-30)

Scope: REPORT.md, reports/briefing/BRIEFING.md, reports/technical/TECHNICAL_REPORT.md, reports/manuscript/MANUSCRIPT.md, and reports/figures (CAPTIONS.md, make_figures.py). Every quantitative claim was checked against the primary sources: results/*/metrics.json and NOTES.md, E12 at `7f64a72`, RESEARCH_TREE.md decisions, results/SEALED_LEDGER.md, DESIGN.md, `bench/sealed_run.py` and config/a5.yaml. All numbers remain **development data**. No held-out test-set result exists.

## Fixes applied in all documents

| Issue | Correct statement now used | Primary source |
|---|---|---|
| (a) "40–65× fewer core-s than Sage" for closed search | The search stage uses about 28–48× fewer CPU core-s (CPU only; GPU time excluded): 4.9 core-s vs Sage 237 at 32 threads (48×) and 137 at 4 threads (28×). Open search is quoted separately (about 65–66×). RAW→proteins is 9.1 core-s vs Sage 70 (4 thr) to 120 (32 thr) from mzML, plus 250 core-s for RAW→mzML | E2 NOTES speed table; M4 NOTES; E12 NOTES headline |
| (a′) E2 erratum | E2's 4.9 core-s was derived. The M3 erratum gives a measured 8.59 s / 6.0 core-s, i.e. 40× vs Sage at 32 threads and about 23× vs Sage at 4 threads. This is noted next to the 28–48× figure in REPORT §2.1, the technical report §6.1.2 and the manuscript Results | results/2026-10-27_M3_scoring/NOTES.md l.210 |
| (b) Retries | A5 allows at most two retries. Only attempt 1 was used and it failed, so **up to two retries remain** | `bench/sealed_run.py` MAX_RETRIES = 2; D29.4; D36; D38 |
| (c) A4 power | P(v1 minimum) 0.08–0.43 is kept distinct from P(clean v1 passes) about 0.06–0.4 and P(leak caught) about 0.07–0.4. The 0.06–0.43 "power" range is gone | D28; D30 |
| (d) Look start | 14:23 CET (ledger 2026-10-28 13:23:46Z; CET = UTC+1, since DST ended 2026-10-25). D36's 14:25 is the timestamp of the decision entry | results/SEALED_LEDGER.md |
| GPU time | Wherever core-s are compared with CPU-only tools, the text now says core-s exclude GPU time | method convention |
| Dev vs held-out | D8h is "10 dev runs held out from scorer selection; not a sealed test set". Figure labels now read "D8h dev holdout" | DATA.md; M7/E7b NOTES |

## REPORT.md

| Section | Old → new | Source |
|---|---|---|
| top | (none) → AI-team disclaimer added as the first line | task requirement |
| header | "RESEARCH_TREE D0-D36" → "D0–D39". "Node allocation ended about 22:00" → planned ~22:00 (D31), actually wiped ~00:30 on 10-29 (D38). Clock times are CET; the ledger is in UTC | RESEARCH_TREE D31, D38 |
| Conventions | "Verified FDP … (library design, spike-in, entrapment)" → entrapment gives an estimate (FDP̂), not a verified FDP. Core-s exclude GPU time | DESIGN §6 terms |
| §0 bullets | "One A5 retry … is still allowed" → "at most two retries; one attempt failed, so up to two remain"; source cited as D38/D39 | sealed_run.py; D29.4, D38 |
| §0 table, DDA | "search 40–65×" → "search stage about 28–48× … at 4–32 threads (CPU only; excludes GPU time)"; Sage 70–120 at 4–32 threads | E2; E12 |
| §0 table, chimeric | "10 held-out dev runs" → "D8h (10 dev runs held out from scorer selection; not a sealed test set)" | E7b |
| §0 table, IT-CID | "about 1 s per fraction vs Sage about 54 s per file" (mixed units) → "about 0.9 core-s per fraction vs Sage about 54 core-s per file (16 threads, unlocked, indicative)" | 3c_itcid NOTES l.125–130 (643 core-s / 12 files) |
| §0 table, open | "about 65× fewer core-s" → "search stage about 65× fewer core-s than Sage-open (32 threads)" | M4 NOTES l.113 |
| §0 table, XL | "Cross-linking (DSSO) … 460–8,900× faster" → accuracy is from DSSO (D3b); the speed is from DSS (D3): about 460× whole-pipeline wall, about 8,900× search-stage wall | M5_xl NOTES l.97 |
| §0 table, L6 | "about 1% global FDP" → 1.16% (P) / 1.21% (C) peptide and 1.00% protein at N = 1,000 (982 runs completed); union 9.3% / 19.1% | E9 NOTES l.116–120 |
| §1 A4 | "simulated power 0.06-0.43" → "P(v1 minimum) 0.08–0.43 (clean v1 passes about 0.06–0.4)" | D28 |
| §1, §3, §4 | look time "14:23" → "14:23 CET (ledger 13:23:46Z)" | ledger |
| §2.1 v1 | "D1 0.98-1.004" (the range over all comparators) → "D1 0.998 [0.997, 1.000]" vs the strongest comparator | T1driver_dryrun/D1/tables.md l.53 |
| §2.1 v1 | "**40-65x fewer core-s**" → "about 48x (32 thr) and 28x (4 thr) fewer CPU core-s (CPU only; excludes GPU time)" + M3 erratum | E2; M3 erratum |
| §2.1 IT-CID, open | Sage IT-CID marked as 16 threads; open-search ratio marked "GPU time excluded" | 3c NOTES |
| §2.1b | directLFQ "22–26 → 5 s" → "22–28 → about 5 s" | E12 NOTES headline |
| §2.4 | PXD009449 inconclusive at 1% → also at 5% (4.44%, UCB 6.93%), "183 single-Y peptides" | T2_driver_dryrun NOTES l.39 |
| §2.7 | TRFP 531 core-s/1e5 MS2 → caveat added: E1 cold-Lustre inflated; 235 on the node of record | hw_adapt NOTES l.28 |
| §3 | Retry wording → "none used, up to two remain"; "about 11 h of look time" → D39's total of about 13 h | D39 |
| §6 | Duplicate item numbers 3 and 4 renumbered (Appendix E #5); "Retry the T1+T5 look" added as the first follow-up | RESEARCH_TREE status summary |
| lint | A5 wording lint: 0 violations before and after | bench/t1/wording_lint.py |

## BRIEFING.md

| Section | Old → new | Source |
|---|---|---|
| Read-this-first box | Look time added (14:23 CET). "One retry … is still allowed" → "A5 allows at most two retries …; none has been used, so up to two retries remain" | ledger; sealed_run.py |
| §3 table | DDA, chimeric, IT-CID, XL and L6 rows corrected as in REPORT §0 (28–48×; D8h wording; core-s units; DSS/DSSO split; 1.16–1.21% / 1.00%) | as above |
| §4 A4 | "simulated power was only 0.06–0.43" → "chance of reaching its minimum sample only 0.08–0.43 (clean scorer would pass with probability about 0.06–0.4)" | D28 |
| §5 | "held-out dev subset" → "D8h, a dev subset held out from selection" | M7 NOTES |

## Figures (make_figures.py, CAPTIONS.md; all six figures regenerated)

| Item | Old → new | Source |
|---|---|---|
| Fig 2a x-axis | "CPU core-seconds (user + sys)" → "… Vulcan GPU time not included" | convention |
| Fig 3 / Fig 4a labels | "D8h holdout", "D8 holdout" → "D8h dev holdout" | M7, E7b |
| Fig 6 A4 label | "not run (power 0.06–0.43)" → "not run (P(v1 min) 0.08–0.43)" | D28 |
| Fig 6 look label | "opened 14:23" → "opened 14:23 CET" | ledger |
| CAPTIONS general | Core-s now noted as excluding GPU time (comparators are CPU-only) | — |
| CAPTIONS Fig 2a | 1.26 × 10⁵× XL ratio: now says it sets Vulcan's search stage (0.43 core-s) against the full xiSEARCH run (53,671 core-s). Source wording added (about 50,000× incl. sites; about 460× whole-pipeline wall). xiSEARCH D3 was run under the bench lock | M5_xl NOTES l.90, l.97; comparator_runtime.json |
| CAPTIONS Fig 3/4 | D8h described as a dev holdout, not a sealed test set | — |
| CAPTIONS Fig 6 | A4 "power 0.06–0.43" → P(v1 minimum) 0.08–0.43 / clean pass 0.06–0.4. Look 14:23 CET (13:23:46Z). "One retry" → up to two retries remain | D28; ledger; sealed_run.py |

## Manuscript, follow-up by the editor
- The Results XL sentence said the xiSEARCH timing was "taken without the bench lock". It now says the D3 xiSEARCH run was timed under the bench lock, and the DSSO figure is a subset extrapolation on a contended node (M5_xl NOTES l.90, l.105). Methods and the Fig. 2 legend were adjusted to match.
- The "⚠" symbol was removed from the disclaimer for font portability. The wording is unchanged.

## TECHNICAL_REPORT.md

| Location | Old → new | Source |
|---|---|---|
| TECHNICAL_REPORT Exec summary table, standard DDA | "search 40–65× fewer CPU core-s than Sage" → "search stage about 28–48× fewer CPU core-s than Sage (4 and 32 threads; CPU only, excludes GPU time)"; Sage "70–120" → "70 (4 thr)–120 (32 thr)" | results/2026-10-27_E2_gpu_search/NOTES.md (4.9 vs 237/137 core-s); E12 NOTES @7f64a72 |
| TECHNICAL_REPORT Exec summary table, chimeric | "10-run D8h holdout" → "10-run D8h dev holdout (held out from selection, not a sealed test set)" | DATA.md / §5.1 roles |
| TECHNICAL_REPORT Exec summary table, IT-CID | "about 1 s per fraction vs Sage about 54 s per file (unlocked)" → "about 0.9 core-s per fraction vs Sage about 54 core-s per file (16 threads, unlocked, indicative)" | results/2026-10-27_3c_itcid/NOTES.md l.125–129 (0.91 core-s; 643 core-s / 12 files, 16 thr) |
| TECHNICAL_REPORT Exec summary table, open search | "about 65× fewer core-s" → "about 65–67× fewer core-s (search stage; 17.2 vs 1,153)" | results/2026-10-27_M4_ptm/NOTES.md l.108–113 |
| TECHNICAL_REPORT Exec summary table, XL | "Cross-linking (DSSO) … 460–8,900× faster" → accuracy on DSSO D3b; speed on DSS D3: ~8,900× less search-stage wall, ~460× whole-pipeline wall, ~50,000× fewer core-s (search + sites) | results/2026-10-27_M5_xl/NOTES.md l.97 |
| TECHNICAL_REPORT Exec summary, note under table | added: core-s excludes GPU time; comparators CPU-only on same 32 cores | E2/E12 NOTES (core-s = user+sys) |
| TECHNICAL_REPORT Exec summary "What is missing" | "14:23" → "14:23 CET (ledger 13:23:46Z)"; "One A5 retry … remains allowed" → "A5 allows at most two retries and only one attempt has been used, so up to two retries … remain" | results/SEALED_LEDGER.md; bench/sealed_run.py MAX_RETRIES=2; D29.4, D38 |
| TECHNICAL_REPORT Exec summary "Speed in one sentence" | "about 9 core-s" → "about 9 CPU core-s (GPU time not counted)", "Sage's search" → "Sage's CPU-only search" | E12 NOTES |
| TECHNICAL_REPORT How-to-read / Glossary | "Core-s = user + sys CPU seconds" → adds "excludes GPU time (reported separately as GPU-s)" | E2 NOTES |
| TECHNICAL_REPORT §3.6 scorer table | v1 "trained on 13 runs from 13 projects" → "12 runs from 12 projects (13 selected; one IT-CID run excluded automatically)" | results/2026-10-27_M3_scoring/NOTES.md l.27–28 |
| TECHNICAL_REPORT §6.1.3 | "12 usable runs from 13 PRIDE projects" → "13 runs from 13 projects selected, one IT-CID run excluded, leaving 12 runs / 12 projects" | M3_scoring NOTES l.27–28 |
| TECHNICAL_REPORT §4.4 | cross-reference "Section 6.13" → "Section 6.12" (decoy realism) | internal |
| TECHNICAL_REPORT §4.10 Fig. 6 caption | "A4/DN3 not run (simulated power 0.06–0.43 …)" → "simulated P(v1 minimum) 0.08–0.43; P(clean v1 passes) about 0.06–0.4"; "14:23" → "14:23 CET"; "one retry … still allowed" → "up to two retries … remain allowed under A5" | RESEARCH_TREE D28/D30; ledger; sealed_run.py |
| TECHNICAL_REPORT §6.1.2 | "40–65× … spans these two readings" paragraph → "about 48× and 28× fewer CPU core-s (237/4.9, 137/4.9) for the search stage; CPU only"; M3 erratum kept (measured 6.0 core-s = 40× vs Sage-32, ~23× vs Sage-4); notes REPORT's 40–65× mixed in the open-search ratio and is corrected | E2 NOTES l.61–72; reviews/m3_independent_review.md l.161; M3_scoring NOTES l.210 |
| TECHNICAL_REPORT §6.1.4 / App. A D17 | "no loss at much higher speed" → qualified: "non-inferior identifications at matched entrapment-estimated FDP; dev only / untested on held-out data" | D17; T1driver_dryrun NOTES |
| TECHNICAL_REPORT §6.1.5, §5.3, Figs 3/4 captions, §6.4, App. D | "D8h holdout"/"D8 holdout" → "D8h dev holdout" (held out from selection; not a sealed test set) | DATA.md roles; §5.1 |
| TECHNICAL_REPORT §6.1.6 | Sage "70–120 core-s/run" → "70 (4 threads)–120 (32 threads) core-s/run … (CPU only; Vulcan core-s exclude GPU time)" | E12 NOTES headline table |
| TECHNICAL_REPORT §6.5 open search | "about 65× fewer core-s" → "about 65–67× fewer CPU core-s (1,153 / 17.2; GPU time excluded)" | M4_ptm NOTES l.108–113 |
| TECHNICAL_REPORT §6.10 Fig. 2 caption | added "Vulcan's GPU time is not included (comparators CPU-only, ratios favour the GPU node)" | CAPTIONS.md (revised) |
| TECHNICAL_REPORT §6.10.1 ingest | "TRFP → mzML 531" → adds "(cold-Lustre inflated; 235 on the node of record)" | results/2026-10-27_hw_adapt/NOTES.md l.28 |
| TECHNICAL_REPORT §6.10.3 | "directLFQ 22–26 s → 5 s" → "22–28 s → 4.7–5.0 s per experiment" | E12 NOTES @7f64a72 "Biggest wins" |
| TECHNICAL_REPORT §6.11 1 TB table | "about 15 core-s" → "about 15 core-s (GPU time not counted)" | E12 NOTES |
| TECHNICAL_REPORT §7.2 | "the node allocation ended about 22:00" → "was planned to end about 22:00 (D31/D36); in the event the node was wiped about 00:30 on 10-29 (D38)" | RESEARCH_TREE D31, D36, D38 |
| TECHNICAL_REPORT §7.3 table header | "Time (local)" → "Time (CET, UTC+1)" | ledger 13:23:46Z = 14:23 CET |
| TECHNICAL_REPORT §7.4 / §11 | retry text → "MAX_RETRIES = 2 … only attempt 1 used, so up to two retries remain"; follow-up 1 "(up to two retries remain)" | bench/sealed_run.py l.68; D29.4; D38 |
| TECHNICAL_REPORT App. E items 1, 3, 4, 5, 10 | rewritten as "corrected in this revision"/"resolved" with resolution stated (item 3 kept as a primary-source inconsistency: D39/status summary "one retry" vs code+D38 "≤2"); item 19 updated with commits 7c44c12/f466a54 (R1 D7 arm-M end-to-end test, contained in prereg-v0.4.4); new item 20 (T2 dry-run PXD000138 Ala 5% band [3.38–3.62] vs [3.43–3.71] in the same NOTES) | git log; T2_driver_dryrun NOTES l.38 vs l.78 |

Unresolved (TECHNICAL_REPORT):
- E2 closed-search ratio: the primary sources disagree (E2 derived 4.9 core-s → 28–48×; M3 erratum measured 6.0 core-s → ~23–40×). Directive wording (28–48×) used; the erratum is stated next to it in §6.1.2 and Appendix E item 10.
- Retry count: RESEARCH_TREE D39 and status summary still say "one retry … remains" (not edited; outside scope). Reports follow code + D38 (up to two).
- Appendix E items 2, 6–9, 11–18, 20 remain primary-source inconsistencies (not fixable by rewording).
- Figure 4a values (e.g. v1 D0 arm P 1.007%, comparator 1.045%) were taken as given from SOURCES.md and not re-derived digit by digit.

## MANUSCRIPT.md

| Location | Old → new | Source |
|---|---|---|
| MANUSCRIPT Abstract | "about 12-fold fewer than the Sage search alone" → "about 7–12-fold fewer than a Sage search from mzML alone (Sage at 4–32 threads, CPU only)" | E12 NOTES @7f64a72 headline: Sage 70 (4 thr)–120 (32 thr) core-s/run; 9.9 core-s default |
| MANUSCRIPT Abstract | "within 1% of the strongest comparator" → "about 1% below (0.990-fold)" | T1driver_dryrun D0 0.990 [0.989, 0.992] |
| MANUSCRIPT Results (E2) | "48-fold and 28-fold fewer core-s" → "about 48-fold and 28-fold fewer CPU core-s for the search stage (CPU only; GPU time excluded)" | E2 NOTES: 4.9 vs 237 (32 thr) / 137 (4 thr) |
| MANUSCRIPT Results (workflows) | "against 120 core-s for a Sage search from mzML (13-fold)" → "70–120 core-s ... at 4–32 threads (8–13-fold)"; 9.1 flagged as opt-in | E12 headline |
| MANUSCRIPT Results (workflows) | "Open search ... 66-fold cheaper" → "about 66-fold fewer core-s than Sage's open mode at 32 threads (17.5 vs 1,153 core-s, two runs)" | M4 metrics speed_ours.d0_open.warm_total_core_s 17.54; Sage 1,153.4 (NOTES rounds to ~65×) |
| MANUSCRIPT Results (workflows) | IT-CID "59-fold" → adds "0.91 vs about 54 core-s; Sage 16 threads, unlocked, indicative" | 3c_itcid NOTES l.125–130 |
| MANUSCRIPT Results (workflows) | XL "1.26e5-fold cheaper in search stage" → adds source wording: ~5e4× core-s for search+sites, ~8,900× search-stage wall, ~460× whole-pipeline wall; lock status later corrected by the editor to "under the bench lock" (see above) | M5_xl NOTES l.97; comparator_runtime.json 53,671 core-s vs 0.427 |
| MANUSCRIPT Results (ingest) | "against 531 for ThermoRawFileParser" → adds caveat: cold-Lustre inflated; 235 per 1e5 MS2 on node of record | hw_adapt NOTES l.28 |
| MANUSCRIPT Results (cohort) | "all 288 output tables bitwise identical in each comparison" → "all compared output tables (288 or 480 per comparison)" | E12 NOTES (closed 288/288; LFQ, chimeric 480/480) |
| MANUSCRIPT Results (sensitivity) | "ten runs held out from selection (D8h) ... 1.55 times ... per run" → "ten development runs held out from scorer selection (D8h; not a sealed test set) ... at each method's own 1% q (geo-mean over 9 evaluable runs)" | M3e NOTES D8 confirmation table |
| MANUSCRIPT Results (L6) | "360,392 peptidoforms were accepted" → "about 358,000–360,000 (358,164 arm-P table; 360,392 master-list example)" | E9 NOTES l.25 vs l.116 (source inconsistency, App. E #8) |
| MANUSCRIPT Results (FLR) | "On PXD009449 the upper bound was 1.34% at 1%" → "inconclusive at both levels (UCB 1.34% at 1%, 6.93% at 5%)" | T2_driver_dryrun NOTES l.39 |
| MANUSCRIPT Results (CF) | "no off-peak values"; "0.088 vs 0.082" → defined off-peak (≥3 fractions from MS2 evidence); "on the common protein set; no CI" | M5_lip_cf NOTES E10 table |
| MANUSCRIPT Results (A1 G1) | "ρ = 0.048 (99% CI 0.019–0.076), outside the ±0.05 band" → "its 99% CI was not contained in the ±0.05 band" (point estimate is inside) | M3e NOTES l.28 |
| MANUSCRIPT Results (A4) | "simulated power of 0.06–0.43" → "P(reaching v1 minimum) 0.08–0.43 (P(clean v1 passes) about 0.06–0.4)" | RESEARCH_TREE D28 |
| MANUSCRIPT Results (look) | "opened on 28 October 2026" → "... at 14:23 CET (13:23 UTC in the ledger)" | SEALED_LEDGER 13:23:46Z |
| MANUSCRIPT Results (look) | "One retry from the same frozen commit is allowed under A5" → "A5 allows at most two retries ...; one attempt used, so up to two retries remain; needs mokapot env rebuild + dev check" | bench/sealed_run.py MAX_RETRIES=2; D29.4; D38 |
| MANUSCRIPT Discussion | "about 10 core-s ... one to two orders of magnitude below" → "(GPU time excluded), roughly 30–40-fold below a CPU-only search engine plus file conversion" | (70–120+250)/9.9 = 32–37 |
| MANUSCRIPT Discussion | "Nineteen additional label-blind features gave no held-out gain" → "no gain on held-out training folds" (avoid dev/held-out confusion) | E7_chimeric_v2 NOTES l.46 (leave-one-fold-out) |
| MANUSCRIPT Methods (scoring) | "tryptic human library of 7.06 M peptidoforms at charges 2 and 3" → "human prod index (7.06 M) + P. furiosus index, charges 2–4, all six NCE grid values (24.9 GB, on GPU)" | scoring_speed NOTES l.72–75; M2b NOTES l.4 |
| MANUSCRIPT Methods (scoring) | "chosen per run from 20 to 35" → "from a grid of six values between 20 and 35" | vulcan/rescore/features.py NCE_GRID |
| MANUSCRIPT Fig 3 legend | "the D8h holdout" → "the D8h development holdout" | dev/held-out wording |
| MANUSCRIPT Fig 5 legend | "PXD009449 at 1% is inconclusive" → "inconclusive at 1% and 5%" | T2_driver_dryrun NOTES l.39 |
| MANUSCRIPT Fig 6 legend | "A4: not run (power 0.06–0.43)" → "(simulated P(v1 minimum) 0.08–0.43)"; "opened 28 October 2026" → "+ 14:23 CET"; "One retry ... remains allowed" → "Up to two retries ... remain allowed under A5" | D28; ledger; sealed_run.py |

Verified unchanged (match source): E2 recall 99.90/99.78 (D0), 99.99/99.91 (D1); 19,601/19,605 and 5,758; 175,429,652 postings, 1.5 GB, 0.53 s; D1 v1 ratio 0.998 [0.997, 1.000]; v2_sep 1.149/1.137/1.150 and protein 0.999; IT-CID 1.034/1.093; ce_s 0.983 [0.968, 0.999]; 12 runs / 829k PSMs; 894,579 nulls, 45 mixtures, 81,390 clean; tail 751, 0.715, UCB 0.742, diff UCB 0.124; 9.42% / 5.34% Sage rank≥2; GW stratum ≥3 3.09% (2.38–3.79); v1 arm-P 1.01/0.87/1.02% with CIs; 1.16%/1.00%/9.3%/19.1%/arm C 1.21%; FLR 0.29/0.86/3.43/4.44%; pyAscore +2.0 (1.1–2.8), −1.6 (−3.2 to −0.3); 313/45,845; 1.7%; 13,674 matched, 0.118 vs 0.1145, E. coli 0.160 vs 0.121; MBR 0.63–1.23% / 0.18–0.79%, 2.4×; r = 0.9999996; XL 334 @1.2%, 459 @2.0%, D3 119 @9.2%; LiP 26,976, 20/20, 37–99; 24 decode + 3 post workers; 0.11–0.20 GPU-s; E12 4,570→5,753, 9.9→9.1, 561→1,280, 572→1,194, 2,567→3,867, 2,362→3,172; folds 17 (266/15.3), 158 (1,543.5/9.8), 492 (446/0.906); A6 0.63/0.8, 0.923–0.944; 1.04%→1.24% simulation (DESIGN §6.3, E9 NOTES). Disclaimer is first; no author names/affiliations/placeholders; no "validated FDR control"/"FDR-validated"/"no loss".

Unresolved:
- 360,392 vs 358,164 accepted at N=1000 are both in E9 NOTES (source inconsistency, App. E #8); manuscript now gives both.
- Open-search ratio: M4 NOTES says "~65×", M4 metrics give 1,153.4/17.54 = 65.8× (figure label 66); manuscript keeps "about 66-fold".
- XL 1.26e5× compares our search-stage core-s with a complete xiSEARCH run (not stage-matched); kept because Fig. 2a uses it, now with the source's own ~5e4× (search+sites) wording alongside.
- Methods says the first ingest benchmark ran on a 4-core node "with an RTX 3090 GPU": E1 summary confirms 4 pinned cores; GPU model only indirectly supported (hw_adapt NOTES mentions replaced "3090" assumptions).

Addendum:

| Location | Old → new | Source |
|---|---|---|
| MANUSCRIPT Results (E2) | added: "These E2 figures were derived (process total minus warm passes); a direct single-process measurement gave 8.59 s / 6.0 core-s (excl. index build and RAW ingest), 40-fold and 23-fold fewer core-s" | M3_scoring NOTES l.210 erratum; 237/6.0 = 40, 137/6.0 = 23 |

## Remaining unresolved items (all documents)
1. **Closed-search CPU ratio.** The primary sources disagree. E2 NOTES give 4.9 core-s, which is derived and yields 28–48×. The M3 erratum gives a measured 6.0 core-s, which yields 40× and about 23×. The reports use 28–48× as directed and state the erratum next to it.
2. **Retry count in RESEARCH_TREE.md.** The D39 entry and the status summary still say "one retry … remains". The code (MAX_RETRIES = 2), D29.4 and D38 say at most two. The decision log was not edited; all reports now say "up to two retries remain".
3. **L6 accepted count.** E9 NOTES give 358,164 in the arm-P table and 360,392 in the master-list example (Appendix E #8).
4. **E7b rank ≥ 2 FDP̂ range.** D25 and E7b NOTES say "4–13%", but the E7b table has 3.7% on D1 (not evaluable). The reports quote 4–13% as the source does (Appendix E #12).
5. **Open-search ratio rounding.** The M4 NOTES say about 65×; the metrics give 65.8× (Fig. 2a label 66×). Both readings appear and are consistent with each other.
6. **XL 1.26 × 10⁵×.** This ratio sets the search stage against a complete xiSEARCH run, so it is not stage-matched. The captions and the manuscript now say so and give the source's ~5 × 10⁴× (search + sites).
7. **Other primary-source disagreements.** Appendix E items 2, 6–9, 11, 13–18 and 20 of the technical report are disagreements inside the primary sources and cannot be fixed by rewording. They remain listed there.
8. **Fig. 4a values.** These were checked through SOURCES.md key paths, not re-derived digit by digit.

## Checks
- **Disclaimer:** each of the four documents now opens with the disclaimer. It says the work was produced by an AI agent team (Claude Opus 5.5 / Sonnet 5.5) and has not been reviewed or verified by human experts. REPORT.md had no disclaimer before; one was added.
- **Names:** no author names, affiliations or placeholders appear.
- **Wording:**
  - None of the documents contains "validated FDR control" or "FDR-validated".
  - "No loss" appears only where it is defined as non-inferiority at matched entrapment-estimated FDP̂ on dev data.
- **Held-out comparison:** every document states that it is missing. D8h is described as a dev holdout, not a sealed test set.
- **A5 wording lint on REPORT.md:** `[]`, i.e. 0 violations.
- **PDFs:**
  - The PDFs are rendered with pandoc 3.13 and tectonic in a scratch env (/tmp/render). The reproducible script is `reports/render/render_pdfs.py`.
  - Layout: A4, 2 cm margins, Libertinus fonts, and column widths set from cell content.
  - Every page footer carries the AI-generated / dev-only notice.
  - The manuscript PDF places each figure above its legend.
  - Sample pages were checked visually: there are no missing glyphs and no tables run past the margins.

## PI resolution (2026-10-30)
- **Closed-search CPU ratio.** The measured single-process figure (M3 erratum: 8.59 s / 6.0 core-s) is now primary in all documents: about **23–40× fewer CPU core-s** than Sage at 4–32 threads. The derived E2 figure (28–48×) is mentioned only as an earlier estimate.
- **Retry count.** RESEARCH_TREE D40 corrects D39: up to two retries remain.
