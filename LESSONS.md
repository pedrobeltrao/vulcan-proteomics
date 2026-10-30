# LESSONS — process rules (pruned at milestone M2, 2026-10-28; retrospective reviews/retro_milestone2.md)

1. Check hardware first and size every pool from `os.sched_getaffinity(0)` via `vulcan.hw`; every benchmark records hwinfo and load average, and speed counts only within one node and session. [D5; afe8d0e; M1 retro; M6 bench lost to 2 h lock hold]
2. Timing needs an exclusive slot with short holds; long comparators run unlocked with load noted; a speedup is quoted only with comparator threads and load. [M6; load 6-55 in locked runs; D31.4 still uses the shared lock]
3. Measure before optimizing (I/O path, warm persistent process), and require bit-identity before speed; new fast paths merge opt-in after the freeze. [E1 30x; M2b 65->14 core-s; ed98d7c 7.2x exact; D27 E11]
4. RAW to local NVMe before decoding (no mmap on Lustre); supervise third-party readers with timeout/no-progress/retry; torch stays out of CPU workers; GIL-bound work goes in processes. [M2b 5 s -> 77 s; 4 hangs/6k runs; E9 OOM]
5. Verify downloads by repository SHA-1, never size; write a manifest-vs-disk verifier before re-fetching. [E9 7% corrupt; data_verify; 8fb2330]
6. Preflight each dataset/index in five minutes: FASTA newline-safe with expected accessions present, MS2 analyzer/activation checked, decoy exchangeability re-checked per instrument. [M5 SNCA loss; D10cf ion trap; M6d 2.4x MBR under-count]
7. Durable state never lives in node-local /tmp: use a persistent cache plus a free-space guard that kills the job below a floor; budget disk for decode running ahead of the GPU. [node move c49a0e7; A2 ENOSPC 2c734fa; E7b /tmp full 7865f8e]
8. Commit code first, then results from a clean tree at the commit named in git_sha; commit WIP at least daily. [7fc78f7, 0165732, 7bad342; 5c90498 clean-worktree run]
9. Wait on background jobs by PID (`kill -0`), never `pgrep -f <pattern>`; sub-agents write `metrics.json` and `NOTES.md`, not summary files. [M7 ~30 min lost; summary.md lost twice]
10. Pin every remote predictor in a content-addressed store and reuse it for all replicates and swaps. [3c G0 TMT failure; b574398; D20]
11. Write the label-free yield/power simulation and its go/no-go threshold before drafting any gate or null-set amendment or fetching data; stop on "not evaluable". [D24 DN2; D20 DNIT; D23; D28 A4 not run vs D30 C9 1.0]
12. Independent review before every adoption and before any sealed look; reuse one regression-tested stats library across amendments instead of re-deriving it. [A3-A6 five "not as written" verdicts; D9; A6 B1-B3; D28 review]
13. Anything chosen by looking at a truth set is tuned there: declare tuning and claim sets in advance, write admission/selection rules before fetching, run bug fixes from a clean worktree with a dated note, and record gate failures as final. [D7; D11; D12; D16; 90b64d8]
14. Compare to the strongest, fairly configured comparator before any headline, label every ID gain "nominal FDR, not gate-validated" until gates pass, and refresh the status block in the same commit as the decision. [D17 0.990x; D21 1.22 -> D25 1.15; stale 10-27 status]
15. Cap every side track with a kill criterion and time box in its first commit, and keep RESEARCH_TREE decisions to one headline line with detail linked and superseded entries marked `[-]`. [D19 E5/E6 NO-GO; D31 re-cut; D0/D2 stale; 60 KB tree]

---
- Time the full sealed look (all comparators × replicate DBs) on a dev set of the same size before the session, and start it with ≥ 2× that time left; the comparator stage, not our engine, dominates. [D36/D37: T1 reached comparators only after 3 h; 6 runs × 25 DBs did not fit the last 8 h of the node]
- Record full manifests (package lists, not just digests) of every env a sealed look checks, in git, and keep the envs on persistent storage; a digest alone cannot be reproduced. [D38: node wipe lost /tmp mokapot env; rebuild matched core versions but not the digest]
