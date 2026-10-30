# Proposed pruned LESSONS.md (12 rules)

1. Check hardware (nproc, affinity, GPU, RAM, network) before planning, and size all thread pools from `os.sched_getaffinity(0)` via `vulcan.hw`, never a constant or `os.cpu_count()`. [HPC review; D5 node change; afe8d0e]
2. Timing needs an exclusive resource: short lock holds only, long comparators unlocked with load noted, and every benchmark records load average and hwinfo. [M6 lost its bench to a 2 h xiSEARCH hold; load 6-55 in locked runs]
3. Measure before optimizing: I/O path and per-run cost in a persistent warm process, not a one-shot. [E1 30x; M2b 17.8 vs 5-8 core-s/1e5 MS2]
4. Copy RAW to local NVMe before decoding (no mmap on Lustre), supervise every third-party reader with timeout, no-CPU-progress check and retry, and keep torch out of CPU workers. [M2b: 5 s -> 77 s; 4 hangs/6k runs; 65 -> 14 core-s]
5. GIL-bound host work belongs in worker processes, not Python threads; long-lived torch processes on a shared GPU must free cache and retry before deleting inputs. [M2b round 2; E9 OOM and ~200 re-downloads]
6. Verify downloads by repository checksum (PRIDE SHA-1), never size alone, and write a manifest-vs-disk verifier before re-fetching. [E9 7% corrupt; data_verify]
7. Preflight every new dataset/index: FASTA concatenated newline-safe with all expected accessions present, MS2 analyzer/activation checked before tolerance and scorer choice, and charge-set/missingness separated in quant comparisons. [M5 SNCA loss; D10cf ion trap; M6d]
8. Never leave durable state in node-local /tmp (envs, caches, builds) unless a script rebuilds it; commit WIP code and results to git at least daily. [node move lost mokapot env, mzML, caches; unversioned A2 gate code]
9. Have a second, independent implementation or reviewer check WIP and every gate/protocol before a sealed look; write admission rules for test/null sets before fetching them. [2c 650 -> 29 core-s; D9 review; D12 contaminated nulls]
10. Anything chosen by looking at a truth set is "tuned" there: declare tuning set and claim set in advance, confirm on a fresh set, and record gate failures as final rather than explaining them away post hoc. [D7 variant selection; M6d/D11; D12]
11. Re-check MBR/entrapment decoy exchangeability on each new instrument or gradient. [M6d: 2.4x under-count on D12 Lumos vs exact on D9 Astral]
12. Ask sub-agents to write `metrics.json` and `NOTES.md` (not summary/report files), and treat headline numbers as unverified until they cite a results path, nominal-vs-verified FDR and comparator load. [stats and E1 agents lost summary.md; status summary vs D12]

Dropped as stale or merged: 4-thread pinning and 4-core lock rules (D5), duplicate RAW/Lustre, reader-supervision, checksum and FASTA rules, `pgrep -f` self-match (minor, shell hygiene), nvcc smoke test (env-specific, kept in ENVIRONMENT.md).
