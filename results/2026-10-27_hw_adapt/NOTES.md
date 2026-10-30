# Hardware adaptation to the A100 / 32-core Zen 4 node (D5), 2026-10-27

Node: gpu-node, 32 affinity cores (os.cpu_count()=64), 2x EPYC 9354 (AVX-512), A100 80GB PCIe (drv 580.178.04), 1.5 TB RAM. Numbers: metrics.json.

## Code
- `vulcan/hw.py` (n_threads, pin_thread_pools, thread_env, fingerprint), `bench/hwinfo.py`, `tests/hw/test_hw.py` (5 tests).
- Replaced hard-coded 4-thread / 3090 / cores-112-115 assumptions: `scripts/baselines/{timed.py,common.sh,run_searches.py,write_env.py,make_summary.py}`,
  `scripts/run_stats_sims.py`, `tests/stats/conftest.py`, `scripts/ingest/{run_e1.sh,e1_bench.py}`, build scripts (CARGO_BUILD_JOBS), `bench/README.md`.
  4-thread runs: `VULCAN_THREADS=4` (or `timed.py --threads 4`). `timed.py` records `hw` and `threads_env` per line; lock convention unchanged
  (`$VULCAN_BENCH_LOCK` can relocate the lock file). Left as history: the E1/M1 result files and `e1_summarize.py` hardware string (old node).
- Not changed: `/tmp` staging defaults for M1e baselines and E1 stores (rebuildable scratch; scripts recreate them). The mokapot env still defaults to
  `/tmp/vulcan_m1e/envs/mokapot` (`setup_mokapot_env.sh` rebuilds it); whoever runs baselines must run it first.
- Persistent: cargo targets `data/cache/cargo_target/<crate>`, conda pkgs `data/cache/conda_pkgs`, `TORCH_EXTENSIONS_DIR=data/cache/torch_extensions`.

## Toolchain (all verified)
- Found: the pe env's nvcc could not compile (no cuda_runtime.h, no cicc). Installed CUDA 12.9 dev packages (listed in env/setup.sh and ENVIRONMENT.md).
- nvcc kernel with sm_80 + sm_86 + compute_86 PTX: correct on the A100, cuobjdump shows both cubins and PTX.
- torch 2.14.1: fp32 matmul 18.6 TFLOPS (TF32 off), bf16 260 TFLOPS, bf16 rel err 0.25% vs fp32; `load_inline` CUDA extension builds and runs
  (first build 220 s on cold Lustre; arch list from activate.sh).
- cargo/rustc 1.99.1; Sage smoke test identifies LQSRPAAPPAPGPGQLTLR (2.1 s); Comet 2026.01, TRFP 2.0.0.0, DIA-NN 2.7.0 `--threads 8` all start.
- Sage znver4 (`tools/bin/sage-znver4`, same v0.14.7 source, `-C target-cpu=znver4`): identical PSMs and hyperscores, no speed gain
  (Velos005137, 15,374 MS2: wall 5.2 vs 5.3 s at 32 thr, 4.19 vs 4.16 s at 4 thr; 3 reps each, locked). Stock binary stays the default.
  Side finding: Sage at 32 threads on this small file burns 56 core-s vs 8.6 at 4 threads for essentially the same wall (mzML parse is serial);
  comparator thread counts must be chosen per workload and core-s reported.
- Direct Thermo reader (1 thread/process, warm page cache, same input in all concurrent processes):
  b1906 (3.4 GB): 19.2 core-s/1e5 MS2 at 1 file (E1 old node: 19.1 same file, 17.8 headline), 25.0 at 8 concurrent (higher sys time from shared mmap of one file; 2152 files/h/node
  extrapolated). Velos005137 (small): 25-27 (fixed startup ~2 core-s dominates), 28.8 at 8. So Zen 4 is on par per core with the 7742 for ingest; no per-core gain.
  TRFP->mzML on Velos: 36 core-s (235 per 1e5 MS2) vs 531 in E1 (E1 was cold-Lustre inflated; not a like-for-like).
- `scripts/m2/build_vulcan_core.sh`: builds as-is with the persistent target dir (the first build was done by another agent invoking the script; a from-scratch
  build in a fresh target dir took 26.7 s on 32 jobs); module imports from tools/vulcan_core. No crate source touched.
  `build_vulcan_tims.sh` got the same plumbing change but was not rebuilt.

## Caveats
- Fingerprint does not import torch (slow); torch fields appear only if torch is already loaded.
- Ingest 8-concurrent figure used one file read 8 times; distinct-file concurrency (E1 methodology) will differ with real I/O.
- Timings were taken under the bench lock but other agents' non-locked work (cargo builds etc.) shared the node.
