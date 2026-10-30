> **Public export note.** This file is the historical record of the two development nodes. The cluster-specific `env/setup.sh` and the original `env/pe.yml` are not included. To install, follow README.md (`environment.yml`, generic `env/activate.sh`); exact package versions are in `env/environment.lock.yml`. Host names and paths have been replaced by placeholders.

# Environment (set up 2026-10-26)

Node: 4 CPU cores, 1x RTX 3090, driver 580 (CUDA 13.0), 501 GB RAM, no root.
Reproduce: `bash env/setup.sh [pe_env rust sage trfp comet diann]`; use: `. env/activate.sh`.
`env/pe.yml` = `mamba env export --no-builds` of the pe env. `tools/` and `conda_envs/` are gitignored.

## pe env (`conda_envs/pe`, python 3.12.16, conda-forge + pip)
numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, polars 2.1.0, pyarrow 25.0.0, numba 0.68.0, scikit-learn 1.9.1,
matplotlib 3.11.2, pytest 9.1.1, maturin 1.15.0, pytorch 2.14.1 (CUDA 12.9 build; `torch.cuda.is_available()` True,
2000x2000 matmul on GPU OK). pip: pyteomics 4.7.5, mokapot 0.10.0, psm-utils 1.5.5.

Optional packages (pip, all installed and import OK):
- alphatims 1.1.1  (Bruker .d reader)
- timsrust_pyo3 0.4.1  (Rust-backed Bruker reader)
- pyopenms 3.6.0
- koinapy 0.0.11  (Koina client; needs internet to use the Koina server)
- psm-utils 1.5.5
None failed.

## Rust
rustup installed to `tools/rustup`, `tools/cargo` (minimal profile, stable): rustc/cargo 1.99.1 (2026-10-13).
Needs `RUSTUP_HOME`/`CARGO_HOME` (set by activate.sh); cargo without them fails by design.

## Engines
- Sage: release tag v0.14.7 (latest release), commit 99407db6e3754b31a9b88b7316a0aee67293c93f, source in `tools/sage`,
  binary `tools/bin/sage` (`sage --version` prints 0.14.6; crate version string not bumped upstream).
  Build gotcha: locked `time 0.3.34` fails on current rustc (E0282); fixed with `cargo update -p time` (local Cargo.lock change only).
  Build ~5 min on 4 cores.
- ThermoRawFileParser 2.0.0.dev (bioconda) in `conda_envs/trfp`; `ThermoRawFileParser --help` works.
- Comet 2026.01 rev.1 (e4f767c) (bioconda comet-ms) in `conda_envs/comet`; binary `comet`. Runs.
- DIA-NN 2.7.0 Academia Linux (github.com/vdemichev/DiaNN releases; compiled Sep 15 2026), in `tools/diann/diann-2.7.0/diann-linux`.
  Runs (prints banner). Licence: the "Academia" edition is free for academic/non-profit use only; commercial use needs a
  different licence. The zip contains no LICENSE file; terms are in the DiaNN repo (LICENSE.md refers to a per-version LICENSE.txt).
  Check terms before any non-academic use or redistribution. Note: it reports 128 logical cores on this 4-core node; pass `--threads 4`.
- MSFragger: NOT downloaded (requires licence click-through at msfragger-upgrader.nesvilab.org). Needs manual action if wanted.

## Smoke test
Sage v0.14.7 on repo test data (tests/config.json adapted to absolute paths, `tools/smoke/config.json`):
1 spectrum mzML, Q99536 FASTA. Wall time 1.40 s (search itself ~0 ms, DB build 46 ms), max RSS 18.5 MB. It identified
LQSRPAAPPAPGPGQLTLR (hyperscore 72.3, charge 3). FDR counts are 0 because there is only one spectrum (expected; not a failure).
Cloud/tooling caveats: first import of heavy python packages (torch, pyopenms) is slow (~1-2 min) on cluster FS.

---

# New-node section (2026-10-27, decision D5). The sections above are history of the 4-core/3090 node.

Node gpu-node: 32 affinity cores (of 64 logical; 2-socket AMD EPYC 9354, Zen 4, SMT off, AVX-512 incl. VNNI/BF16/VPOPCNTDQ),
NVIDIA A100 80GB PCIe (sm_80, driver 580.178.04, CUDA 13.0 driver), 1.5 TB RAM, NVMe `/tmp` 862 GB (node-local, vanishes with the node).
Fingerprint of any node: `python bench/hwinfo.py [out.json]` (`vulcan.hw.fingerprint()`).

## Rules that make this portable
- Threads: `vulcan.hw.n_threads()` = affinity mask size, overridden by `$VULCAN_THREADS`. `VULCAN_THREADS=4` reproduces the M1 4-thread runs
  (`timed.py --threads 4` for one command). `vulcan.hw.pin_thread_pools(n)` / `hw.thread_env(n)` set OMP/MKL/OPENBLAS/NUMBA/RAYON/NUMEXPR/VECLIB/VULCAN_THREADS
  (and torch/pyarrow/numba live if imported). `os.cpu_count()` is wrong here (64 vs 32 usable).
- Every benchmark record carries `hw` (the fingerprint): `scripts/baselines/timed.py` (each timings.jsonl line), `scripts/ingest/e1_bench.py`,
  `scripts/ingest/hw_ingest_check.py`, `scripts/baselines/write_env.py`, `scripts/run_stats_sims.py`.
- Benchmark lock unchanged: `flock /tmp/vulcan_bench.lock` (override path with `$VULCAN_BENCH_LOCK`; the lock file is not state).
- Persistent state only under the project: cargo targets in `data/cache/cargo_target/<crate>` (gitignored), conda pkg cache `data/cache/conda_pkgs`,
  torch extension builds `data/cache/torch_extensions` (`TORCH_EXTENSIONS_DIR`, set by activate.sh). `/tmp` holds only rebuildable scratch
  (`$VULCAN_M1E_WORK` staging, E1 stores, the synthetic .NET root). Rebuild scripts: `env/setup.sh`, `scripts/baselines/01_convert.sh`,
  `scripts/baselines/setup_mokapot_env.sh`, `scripts/m2/build_vulcan_core.sh`, `scripts/ingest/build_vulcan_tims.sh`, `scripts/ingest/build_vulcanraw.sh`.
- CUDA targets: sm_80 + sm_86 + PTX. activate.sh exports `TORCH_CUDA_ARCH_LIST="8.0;8.6+PTX"`, `CUDA_HOME=conda_envs/pe` and
  `VULCAN_NVCC_GENCODE="-gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_86,code=compute_86"`; use it on every nvcc line.

## Toolchain changes on this node
- pe env had nvcc but no CUDA headers or cicc (conda-forge `cuda-nvcc-tools` only), so nvcc could not compile. Added (also in `env/setup.sh`):
  cuda-nvcc-impl, cuda-nvvm-impl, cuda-cudart-dev, cuda-cccl, cuda-driver-dev, cuda-nvrtc-dev, cuda-profiler-api, cuda-nvtx-dev, cuda-cupti-dev,
  libcublas/cusparse/cusolver/curand/cufft-dev, ninja (all CUDA 12.9). Host compiler is the system gcc 13.3 (supported by CUDA 12.9).
  `env/pe.yml` should be regenerated with `mamba env export --no-builds` if the env is rebuilt.
- Sage: stock `tools/bin/sage` (v0.14.7, 99407db6, `--version` prints 0.14.6) is kept. `tools/bin/sage-znver4` is the same source built with
  `RUSTFLAGS="-C target-cpu=znver4"` (separate target dir `tools/sage/target-znver4`, 3m57s on 32 cores). Measured no speed gain and identical output
  (see results/2026-10-27_hw_adapt); stock stays the comparator default.
- Sage at 32 threads on one small file wastes CPU (spinning rayon workers): ~56 core-s vs 8.6 at 4 threads for the same ~5 s wall. Count core-s, not
  just wall, when reporting full-node comparator numbers, and look at scaling per file size before choosing thread counts.
- DIA-NN reports "Logical CPU cores: 64" (the host); always pass `--threads N` with N from `vulcan.hw.n_threads()`.

## Verified on this node (details: results/2026-10-27_hw_adapt/{metrics.json,NOTES.md})
torch 2.14.1 (CUDA 12.9) on the A100: matmul and bf16 OK; nvcc 12.9 compiles and runs a kernel for sm_80+sm_86+PTX; torch `load_inline` CUDA extension
builds and runs; cargo/rustc 1.99.1; Sage smoke test; Comet 2026.01; ThermoRawFileParser 2.0.0.0; DIA-NN 2.7.0; direct Thermo reader; PyO3 crate vulcan_core.
