# Generic environment setup for Vulcan. Source it from the repository root after activating the conda env:
#     conda activate vulcan        # created from environment.yml
#     . env/activate.sh
# Every variable below can be overridden by exporting it before sourcing this file.
_R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export VULCAN_ROOT="${VULCAN_ROOT:-$_R}" PROJECT_ROOT="${PROJECT_ROOT:-$_R}"
# Node-local scratch for staged RAW files, Arrow stores and work dirs (fast local disk recommended).
export VULCAN_SCRATCH="${VULCAN_SCRATCH:-${TMPDIR:-/tmp}}"
# Optional project-local Rust toolchain (otherwise the rustup/cargo on PATH is used).
if [ -x "$_R/tools/cargo/bin/cargo" ]; then export RUSTUP_HOME="$_R/tools/rustup" CARGO_HOME="$_R/tools/cargo"; export PATH="$CARGO_HOME/bin:$PATH"; fi
# Optional comparator binaries (Sage, DIA-NN, ...) placed in tools/bin.
if [ -d "$_R/tools/bin" ]; then export PATH="$_R/tools/bin:$PATH"; fi
# CUDA: conda's nvcc lives in the active env; builds target sm_80 (A100) + sm_86 (RTX 30xx / A40) + PTX.
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/nvcc" ]; then export CUDA_HOME="${CUDA_HOME:-$CONDA_PREFIX}"; fi
export TORCH_CUDA_ARCH_LIST="${TORCH_CUDA_ARCH_LIST:-8.0;8.6+PTX}"
export VULCAN_NVCC_GENCODE="${VULCAN_NVCC_GENCODE:--gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_86,code=compute_86}"
# Persistent build cache for torch CUDA extensions.
export TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$_R/data/cache/torch_extensions}"
# Thermo RAW reading (vulcan/io/thermo_direct.py): VULCAN_TRFP_BIN = directory holding ThermoFisher.CommonCore.*.dll
# (obtained separately under Thermo's licence, e.g. a ThermoRawFileParser install; default conda_envs/trfp/bin),
# VULCAN_RAW_DLL (default tools/vulcanraw/VulcanRaw.dll), VULCAN_DOTNET_ROOT (default tools/dotnet-sdk).
# Threads: vulcan.hw uses the CPU affinity mask; set VULCAN_THREADS to override.
