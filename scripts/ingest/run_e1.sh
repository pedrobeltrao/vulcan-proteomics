#!/usr/bin/env bash
# E1 ingest benchmark: env capture + driver under the node-wide bench lock (bench/README.md).
# usage: bash scripts/ingest/run_e1.sh [results_dir] [extra driver args]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${1:-$ROOT/results/2026-10-26_E1_ingest}"; shift || true
mkdir -p "$OUT"
set +u; . "$ROOT/env/activate.sh"; set -u
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 RAYON_NUM_THREADS=1 VULCAN_THREADS=1
export DOTNET_PROCESSOR_COUNT=${DOTNET_PROCESSOR_COUNT:-4}  # per-process .NET runtime cap; the bench itself runs 1 thread per process
( cd "$ROOT" && git rev-parse HEAD && (git diff --quiet HEAD -- . 2>/dev/null && echo clean || echo "dirty (uncommitted vulcan/io, scripts/ingest)") ) > "$OUT/git_sha"
{
  echo "date: $(date -Is)"; echo "host: $(hostname)"; echo "kernel: $(uname -r)"
  echo "nproc: $(nproc)  affinity: $(taskset -pc $$ | sed 's/.*: //')"
  lscpu | grep -E "Model name|Socket|NUMA node\(s\)" ; free -g | head -2
  df -h /tmp "$ROOT" | sed 1d
  echo "cgroup: $(cat /proc/self/cgroup)"
  echo "python: $(python -V 2>&1)"; python -c "import numpy, pyarrow, pythonnet, clr_loader, pyteomics; print('numpy', numpy.__version__, 'pyarrow', pyarrow.__version__, 'pythonnet', pythonnet.__version__ if hasattr(pythonnet,'__version__') else '3.2.1', 'pyteomics', pyteomics.__version__ if hasattr(pyteomics,'__version__') else '4.7.5')" 2>/dev/null || true
  pip list 2>/dev/null | grep -i -E "^(timsrust|pythonnet|clr.loader|pyteomics|pyarrow) " || true
  echo "ThermoRawFileParser: $(ThermoRawFileParser --version 2>&1 | tail -1)"
  echo "dotnet runtime (pythonnet host): $(ls "$ROOT/tools/dotnet-sdk/shared/Microsoft.NETCore.App" 2>/dev/null | tr '\n' ' ') (SDK $( "$ROOT/tools/dotnet-sdk/dotnet" --version 2>/dev/null || echo n/a))"
  echo "RawFileReader: $(strings "$ROOT/conda_envs/trfp/bin/ThermoRawFileParser.deps.json" | grep -o 'ThermoFisher.CommonCore.RawFileReader/[0-9.]*' | head -1)"
  echo "rustc: $(rustc --version 2>/dev/null); vulcan_tims: timsrust 0.4.2 + pyo3 0.22 (vulcan/io/rust/vulcan_tims)"
  echo "threads: OMP=$OMP_NUM_THREADS RAYON=$RAYON_NUM_THREADS VULCAN_THREADS=$VULCAN_THREADS DOTNET_PROCESSOR_COUNT=$DOTNET_PROCESSOR_COUNT"
  echo "lock: flock ${VULCAN_BENCH_LOCK:-/tmp/vulcan_bench.lock}"
  python "$ROOT/bench/hwinfo.py" | tr -d "\n"; echo
} > "$OUT/env.txt" 2>&1
echo "[run_e1] waiting for /tmp/vulcan_bench.lock ..."
flock "${VULCAN_BENCH_LOCK:-/tmp/vulcan_bench.lock}" python "$ROOT/scripts/ingest/e1_bench.py" --out "$OUT" "$@"
