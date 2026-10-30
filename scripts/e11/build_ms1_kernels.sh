#!/usr/bin/env bash
# Build the E11 GPU MS1 kernels (vulcan/ms1/ms1_kernels.cu) into tools/vulcan_ms1_kernels/libvulcan_ms1.so.
# sm_80 (A100) + sm_86 (RTX 3090) SASS + compute_86 PTX (DESIGN §2 portability rule).
# -fmad=false: no FMA contraction, so f32/f64 arithmetic rounds exactly like the numba CPU path (bitwise equality).
# Static cudart: the library needs only the driver (no torch), so CPU decode workers can load it.
# Idempotent; vulcan/ms1/gpu.py calls this automatically when the recorded source hash differs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if ! command -v nvcc >/dev/null 2>&1; then
  set +u; . "$ROOT/env/activate.sh" >/dev/null 2>&1 || true; set -u
fi
if ! command -v nvcc >/dev/null 2>&1 && [ -n "${CUDA_HOME:-}" ]; then export PATH="$CUDA_HOME/bin:$PATH"; fi
SRC="$ROOT/vulcan/ms1/ms1_kernels.cu"
OUT="$ROOT/tools/vulcan_ms1_kernels"
GENCODE="${VULCAN_NVCC_GENCODE:--gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_86,code=compute_86}"
mkdir -p "$OUT"
# shellcheck disable=SC2086
nvcc -O3 -std=c++17 -shared -Xcompiler -fPIC -fmad=false -cudart static $GENCODE \
  -o "$OUT/libvulcan_ms1.so.partial.$$" "$SRC"
mv "$OUT/libvulcan_ms1.so.partial.$$" "$OUT/libvulcan_ms1.so"
sha256sum "$SRC" | cut -d' ' -f1 > "$OUT/source.sha256"
echo "built $OUT/libvulcan_ms1.so ($(nvcc --version | tail -1))"
