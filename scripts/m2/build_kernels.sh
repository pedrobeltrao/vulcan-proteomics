#!/usr/bin/env bash
# Build the Vulcan CUDA kernels (vulcan/gpu/kernels.cu) into tools/vulcan_kernels/libvulcan_kernels.so.
# sm_80 (A100) + sm_86 (RTX 3090) SASS + compute_86 PTX (DESIGN §2 portability rule).
# -fmad=false: no FMA contraction, so f32 fragment/tolerance arithmetic rounds exactly like Sage (Rust).
# Idempotent; vulcan/gpu/_lib.py calls this automatically when the .so is older than the source.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh" >/dev/null 2>&1 || true; set -u
SRC="$ROOT/vulcan/gpu/kernels.cu"
OUT="$ROOT/tools/vulcan_kernels"
mkdir -p "$OUT"
nvcc -O3 -std=c++17 -shared -Xcompiler -fPIC -fmad=false \
  -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_86,code=compute_86 \
  -o "$OUT/libvulcan_kernels.so.partial" "$SRC"
mv "$OUT/libvulcan_kernels.so.partial" "$OUT/libvulcan_kernels.so"
sha256sum "$SRC" | cut -d' ' -f1 > "$OUT/source.sha256"
echo "built $OUT/libvulcan_kernels.so ($(nvcc --version | tail -1))"
