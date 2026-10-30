#!/usr/bin/env bash
# Build the numpy-returning Bruker TDF reader (vulcan/io/rust/vulcan_tims; timsrust 0.4.2 + pyo3 0.22)
# into tools/vulcan_tims/vulcan_tims.abi3.so (imported by vulcan/io/bruker.py). 
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
# persistent target dir (survives node changes); jobs from the affinity mask or $VULCAN_THREADS
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$ROOT/data/cache/cargo_target/vulcan_tims}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-$(python -c "import sys; sys.path.insert(0,'$ROOT'); from vulcan import hw; print(hw.n_threads())")}"
mkdir -p "$CARGO_TARGET_DIR"
cd "$ROOT/vulcan/io/rust/vulcan_tims"
cargo build --release
mkdir -p "$ROOT/tools/vulcan_tims"
cp "$CARGO_TARGET_DIR/release/libvulcan_tims.so" "$ROOT/tools/vulcan_tims/vulcan_tims.abi3.so"
echo "built $ROOT/tools/vulcan_tims/vulcan_tims.abi3.so"
