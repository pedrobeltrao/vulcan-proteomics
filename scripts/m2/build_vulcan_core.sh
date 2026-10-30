#!/usr/bin/env bash
# Build the PyO3 crate vulcan/rust/vulcan_core (sage-core 99407db6 + pyo3 0.22) into
# tools/vulcan_core/vulcan_core.abi3.so (imported by vulcan/index/_core.py). 
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
# persistent target dir (survives node changes); jobs from the affinity mask or $VULCAN_THREADS
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$ROOT/data/cache/cargo_target/vulcan_core}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-$(python -c "import sys; sys.path.insert(0,'$ROOT'); from vulcan import hw; print(hw.n_threads())")}"
mkdir -p "$CARGO_TARGET_DIR"
cd "$ROOT/vulcan/rust/vulcan_core"
cargo build --release
mkdir -p "$ROOT/tools/vulcan_core"
cp "$CARGO_TARGET_DIR/release/libvulcan_core.so" "$ROOT/tools/vulcan_core/vulcan_core.abi3.so"
echo "built $ROOT/tools/vulcan_core/vulcan_core.abi3.so"
