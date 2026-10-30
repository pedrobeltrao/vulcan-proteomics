# Shared paths for the M1e baseline scripts.  Source me.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
export WORK=${VULCAN_M1E_WORK:-/tmp/vulcan_m1e}          # local-disk staging (DESIGN §2: never run off Lustre)
export CACHE="$WORK/cache"                                # search outputs (local disk; Lustre quota nearly full)
export RES="$ROOT/results/2026-10-27_M1e_baselines"
export TIMINGS="$RES/timings.jsonl"
export TIMED="python $ROOT/scripts/baselines/timed.py"
# Thread budget: $VULCAN_THREADS if set (VULCAN_THREADS=4 reproduces the M1 4-core runs), else the affinity mask.
export VULCAN_THREADS=${VULCAN_THREADS:-$(python -c "import sys; sys.path.insert(0,'$ROOT'); from vulcan import hw; print(hw.n_threads())")}
export OMP_NUM_THREADS=$VULCAN_THREADS MKL_NUM_THREADS=$VULCAN_THREADS NUMBA_NUM_THREADS=$VULCAN_THREADS
export RAYON_NUM_THREADS=$VULCAN_THREADS OPENBLAS_NUM_THREADS=$VULCAN_THREADS
mkdir -p "$WORK" "$CACHE" "$RES/configs" "$RES/logs"
# dataset -> raw files
D0_RAWS=("$ROOT/data/raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw" "$ROOT/data/raw/PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw")
D1_RAWS=("$ROOT/data/raw/PXD001077/Velos005137.raw")
