#!/usr/bin/env bash
# M4 timed GPU searches (each under the node benchmark lock), outputs to $W/bench/<name>/ (PSMs + stats.json).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W=${VULCAN_M4_WORK:-/tmp/vulcan_m4}
PY="$ROOT/scripts/m2b/vpy.sh"
D0="/tmp/vulcan_e2/store/b1906_293T_proteinID_01A_QE3_122212.vspec.arrow /tmp/vulcan_e2/store/b1937_293T_proteinID_01B_QE3_122212.vspec.arrow"
IDX0=/tmp/vulcan_e2/index/human_contam_prod
PIDX=$W/index/d2_phospho_prod
PV='{"S": 79.966331, "T": 79.966331, "Y": 79.966331, "M": 15.994915, "^": 42.010565}'
PP='{"chunk_size": 32768, "max_slice_entries": 536870912}'
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
run() { name=$1; shift; flock /tmp/vulcan_bench.lock "$PY" "$ROOT/scripts/m4/search_runs.py" --out "$W/bench/$name" --repeat 2 "$@" > "$W/logs/bench_$name.log" 2>&1; echo "done $name"; }
run d0_closed --index $IDX0 --var '{"M": 15.9949}' --stores $D0 --mode closed --report 5
run d0_offset --index $IDX0 --var '{"M": 15.9949}' --stores $D0 --mode offset --offsets-json $W/offsets_D0.json --report 5
run d0_open --index $IDX0 --var '{"M": 15.9949}' --stores $D0 --mode open --report 5
run d2_phospho --index $PIDX --var "$PV" --params "$PP" --stores $(cat $W/d2_stores.txt) --mode closed --report 1
run d0_phospho --index $PIDX --var "$PV" --params "$PP" --stores $D0 --mode closed --report 1
