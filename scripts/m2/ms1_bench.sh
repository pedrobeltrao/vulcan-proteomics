#!/usr/bin/env bash
# Timed MS1 feature finder + hypothesis table runs (benchmark lock; core-s from /usr/bin/time -v via timed.py).
# Calibrated hypotheses use the pilot calibration JSONs if present in $RES.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
WORK=${VULCAN_MS1_WORK:-/tmp/vulcan_ms1}
RES=${RES:-$ROOT/results/2026-10-27_M2_ms1_pilot}
cd "$ROOT"
for run in b1906_293T_proteinID_01A_QE3_122212 b1937_293T_proteinID_01B_QE3_122212 Velos005137; do
  cal=""; [ -f "$RES/$run.calib.json" ] && cal="--calib $RES/$run.calib.json"
  for thr in ${THREADS_LIST:-32 4}; do
    for rep in 1 2; do
      python scripts/baselines/timed.py --lock --threads $thr --label "vulcan_ms1:$run:t$thr:r$rep" --log "$RES/timings.jsonl" \
        --stdout "$WORK/bench_${run}_t${thr}_r${rep}.log" -- env PYTHONPATH="$ROOT" python -m vulcan.ms1 \
        "$WORK/store/$run.vspec.arrow" "$WORK/ms1_bench/t$thr" --threads $thr $cal
    done
  done
done
