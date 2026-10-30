#!/usr/bin/env bash
# M2b cohort scaling: one fresh cohort process per config, each under the bench lock (one lock hold per config).
#   scripts/m2b/cohort_bench.sh TAG "WORKERS..." [extra cohort args]
# Stats -> results/2026-10-27_M2b_cohort/runs/<TAG>_w<N>.json ; PSMs/stage on /tmp (deleted after each config).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TAG=$1; WS=$2; shift 2
RES="$ROOT/results/2026-10-27_M2b_cohort/runs"; mkdir -p "$RES"
WORK="${VULCAN_M2B_WORK:-/tmp/vulcan_m2b_$(id -u)}"
INDEX="${VULCAN_M2B_INDEX:-/tmp/vulcan_e2/index/human_contam_prod}"
for W in $WS; do
  OUT="$WORK/cohort_${TAG}_w$W"; rm -rf "$OUT" "$OUT.stage"
  echo "[$(date +%T)] $TAG workers=$W (waiting for lock)"
  flock "${VULCAN_BENCH_LOCK:-/tmp/vulcan_bench.lock}" \
    "$ROOT/scripts/m2b/vpy.sh" -m vulcan.cohort --index "$INDEX" --raw-list "$ROOT/scripts/m2b/dev_dda_raws.txt" \
      --out "$OUT" --stage "$OUT.stage" --workers "$W" --stats "$RES/${TAG}_w$W.json" "$@" > "$WORK/${TAG}_w$W.log" 2>&1 \
    || { echo "FAILED (see $WORK/${TAG}_w$W.log)"; tail -5 "$WORK/${TAG}_w$W.log"; }
  tail -1 "$WORK/${TAG}_w$W.log" | cut -c1-400
  rm -rf "$OUT" "$OUT.stage"
done
