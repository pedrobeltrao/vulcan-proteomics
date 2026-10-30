#!/usr/bin/env bash
# M2b side measurements, one bench-lock hold each: startup, H2D, Lustre read, ingest profile.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RES="$ROOT/results/2026-10-27_M2b_cohort"
L="${VULCAN_BENCH_LOCK:-/tmp/vulcan_bench.lock}"
VPY="$ROOT/scripts/m2b/vpy.sh"
S=/tmp/vulcan_e2/store
flock "$L" "$VPY" "$ROOT/scripts/m2b/startup_measure.py" --out "$RES/startup.json" --reps 3
flock "$L" "$VPY" "$ROOT/scripts/m2b/h2d_measure.py" --out "$RES/h2d.json" \
  --stores $S/b1906_293T_proteinID_01A_QE3_122212.vspec.arrow $S/b1937_293T_proteinID_01B_QE3_122212.vspec.arrow $S/Velos005137.vspec.arrow
flock "$L" bash "$ROOT/scripts/m2b/lustre_read.sh" "$RES/lustre_read.json"
flock "$L" "$VPY" "$ROOT/scripts/m2b/ingest_profile.py" --out "$RES/ingest_profile.json" --reps 3 \
  "$ROOT/data/raw/PXD001077/Velos005137.raw" "$ROOT/data/raw/PXD016766/dq_00095_11cell_90min_sps_A1.raw" \
  "$ROOT/data/raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw" "$ROOT/data/raw/PXD031616/LM2480_DDA.raw"
