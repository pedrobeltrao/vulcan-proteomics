#!/usr/bin/env bash
# Idempotent staging for the M2 MS1/pilot work: Vulcan stores (direct Thermo) + TRFP mzML (for biosaur2 / Sage parity)
# on node-local disk.  Rebuilds whatever is missing; never touches data/sealed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
WORK=${VULCAN_MS1_WORK:-/tmp/vulcan_ms1}
mkdir -p "$WORK/store" "$WORK/mzml"
RAWS=("$ROOT/data/raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw"
      "$ROOT/data/raw/PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw"
      "$ROOT/data/raw/PXD001077/Velos005137.raw")
for r in "${RAWS[@]}"; do
  b=$(basename "$r" .raw)
  if [ ! -s "$WORK/store/$b.vspec.arrow" ]; then
    python -m vulcan.io.convert "$r" "$WORK/store/$b.vspec.arrow" > "$WORK/store/$b.log" 2>&1
  fi
  if [ "${WITH_MZML:-1}" = 1 ] && [ ! -s "$WORK/mzml/$b.mzML" ]; then
    ThermoRawFileParser -i "$r" -o "$WORK/mzml" -f 1 > "$WORK/mzml/$b.trfp.log" 2>&1
  fi
done
ls -la "$WORK/store" "$WORK/mzml"
