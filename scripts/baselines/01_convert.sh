#!/usr/bin/env bash
# Stage raw files to local /tmp and convert RAW -> indexed mzML with ThermoRawFileParser (timed, under the lock).
# Both Sage 0.14.7 and Comet need mzML, so this conversion is part of both tools' end-to-end cost.
# TRFP flags: -f 2 (indexed mzML), vendor peak picking ON (default), zlib ON (default), -l 3 (warnings).
. "$(dirname "$0")/common.sh"
mkdir -p "$WORK/raw" "$WORK/mzml" "$WORK/meta"
convert() {  # ds rawpath
  local ds=$1 raw=$2 b; b=$(basename "$raw" .raw)
  [ -f "$WORK/raw/$b.raw" ] || cp "$raw" "$WORK/raw/$b.raw"
  [ -f "$WORK/meta/$b-metadata.json" ] || ThermoRawFileParser -i "$WORK/raw/$b.raw" -f 4 -m 0 -o "$WORK/meta" -l 4
  if [ ! -f "$WORK/mzml/$b.mzML" ] || [ "${FORCE:-0}" = 1 ]; then
    rm -f "$WORK/mzml/$b.mzML"
    $TIMED --lock --label "trfp_mzml:$ds:$b" --log "$TIMINGS" --stdout "$RES/logs/trfp_$b.log" \
      --meta "{\"stage\":\"conversion\",\"tool\":\"trfp\",\"dataset\":\"$ds\",\"file\":\"$b\"}" \
      -- ThermoRawFileParser -i "$WORK/raw/$b.raw" -o "$WORK/mzml" -f 2 -l 3
  fi
}
for r in "${D0_RAWS[@]}"; do convert D0 "$r"; done
for r in "${D1_RAWS[@]}"; do convert D1 "$r"; done
if [ -d "$ROOT/data/raw/PXD016766" ]; then
  for r in "$ROOT"/data/raw/PXD016766/*.raw; do [ -f "$r" ] && convert D4 "$r"; done
fi
ls -la "$WORK/mzml"
