#!/usr/bin/env bash
# Uncached Lustre read bandwidth for RAW files (O_DIRECT, bypasses the page cache), P concurrent distinct files.
#   scripts/m2b/lustre_read.sh OUT.json
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT=$1
mapfile -t F < <(grep -v '^#' "$ROOT/scripts/m2b/dev_dda_raws.txt" | sed "s|^|$ROOT/|" | xargs ls -S)   # largest first
echo "{" > "$OUT.partial"
first=1
for P in 1 4 8 16; do
  sel=("${F[@]:0:$P}")
  bytes=$(stat -c %s "${sel[@]}" | awk '{s+=$1} END {print s}')
  t0=$(date +%s.%N)
  printf '%s\n' "${sel[@]}" | xargs -P "$P" -I{} dd if={} of=/dev/null bs=64M iflag=direct status=none
  t1=$(date +%s.%N)
  sec=$(python3 -c "print($t1 - $t0)")
  gbps=$(python3 -c "print($bytes / ($t1 - $t0) / 1e9)")
  [ $first = 1 ] || echo "," >> "$OUT.partial"; first=0
  printf '"P%d": {"files": %d, "bytes": %d, "seconds": %.3f, "GBps": %.3f}' "$P" "$P" "$bytes" "$sec" "$gbps" >> "$OUT.partial"
  echo "P=$P ${gbps} GB/s"
done
echo "}" >> "$OUT.partial"
python3 - "$OUT.partial" "$OUT" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
d["aggregate_GBps_best"] = max(v["GBps"] for v in d.values())
d["note"] = "dd iflag=direct (no page cache), distinct dev RAW files, largest first; shared Lustre (shared-fs), other users active"
json.dump(d, open(sys.argv[2], "w"), indent=1)
EOF
rm -f "$OUT.partial"
