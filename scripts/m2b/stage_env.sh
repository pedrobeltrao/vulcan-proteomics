#!/usr/bin/env bash
# Stage the `pe` conda env (python, torch, pyarrow, pythonnet) onto node-local NVMe so that process startup
# does not pay Lustre metadata/latency (M2b: torch import 4.6 s warm / 183 s cold from Lustre).
# Idempotent: re-copies only if conda-meta/history changed. Usage: scripts/m2b/stage_env.sh [DEST]
# Then run   $DEST/pe/bin/python -m vulcan.cohort ...   (PYTHONPATH=<repo>). The Lustre env stays the source of truth.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$ROOT/conda_envs/pe"
DEST="${1:-${VULCAN_LOCAL_ENV:-/tmp/vulcan_env_$(id -u)}}"
STAMP="$DEST/pe.stamp"
want="$(stat -c '%Y %s' "$SRC/conda-meta/history")"
if [ -x "$DEST/pe/bin/python" ] && [ -f "$STAMP" ] && [ "$(cat "$STAMP")" = "$want" ]; then
  echo "$DEST/pe up to date"; SKIP_PE=1
fi
mkdir -p "$DEST"
if [ -z "${SKIP_PE:-}" ]; then
rm -rf "$DEST/pe.partial"
mkdir -p "$DEST/pe.partial"
t0=$(date +%s.%N)
# parallel tar pipes per lib/ and site-packages entry (no rsync on the node; one serial stream crawled on a loaded
# Lustre); skip the package cache and static libs, which are not needed at runtime.
copy() { mkdir -p "$DEST/pe.partial/$(dirname "$1")"; tar -C "$SRC" --exclude='*.a' -cf - "$1" | tar -C "$DEST/pe.partial" -xf -; }
export -f copy; export SRC DEST
SP=lib/python3.12/site-packages
( cd "$SRC"
  { for e in *; do [ "$e" = lib ] || [ "$e" = pkgs ] || echo "$e"; done
    for e in lib/*; do [ "$e" = lib/python3.12 ] || echo "$e"; done
    for e in lib/python3.12/*; do [ "$e" = "$SP" ] || echo "$e"; done
    for e in "$SP"/*; do echo "$e"; done; } | xargs -P 16 -I{} bash -c 'copy "$@"' _ {} )
rm -rf "$DEST/pe"
mv "$DEST/pe.partial" "$DEST/pe"
echo "$want" > "$STAMP"
echo "staged $SRC -> $DEST/pe in $(python3 -c "print(round($(date +%s.%N) - $t0, 1))") s ($(du -sh "$DEST/pe" | cut -f1))"
fi

# .NET runtime + RawFileReader + VulcanRaw helper (small; always refreshed). Use with:
#   VULCAN_DOTNET_ROOT=$DEST/dotnet VULCAN_TRFP_BIN=$DEST/trfp VULCAN_RAW_DLL=$DEST/vulcanraw/VulcanRaw.dll
mkdir -p "$DEST/dotnet/shared" "$DEST/trfp" "$DEST/vulcanraw"
cp -a "$ROOT/tools/dotnet-sdk/host" "$ROOT/tools/dotnet-sdk/dotnet" "$DEST/dotnet/"
cp -a "$ROOT/tools/dotnet-sdk/shared/Microsoft.NETCore.App" "$DEST/dotnet/shared/"
cp -a "$ROOT"/conda_envs/trfp/bin/ThermoFisher.CommonCore.*.dll "$DEST/trfp/"
cp -a "$ROOT"/tools/vulcanraw/* "$DEST/vulcanraw/"
echo "staged .NET runtime + RawFileReader + VulcanRaw -> $DEST/{dotnet,trfp,vulcanraw}"
