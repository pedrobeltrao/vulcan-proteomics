#!/usr/bin/env bash
# Run python from the node-local staged env (scripts/m2b/stage_env.sh) with the staged .NET/RawFileReader; falls back
# to the Lustre env if nothing is staged.   scripts/m2b/vpy.sh -m vulcan.cohort ...
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="${VULCAN_LOCAL_ENV:-/tmp/vulcan_env_$(id -u)}"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
if [ -x "$DEST/pe/bin/python" ]; then
  PY="$DEST/pe/bin/python"
  [ -f "$DEST/vulcanraw/VulcanRaw.dll" ] && export VULCAN_DOTNET_ROOT="$DEST/dotnet" VULCAN_TRFP_BIN="$DEST/trfp" \
     VULCAN_RAW_DLL="$DEST/vulcanraw/VulcanRaw.dll"
else
  PY="$ROOT/conda_envs/pe/bin/python"
fi
exec "$PY" "$@"
