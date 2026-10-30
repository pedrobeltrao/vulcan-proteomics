#!/usr/bin/env bash
# Stage the pred venv (pe + peptdeep) onto node-local NVMe next to the staged pe env (scripts/m2b/stage_env.sh), so
# score workers never import from Lustre (M2b lesson: Lustre imports took 13-180 s, minutes under load).  Idempotent.
#   scripts/e9/stage_pred_env.sh  -> prints the staged python path
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="${VULCAN_LOCAL_ENV:-/tmp/vulcan_env_$(id -u)}"
[ -x "$DEST/pe/bin/python" ] || "$ROOT/scripts/m2b/stage_env.sh" >&2
if [ ! -f "$DEST/pred.stamp" ]; then
  rm -rf "$DEST/pred.partial"
  cp -a "$ROOT/conda_envs/pred" "$DEST/pred.partial"
  sed -i "s#^home = .*#home = $DEST/pe/bin#; s#^executable = .*#executable = $DEST/pe/bin/python3.12#" "$DEST/pred.partial/pyvenv.cfg"
  for f in python python3 python3.12; do ln -sfn "$DEST/pe/bin/python3.12" "$DEST/pred.partial/bin/$f"; done
  mv "$DEST/pred.partial" "$DEST/pred"
  date > "$DEST/pred.stamp"
fi
echo "$DEST/pred/bin/python"
