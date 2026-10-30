#!/usr/bin/env bash
# 3c dev L3 features (Koina Prosit CID / TMT-CID; label-blind) for the dev PSM directories, after dev_search.sh.
set -euo pipefail
cd "$(dirname "$0")/../.."
. env/activate.sh >/dev/null 2>&1
for ds in ${DSETS:-D10cf D4}; do
  fam=$([ $ds = D4 ] && echo human_tmt || echo human)
  for k in ${KINDS:-prod E2001 E2002 E2003 C}; do
    while [ ! -f /tmp/vulcan_3c/psms/${ds}_$k/stats.json ]; do sleep 30; done
    python scripts/m3c/features.py --psms /tmp/vulcan_3c/psms/${ds}_$k --family $fam --index-kind $k | tail -1
    echo "$ds $k features done"
  done
done
