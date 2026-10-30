#!/usr/bin/env bash
# 3c dev searches (data generation; not timed): D10cf (human) and D4 (human_tmt) x {prod, E2001-3, C}
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."   # repository root
. env/activate.sh >/dev/null 2>&1
for ds in ${DSETS:-D10cf D4}; do
  fam=$([ $ds = D4 ] && echo human_tmt || echo human)
  stores=$(python -c "import sys; sys.path.insert(0,'scripts/m3c'); import common as C; print(' '.join(C.DEV['$ds']['stores']()))")
  for k in ${KINDS:-prod E2001 E2002 E2003 C}; do
    out=/tmp/vulcan_3c/psms/${ds}_$k
    [ -f $out/stats.json ] && continue
    python scripts/m3c/search.py --family $fam --kind $k --stores $stores --out $out | tail -1
    echo "$ds $k searched"
  done
done
