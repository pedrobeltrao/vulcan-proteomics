#!/usr/bin/env bash
# M3 searches (data generation, not timed): dev D0/D1 x {prod, E2001-3, C} and training runs, top-5 PSMs per spectrum.
#   KINDS="prod E2001" scripts/m3/search_dev.sh      TRAIN=1 scripts/m3/search_dev.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
W=${VULCAN_M3_WORK:-/tmp/vulcan_m3}
E2=${VULCAN_E2_WORK:-/tmp/vulcan_e2}
P='{"report_psms": 5}'
D0="$E2/store/b1906_293T_proteinID_01A_QE3_122212.vspec.arrow $E2/store/b1937_293T_proteinID_01B_QE3_122212.vspec.arrow"
D1="$E2/store/Velos005137.vspec.arrow"
idx() { case $2 in prod|E2001) if [ "$1" = human_contam ] || [ "$1" = pfu_contam ]; then echo "$E2/index/$1_$2"; else echo "$W/index/$1_$2"; fi;; *) echo "$W/index/$1_$2";; esac; }
if [ "${TRAIN:-0}" = 1 ]; then
  for sp in human mouse ecoli; do
    stores=$(python - <<EOF
import sys, os; sys.path.insert(0, "$ROOT/scripts/m3"); import m3_common as C
print(" ".join(C.train_store(r) for r in C.train_rows() if r["species"] == "$sp" and os.path.exists(C.train_store(r))))
EOF
)
    out="$W/psms/train_$sp"
    [ -f "$out/stats.json" ] && continue
    python -m vulcan.gpu.run --index "$(idx ${sp}_contam prod)" --stores $stores --out "$out" --params "$P" | cut -c1-200
  done
  exit 0
fi
for kind in ${KINDS:-prod E2001 E2002 E2003 C}; do
  for ds in D0 D1; do
    fam=$([ $ds = D0 ] && echo human_contam || echo pfu_contam)
    out="$W/psms/${ds}_${kind}"
    [ -f "$out/stats.json" ] && continue
    stores=$([ $ds = D0 ] && echo $D0 || echo $D1)
    python -m vulcan.gpu.run --index "$(idx $fam $kind)" --stores $stores --out "$out" --params "$P" | cut -c1-200
  done
done
