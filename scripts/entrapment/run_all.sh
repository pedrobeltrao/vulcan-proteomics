#!/usr/bin/env bash
# Build every M1e search database (D0: human+contaminants; D1: Pfu+contaminants) and record checksums.
# Seeds: ENTRAP_SEEDS_DEV = 2001 2002 2003 (DESIGN §6.1). Decoy seed for collision shuffles: 1001.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
set +u; . "$ROOT/env/activate.sh"; set -u
# Heavy files live on local disk (the project Lustre quota is nearly full); data/cache/M1e/entrapment/{db,inputs}
# are symlinks to them. Manifests (SHA256SUMS, build.log) stay in data/cache/M1e/entrapment.
W=${VULCAN_M1E_WORK:-/tmp/vulcan_m1e}/entrapment; mkdir -p "$W/inputs" "$W/db"
E="$ROOT/data/cache/M1e/entrapment"; mkdir -p "$E"
for d in inputs db; do [ -e "$E/$d" ] || ln -s "$W/$d" "$E/$d"; done
F="$ROOT/data/fasta"
cat "$F/human_sp_canonical.fasta" "$F/contaminants_frankenfield2022.fasta" > "$E/inputs/human_contam.fasta"
cat "$F/pfu_UP000001013.fasta"    "$F/contaminants_frankenfield2022.fasta" > "$E/inputs/pfu_contam.fasta"
B="python $ROOT/scripts/entrapment/build_dbs.py"
# Independent builds run in parallel (each repair pass is single-threaded Python; a sequential build takes ~1.5 h).
# D0 (HEK293): P = FDRBench shuffles x3 seeds; C = Arabidopsis foreign proteome
for s in 2001 2002 2003; do
  $B --name human_contam --target "$E/inputs/human_contam.fasta" --seeds $s --arms P --out "$E/db" > "$E/build_human_P$s.log" 2>&1 &
done
$B --name human_contam --target "$E/inputs/human_contam.fasta" --foreign "$F/arath_UP000006548.fasta" --arms target C S --out "$E/db" \
   > "$E/build_human_C.log" 2>&1 &
# D1 (P. furiosus): C = human proteome as entrapment (r = |human| / |Pfu|); P = Pfu shuffles x3 seeds
$B --name pfu_contam --target "$E/inputs/pfu_contam.fasta" --foreign "$F/human_sp_canonical.fasta" \
   --seeds 2001 2002 2003 --arms target P C S --out "$E/db" > "$E/build_pfu.log" 2>&1 &
wait
python "$ROOT/scripts/entrapment/verify_dbs.py" "$E/db" > "$E/verify.log"
(cd "$ROOT" && sha256sum data/cache/M1e/entrapment/inputs/*.fasta data/cache/M1e/entrapment/db/*.fasta \
   data/cache/M1e/entrapment/db/*.peptides.tsv data/cache/M1e/entrapment/db/fdrbench/*.txt) > "$E/SHA256SUMS"
cat "$E/verify.log"
