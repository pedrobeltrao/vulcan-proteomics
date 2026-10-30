#!/usr/bin/env bash
# M1e baselines end to end (re-runnable; finished steps are skipped, FORCE=1 redoes conversions).
#   bash scripts/baselines/run_all.sh
# Node-local state under /tmp/vulcan_m1e (mokapot env, entrapment DBs, mzML, search caches) is rebuilt by this script.
. "$(dirname "$0")/common.sh"
bash "$ROOT/scripts/baselines/setup_mokapot_env.sh"
[ -f "$ROOT/data/cache/M1e/entrapment/db/human_contam.P2003.json" ] && [ -f "$ROOT/data/cache/M1e/entrapment/db/pfu_contam.C.json" ] \
  || bash "$ROOT/scripts/entrapment/run_all.sh" > "$ROOT/data/cache/M1e/entrapment/build.log" 2>&1
bash "$ROOT/scripts/baselines/01_convert.sh"
# Pass 1: full node budget (32 affinity threads): all DBs. prod and P2001 are timed under the benchmark lock
# (sensitivity and FDP come from this pass); the other DBs are correctness replicates (timed, flagged unlocked).
python "$ROOT/scripts/baselines/run_searches.py" --datasets D1 D0 D4 --tools sage comet --lock-dbs prod P2001
# Pass 2: 4-thread continuity runs. Sage: prod and P2001 on all datasets. Comet: prod on all datasets, P2001 only on D1
# (Comet on a 2.4M-peptide DB at 4 threads is 80 min per file; Comet core-s is thread-independent to first order).
for t in 4; do
  VULCAN_THREADS=$t python "$ROOT/scripts/baselines/run_searches.py" --datasets D1 D0 D4 --tools sage --dbs prod P2001 --lock-dbs prod P2001
  VULCAN_THREADS=$t python "$ROOT/scripts/baselines/run_searches.py" --datasets D1 D0 D4 --tools comet --dbs prod --lock-dbs prod
  VULCAN_THREADS=$t python "$ROOT/scripts/baselines/run_searches.py" --datasets D1 --tools comet --dbs P2001 --lock-dbs P2001
done
# Decoy-variability arm: external shuffled decoys (DECOY_SEEDS), Sage on D1/D0/D4, Comet on D1 (peptide-level DBs; unlocked)
python "$ROOT/scripts/baselines/run_searches.py" --datasets D1 D0 D4 --tools sage --dbs S1001 S1002 S1003 S1004 S1005 --lock-dbs
python "$ROOT/scripts/baselines/run_searches.py" --datasets D1 --tools comet --dbs S1001 S1002 S1003 S1004 S1005 --lock-dbs
mkdir -p "$RES/entrapment"
cp "$ROOT"/data/cache/M1e/entrapment/db/*.json "$ROOT/data/cache/M1e/entrapment/SHA256SUMS" "$RES/entrapment/"
cp "$ROOT"/data/cache/M1e/entrapment/db/fdrbench/*.log "$ROOT/data/cache/M1e/entrapment/verify.log" "$RES/entrapment/" 2>/dev/null || true
python "$ROOT/scripts/baselines/write_env.py"
python "$ROOT/scripts/baselines/analyze.py"
python "$ROOT/scripts/baselines/make_summary.py" > /dev/null
