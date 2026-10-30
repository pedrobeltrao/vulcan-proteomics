#!/usr/bin/env bash
# M6 runtime per stage (each stage timed under flock /tmp/vulcan_bench.lock by scripts/baselines/timed.py, records
# carry the vulcan.hw fingerprint).  D9 = 6 HYE Astral runs; D4 = 12 TMT fractions.  One-shot processes (setup
# included); the cohort driver amortises setup (M2b).   scripts/m6/bench.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
W=${VULCAN_M6_WORK:-/tmp/vulcan_m6}
B=$W/bench
RES=$ROOT/results/2026-10-27_M6_quant
LOG=$RES/timings.jsonl
V="$ROOT/scripts/m2b/vpy.sh"
T() { local label=$1; shift; "$V" "$ROOT/scripts/baselines/timed.py" --label "$label" --log "$LOG" --lock --stdout "$B/$label.log" -- "$@"; }
mkdir -p "$B"/{store,ms1,psms}
R1=LFQ_Astral_DDA_15min_50ng_Condition_A_REP1
HYE=$(for c in A B; do for r in 1 2 3; do echo -n "LFQ_Astral_DDA_15min_50ng_Condition_${c}_REP$r "; done; done)
STORES=$(for r in $HYE; do echo -n "$W/store/$r.vspec.arrow "; done)
cp -n "$W/raw_local/$R1.raw" "$B/" 2>/dev/null || true
# 1. ingest (one-shot, 4 threads) and 2. MS1 features + hypotheses (4 threads) for one Astral run
T m6_ingest_${R1}_t4 "$V" -m vulcan.io.convert "$B/$R1.raw" "$B/store/" --threads 4
T m6_ms1_${R1}_t4 "$V" -m vulcan.ms1 "$W/store/$R1.vspec.arrow" "$B/ms1" --threads 4
# 3. GPU search of the 6 HYE runs (one process: setup + 6 searches)
T m6_search_d9hye "$V" -m vulcan.gpu.run --index "$W/index/hye_prod" --stores $STORES --out "$B/psms" --params '{"report_psms": 5}'
# 4. L3 features + frozen static_gbm_pred (staged env + pred site-packages; prediction cache warm from the data run)
T m6_score_d9hye env PYTHONPATH="$ROOT:$W/pred_site" /tmp/vulcan_env_$(id -u)/pe/bin/python "$ROOT/scripts/m6/score.py" --psms "$B/psms" --index "$W/index/hye_prod"
# 5. pooled paired TDC
T m6_ids_d9hye "$V" "$ROOT/scripts/m6/ids.py" --psms "$B/psms" --index "$W/index/hye_prod" --out "$B/ids" --runs $HYE
# 6. quant: link + align (+ MBR) without evaluation
T m6_lfq_link_align_d9hye "$V" "$ROOT/scripts/m6/lfq_d9.py" --ids "$B/ids" --out "$B/lfq_nombr" --no-eval
T m6_lfq_link_align_mbr_d9hye "$V" "$ROOT/scripts/m6/lfq_d9.py" --ids "$B/ids" --out "$B/lfq_mbr" --mbr --no-eval
# 7. TMT reporter quant (D4, 12 fractions; extraction + linkage + purity + correction)
T m6_tmt_quant_d4 "$V" "$ROOT/scripts/m6/tmt_d4.py" quant
