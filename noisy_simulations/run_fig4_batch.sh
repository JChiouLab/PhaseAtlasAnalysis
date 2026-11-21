#!/bin/sh
# Fig. 4 batch: 50 noisy replicates (A = 0.001, rep 0-49) and the noise-free control (A = 0, rep 0),
# all with the driver's defaults (mode massfix, T = 4000 s, dt = 0.0025 s).
# Fig. 4C is a separate long run (rep 100, --T 14400); see README.md.
# Usage: sh run_fig4_batch.sh [output folder] [parallel jobs]
# Each run takes about 20-30 min on one core; finished runs are skipped, so the batch can be restarted.
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-"$HERE/runs"}
JOBS=${2:-4}
mkdir -p "$OUT"
{ echo "0 0"; for r in $(seq 0 49); do echo "0.001 $r"; done; } | \
  xargs -P "$JOBS" -n 2 sh -c 'python3 "$0" --out "$1" --A "$2" --rep "$3"; echo "A$2 rep$3 exit $?"' \
  "$HERE/noise_driver.py" "$OUT"
# Classify each run into oscillatory / monopolar epochs (writes <OUT>/fig4_classification.{json,csv}):
python3 "$HERE/analyze_runs.py" "$OUT" "$OUT/fig4"
