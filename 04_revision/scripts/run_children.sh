#!/usr/bin/env bash
set -uo pipefail
SCRIPT=/Volumes/ExtremeSSD/ibmfs/revision_analysis/scripts/run_pipeline.sh
LOG=/Volumes/ExtremeSSD/ibmfs/revision_analysis/logs/children_pipeline.log
echo "[$(date)] start Child1/2/3" >> "$LOG"
for s in Child1 Child2 Child3; do
  echo "[$(date)] >>> $s" >> "$LOG"
  bash "$SCRIPT" sample "$s" 2>&1 | tee -a "$LOG"
  echo "[$(date)] <<< $s" >> "$LOG"
done
echo "[$(date)] DONE all Child" >> "$LOG"
