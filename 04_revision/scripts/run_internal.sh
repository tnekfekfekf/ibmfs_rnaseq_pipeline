#!/usr/bin/env bash
# Run trim+align on the 14 internal samples sequentially.
set -euo pipefail
SCRIPT=/Volumes/ExtremeSSD/ibmfs/revision_analysis/scripts/run_pipeline.sh
LOG=/Volumes/ExtremeSSD/ibmfs/revision_analysis/logs/internal_pipeline.log
INTERNAL=(AA-KEW AA-PRO AA-HMH AA-PJH AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18 AA-RNA-DKC AA-RNA-FA AA-RNA-FA2 AA-RNA-FA3)
echo "[$(date)] internal pipeline start (${#INTERNAL[@]} samples)" >>"$LOG"
for s in "${INTERNAL[@]}"; do
  echo "[$(date)] >>> $s" >>"$LOG"
  "$SCRIPT" sample "$s" 2>&1 | tee -a "$LOG" | grep -E "^\[(trim|align|count)\]" || true
  echo "[$(date)] <<< $s" >>"$LOG"
done
echo "[$(date)] internal pipeline DONE" >>"$LOG"
