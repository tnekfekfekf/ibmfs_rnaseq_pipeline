#!/usr/bin/env bash
set -uo pipefail
SCRIPT=/Volumes/ExtremeSSD/ibmfs/revision_analysis/scripts/run_pipeline.sh
LOG=/Volumes/ExtremeSSD/ibmfs/revision_analysis/logs/internal_remaining.log
# All 13 remaining internal samples (AA-KEW already done)
INTERNAL=(AA-PRO AA-HMH AA-PJH AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18 AA-RNA-DKC AA-RNA-FA AA-RNA-FA2 AA-RNA-FA3)
echo "[$(date)] start ${#INTERNAL[@]} samples" >> "$LOG"
for s in "${INTERNAL[@]}"; do
  echo "[$(date)] >>> $s" >> "$LOG"
  bash "$SCRIPT" sample "$s" 2>&1 | tee -a "$LOG"
  echo "[$(date)] <<< $s" >> "$LOG"
done
echo "[$(date)] DONE all internal" >> "$LOG"
