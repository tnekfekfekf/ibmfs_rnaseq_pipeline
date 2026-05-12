#!/bin/bash
# Manuscript-validated featureCounts (per-sample, sequential)
#   Validated 2026-05-11: Pearson r=0.9998 vs manuscript CPM for AA-PRO/KEW/HMH/PJH
# Options: -T 8 -p -s 2 -t exon -g gene_id  (no -B, no -M, no --primary)
# GTF: gencode.v44.annotation.no_rRNA.gtf

set -euo pipefail

FEATURECOUNTS=${FEATURECOUNTS:-/opt/miniconda3/bin/featureCounts}
GTF=${GTF:-/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf}
BAM_DIR=${BAM_DIR:-/Users/jaeeunyoo/Desktop/star_workdir/local_bams}
OUT_DIR=${OUT_DIR:-/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3}
SAMPLES=${SAMPLES:-"AA-RNA-FA AA-RNA-DKC AA-RNA-FA2 AA-RNA-FA3 AA-PRO AA-KEW AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18 AA-HMH AA-PJH"}

mkdir -p "$OUT_DIR"
SAMPLES_ARR=($SAMPLES)
TOTAL=${#SAMPLES_ARR[@]}

echo "=========================================================="
echo "featureCounts (v3 validated) — $TOTAL samples sequential"
echo "Started: $(date)"
echo "Options: -T 8 -p -s 2 -t exon -g gene_id"
echo "GTF:    $GTF"
echo "BAM:    $BAM_DIR"
echo "Output: $OUT_DIR"
echo "=========================================================="

START_ALL=$(date +%s)
for i in "${!SAMPLES_ARR[@]}"; do
  s=${SAMPLES_ARR[$i]}
  n=$((i+1))
  out="$OUT_DIR/${s}.counts.txt"
  bam="$BAM_DIR/${s}_sorted.bam"
  [ -f "$bam" ] || bam="$BAM_DIR/${s}.bam"
  [ -f "$bam" ] || { echo "[$(date +%T)] SKIP: $s — BAM not found"; continue; }

  if [ -s "$out" ]; then
    echo "[$(date +%T)] ($n/$TOTAL) $s — already done, skip"
    continue
  fi

  echo "[$(date +%T)] ($n/$TOTAL) START: $s"
  START=$(date +%s)
  "$FEATURECOUNTS" -T 8 -p -s 2 -t exon -g gene_id \
    -a "$GTF" -o "$out" "$bam" > /dev/null 2>&1
  ELAPSED=$(($(date +%s) - START))

  if [ -f "${out}.summary" ]; then
    assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "${out}.summary")
    total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "${out}.summary")
    pct=$(awk -v a=$assigned -v t=$total 'BEGIN{if(t>0) printf "%.1f", a/t*100; else print "?"}')
    echo "[$(date +%T)] DONE:  $s — ${ELAPSED}s, ${assigned} assigned (${pct}%)"
  fi
done

ELAPSED_ALL=$(($(date +%s) - START_ALL))
echo "=========================================================="
echo "ALL DONE in $((ELAPSED_ALL/60))min $((ELAPSED_ALL%60))sec"
echo "Output files in: $OUT_DIR"
echo "=========================================================="
