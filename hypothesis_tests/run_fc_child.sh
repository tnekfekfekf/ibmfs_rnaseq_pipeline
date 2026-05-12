#!/bin/bash
# Run featureCounts on Child1/2/3 with v3 options (same as 14 manuscript samples)
# GTF: Ensembl-chr version (chr prefix stripped to match Child BAM headers)

set -euo pipefail
FEATURECOUNTS=/opt/miniconda3/bin/featureCounts
GTF=/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.ensembl_chr.gtf
BAM_DIR=/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/aligned
OUT_DIR=/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3
mkdir -p "$OUT_DIR"

SAMPLES=(Child1 Child2 Child3)
TOTAL=${#SAMPLES[@]}

echo "=========================================================="
echo "featureCounts on $TOTAL Child samples (v3 options + Ensembl-chr GTF)"
echo "Started: $(date)"
echo "=========================================================="

START_ALL=$(date +%s)
for i in "${!SAMPLES[@]}"; do
  s=${SAMPLES[$i]}
  n=$((i+1))
  out="$OUT_DIR/${s}.counts.txt"

  if [ -s "$out" ]; then
    echo "[$(date +%T)] ($n/$TOTAL) $s — already done, skip"
    continue
  fi

  echo "[$(date +%T)] ($n/$TOTAL) START: $s"
  START=$(date +%s)
  "$FEATURECOUNTS" -T 8 -p -s 2 -t exon -g gene_id \
    -a "$GTF" -o "$out" "$BAM_DIR/${s}.bam" > /dev/null 2>&1
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
echo "=========================================================="
