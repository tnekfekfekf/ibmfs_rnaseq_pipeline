#!/bin/bash
# v4 hypothesis test: add -M --primary for the 10 AA-RNA patient samples
# (only those that had Pearson < 0.999 with manuscript CPM)
# Options: -T 8 -p -M --primary -s 2 -t exon -g gene_id
# All else same as v3 (gencode v44 no_rRNA GTF, Macrogen BAMs)

set -euo pipefail
FEATURECOUNTS=/opt/miniconda3/bin/featureCounts
GTF=/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf
BAM_DIR=/Users/jaeeunyoo/Desktop/star_workdir/local_bams
OUT_DIR=/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v4_Mprimary
mkdir -p "$OUT_DIR"

# Only the 10 patient samples (4 controls already match v3 perfectly)
SAMPLES=(AA-RNA-FA AA-RNA-DKC AA-RNA-FA2 AA-RNA-FA3 AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18)
TOTAL=${#SAMPLES[@]}

echo "================================================================="
echo "v4 hypothesis: featureCounts with -M --primary on $TOTAL patients"
echo "Started: $(date)"
echo "Options: -T 8 -p -M --primary -s 2 -t exon -g gene_id"
echo "================================================================="

START_ALL=$(date +%s)
for i in "${!SAMPLES[@]}"; do
  s=${SAMPLES[$i]}
  n=$((i+1))
  out="$OUT_DIR/${s}.counts.txt"
  bam="$BAM_DIR/${s}_sorted.bam"

  if [ -s "$out" ]; then
    echo "[$(date +%T)] ($n/$TOTAL) $s — already done, skip"
    continue
  fi

  echo "[$(date +%T)] ($n/$TOTAL) START: $s"
  START=$(date +%s)
  "$FEATURECOUNTS" -T 8 -p -M --primary -s 2 -t exon -g gene_id \
    -a "$GTF" -o "$out" "$bam" > /dev/null 2>&1
  ELAPSED=$(($(date +%s) - START))

  if [ -f "${out}.summary" ]; then
    assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "${out}.summary")
    total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "${out}.summary")
    pct=$(awk -v a=$assigned -v t=$total 'BEGIN{if(t>0) printf "%.1f", a/t*100; else print "?"}')
    # Compare to v3 assignment rate
    v3_sum="/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3/${s}.counts.txt.summary"
    if [ -f "$v3_sum" ]; then
      v3_assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "$v3_sum")
      v3_total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "$v3_sum")
      v3_pct=$(awk -v a=$v3_assigned -v t=$v3_total 'BEGIN{printf "%.1f", a/t*100}')
      delta=$(awk -v new=$pct -v old=$v3_pct 'BEGIN{printf "%+.1f", new-old}')
      echo "[$(date +%T)] DONE:  $s — ${ELAPSED}s, ${assigned} (${pct}%) | v3: ${v3_pct}% | Δ=${delta}pp"
    else
      echo "[$(date +%T)] DONE:  $s — ${ELAPSED}s, ${assigned} (${pct}%)"
    fi
  fi
done

ELAPSED_ALL=$(($(date +%s) - START_ALL))
echo "================================================================="
echo "ALL DONE in $((ELAPSED_ALL/60))min $((ELAPSED_ALL%60))sec"
echo "================================================================="
