#!/bin/bash
# featureCounts per-sample (Sept 12 KPbe.sh options)
# Manuscript options: -T 8 -p -B -M -s 2 --primary -t exon -g gene_id -F GTF
# GTF: gencode.v44.annotation.no_rRNA.gtf

FEATURECOUNTS=/opt/miniconda3/bin/featureCounts
GTF=/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf
BAM_DIR=/Users/jaeeunyoo/Desktop/star_workdir/local_bams
OUT_DIR=/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample
mkdir -p "$OUT_DIR"

SAMPLES=(AA-RNA-FA AA-RNA-DKC AA-RNA-FA2 AA-RNA-FA3 AA-PRO AA-KEW \
         AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18 AA-HMH AA-PJH)

TOTAL=${#SAMPLES[@]}
echo "===== featureCounts PER-SAMPLE ($TOTAL samples) ====="
echo "Started: $(date)"
echo "Options: -T 8 -p -B -M -s 2 --primary -t exon -g gene_id -F GTF"
echo ""

START_ALL=$(date +%s)
for i in "${!SAMPLES[@]}"; do
  s=${SAMPLES[$i]}
  n=$((i+1))
  out="$OUT_DIR/${s}.counts.txt"

  if [ -s "$out" ]; then
    echo "[$(date +%T)] ($n/$TOTAL) $s — already done, skipping"
    continue
  fi

  echo "[$(date +%T)] ($n/$TOTAL) START: $s"
  START=$(date +%s)
  "$FEATURECOUNTS" \
    -T 8 -p -B -M -s 2 --primary \
    -g gene_id -t exon -F GTF \
    -a "$GTF" \
    -o "$out" \
    "$BAM_DIR/${s}_sorted.bam" 2>&1 | tail -3
  ELAPSED=$(($(date +%s) - START))

  # Extract assignment rate from summary
  if [ -f "${out}.summary" ]; then
    assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "${out}.summary")
    total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "${out}.summary")
    pct=$(awk -v a=$assigned -v t=$total 'BEGIN{if(t>0) printf "%.1f", a/t*100; else print "?"}')
    echo "[$(date +%T)] DONE: $s — ${ELAPSED}s, $assigned reads assigned (${pct}%)"
  else
    echo "[$(date +%T)] FAIL: $s — summary missing"
  fi
  echo ""
done

ELAPSED_ALL=$(($(date +%s) - START_ALL))
echo "===== ALL DONE in $((ELAPSED_ALL/60))분 $((ELAPSED_ALL%60))초 ====="
ls -lh "$OUT_DIR/"
