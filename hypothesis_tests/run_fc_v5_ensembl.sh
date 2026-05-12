#!/bin/bash
# v5 hypothesis: featureCounts with Ensembl GRCh38.110 GTF (instead of Gencode v44 no_rRNA)
# Options: -T 8 -p -s 2 -t exon -g gene_id  (same as v3, just different GTF)
# GTF: Homo_sapiens.GRCh38.110.ucsc_chr.gtf (chr-prefixed for Macrogen BAMs)
# Sample: 10 patient samples (4 controls already perfect-match with v3)

set -euo pipefail
FEATURECOUNTS=/opt/miniconda3/bin/featureCounts
GTF=/Users/jaeeunyoo/Desktop/star_workdir/Homo_sapiens.GRCh38.110.ucsc_chr.gtf
BAM_DIR=/Users/jaeeunyoo/Desktop/star_workdir/local_bams
OUT_DIR=/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v5_ensembl
mkdir -p "$OUT_DIR"

SAMPLES=(AA-RNA-FA AA-RNA-DKC AA-RNA-FA2 AA-RNA-FA3 AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18)
TOTAL=${#SAMPLES[@]}

echo "================================================================="
echo "v5 hypothesis: featureCounts with ENSEMBL GRCh38.110 GTF"
echo "Started: $(date)"
echo "Options: -T 8 -p -s 2 -t exon -g gene_id"
echo "GTF: Homo_sapiens.GRCh38.110.ucsc_chr.gtf"
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
  "$FEATURECOUNTS" -T 8 -p -s 2 -t exon -g gene_id \
    -a "$GTF" -o "$out" "$bam" > /dev/null 2>&1
  ELAPSED=$(($(date +%s) - START))

  if [ -f "${out}.summary" ]; then
    assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "${out}.summary")
    total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "${out}.summary")
    pct=$(awk -v a=$assigned -v t=$total 'BEGIN{if(t>0) printf "%.1f", a/t*100; else print "?"}')
    v3_sum="/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3/${s}.counts.txt.summary"
    if [ -f "$v3_sum" ]; then
      v3_assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "$v3_sum")
      v3_total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "$v3_sum")
      v3_pct=$(awk -v a=$v3_assigned -v t=$v3_total 'BEGIN{printf "%.1f", a/t*100}')
      delta=$(awk -v new=$pct -v old=$v3_pct 'BEGIN{printf "%+.1f", new-old}')
      echo "[$(date +%T)] DONE:  $s — ${ELAPSED}s, ${assigned} (${pct}%) | v3 Gencode: ${v3_pct}% | Δ=${delta}pp"
    fi
  fi
done

ELAPSED_ALL=$(($(date +%s) - START_ALL))
echo "================================================================="
echo "ALL DONE in $((ELAPSED_ALL/60))min $((ELAPSED_ALL%60))sec"
echo "================================================================="
