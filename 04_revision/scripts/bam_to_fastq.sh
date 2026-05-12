#!/usr/bin/env bash
# Convert existing Child*_sorted.bam back to paired FASTQ.gz for our pipeline.
# These BAMs are coordinate-sorted; samtools fastq needs name-sorted, so we sort first.
set -euo pipefail

SRC=/Volumes/ExtremeSSD/ibmfs/ibmfs_RNA_raw_data
OUT=/Volumes/ExtremeSSD/ibmfs/revision_analysis/fastq
LOG=/Volumes/ExtremeSSD/ibmfs/revision_analysis/logs
TMP=/tmp/bam2fq
mkdir -p "$OUT" "$LOG" "$TMP"

SAMTOOLS=/opt/miniconda3/bin/samtools

for s in Child1 Child2 Child3; do
  bam="$SRC/${s}_sorted.bam"
  out1="$OUT/${s}_R1.fastq.gz"
  out2="$OUT/${s}_R2.fastq.gz"
  if [ -s "$out1" ] && [ -s "$out2" ]; then
    echo "[skip] $s already done ($(du -h "$out1" | cut -f1) / $(du -h "$out2" | cut -f1))"
    continue
  fi
  nsbam="$TMP/${s}.namesorted.bam"
  if [ -s "$nsbam" ]; then
    echo "[$(date +%T)] === $s : namesort already done ($(du -h "$nsbam" | cut -f1)) ==="
  else
    echo "[$(date +%T)] === $s : namesort ==="
    "$SAMTOOLS" sort -n -@ 8 -o "$nsbam" "$bam"
  fi
  echo "[$(date +%T)] === $s : fastq extract ==="
  # samtools 1.3.1 doesn't have -@ for fastq subcommand
  "$SAMTOOLS" fastq -1 "$TMP/${s}_R1.fq" -2 "$TMP/${s}_R2.fq" -0 /dev/null -s /dev/null -n "$nsbam"
  echo "[$(date +%T)] === $s : gzip ==="
  gzip -c "$TMP/${s}_R1.fq" > "$out1" &
  gzip -c "$TMP/${s}_R2.fq" > "$out2" &
  wait
  rm -f "$nsbam" "$TMP/${s}_R1.fq" "$TMP/${s}_R2.fq"
  echo "[$(date +%T)] === $s done : $(du -h "$out1" | cut -f1) / $(du -h "$out2" | cut -f1) ==="
done
echo "[$(date +%T)] ALL BAM2FASTQ DONE"
ls -la "$OUT"/Child*.fastq.gz
