#!/bin/bash
# Add NEW samples to manuscript count matrix (Option 1 strategy)
#
# 1. Quantify new samples with v3 pipeline (consistent featureCounts options)
# 2. Append columns to manuscript_count_matrix_19samples.txt
# 3. Re-run DESeq2 with combined matrix
#
# Usage:
#   bash 06_add_new_samples.sh <BAM_DIR> <SAMPLE1> [SAMPLE2 ...]
#
# Inputs (per sample SAMPLE_NAME):
#   - BAM_DIR/SAMPLE_NAME_sorted.bam OR BAM_DIR/SAMPLE_NAME.bam
#
# Outputs:
#   - /Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/added_samples/SAMPLE_NAME.counts.txt
#   - /Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/combined_matrix_with_new.txt

set -euo pipefail

FEATURECOUNTS=${FEATURECOUNTS:-/opt/miniconda3/bin/featureCounts}
GTF=${GTF:-/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf}
MANUSCRIPT_DIR=/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS
NEW_DIR="$MANUSCRIPT_DIR/added_samples"
mkdir -p "$NEW_DIR"

if [ $# -lt 2 ]; then
  echo "Usage: $0 <BAM_DIR> <SAMPLE_NAME> [SAMPLE_NAME ...]"
  echo "Example: $0 /path/to/bams sample01 sample02"
  exit 1
fi

BAM_DIR="$1"
shift
SAMPLES=("$@")

echo "================================================="
echo "Add new samples to manuscript matrix (Option 1)"
echo "================================================="
echo "BAM directory: $BAM_DIR"
echo "Samples: ${SAMPLES[@]}"
echo "GTF:    $GTF"
echo "Output: $NEW_DIR"
echo ""

# Step 1: Quantify each new sample with v3 options (consistent with what we use going forward)
for s in "${SAMPLES[@]}"; do
  bam="$BAM_DIR/${s}_sorted.bam"
  [ -f "$bam" ] || bam="$BAM_DIR/${s}.bam"
  [ -f "$bam" ] || { echo "MISSING BAM: $s"; continue; }

  out="$NEW_DIR/${s}.counts.txt"
  if [ -s "$out" ]; then echo "Skip $s (already done)"; continue; fi

  echo "[$(date +%T)] Quantify: $s"
  "$FEATURECOUNTS" -T 8 -p -s 2 -t exon -g gene_id \
    -a "$GTF" -o "$out" "$bam" > /dev/null 2>&1

  # Print assignment rate
  if [ -f "${out}.summary" ]; then
    assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "${out}.summary")
    total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "${out}.summary")
    pct=$(awk -v a=$assigned -v t=$total 'BEGIN{printf "%.1f", a/t*100}')
    echo "[$(date +%T)] DONE: $s — $assigned reads (${pct}%)"
  fi
done

# Step 2: Append to manuscript matrix
echo ""
echo "Merging new samples into manuscript matrix..."
Rscript --vanilla -e "
suppressPackageStartupMessages({library(dplyr); library(readr)})

# Load manuscript matrix
ms <- read_tsv('$MANUSCRIPT_DIR/manuscript_count_matrix_19samples.txt', show_col_types=FALSE)
cat(sprintf('Manuscript matrix: %d genes × %d cols\n', nrow(ms), ncol(ms)))

# Add new sample columns
samps <- c('${SAMPLES[@]/ /\',\'}')
for (s in samps) {
  f <- paste0('$NEW_DIR/', s, '.counts.txt')
  if (!file.exists(f)) next
  d <- read_tsv(f, comment='#', show_col_types=FALSE)
  # Match by EnsemblID
  d_matched <- d[[7]][match(ms\$EnsemblID, d\$Geneid)]
  ms[[s]] <- d_matched
}

write_tsv(ms, '$MANUSCRIPT_DIR/combined_matrix_with_new.txt')
cat(sprintf('Output: %d genes × %d cols (with %d new samples)\n', nrow(ms), ncol(ms), length(samps)))
"

echo ""
echo "================================================="
echo "DONE. Combined matrix: $MANUSCRIPT_DIR/combined_matrix_with_new.txt"
echo ""
echo "Next step: re-run DESeq2 with this matrix"
echo "  Rscript ../PIPELINE/05_use_manuscript_matrix.R"
echo "  (or modify to use combined_matrix_with_new.txt)"
echo "================================================="
