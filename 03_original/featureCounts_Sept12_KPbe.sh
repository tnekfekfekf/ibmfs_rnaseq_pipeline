#!/bin/bash

# Run gene quantification for remaining 6 external SSD AA samples
# Using -M flag to include multi-mapping reads for better assignment rates

set -e

# Configuration
BAM_DIR="/Volumes/ExtremeSSD/RNAseq/alignment_ucsc"
GTF_FILE="/Users/jaeeunyoo/Desktop/researches/rna_analysis_desktop/reference_data/annotations/gencode.v44.annotation.no_rRNA.gtf"
OUTPUT_DIR="/Users/jaeeunyoo/Desktop/researches/rna_analysis_desktop/processed_data/quantification"
THREADS=8

# Remaining samples to quantify
REMAINING_SAMPLES=("AA-HMH" "AA-KEW" "AA-LES" "AA-NMS" "AA-PJH" "AA-PRO")

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "🧬 GENE QUANTIFICATION FOR REMAINING SAMPLES"
echo "=" * 60
echo "Samples: ${REMAINING_SAMPLES[*]}"
echo "BAM Directory: $BAM_DIR"
echo "GTF File: $GTF_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

# Check if BAM files exist
echo "🔍 Checking BAM files..."
for sample in "${REMAINING_SAMPLES[@]}"; do
    bam_file="${BAM_DIR}/${sample}_sorted.bam"
    if [[ -f "$bam_file" ]]; then
        echo "  ✅ $sample: $bam_file"
    else
        echo "  ❌ $sample: $bam_file (NOT FOUND)"
        exit 1
    fi
done

echo ""
echo "🚀 Starting quantification with -M flag (include multi-mapping reads)..."
echo ""

# Run featureCounts on all remaining samples
featureCounts \
    -a "$GTF_FILE" \
    -o "$OUTPUT_DIR/gene_counts_remaining_samples.txt" \
    -T "$THREADS" \
    -p \
    -B \
    -M \
    -s 2 \
    --primary \
    -g gene_id \
    -t exon \
    -F GTF \
    "${BAM_DIR}/AA-HMH_sorted.bam" \
    "${BAM_DIR}/AA-KEW_sorted.bam" \
    "${BAM_DIR}/AA-LES_sorted.bam" \
    "${BAM_DIR}/AA-NMS_sorted.bam" \
    "${BAM_DIR}/AA-PJH_sorted.bam" \
    "${BAM_DIR}/AA-PRO_sorted.bam" \
    > "$OUTPUT_DIR/quantification_remaining_samples.log" 2>&1

echo "✅ Quantification completed!"
echo ""

# Calculate assignment percentages
echo "📊 CALCULATING ASSIGNMENT PERCENTAGES..."
echo "=" * 50

if [[ -f "$OUTPUT_DIR/gene_counts_remaining_samples.txt.summary" ]]; then
    python3 -c "
import pandas as pd
import os

# Read the summary file
df = pd.read_csv('$OUTPUT_DIR/gene_counts_remaining_samples.txt.summary', sep='\t')

print('📊 REMAINING SAMPLES QUANTIFICATION RESULTS:')
print('-' * 60)
print(f\"{'Sample':<15} {'Assigned Reads':<15} {'Total Reads':<15} {'Assignment Rate':<15}\")
print('-' * 60)

for col in df.columns[1:]:
    sample_name = os.path.basename(col).replace('_sorted.bam', '')
    assigned_count = df[df['Status'] == 'Assigned'][col].iloc[0]
    total_count = df[col].sum()
    percentage = (assigned_count / total_count) * 100 if total_count > 0 else 0
    print(f\"{sample_name:<15} {assigned_count:>12,} {total_count:>12,} {percentage:>12.1f}%\")

# Calculate average
assigned_total = df[df['Status'] == 'Assigned'].iloc[0, 1:].sum()
total_total = df.iloc[:, 1:].sum().sum()
avg_percentage = (assigned_total / total_total) * 100 if total_total > 0 else 0

print('-' * 60)
print(f\"{'AVERAGE':<15} {assigned_total:>12,} {total_total:>12,} {avg_percentage:>12.1f}%\")
print()
print('🎯 COMPARISON WITH PREVIOUS RESULTS:')
print('• AA-CSB (with -M flag): 35.6%')
print('• Local AA-RNA samples: 41.5% average')
print('• Expected: 60-80%+')
"
else
    echo "❌ Summary file not found. Check log for errors:"
    echo "   $OUTPUT_DIR/quantification_remaining_samples.log"
fi

echo ""
echo "📁 Output files:"
echo "• Gene counts: $OUTPUT_DIR/gene_counts_remaining_samples.txt"
echo "• Summary: $OUTPUT_DIR/gene_counts_remaining_samples.txt.summary"
echo "• Log: $OUTPUT_DIR/quantification_remaining_samples.log"
echo ""
echo "✅ Quantification of remaining samples completed!"
