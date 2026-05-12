#!/bin/bash
# Convert sorted BAM files to paired FASTQ for SRA submission
# Compatible with samtools 1.3.1

set -e

WORKDIR="/Volumes/ExtremeSSD/ibmfs_RNA_raw_data"
cd "$WORKDIR"

# Get all sorted BAM files (excluding mini files)
BAM_FILES=$(ls *_sorted.bam 2>/dev/null)

total=$(echo "$BAM_FILES" | wc -w | tr -d ' ')
count=0

echo "=========================================="
echo "BAM to FASTQ Conversion for SRA Submission"
echo "Found $total BAM files to process"
echo "=========================================="

for bam in $BAM_FILES; do
    count=$((count + 1))
    base="${bam%_sorted.bam}"
    
    echo ""
    echo "[$count/$total] Processing: $bam"
    echo "  Output: ${base}_R1.fastq.gz, ${base}_R2.fastq.gz"
    
    # Check if output already exists
    if [[ -f "${base}_R1.fastq.gz" && -f "${base}_R2.fastq.gz" ]]; then
        echo "  SKIPPED: Output files already exist"
        continue
    fi
    
    # Sort by read name and convert to paired FASTQ
    echo "  Step 1: Sorting by read name..."
    samtools sort -n -@ 4 -o "${base}_namesorted.bam" "$bam"
    
    echo "  Step 2: Converting to FASTQ..."
    # Use bedtools for conversion (more compatible) or samtools without threading
    samtools fastq \
        -1 "${base}_R1.fastq" \
        -2 "${base}_R2.fastq" \
        -0 /dev/null \
        -s /dev/null \
        -n \
        "${base}_namesorted.bam"
    
    echo "  Step 3: Compressing FASTQ files..."
    gzip "${base}_R1.fastq" &
    gzip "${base}_R2.fastq" &
    wait
    
    # Remove temporary name-sorted BAM
    rm -f "${base}_namesorted.bam"
    
    echo "  DONE: $(ls -lh ${base}_R1.fastq.gz ${base}_R2.fastq.gz | awk '{print $5}' | paste -sd ', ')"
done

echo ""
echo "=========================================="
echo "Conversion complete!"
echo "=========================================="
