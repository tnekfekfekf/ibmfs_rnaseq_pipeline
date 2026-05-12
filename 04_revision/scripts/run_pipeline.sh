#!/usr/bin/env bash
# Unified RNA-seq pipeline (HISAT2 edition — STAR was unstable on 16GB Mac).
# Per-sample: Trimmomatic -> HISAT2 -> samtools sort/index/flagstat
# Cohort:    featureCounts (-s 2 for TruSeq Stranded Total RNA + Ribo-Zero)

set -euo pipefail

ROOT="/Volumes/ExtremeSSD/ibmfs/revision_analysis"
REF_GTF="/Volumes/ExtremeSSD/ibmfs/revision_analysis/reference/gencode.v45.ensembl_chr.gtf"  # Ensembl chr names match HISAT2 grch38_tran BAM headers
HISAT2_INDEX="/Users/jaeeunyoo/hisat2_index_local/genome_tran"

FASTQ_DIR="$ROOT/fastq"
TRIM_DIR="$ROOT/trimmed"
ALIGN_DIR="$ROOT/aligned"
COUNTS_DIR="$ROOT/counts"
QC_TRIM="$ROOT/qc/trimmed"
LOGS="$ROOT/logs"

THREADS="${THREADS:-8}"

TRIMMOMATIC=/opt/miniconda3/bin/trimmomatic
FASTQC=/opt/miniconda3/bin/fastqc
HISAT2=/opt/miniconda3/bin/hisat2
SAMTOOLS=/opt/homebrew/bin/samtools  # 1.22, much faster sort than miniconda's 1.3.1
FEATURECOUNTS=/opt/miniconda3/bin/featureCounts
ADAPTERS=/opt/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa

mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$COUNTS_DIR" "$QC_TRIM" "$LOGS"

SAMPLES=(
  AA-KEW AA-PRO
  AA-HMH AA-PJH
  AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18
  AA-RNA-DKC AA-RNA-FA AA-RNA-FA2 AA-RNA-FA3
  Child1 Child2 Child3
)

stage_inputs() {
  local SRC=/Volumes/ExtremeSSD/ibmfs/ibmfs_fastq_raw_data
  for s in AA-KEW AA-PRO AA-HMH AA-PJH \
           AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18 \
           AA-RNA-DKC AA-RNA-FA AA-RNA-FA2 AA-RNA-FA3; do
    for r in R1 R2; do
      [ -e "$FASTQ_DIR/${s}_${r}.fastq.gz" ] || ln -sf "$SRC/${s}_${r}.fastq.gz" "$FASTQ_DIR/${s}_${r}.fastq.gz"
    done
  done
  echo "[stage] internal FASTQs symlinked into $FASTQ_DIR"
}

trim_one() {
  local s="$1"
  local r1="$FASTQ_DIR/${s}_R1.fastq.gz"
  local r2="$FASTQ_DIR/${s}_R2.fastq.gz"
  local p1="$TRIM_DIR/${s}_R1.paired.fq.gz"
  local p2="$TRIM_DIR/${s}_R2.paired.fq.gz"
  local u1="$TRIM_DIR/${s}_R1.unpaired.fq.gz"
  local u2="$TRIM_DIR/${s}_R2.unpaired.fq.gz"
  if [ ! -f "$r1" ]; then echo "[trim] $s: missing R1, skipping"; return; fi
  if [ -s "$p1" ] && [ -s "$p2" ]; then
    echo "[trim] $s: already trimmed, skipping"
    return
  fi
  echo "[trim] $s start $(date +%T)"
  "$TRIMMOMATIC" PE -threads "$THREADS" -phred33 \
    "$r1" "$r2" "$p1" "$u1" "$p2" "$u2" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
    2>&1 | tee "$LOGS/trimmomatic.${s}.log" | tail -3
  "$FASTQC" -t "$THREADS" -o "$QC_TRIM" "$p1" "$p2" 2>&1 | tee "$LOGS/fastqc_trim.${s}.log" | tail -1
  echo "[trim] $s done $(date +%T)"
}

align_one() {
  local s="$1"
  local p1="$TRIM_DIR/${s}_R1.paired.fq.gz"
  local p2="$TRIM_DIR/${s}_R2.paired.fq.gz"
  local sam="$ALIGN_DIR/${s}.sam"
  local bam="$ALIGN_DIR/${s}.bam"
  if [ -s "$bam" ]; then echo "[align] $s: BAM present, skip"; return; fi
  if [ ! -s "$p1" ] || [ ! -s "$p2" ]; then echo "[align] $s: trim missing, skip"; return; fi

  echo "[align] $s start $(date +%T)"
  # Pipe HISAT2 -> samtools sort directly to avoid huge intermediate SAM file.
  # --rna-strandness RF for TruSeq Stranded (R2 is sense). --dta for downstream tools.
  # samtools sort temp on /tmp (POSIX, FIFO/locking-friendly)
  local sort_tmp="/tmp/sort_${s}_$$"
  rm -rf "$sort_tmp"; mkdir -p "$sort_tmp"
  set -o pipefail
  "$HISAT2" -p "$THREADS" -x "$HISAT2_INDEX" \
    -1 "$p1" -2 "$p2" \
    --dta --rna-strandness RF --no-unal \
    --summary-file "$ALIGN_DIR/${s}.hisat2_summary.txt" 2> "$LOGS/hisat2.${s}.log" \
    | "$SAMTOOLS" sort -@ 4 -m 512M -T "$sort_tmp/x" -O bam -o "$bam" -
  set +o pipefail
  rm -rf "$sort_tmp"
  "$SAMTOOLS" index -@ 4 "$bam"
  "$SAMTOOLS" flagstat -@ 4 "$bam" > "$ALIGN_DIR/${s}.flagstat.txt"
  cat "$LOGS/hisat2.${s}.log"   # echo summary to main log
  echo "[align] $s done $(date +%T)"
}

count_all() {
  local bams=()
  for s in "${SAMPLES[@]}"; do
    local b="$ALIGN_DIR/${s}.bam"
    if [ -s "$b" ]; then bams+=("$b"); else echo "[count] WARN missing $b"; fi
  done
  echo "[count] featureCounts on ${#bams[@]} BAMs (-s 2 reverse-stranded)"
  # -p : paired; -C : exclude chimeric (mates on diff chr).
  # NO -B : allow singletons (HISAT2 --no-unal removes the unmapped mate, otherwise many "Unassigned_Singleton")
  # -s 2 : reverse-stranded TruSeq Stranded
  # -t gene (not exon): captures pre-mRNA / intronic reads typical of total-RNA
  # Ribo-Zero libraries. Internal samples have ~55% reads in introns vs exons.
  "$FEATURECOUNTS" -T "$THREADS" \
    -p \
    -C -s 2 \
    -a "$REF_GTF" -t gene -g gene_id \
    -o "$COUNTS_DIR/featureCounts.txt" \
    "${bams[@]}" 2>&1 | tee "$LOGS/featureCounts.log"
  awk 'BEGIN{FS=OFS="\t"} /^#/{print; next} NR==2 || $1=="Geneid"{
        for(i=7;i<=NF;i++){
          n=$i; gsub(/^.*\//,"",n); gsub(/\.bam$/,"",n); $i=n
        } print; next
      } {print}' "$COUNTS_DIR/featureCounts.txt" > "$COUNTS_DIR/featureCounts.cleaned.txt"
  echo "[count] done"
}

usage() { echo "Usage: $0 {stage|trim|align|count|all|sample SAMPLE_ID}"; }

case "${1:-}" in
  stage)  stage_inputs ;;
  trim)   for s in "${SAMPLES[@]}"; do trim_one "$s"; done ;;
  align)  for s in "${SAMPLES[@]}"; do align_one "$s"; done ;;
  count)  count_all ;;
  sample)
    shift; s="$1"
    trim_one "$s"; align_one "$s" ;;
  all)
    stage_inputs
    for s in "${SAMPLES[@]}"; do trim_one "$s"; align_one "$s"; done
    count_all ;;
  *) usage; exit 1 ;;
esac
