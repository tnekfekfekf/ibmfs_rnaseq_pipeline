#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/Users/jaeeunyoo/Desktop/researches/ibmfs_fastq_raw_data"
FASTQ_DIR="${FASTQ_DIR:-$ROOT_DIR}"

REF_FASTA="${REF_FASTA:-$ROOT_DIR/reference/gencode/GRCh38.primary_assembly.genome.fa}"
REF_GTF="${REF_GTF:-$ROOT_DIR/reference/gencode/gencode.v45.annotation.gtf}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-$ROOT_DIR/reference/star_index}"

THREADS="${THREADS:-8}"
TRIMMOMATIC_BIN="${TRIMMOMATIC_BIN:-/opt/miniconda3/bin/trimmomatic}"
FASTQC_BIN="${FASTQC_BIN:-/opt/homebrew/bin/fastqc}"
STAR_BIN="${STAR_BIN:-/opt/homebrew/bin/star}"
FEATURECOUNTS_BIN="${FEATURECOUNTS_BIN:-/opt/miniconda3/bin/featureCounts}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-/opt/homebrew/bin/samtools}"

OUT_QC_RAW="$ROOT_DIR/qc/raw"
OUT_QC_TRIM="$ROOT_DIR/qc/trimmed"
OUT_TRIM="$ROOT_DIR/trimmed"
OUT_ALIGN="$ROOT_DIR/aligned"
OUT_COUNTS="$ROOT_DIR/counts"
OUT_LOGS="$ROOT_DIR/logs"

SAMPLES_DEFAULT=("AA-RNA-1" "AA-RNA-4" "AA-RNA-5")

usage() {
  cat <<'EOF'
Usage:
  rnaseq_pipeline.sh [--config /path/to/env] [--threads N] [--samples "S1,S2,S3"] <command>

Commands:
  qc_raw        FastQC on raw FASTQs
  trim          Trimmomatic PE trimming (+ FastQC on trimmed)
  align         STAR alignment (sorted BAM)
  count         featureCounts gene-level counts
  all           qc_raw + trim + align + count

Notes:
  Samples assume files like <SAMPLE>_R1.fastq.gz and <SAMPLE>_R2.fastq.gz in FASTQ_DIR.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

need() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required executable: $1"
}

ensure_dirs() {
  mkdir -p "$OUT_QC_RAW" "$OUT_QC_TRIM" "$OUT_TRIM" "$OUT_ALIGN" "$OUT_COUNTS" "$OUT_LOGS"
}

detect_adapters() {
  if [ "${ADAPTERS_FA:-}" != "" ] && [ -f "$ADAPTERS_FA" ]; then
    echo "$ADAPTERS_FA"
    return
  fi
  for p in \
    "/opt/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa" \
    "/opt/miniconda3/share/trimmomatic/adapters/TruSeq3-PE-2.fa" \
    "/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa" \
    "/usr/share/trimmomatic/adapters/TruSeq3-PE.fa"
  do
    if [ -f "$p" ]; then
      echo "$p"
      return
    fi
  done
  die "Could not find adapters FASTA for Trimmomatic. Set ADAPTERS_FA=/path/to/adapters.fa"
}

samples_from_arg() {
  local s="${1:-}"
  if [ "$s" = "" ]; then
    printf "%s\n" "${SAMPLES_DEFAULT[@]}"
    return
  fi
  IFS=',' read -r -a arr <<<"$s"
  printf "%s\n" "${arr[@]}"
}

fastq_paths() {
  local sample="$1"
  local r1="$FASTQ_DIR/${sample}_R1.fastq.gz"
  local r2="$FASTQ_DIR/${sample}_R2.fastq.gz"
  [ -f "$r1" ] || die "Missing: $r1"
  [ -f "$r2" ] || die "Missing: $r2"
  echo "$r1|$r2"
}

trimmed_paths() {
  local sample="$1"
  local r1p="$OUT_TRIM/${sample}_R1.trimmed.paired.fastq.gz"
  local r2p="$OUT_TRIM/${sample}_R2.trimmed.paired.fastq.gz"
  echo "$r1p|$r2p"
}

cmd_qc_raw() {
  ensure_dirs
  need "$FASTQC_BIN"
  while read -r sample; do
    local pair; pair="$(fastq_paths "$sample")"
    local r1="${pair%%|*}"
    local r2="${pair##*|}"
    "$FASTQC_BIN" -t "$THREADS" -o "$OUT_QC_RAW" "$r1" "$r2" |& tee "$OUT_LOGS/fastqc_raw.${sample}.log"
  done < <(samples_from_arg "$SAMPLES_ARG")
}

cmd_trim() {
  ensure_dirs
  need "$TRIMMOMATIC_BIN"
  need "$FASTQC_BIN"
  local adapters; adapters="$(detect_adapters)"
  while read -r sample; do
    local pair; pair="$(fastq_paths "$sample")"
    local r1="${pair%%|*}"
    local r2="${pair##*|}"

    local r1p="$OUT_TRIM/${sample}_R1.trimmed.paired.fastq.gz"
    local r1u="$OUT_TRIM/${sample}_R1.trimmed.unpaired.fastq.gz"
    local r2p="$OUT_TRIM/${sample}_R2.trimmed.paired.fastq.gz"
    local r2u="$OUT_TRIM/${sample}_R2.trimmed.unpaired.fastq.gz"

    "$TRIMMOMATIC_BIN" PE -threads "$THREADS" -phred33 \
      "$r1" "$r2" \
      "$r1p" "$r1u" \
      "$r2p" "$r2u" \
      ILLUMINACLIP:"$adapters":2:30:10 \
      SLIDINGWINDOW:4:20 \
      MINLEN:36 |& tee "$OUT_LOGS/trimmomatic.${sample}.log"

    "$FASTQC_BIN" -t "$THREADS" -o "$OUT_QC_TRIM" "$r1p" "$r2p" |& tee "$OUT_LOGS/fastqc_trimmed.${sample}.log"
  done < <(samples_from_arg "$SAMPLES_ARG")
}

cmd_align() {
  ensure_dirs
  need "$STAR_BIN"
  need "$SAMTOOLS_BIN"
  [ -d "$STAR_INDEX_DIR" ] || die "Missing STAR index dir: $STAR_INDEX_DIR"
  [ -f "$REF_GTF" ] || die "Missing GTF: $REF_GTF"

  while read -r sample; do
    local tp; tp="$(trimmed_paths "$sample")"
    local r1p="${tp%%|*}"
    local r2p="${tp##*|}"
    [ -f "$r1p" ] || die "Missing trimmed FASTQ: $r1p (run trim first)"
    [ -f "$r2p" ] || die "Missing trimmed FASTQ: $r2p (run trim first)"

    local prefix="$OUT_ALIGN/${sample}."
    "$STAR_BIN" \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX_DIR" \
      --readFilesIn "$r1p" "$r2p" \
      --readFilesCommand "gunzip -c" \
      --sjdbGTFfile "$REF_GTF" \
      --outFileNamePrefix "$prefix" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts |& tee "$OUT_LOGS/star.${sample}.log"

    local bam="${prefix}Aligned.sortedByCoord.out.bam"
    [ -f "$bam" ] || die "STAR did not produce BAM: $bam"
    "$SAMTOOLS_BIN" index -@ "$THREADS" "$bam"
    "$SAMTOOLS_BIN" flagstat -@ "$THREADS" "$bam" >"$OUT_ALIGN/${sample}.flagstat.txt"
  done < <(samples_from_arg "$SAMPLES_ARG")
}

cmd_count() {
  ensure_dirs
  need "$FEATURECOUNTS_BIN"
  [ -f "$REF_GTF" ] || die "Missing GTF: $REF_GTF"

  local bam_list=()
  while read -r sample; do
    local bam="$OUT_ALIGN/${sample}.Aligned.sortedByCoord.out.bam"
    [ -f "$bam" ] || die "Missing BAM: $bam (run align first)"
    bam_list+=("$bam")
  done < <(samples_from_arg "$SAMPLES_ARG")

  "$FEATURECOUNTS_BIN" \
    -T "$THREADS" \
    -p -B -C \
    -a "$REF_GTF" \
    -o "$OUT_COUNTS/featureCounts.txt" \
    "${bam_list[@]}" |& tee "$OUT_LOGS/featureCounts.log"

  awk 'BEGIN{FS=OFS="\t"} NR==1{next} NR==2{for(i=7;i<=NF;i++){gsub(/.*\\//,"",$i); gsub(/\\.Aligned\\.sortedByCoord\\.out\\.bam/,"",$i)} print; next} {print}' \
    "$OUT_COUNTS/featureCounts.txt" >"$OUT_COUNTS/featureCounts.cleaned.txt"
}

CONFIG_PATH=""
SAMPLES_ARG=""
THREADS_ARG=""

if [ "${1:-}" = "" ]; then
  usage
  exit 1
fi

while [ $# -gt 0 ]; do
  case "$1" in
    --config)
      CONFIG_PATH="${2:-}"; shift 2;;
    --samples)
      SAMPLES_ARG="${2:-}"; shift 2;;
    --threads)
      THREADS_ARG="${2:-}"; shift 2;;
    -h|--help)
      usage; exit 0;;
    qc_raw|trim|align|count|all)
      COMMAND="$1"; shift; break;;
    *)
      die "Unknown arg: $1";;
  esac
done

if [ "$CONFIG_PATH" != "" ]; then
  [ -f "$CONFIG_PATH" ] || die "Config not found: $CONFIG_PATH"
  set -a
  # shellcheck disable=SC1090
  . "$CONFIG_PATH"
  set +a
fi

if [ "$THREADS_ARG" != "" ]; then
  THREADS="$THREADS_ARG"
fi

case "${COMMAND:-}" in
  qc_raw) cmd_qc_raw;;
  trim) cmd_trim;;
  align) cmd_align;;
  count) cmd_count;;
  all) cmd_qc_raw; cmd_trim; cmd_align; cmd_count;;
  *) usage; exit 1;;
esac

