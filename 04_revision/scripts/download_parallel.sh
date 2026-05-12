#!/usr/bin/env bash
# Parallel ENA download — 6 curl streams concurrently to maximize throughput.
set -uo pipefail
OUT="/Volumes/ExtremeSSD/ibmfs/revision_analysis/fastq"
LOG_DIR="/Volumes/ExtremeSSD/ibmfs/revision_analysis/logs"
mkdir -p "$OUT" "$LOG_DIR"

# Hard-coded URLs (already verified via filereport)
declare -A URLS=(
  [Child1_R1]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/021/SRR11414221/SRR11414221_1.fastq.gz"
  [Child1_R2]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/021/SRR11414221/SRR11414221_2.fastq.gz"
  [Child2_R1]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/025/SRR11414225/SRR11414225_1.fastq.gz"
  [Child2_R2]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/025/SRR11414225/SRR11414225_2.fastq.gz"
  [Child3_R1]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/027/SRR11414227/SRR11414227_1.fastq.gz"
  [Child3_R2]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/027/SRR11414227/SRR11414227_2.fastq.gz"
)

dl_one() {
  local key="$1"
  local url="${URLS[$key]}"
  local dst="$OUT/${key}.fastq.gz"
  local lg="$LOG_DIR/dl_${key}.log"
  if [ -s "$dst" ]; then
    echo "[skip] $key ($(du -h "$dst" | cut -f1))"
    return
  fi
  echo "[start $(date +%T)] $key"
  curl -fL --retry 10 --retry-delay 10 --retry-all-errors \
    -C - -o "$dst" "$url" \
    --silent --show-error \
    -w "[done $(date +%T)] $key  size=%{size_download}  time=%{time_total}s  speed=%{speed_download}\n" \
    >"$lg" 2>&1
  local rc=$?
  if [ $rc -eq 0 ]; then
    echo "[ok $(date +%T)] $key  $(du -h "$dst" | cut -f1)"
  else
    echo "[FAIL $(date +%T)] $key  curl rc=$rc  see $lg"
  fi
}

export -f dl_one
export OUT LOG_DIR
# shellcheck disable=SC2034
declare -p URLS >/tmp/urls_export
# can't easily export assoc arrays — call sequentially via & but in same shell

pids=()
for key in Child1_R1 Child1_R2 Child2_R1 Child2_R2 Child3_R1 Child3_R2; do
  dl_one "$key" &
  pids+=($!)
done

echo "[$(date +%T)] launched ${#pids[@]} parallel downloads, pids: ${pids[*]}"
wait "${pids[@]}"
echo "[$(date +%T)] ALL DOWNLOADS DONE"
ls -la "$OUT"/Child*.fastq.gz
