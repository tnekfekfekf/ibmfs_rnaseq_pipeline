#!/usr/bin/env bash
# Direct FASTQ download from EBI ENA mirror using filereport API for canonical URLs.
set -euo pipefail

OUT="/Volumes/ExtremeSSD/ibmfs/revision_analysis/fastq"
mkdir -p "$OUT"

declare -A MAP=(
  [SRR11414221]=Child1
  [SRR11414225]=Child2
  [SRR11414227]=Child3
)

for srr in SRR11414221 SRR11414225 SRR11414227; do
  alias_name="${MAP[$srr]}"
  echo "[$(date +%T)] === $srr -> $alias_name ==="
  urls=$(curl -sL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$srr&result=read_run&fields=fastq_ftp" | tail -1 | awk -F'\t' '{print $2}')
  IFS=';' read -r u1 u2 <<< "$urls"
  for pair in "1:$u1" "2:$u2"; do
    r="${pair%%:*}"
    u="${pair#*:}"
    src="https://${u}"
    dst="$OUT/${alias_name}_R${r}.fastq.gz"
    if [ -s "$dst" ]; then
      echo "[skip] $dst exists ($(du -h "$dst"|cut -f1))"
      continue
    fi
    echo "[$(date +%T)] curl  $src  ->  $dst"
    curl -fL --retry 5 --retry-delay 5 -C - -o "$dst" "$src"
    echo "[$(date +%T)] size: $(du -h "$dst" | cut -f1)"
  done
done

echo "[$(date +%T)] DONE."
ls -la "$OUT"/Child*.fastq.gz
