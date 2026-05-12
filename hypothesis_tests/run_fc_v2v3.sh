#!/bin/bash
# Run v2 and v3 featureCounts per-sample sequentially
#   v2: -T 8 -p -B -C -s 2 -t exon -g gene_id    (Oct 13-15 pattern)
#   v3: -T 8 -p    -s 2 -t exon -g gene_id       (Aug 16 minimal)

FEATURECOUNTS=/opt/miniconda3/bin/featureCounts
GTF=/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf
BAM_DIR=/Users/jaeeunyoo/Desktop/star_workdir/local_bams

SAMPLES=(AA-RNA-FA AA-RNA-DKC AA-RNA-FA2 AA-RNA-FA3 AA-PRO AA-KEW \
         AA-RNA-1 AA-RNA-4 AA-RNA-5 AA-RNA-13 AA-RNA-16 AA-RNA-18 AA-HMH AA-PJH)

run_version() {
  local VERSION="$1"; shift
  local OUT_DIR="$1"; shift
  local OPTS=("$@")
  mkdir -p "$OUT_DIR"

  echo ""
  echo "================================================"
  echo "===== $VERSION ====="
  echo "Options: ${OPTS[@]}"
  echo "================================================"

  local TOTAL=${#SAMPLES[@]}
  for i in "${!SAMPLES[@]}"; do
    local s=${SAMPLES[$i]}
    local n=$((i+1))
    local out="$OUT_DIR/${s}.counts.txt"
    if [ -s "$out" ]; then
      echo "[$(date +%T)] [$VERSION $n/$TOTAL] $s — already done, skip"
      continue
    fi

    echo "[$(date +%T)] [$VERSION $n/$TOTAL] START: $s"
    local START=$(date +%s)
    "$FEATURECOUNTS" "${OPTS[@]}" -a "$GTF" -o "$out" \
      "$BAM_DIR/${s}_sorted.bam" > /dev/null 2>&1
    local ELAPSED=$(($(date +%s) - START))

    if [ -f "${out}.summary" ]; then
      local assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "${out}.summary")
      local total=$(awk -F'\t' 'NR>1{s+=$2} END{print s}' "${out}.summary")
      local pct=$(awk -v a=$assigned -v t=$total 'BEGIN{if(t>0) printf "%.1f", a/t*100; else print "?"}')
      echo "[$(date +%T)] [$VERSION] DONE: $s — ${ELAPSED}s, ${assigned} assigned (${pct}%)"
    else
      echo "[$(date +%T)] [$VERSION] FAIL: $s"
    fi
  done
}

echo "Started: $(date)"

run_version "v2" /Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v2 \
  -T 8 -p -B -C -s 2 -t exon -g gene_id

run_version "v3" /Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3 \
  -T 8 -p -s 2 -t exon -g gene_id

echo ""
echo "===== ALL DONE at $(date) ====="
