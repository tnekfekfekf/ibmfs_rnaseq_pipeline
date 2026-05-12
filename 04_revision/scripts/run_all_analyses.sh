#!/usr/bin/env bash
# Run all R analyses end-to-end once the count matrix is ready.
# Each step is idempotent — re-running picks up where it left off.
set -euo pipefail
ROOT=/Volumes/ExtremeSSD/ibmfs/revision_analysis
LOG="$ROOT/logs/analyses.log"
SCRIPTS="$ROOT/scripts"

run() {
  local label="$1"; shift
  echo ""
  echo "================================================================"
  echo "[$(date)]  $label"
  echo "================================================================"
  "$@" 2>&1 | tee -a "$LOG"
}

cd "$ROOT"
: > "$LOG"

# Fail early if count matrix isn't there
if [ ! -s "$ROOT/counts/featureCounts.cleaned.txt" ]; then
  echo "ERROR: $ROOT/counts/featureCounts.cleaned.txt not found. Run pipeline first."
  exit 1
fi

run "01 build_dds (load counts + metadata + GTF biotype)"     Rscript "$SCRIPTS/01_build_dds.R"
run "02 qc_pca (PCA, ComBat-seq, hierarchical clustering)"    Rscript "$SCRIPTS/02_qc_pca.R"
run "03 main_contrasts (~ cohort + group_combined)"           Rscript "$SCRIPTS/03_main_contrasts.R"
run "04 sensitivity (internal-only / public-only / drop-Child1 / FA-only / DKC-only)" \
                                                              Rscript "$SCRIPTS/04_sensitivity.R"
run "05 lncRNA_filter (4 strategies, robust intersect)"       Rscript "$SCRIPTS/05_lncRNA_filter.R"
run "06 qc_table (per-sample QC summary)"                     Rscript "$SCRIPTS/06_qc_table.R"
run "07 deconvolution (xCell + MCP + EPIC)"                   Rscript "$SCRIPTS/07_deconvolution.R"
run "08 gsea (multi-threshold)"                               Rscript "$SCRIPTS/08_gsea.R"

echo ""
echo "[$(date)] ALL ANALYSES DONE.  Outputs in:"
echo "  $ROOT/deseq2/        DESeq2 results, sensitivity, lncRNA filters, GSEA"
echo "  $ROOT/figures/       PCA, volcanoes, clustering heatmap"
echo "  $ROOT/qc/            QC summary table"
echo "  $ROOT/deconv/        immunedeconv outputs"
