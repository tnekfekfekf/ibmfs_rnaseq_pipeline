#!/bin/bash

# ==============================================================================
# FLEXIBLE RNA-seq ANALYSIS PIPELINE - MASTER SCRIPT
# ==============================================================================
# Updated: 2025-08-05
# Description: Complete RNA-seq pipeline with flexible execution modes
# Features: Individual sample processing OR group DEG analysis OR both
# NEW: Added child sample processing and comprehensive QC
# ==============================================================================

# Configuration
SCRIPT_VERSION="5.0-CLEAN-ORGANIZED-HYBRID"
# Location of external data workspace (raw_data, reference, tools)
# Prefer environment variable if provided to support code/data separation
WORKDIR="${RNASEQ_DATA_DIR:-/Volumes/Extreme SSD/RNAseq}"
# Location of this script repository (internal Mac storage). By default, resolve from this script's path
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Results will be saved on Mac (same directory as scripts)
RESULTS_DIR="$SCRIPTS_DIR/results"
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

# Performance optimizations
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
MEMORY_LIMIT="8G"
PARALLEL_JOBS=2

# Create results directory on Mac
mkdir -p "$RESULTS_DIR"

# Navigate to working directory for raw data access
cd "$WORKDIR" || { echo "❌ Cannot access $WORKDIR"; exit 1; }

# Pipeline header
echo ""
echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║              RNA-seq Analysis Pipeline v${SCRIPT_VERSION}              ║"
echo "║           🔥 CLEAN & ORGANIZED + CHILD SAMPLE SUPPORT 🔥                   ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Performance monitoring
monitor_resources() {
    echo "📊 System Resources:"
    echo "   CPU Cores: $THREADS"
    echo "   Memory: $(free -h | grep Mem | awk '{print $2}')"
    echo "   Available: $(free -h | grep Mem | awk '{print $7}')"
    echo "   Disk Space: $(df -h . | tail -1 | awk '{print $4}') available"
    echo ""
}

# Progress tracking
PROGRESS_FILE="$RESULTS_DIR/logs/pipeline_progress_${TIMESTAMP}.log"
log_progress() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$PROGRESS_FILE"
}

# Enhanced error handling with recovery
handle_error_with_recovery() {
    local error_msg="$1"
    local recovery_cmd="$2"
    
    log_message "ERROR: $error_msg"
    echo "🔄 Attempting recovery..."
    
    if [ -n "$recovery_cmd" ]; then
        eval "$recovery_cmd" && {
            log_message "Recovery successful"
            return 0
        } || {
            log_message "Recovery failed"
            exit 1
        }
    else
        exit 1
    fi
}

# ==============================================================================
# EXECUTION MODE DETECTION
# ==============================================================================

# Check execution mode BEFORE setting SAMPLE_NAME
CHILD_SUBSET=""
if [ "$1" = "--deg-analysis" ] || [ "$1" = "-d" ]; then
    MODE="deg_analysis"
# lncRNA mode removed
elif [ "$1" = "--copy-child-data" ] || [ "$1" = "-cc" ]; then
    MODE="copy_child_data"
elif [ "$1" = "--child-processing" ] || [ "$1" = "-cp" ]; then
    MODE="child_processing"
elif [ "$1" = "--child-subset" ] || [ "$1" = "-cs" ]; then
    MODE="child_processing_subset"
    CHILD_SUBSET="$2"
elif [ "$1" = "--child-qc" ] || [ "$1" = "-cq" ]; then
    MODE="child_qc"
elif [ "$1" = "--comprehensive-qc" ] || [ "$1" = "-qc" ]; then
    MODE="comprehensive_qc"
elif [ "$1" = "--deg-aa-vs-child" ] || [ "$1" = "-dac" ]; then
    MODE="deg_aa_vs_child"
elif [ "$1" = "--clean-results" ] || [ "$1" = "-c" ]; then
    MODE="clean_results"
elif [ -z "$1" ]; then
    MODE="help"
else
    MODE="single_sample"
    SAMPLE_NAME="$1"
fi

case $MODE in
    "help")
        echo "🎯 USAGE OPTIONS:"
        echo "================="
        echo ""
        echo "🗂️ Code/Data separation:"
        echo "   export RNASEQ_DATA_DIR=\"/Volumes/Extreme SSD/RNAseq\""
        echo "   # Scripts run from Mac; raw data read from external SSD; results saved on Mac"
        echo ""
        echo "📊 Individual Sample Processing:"
        echo "   $0 SAMPLE_NAME"
        echo "   Example: $0 AA-RNA-1"
        echo "   Note: HD and AA samples are already finished (BAM files exist)"
        echo ""
        echo "🧬 mRNA DEG Analysis (protein-coding genes with top 100 gene lists):"
        echo "   $0 --deg-analysis    (or -d)"
        echo "   → Runs complete mRNA DEG analysis + automatic interpretation"
        echo "   → Generates top 100 up/down genes for pathway analysis"
        echo ""
        echo "🧬 DEG: General AA vs CHILD (featureCounts + DESeq2):"
        echo "   $0 --deg-aa-vs-child   (or -dac)"
        echo "   → Runs featureCounts on AA + Child BAMs, then DESeq2 (~ batch + group if batch available)"
        echo "   → Outputs: tables (all/significant/top100) and plots (PCA/Volcano/Heatmap) in analysis_results_child_vs_aa/"
        echo ""
        echo "👶 Child Sample Processing:"
        echo "   $0 --child-processing (or -cp)"
        echo "   → Processes Child2, Child3, Child1 samples (fresh alignment + quantification)"
        echo "   → HD and AA samples are already finished (BAM files exist)"
        echo "   → RECOMMENDED: Uses already-trimmed FASTQ files + fresh alignment"
        echo "   → Input: Child2_1_paired.fastq.gz, Child2_2_paired.fastq.gz, etc."
        echo "   → Order: Child2 → Child3 → Child1"
        echo ""
        echo "📋 Copy Existing Child Data:"
        echo "   $0 --copy-child-data    (or -cc)"
        echo "   → Copies existing Child2/3 data from external SSD to Mac"
        echo "   → Completes BAM conversion and quantification"
        echo "   → NOTE: Existing SAM files may be incomplete - use with caution"
        echo ""
        echo "👶 Child Sample Processing (subset):"
        echo "   $0 --child-subset \"Child2,Child3\"   (or -cs \"Child2,Child3\")"
        echo "   → Processes only the specified child samples (comma or space separated)"
        echo ""
        echo "🔍 Child Sample QC:"
        echo "   $0 --child-qc        (or -cq)"
        echo "   → Runs comprehensive rRNA QC on child samples"
        echo ""
        echo "🔍 Comprehensive QC (All samples):"
        echo "   $0 --comprehensive-qc (or -qc)"
        echo "   → Runs comprehensive rRNA QC on all samples (HD + AA + Child)"
        echo ""
        echo "🧹 Clean Results:"
        echo "   $0 --clean-results   (or -c)"
        echo "   → Cleans and organizes all results"
        echo ""
        echo "📋 Examples:"
        echo "   $0 --deg-analysis     # Run mRNA DEG analysis"
        echo "   $0 --deg-aa-vs-child  # Run AA vs Child DEG (featureCounts + DESeq2)"
        echo "   $0 --child-processing # Process child samples"
        echo "   $0 --child-qc         # QC child samples"
        echo "   $0 --comprehensive-qc # QC all samples"
        echo ""
        exit 0
        ;;
    "deg_analysis")
        echo "🧬 MODE: Complete mRNA DEG Analysis + Results Interpretation"
        echo "==========================================================="
        echo "⚡ Running mRNA DEG analysis with top 100 up/down gene lists"
        echo "🔬 Ready for pathway analysis (Metascape/Enrichr) and biological report"
        echo "📂 Results will be saved in: analysis_results_mrna/"
        echo ""
        
        # Check if quantification files exist
        EXPECTED_FILES=20
        EXISTING_FILES=$(ls "$RESULTS_DIR"/*/quantification/*_abundance.tab 2>/dev/null | wc -l)
        
        if [ $EXISTING_FILES -lt $EXPECTED_FILES ]; then
            echo "⚠️  Warning: Only found $EXISTING_FILES of $EXPECTED_FILES expected quantification files"
            echo "📁 Available files:"
            ls "$RESULTS_DIR"/*/quantification/*_abundance.tab 2>/dev/null | head -10
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ DEG analysis cancelled"
                exit 1
            fi
        else
            echo "✅ Found $EXISTING_FILES quantification files - ready for analysis!"
        fi
        
        # Set up conda environment
        source ~/.bash_profile
        conda activate rnaseq
        
        # Run mRNA DEG analysis
        echo "🚀 Starting complete mRNA DEG analysis..."
        Rscript "$SCRIPTS_DIR/three_comparisons_mrna_only.R"
        
        if [ $? -eq 0 ]; then
            echo "🎉 Complete mRNA DEG analysis finished successfully!"
            echo "📂 Results in: analysis_results_mrna/"
            echo "🧬 Enrichment files created (6 files per comparison):"
            echo "   • top100_upregulated_genes_only.txt"
            echo "   • top100_downregulated_genes_only.txt"
            echo "   • top100_upregulated_for_enrichment.csv (with p-values)"
            echo "   • top100_downregulated_for_enrichment.csv (with p-values)"
            echo ""
            echo "🌍 NEXT: Upload *_genes_only.txt to https://maayanlab.cloud/Enrichr/"
            echo ""
            echo "🎯 mRNA ANALYSIS COMPLETE - Ready for publication!"
        else
            echo "❌ mRNA analysis failed"
            exit 1
        fi
        
        exit 0
        ;;
    # lncRNA analysis case removed
    "deg_aa_vs_child")
        echo "🧬 MODE: DEG (AA_GENERAL vs CHILD_HEALTHY) via featureCounts + DESeq2"
        echo "========================================================================"
        echo "⚙️  Building counts with featureCounts from existing BAMs, then running DESeq2"
        echo "📂 Output: analysis_results_child_vs_aa/"
        echo ""

        # Activate environment
        source ~/.bash_profile
        conda activate rnaseq || { echo "❌ Cannot activate rnaseq environment"; exit 1; }

        # Check required tools
        which featureCounts >/dev/null 2>&1 || { echo "❌ featureCounts not found (Subread). Install 'subread' in rnaseq env."; exit 1; }
        which Rscript >/dev/null 2>&1 || { echo "❌ Rscript not found."; exit 1; }

        # Paths
        GTF="$WORKDIR/reference/annotation/Homo_sapiens.GRCh38.110.gtf"
        OUTDIR="$WORKDIR/analysis_results_child_vs_aa"
        mkdir -p "$OUTDIR/tables" "$OUTDIR/plots"

        # Build/Load sample table
        SAMPLE_TABLE="$WORKDIR/results/deg_input/sample_table.csv"
        if [ -f "$SAMPLE_TABLE" ]; then
            echo "📄 Using provided sample table: $SAMPLE_TABLE"
        else
            echo "🔎 No sample_table.csv found – auto-discovering AA + Child BAMs under results/ ..."
            mkdir -p "$WORKDIR/results/deg_input"
            TMP_ST="$SAMPLE_TABLE"
            : > "$TMP_ST"
            echo "sample,group,batch,bam" >> "$TMP_ST"

            # Discover AA BAMs (master single-sample outputs)
            for b in "$RESULTS_DIR"/alignment/AA-RNA-*_sorted.bam; do
                [ -f "$b" ] || continue
                s=$(basename "$b" | sed 's/_sorted\.bam$//')
                echo "$s,AA_GENERAL,,$b" >> "$TMP_ST"
            done

            # Discover child BAMs
            for b in "$RESULTS_DIR"/child_samples/*/alignment/*_sorted.bam; do
                [ -f "$b" ] || continue
                s=$(basename "$(dirname "$(dirname "$b")")")
                echo "$s,CHILD_HEALTHY,,$b" >> "$TMP_ST"
            done

            echo "📋 Auto-generated sample table:"
            cat "$TMP_ST"
        fi

        # Count AA + Child BAMs
        COUNT_TXT="$OUTDIR/tables/featureCounts_gene_counts.txt"
        if [ ! -f "$GTF" ]; then
            echo "❌ GTF not found: $GTF"; exit 1;
        fi

        # Collect BAM list from sample table
        BAM_LIST=$(tail -n +2 "$SAMPLE_TABLE" | cut -d',' -f4 | tr '\n' ' ')
        if [ -z "$BAM_LIST" ]; then
            echo "❌ No BAMs found for counting."; exit 1;
        fi

        echo "🧮 Running featureCounts..."
        featureCounts -T 6 -p -s 2 -t exon -g gene_id \
            -a "$GTF" \
            -o "$COUNT_TXT" \
            $BAM_LIST
        FC_STATUS=$?
        if [ $FC_STATUS -ne 0 ]; then
            echo "❌ featureCounts failed"; exit $FC_STATUS;
        fi

        echo "✅ featureCounts completed: $COUNT_TXT"

        # Create inline R script for DESeq2
        RUN_R="$SCRIPTS_DIR/run_deg_aa_vs_child_inline.R"
        cat > "$RUN_R" <<'RSCRIPT'
        suppressPackageStartupMessages({
          library(DESeq2); library(data.table); library(ggplot2); library(pheatmap); library(ggrepel)
        })

        args <- commandArgs(trailingOnly=TRUE)
        count_txt <- args[1]
        sample_csv <- args[2]
        outdir <- args[3]

        dir.create(file.path(outdir,"tables"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(outdir,"plots"),  showWarnings = FALSE, recursive = TRUE)

        # Read featureCounts txt
        fc <- fread(count_txt, skip=1)
        # featureCounts header: Geneid Chr Start End Strand Length <samples...>
        count_mat <- as.matrix(fc[, -(1:6)])
        rownames(count_mat) <- fc[[1]]
        colnames(count_mat) <- gsub("^.*\\/", "", colnames(count_mat)) # keep filenames at first; will remap

        # Read sample table
        st <- fread(sample_csv)
        # Map sample names to BAM basenames
        st[, bam_base := basename(bam)]
        # Reorder to match columns
        match_idx <- match(colnames(count_mat), st$bam_base)
        if (anyNA(match_idx)) {
          stop("Some BAMs in counts not present in sample_table.csv")
        }
        colnames(count_mat) <- st$sample[match_idx]

        # Build colData
        coldata <- data.frame(sample = st$sample[match_idx],
                              group  = factor(st$group[match_idx], levels=c("CHILD_HEALTHY","AA_GENERAL")),
                              stringsAsFactors = FALSE)
        rownames(coldata) <- coldata$sample

        use_batch <- FALSE
        if ("batch" %in% names(st)) {
          bvec <- st$batch[match_idx]
          if (all(!is.na(bvec)) && length(unique(bvec))>1) {
            coldata$batch <- factor(bvec)
            tab <- table(coldata$batch, coldata$group)
            confounded <- all(rowSums(tab>0)==1 | colSums(tab>0)==1)
            if (!confounded) use_batch <- TRUE
          }
        }

        dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                                      colData = coldata,
                                      design = as.formula(if (use_batch) "~ batch + group" else "~ group"))
        keep <- rowSums(counts(dds) >= 10) >= 3
        dds <- dds[keep, ]
        dds <- DESeq(dds)

        res <- results(dds, contrast = c("group","AA_GENERAL","CHILD_HEALTHY"))
        res <- res[order(res$padj), ]
        res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")

        fwrite(res_dt, file.path(outdir,"tables","AA_GENERAL_vs_CHILD_HEALTHY_all_genes.csv"))

        sig_strict   <- res_dt[!is.na(padj) & padj<0.01 & abs(log2FoldChange)>1]
        sig_standard <- res_dt[!is.na(padj) & padj<0.05 & abs(log2FoldChange)>log2(1.5)]
        fwrite(sig_strict,   file.path(outdir,"tables","significant_strict_padj0.01_fc2.csv"))
        fwrite(sig_standard, file.path(outdir,"tables","significant_standard_padj0.05_fc1.5.csv"))

        up_all   <- res_dt[!is.na(padj) & log2FoldChange>0][order(padj)]
        down_all <- res_dt[!is.na(padj) & log2FoldChange<0][order(padj)]
        fwrite(head(up_all,100)[,.(gene_id,log2FoldChange,padj)],   file.path(outdir,"tables","top100_upregulated_for_enrichment.csv"))
        fwrite(head(down_all,100)[,.(gene_id,log2FoldChange,padj)], file.path(outdir,"tables","top100_downregulated_for_enrichment.csv"))
        write.table(head(up_all$gene_id,100),   file.path(outdir,"tables","top100_upregulated_genes_only.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
        write.table(head(down_all$gene_id,100), file.path(outdir,"tables","top100_downregulated_genes_only.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

        vsd <- vst(dds, blind=FALSE)
        pca_dat <- plotPCA(vsd, intgroup = if ("batch" %in% names(coldata)) c("group","batch") else "group", returnData=TRUE)
        pv <- round(100*attr(pca_dat,"percentVar"))
        p <- ggplot(pca_dat, aes(PC1,PC2,color=group)) +
             { if ("batch" %in% names(coldata)) geom_point(aes(shape=coldata$batch), size=4, alpha=.85) else geom_point(size=4, alpha=.85) } +
             ggrepel::geom_text_repel(aes(label=name), size=3, max.overlaps = 50) +
             xlab(paste0("PC1: ",pv[1],"%")) + ylab(paste0("PC2: ",pv[2],"%")) +
             ggtitle("PCA: AA_GENERAL vs CHILD_HEALTHY") + theme_minimal()
        ggsave(file.path(outdir,"plots","01_PCA.png"), p, width=8, height=6, dpi=300, bg="white")

        volc <- ggplot(res_dt[!is.na(padj)], aes(log2FoldChange, -log10(padj))) +
                geom_point(aes(color = padj<0.01 & abs(log2FoldChange)>1), alpha=.6, size=1.2) +
                scale_color_manual(values=c("grey60","red")) +
                geom_vline(xintercept=c(-1,1), linetype="dashed") +
                geom_hline(yintercept=-log10(0.01), linetype="dashed") +
                labs(title="Volcano: AA_GENERAL vs CHILD_HEALTHY", x="log2FC", y="-log10(FDR)") +
                theme_minimal()
        ggsave(file.path(outdir,"plots","02_Volcano.png"), volc, width=9, height=7, dpi=300, bg="white")

        if (nrow(sig_strict) >= 6) {
          top_ids <- head(sig_strict$gene_id, 30)
          m <- assay(vsd)[top_ids, , drop=FALSE]
          ann <- data.frame(Group = coldata$group); rownames(ann) <- rownames(coldata)
          png(file.path(outdir,"plots","03_Heatmap_top30_strict.png"), width=12, height=10, units="in", res=300, bg="white")
          pheatmap(m, scale="row", cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=ann,
                   show_rownames=TRUE, show_colnames=TRUE)
          dev.off()
        }

        sink(file.path(outdir,"AA_GENERAL_vs_CHILD_HEALTHY_summary.txt"))
        cat("DEG: AA_GENERAL vs CHILD_HEALTHY (featureCounts + DESeq2)\n")
        cat("Date:", format(Sys.time()), "\n\n")
        cat("Samples:\n"); print(with(coldata, table(group)))
        if ("batch" %in% names(coldata)) { cat("\nBatch:\n"); print(table(coldata$batch)) }
        cat("\nTotal genes tested:", nrow(res_dt), "\n")
        cat("Significant (strict FDR<0.01 & |log2FC|>1):", nrow(sig_strict), "\n")
        cat("Significant (standard FDR<0.05 & |log2FC|>log2(1.5)):", nrow(sig_standard), "\n\n")
        cat("Top 10 strict DEGs:\n")
        print(sig_strict[1:min(10, nrow(sig_strict)), .(gene_id, log2FoldChange, padj)])
        sink()
        RSCRIPT

        echo "🚀 Running DESeq2..."
        Rscript "$RUN_R" "$COUNT_TXT" "$SAMPLE_TABLE" "$OUTDIR"
        RSTAT=$?
        if [ $RSTAT -ne 0 ]; then
            echo "❌ R DEG script failed"; exit $RSTAT;
        fi

        echo "🎉 DEG complete! See: $OUTDIR"
        exit 0
        ;;
    "clean_results")
        echo "🧹 MODE: Clean Old Results"
        echo "==========================="
        echo "🔍 Identifying old DEG analysis results..."
        echo ""
        
        # List what will be removed
        echo "🗂️ Files/folders to be removed:"
        
        # Check what exists
        ITEMS_TO_REMOVE=()
        
        # mRNA results
        if [ -d "analysis_results_mrna/01_comparison_general_AA_vs_HD" ]; then
            echo "  • analysis_results_mrna/01_comparison_general_AA_vs_HD/"
            ITEMS_TO_REMOVE+=("analysis_results_mrna/01_comparison_general_AA_vs_HD")
        fi
        
        if [ -d "analysis_results_mrna/02_comparison_genetic_AA_vs_HD" ]; then
            echo "  • analysis_results_mrna/02_comparison_genetic_AA_vs_HD/"
            ITEMS_TO_REMOVE+=("analysis_results_mrna/02_comparison_genetic_AA_vs_HD")
        fi
        
        if [ -d "analysis_results_mrna/03_comparison_general_vs_genetic_AA" ]; then
            echo "  • analysis_results_mrna/03_comparison_general_vs_genetic_AA/"
            ITEMS_TO_REMOVE+=("analysis_results_mrna/03_comparison_general_vs_genetic_AA")
        fi
        
        # lncRNA results removed
        
        # Also clean old non-numbered folders if they exist
        if [ -d "analysis_results/comparison_1_general_AA_vs_HD" ]; then
            echo "  • analysis_results/comparison_1_general_AA_vs_HD/ (old format)"
            ITEMS_TO_REMOVE+=("analysis_results/comparison_1_general_AA_vs_HD")
        fi
        
        if [ -d "analysis_results/comparison_2_genetic_AA_vs_HD" ]; then
            echo "  • analysis_results/comparison_2_genetic_AA_vs_HD/ (old format)"
            ITEMS_TO_REMOVE+=("analysis_results/comparison_2_genetic_AA_vs_HD")
        fi
        
        if [ -d "analysis_results/comparison_3_general_vs_genetic_AA" ]; then
            echo "  • analysis_results/comparison_3_general_vs_genetic_AA/ (old format)"
            ITEMS_TO_REMOVE+=("analysis_results/comparison_3_general_vs_genetic_AA")
        fi
        
        if [ -d "analysis_results/deg_analysis" ]; then
            echo "  • analysis_results/deg_analysis/"
            ITEMS_TO_REMOVE+=("analysis_results/deg_analysis")
        fi
        
        if [ -f "analysis_results/Overall_Three_Comparisons_Summary.txt" ]; then
            echo "  • analysis_results/Overall_Three_Comparisons_Summary.txt"
            ITEMS_TO_REMOVE+=("analysis_results/Overall_Three_Comparisons_Summary.txt")
        fi
        
        if [ ${#ITEMS_TO_REMOVE[@]} -eq 0 ]; then
            echo "✅ No old DEG results found - nothing to clean!"
            exit 0
        fi
        
        echo ""
        echo "⚠️  🛡️ PROTECTED (will NOT be removed):"
        echo "  • Raw data (raw_data/)"
        echo "  • Results on Mac ($RESULTS_DIR/)"
        echo "  • Company analysis (company_analysis/)"
        echo "  • Scripts and references"
        echo ""
        
        # Confirmation
        read -p "🗑️ Proceed with cleanup? This cannot be undone! (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "❌ Cleanup cancelled"
            exit 0
        fi
        
        # Perform cleanup
        echo "🧹 Cleaning up old results..."
        for ITEM in "${ITEMS_TO_REMOVE[@]}"; do
            if [ -e "$ITEM" ]; then
                echo "  🗑️ Removing $ITEM"
                rm -rf "$ITEM"
            fi
        done
        
        # Create backup timestamp
        echo "📅 Cleanup completed at: $(date)" > "analysis_results/last_cleanup.log"
        
        echo "✅ Cleanup completed successfully!"
        echo "🎆 Ready for fresh DEG analysis!"
        echo ""
        echo "🚀 Run analysis now with:"
        echo "   • For mRNA: $0 --deg-analysis"
        
        exit 0
        ;;
    
    "child_processing")
        echo "🧬 MODE: Child Sample Processing (FRESH ALIGNMENT FROM TRIMMED FILES)"
        echo "===================================================================="
        echo "⚡ Processing Child2, Child3, Child1 samples (fresh alignment + quantification)"
        echo "📂 Results will be saved in: $RESULTS_DIR/child_samples/"
        echo "📝 Note: HD and AA samples are already finished (BAM files exist)"
        echo "✂️  Using already-trimmed paired FASTQ files (Child2_1_paired.fastq.gz, etc.)"
        echo "🆕 Running fresh HISAT2 alignment for complete, high-quality results"
        echo ""
        
        # Check if child sample files exist (expecting Child2, Child3, Child1 trimmed FASTQs)
        CHILD_FILES=$(find "$WORKDIR/processed_data/child_samples" -name "*_paired.fastq.gz" 2>/dev/null | wc -l)
        if [ $CHILD_FILES -lt 6 ]; then
            echo "⚠️  Warning: Not all child sample files found. Expected 6, found $CHILD_FILES"
            echo "📁 Expected files (6 total):"
            echo "   • Child2: Child2_1_paired.fastq.gz, Child2_2_paired.fastq.gz"
            echo "   • Child3: Child3_1_paired.fastq.gz, Child3_2_paired.fastq.gz"
            echo "   • Child1: Child1_1_paired.fastq.gz, Child1_2_paired.fastq.gz"
            echo "📁 Available child files:"
            find "$WORKDIR/processed_data/child_samples" -name "*_paired.fastq.gz" 2>/dev/null
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ Child processing cancelled"
                exit 1
            fi
        fi
        
        # Clean up old results first - REMOVED to prevent loss of copied files
        # echo "🧹 Cleaning up old child sample results..."
        # rm -rf "$RESULTS_DIR/child_samples/"
        
        # Activate conda environment
        echo "🐍 Activating R environment..."
        source ~/.bash_profile
        conda activate rnaseq
        
        # Run child sample processing
        echo "🚀 Starting Child Sample Processing..."
        
        # Child sample processing logic (directly in master pipeline)
        echo "👶 CHILD SAMPLE PROCESSING (INTEGRATED VERSION)"
        echo "=============================================="
        
        # Configuration for child processing
        DATA_DIR="$HOME/extreme_ssd_rnaseq"
        GENOME_INDEX="$HOME/extreme_ssd_rnaseq/reference/genome/working_index_clean/genome"
        CHILD_SAMPLES=("Child2" "Child3" "Child1")
        
        # Check if genome index exists
        if [ ! -d "$(dirname "$GENOME_INDEX")" ]; then
            echo "❌ Genome index directory not found: $(dirname "$GENOME_INDEX")"
            exit 1
        fi
        
        echo "✅ Genome index found"
        echo ""
        
        # Function to check if step is completed
        check_step_completed() {
            local sample=$1
            local step=$2
            local sample_dir="$RESULTS_DIR/child_samples/$sample"
            
            case $step in
                "fastqc")
                    [ -f "$sample_dir/alignment/${sample}_1_paired_fastqc.html" ] && \
                    [ -f "$sample_dir/alignment/${sample}_2_paired_fastqc.html" ]
                    ;;
                "alignment")
                    [ -f "$sample_dir/alignment/${sample}_sorted.bam" ] && \
                    [ -f "$sample_dir/alignment/${sample}_sorted.bam.bai" ]
                    ;;
                "quantification")
                    [ -f "$sample_dir/quantification/${sample}.gtf" ] && \
                    [ -f "$sample_dir/quantification/${sample}_abundance.tab" ]
                    ;;
                *)
                    false
                    ;;
            esac
        }
        
        # Function to process a single child sample
        process_child_sample() {
            local sample=$1
            echo "🚀 Processing $sample ..."
            
            # Create directories
            local sample_dir="$RESULTS_DIR/child_samples/$sample"
            local alignment_dir="$sample_dir/alignment"
            local quantification_dir="$sample_dir/quantification"
            
            mkdir -p "$alignment_dir" "$quantification_dir"
            
            # Find paired FASTQ files on external SSD
            local r1_file="$DATA_DIR/processed_data/child_samples/$sample/alignment/${sample}_1_paired.fastq.gz"
            local r2_file="$DATA_DIR/processed_data/child_samples/$sample/alignment/${sample}_2_paired.fastq.gz"
            
            if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
                echo "❌ Paired FASTQ files not found for $sample"
                echo "   Expected: $r1_file"
                echo "   Expected: $r2_file"
                return 1
            fi
            
            echo "    📁 Sample files found: 2"
            echo "    📊 R1 file: $(basename "$r1_file")"
            echo "    📊 R2 file: $(basename "$r2_file")"
            echo "    📊 Using trimmed files: $(basename "$r1_file") and $(basename "$r2_file")"
            
            # Step 1: Quality Control (FastQC)
            if ! check_step_completed "$sample" "fastqc"; then
                echo "  🔍 Running quality control..."
                fastqc -o "$alignment_dir" -t 4 "$r1_file" "$r2_file"
                if [ $? -eq 0 ]; then
                    echo "  ✅ FastQC completed for $sample"
                else
                    echo "  ❌ FastQC failed for $sample"
                    return 1
                fi
            else
                echo "  ✅ FastQC already completed for $sample - skipping"
            fi
            
            # Step 2: Trimming (already completed on external SSD)
            echo "  ✅ Trimming already completed for $sample - using existing paired files"
            
            # Step 3: Alignment with HISAT2
            local sam_file="$alignment_dir/${sample}.sam"
            local bam_file="$alignment_dir/${sample}_sorted.bam"
            
            if ! check_step_completed "$sample" "alignment"; then
                echo "  🧬 Running fresh alignment with HISAT2..."
                
                echo "    📍 Genome index: $GENOME_INDEX"
                echo "    📍 Output SAM: $sam_file"
                
                # Run HISAT2 with proper path handling
                echo "    📍 Running HISAT2 with paths:"
                echo "       Genome index: $GENOME_INDEX"
                echo "       R1 file: $r1_file"
                echo "       R2 file: $r2_file"
                echo "       Output: $sam_file"
                
                # Use eval to handle paths with spaces properly
                hisat2_cmd="hisat2 -x \"$GENOME_INDEX\" -1 \"$r1_file\" -2 \"$r2_file\" -S \"$sam_file\" -p 4 --dta --rna-strandness RF"
                echo "    📍 Command: $hisat2_cmd"
                
                eval "$hisat2_cmd"
                
                if [ $? -eq 0 ]; then
                    echo "    ✅ HISAT2 alignment completed for $sample"
                    
                    # Convert SAM to sorted BAM
                    echo "    🔄 Converting to sorted BAM..."
                    
                    # Remove pre-existing outputs
                    rm -f "$bam_file" "${bam_file}.bai"
                    
                    # Convert to sorted BAM
                    samtools sort -@ 4 -T "$alignment_dir/tmp_$sample" "$sam_file" -o "$bam_file"
                    
                    if [ $? -eq 0 ]; then
                        # Index BAM file
                        echo "    📍 Indexing BAM file..."
                        samtools index "$bam_file"
                        
                        if [ $? -eq 0 ]; then
                            echo "    ✅ BAM indexing completed for $sample"
                            
                            # Clean up SAM file to save space
                            echo "    🧹 Cleaning up SAM file..."
                            rm -f "$sam_file"
                            
                            echo "    ✅ $sample alignment complete!"
                        else
                            echo "    ❌ BAM indexing failed for $sample"
                            return 1
                        fi
                    else
                        echo "    ❌ SAM to BAM conversion failed for $sample"
                        return 1
                    fi
                else
                    echo "    ❌ HISAT2 alignment failed for $sample"
                    return 1
                fi
            else
                echo "  ✅ Alignment already completed for $sample - skipping"
            fi
            
            # Step 4: Quantification with StringTie
            if ! check_step_completed "$sample" "quantification"; then
                echo "  📊 Running quantification with StringTie..."
                
                local gtf_file="$quantification_dir/${sample}.gtf"
                local abundance_file="$quantification_dir/${sample}_abundance.tab"
                
                # Run StringTie
                stringtie "$bam_file" \
                          -o "$gtf_file" \
                          -p 4 \
                          -G "$DATA_DIR/reference/annotation/Homo_sapiens.GRCh38.110.gtf" \
                          -A "$abundance_file"
                
                if [ $? -eq 0 ]; then
                    echo "  ✅ Quantification completed for $sample"
                else
                    echo "  ❌ Quantification failed for $sample"
                    return 1
                fi
            else
                echo "  ✅ Quantification already completed for $sample - skipping"
            fi
            
            echo "  ✅ $sample processing complete!"
            echo ""
        }
        
        # Check child sample files
        echo "🔍 Checking child sample files..."
        for sample in "${CHILD_SAMPLES[@]}"; do
            r1_file="$DATA_DIR/processed_data/child_samples/$sample/alignment/${sample}_1_paired.fastq.gz"
            r2_file="$DATA_DIR/processed_data/child_samples/$sample/alignment/${sample}_2_paired.fastq.gz"
            
            if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
                echo "   ✅ Found: ${sample}_1_paired.fastq.gz and ${sample}_2_paired.fastq.gz"
            else
                echo "   ❌ Missing: ${sample}_1_paired.fastq.gz or ${sample}_2_paired.fastq.gz"
            fi
        done
        echo ""
        
        echo "✅ Processing files directly from external SSD"
        echo ""
        
        # Process each child sample
        for sample in "${CHILD_SAMPLES[@]}"; do
            if process_child_sample "$sample"; then
                echo "🎉 $sample processing completed successfully!"
            else
                echo "❌ $sample processing failed"
                exit 1
            fi
        done
        
        echo "🎉 Child Sample Processing finished successfully!"
        echo "📂 Results in: $RESULTS_DIR/child_samples/"
        echo "🧬 Child samples processed in order: Child2 → Child3 → Child1"
        echo ""
        echo "📊 Each folder contains:"
        echo "   • alignment/ - BAM files"
        echo "   • quantification/ - GTF and abundance files"
        echo ""
        echo "🎯 CHILD SAMPLE PROCESSING COMPLETE!"
        exit 0
        ;;

    "child_processing_subset")
        echo "🧬 MODE: Child Sample Processing (Subset)"
        echo "========================================="
        if [ -z "$CHILD_SUBSET" ]; then
            echo "❌ No child subset provided. Usage: $0 --child-subset \"Child2,Child3\""
            exit 1
        fi
        echo "⚡ Processing subset: $CHILD_SUBSET"
        echo "📂 Results will be saved in: $RESULTS_DIR/child_samples/"
        echo ""

        # Activate environment
        echo "🐍 Activating R environment..."
        source ~/.bash_profile
        conda activate rnaseq

        # Check if child sample files exist
        CHILD_FILES=$(ls raw_data/SRR114142*_1.fastq.gz 2>/dev/null | wc -l)
        if [ $CHILD_FILES -lt 2 ]; then
            echo "⚠️  Warning: Not all child sample files found. Expected 2, found $CHILD_FILES"
            echo "📁 Available child files (R1):"
            ls raw_data/SRR114142*_1.fastq.gz 2>/dev/null
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ Child subset processing cancelled"
                exit 1
            fi
        fi

        # Run child sample processing with subset passed as arg to R script
        echo "🚀 Starting Child Sample Processing (subset)..."
        Rscript "$SCRIPTS_DIR/child_sample_processing.R" "$CHILD_SUBSET"
        
        # Note: No automatic deletion of intermediates. Cleanup must be explicit.

        if [ $? -eq 0 ]; then
            echo "🎉 Child Sample Processing (subset) finished successfully!"
            echo "📂 Results in: $RESULTS_DIR/child_samples/"
            echo "🎯 CHILD SAMPLE PROCESSING (SUBSET) COMPLETE!"
        else
            echo "❌ Child sample processing (subset) failed"
            exit 1
        fi

        exit 0
        ;;

    "child_qc")
        echo "🧬 MODE: Child Sample QC"
        echo "=========================="
        echo "⚡ Running comprehensive rRNA QC on child samples"
        echo "📂 Results will be saved in: $RESULTS_DIR/child_qc/"
        echo ""
        
        # Check if child BAMs exist (expecting at least 2 BAMs: Child2 and Child3)
        CHILD_BAMS=$(ls "$RESULTS_DIR/child_samples"/*/alignment/*_sorted.bam 2>/dev/null | wc -l)
        if [ $CHILD_BAMS -lt 2 ]; then
            echo "⚠️  Warning: Not all child sample alignments found. Expected ≥2, found $CHILD_BAMS"
            echo "📁 Available child BAM files:"
            ls "$RESULTS_DIR/child_samples"/*/alignment/*_sorted.bam 2>/dev/null
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ Child QC cancelled"
                exit 1
            fi
        fi
        
        # Clean up old results first
        echo "🧹 Cleaning up old child QC results..."
        rm -rf "$RESULTS_DIR/child_qc/"
        
        # Activate conda environment
        echo "🐍 Activating R environment..."
        source ~/.bash_profile
        conda activate rnaseq
        
        # Run child QC
        echo "🚀 Starting Child QC..."
        Rscript "$SCRIPTS_DIR/child_sample_qc.R"
        
        if [ $? -eq 0 ]; then
            echo "🎉 Child QC finished successfully!"
            echo "📂 Results in: $RESULTS_DIR/child_qc/"
            echo "🧬 QC folders created:"
            echo "   • alignment/ - QC plots"
            echo "   • quantification/ - QC plots"
            echo ""
            echo "🎯 CHILD QC COMPLETE!"
        else
            echo "❌ Child QC failed"
            exit 1
        fi
        
        exit 0
        ;;

    "comprehensive_qc")
        echo "🧬 MODE: Comprehensive QC"
        echo "=========================="
        echo "⚡ Running comprehensive rRNA QC on all samples (HD + AA + Child)"
        echo "📂 Results will be saved in: $RESULTS_DIR/comprehensive_qc/"
        echo ""
        
        # Check if all sample files exist
        HD_FILES=$(ls raw_data/HD_*_1.fastq.gz 2>/dev/null | wc -l)
        AA_FILES=$(ls raw_data/AA-RNA-*_1.fastq.gz 2>/dev/null | wc -l)
        CHILD_FILES=$(ls processed_data/child_samples_combined/*_1.fastq.gz 2>/dev/null | wc -l)
        
        echo "📊 Sample availability:"
        echo "   HD samples found: $HD_FILES (expected: 10)"
        echo "   AA samples found: $AA_FILES (expected: 9+)"
        echo "   Child samples found: $CHILD_FILES (expected: 2)"
        echo ""
        
        if [ $HD_FILES -lt 10 ]; then
            echo "⚠️  Warning: Not all HD samples found. Expected 10, found $HD_FILES"
            echo "📁 Available HD files:"
            ls raw_data/HD_*_1.fastq.gz 2>/dev/null
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ Comprehensive QC cancelled"
                exit 1
            fi
        fi
        
        if [ $AA_FILES -lt 9 ]; then
            echo "⚠️  Warning: Not all AA samples found. Expected 9, found $AA_FILES"
            echo "📁 Available AA files:"
            ls raw_data/AA-RNA-*_1.fastq.gz 2>/dev/null
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ Comprehensive QC cancelled"
                exit 1
            fi
        fi
        
        if [ $CHILD_FILES -lt 2 ]; then
            echo "⚠️  Warning: Not all child samples found. Expected 2, found $CHILD_FILES"
            echo "📁 Available child files (R1):"
            ls processed_data/child_samples_combined/*_1.fastq.gz 2>/dev/null
            echo ""
            read -p "Continue with available files? (y/n): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "❌ Comprehensive QC cancelled"
                exit 1
            fi
        fi
        
        # Clean up old results first
        echo "🧹 Cleaning up old comprehensive QC results..."
        rm -rf "$RESULTS_DIR/comprehensive_qc/"
        
        # Activate conda environment
        echo "🐍 Activating R environment..."
        source ~/.bash_profile
        conda activate rnaseq
        
        # Run comprehensive QC
        echo "🚀 Starting Comprehensive QC..."
        Rscript "$SCRIPTS_DIR/comprehensive_qc.R"
        
        if [ $? -eq 0 ]; then
            echo "🎉 Comprehensive QC finished successfully!"
            echo "📂 Results in: $RESULTS_DIR/comprehensive_qc/"
            echo "🧬 QC folders created:"
            echo "   • alignment/ - QC plots"
            echo "   • quantification/ - QC plots"
            echo ""
            echo "🎯 COMPREHENSIVE QC COMPLETE!"
        else
            echo "❌ Comprehensive QC failed"
            exit 1
        fi
        
        exit 0
        ;;

    "single_sample")
        echo "📊 MODE: Single Sample Processing"
        echo "================================="
        echo "🧬 Processing sample: $SAMPLE_NAME"
        echo ""
        # Continue with normal single sample processing...
        
        # NOTE: The rest of the single sample processing code follows after esac
        ;;
    "copy_child_data")
        echo "📋 MODE: Copy Existing Child Data + Complete Processing"
        echo "======================================================="
        echo "⚡ Copying existing Child2/3 data from external SSD to Mac"
        echo "🔄 Completing BAM conversion and quantification"
        echo "📂 Results will be saved in: $RESULTS_DIR/child_samples/"
        echo ""
        
        # Check if external SSD child data exists
        if [ ! -d "$WORKDIR/processed_data/child_samples" ]; then
            echo "❌ No existing child sample data found on external SSD"
            echo "   Expected: $WORKDIR/processed_data/child_samples/"
            exit 1
        fi
        
        # Check for Child2 and Child3
        if [ ! -d "$WORKDIR/processed_data/child_samples/Child2" ] || [ ! -d "$WORKDIR/processed_data/child_samples/Child3" ]; then
            echo "❌ Child2 or Child3 directories not found on external SSD"
            echo "   Found:"
            ls "$WORKDIR/processed_data/child_samples/" 2>/dev/null
            exit 1
        fi
        
        echo "✅ Found existing child sample data on external SSD"
        echo "📁 Copying to Mac results directory..."
        
        # Clean up old results first
        rm -rf "$RESULTS_DIR/child_samples"
        
        # Copy existing data to Mac
        cp -r "$WORKDIR/processed_data/child_samples" "$RESULTS_DIR/"
        
        if [ $? -eq 0 ]; then
            echo "✅ Data copied successfully to Mac"
            echo "📂 Location: $RESULTS_DIR/child_samples/"
        else
            echo "❌ Failed to copy data"
            exit 1
        fi
        
        # Activate conda environment
        echo "🐍 Activating R environment..."
        source ~/.bash_profile
        conda activate rnaseq
        
        # Complete the processing (convert SAM to BAM, run quantification)
        echo "🚀 Completing child sample processing..."
        Rscript "$SCRIPTS_DIR/complete_child_processing.R"
        
        if [ $? -eq 0 ]; then
            echo "🎉 Child sample processing completed successfully!"
            echo "📂 Final results in: $RESULTS_DIR/child_samples/"
            echo "🧬 Child samples processed: Child2, Child3"
            echo ""
            echo "📊 Each folder contains:"
            echo "   • alignment/ - BAM files (converted from SAM)"
            echo "   • quantification/ - GTF and abundance files"
            echo ""
            echo "🎯 CHILD SAMPLE PROCESSING COMPLETE!"
        else
            echo "❌ Child sample processing failed"
            exit 1
        fi
        
        exit 0
        ;;
esac

# ==============================================================================
# SINGLE SAMPLE PROCESSING SETUP
# ==============================================================================

# Only run single sample processing if we're in single_sample mode
if [ "$MODE" != "single_sample" ]; then
    exit 0
fi

# Create log directory
mkdir -p "$RESULTS_DIR/logs"
LOGFILE="$RESULTS_DIR/logs/${SAMPLE_NAME}_pipeline_${TIMESTAMP}.log"

# Logging function
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGFILE"
}

# Error handling function
handle_error() {
    log_message "ERROR: $1"
    log_message "Pipeline FAILED for sample: $SAMPLE_NAME"
    exit 1
}

log_message "Starting RNA-seq pipeline for sample: $SAMPLE_NAME"
log_message "Pipeline version: $SCRIPT_VERSION"
log_message "Working directory: $WORKDIR"
START_TIME=$(date +%s)

# Environment Setup
log_message "Setting up conda environment..."
source ~/.bash_profile || handle_error "Cannot source bash profile"
conda activate rnaseq || handle_error "Cannot activate rnaseq conda environment"

# Verify tools are available
log_message "Verifying bioinformatics tools..."
which fastqc >/dev/null || handle_error "FastQC not found in conda environment"
which hisat2 >/dev/null || handle_error "HISAT2 not found in conda environment"
which samtools >/dev/null || handle_error "Samtools not found in conda environment"
which stringtie >/dev/null || handle_error "StringTie not found in conda environment"

# Configuration
THREADS=4
RAW_DATA="raw_data"
PROCESSED="$RESULTS_DIR"
REFERENCE="reference"

log_message "Configuration: $THREADS threads, processing $SAMPLE_NAME"

echo "🧬 Processing sample: $SAMPLE_NAME"
echo "======================================"

# Check input files
R1_FILE="$RAW_DATA/${SAMPLE_NAME}_1.fastq.gz"
R2_FILE="$RAW_DATA/${SAMPLE_NAME}_2.fastq.gz"

if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
    handle_error "FASTQ files not found for $SAMPLE_NAME. Looking for: $R1_FILE and $R2_FILE"
fi

echo "✅ Input files found"
echo "   R1: $R1_FILE"
echo "   R2: $R2_FILE"
log_message "Input files validated for $SAMPLE_NAME"

# Get file sizes for logging
R1_SIZE=$(ls -lh "$R1_FILE" | awk '{print $5}')
R2_SIZE=$(ls -lh "$R2_FILE" | awk '{print $5}')
log_message "File sizes: R1=$R1_SIZE, R2=$R2_SIZE"

# Check if already processed
if [ -f "$PROCESSED/alignment/${SAMPLE_NAME}.sorted.bam" ]; then
    echo "✅ $SAMPLE_NAME already processed, skipping..."
    log_message "Sample $SAMPLE_NAME already processed, skipping pipeline"
    exit 0
fi

# Verify reference files
log_message "Verifying reference files..."
if [ ! -f "$REFERENCE/annotation/Homo_sapiens.GRCh38.110.gtf" ]; then
    handle_error "GTF annotation file not found: $REFERENCE/annotation/Homo_sapiens.GRCh38.110.gtf"
fi

if [ ! -f "$REFERENCE/genome/hisat2_index/grch38/genome.1.ht2" ]; then
    handle_error "HISAT2 index not found: $REFERENCE/genome/hisat2_index/grch38/"
fi

log_message "Reference files validated"

# ==============================================================================
# STEP 1: QUALITY CONTROL
# ==============================================================================
echo ""
echo "🔹 Step 1: Quality Control"
log_message "Starting Step 1: Quality Control"
STEP_START=$(date +%s)

mkdir -p "$PROCESSED/quality_control"

fastqc "$R1_FILE" "$R2_FILE" -o "$PROCESSED/quality_control" -t $THREADS || \
    handle_error "FastQC failed for $SAMPLE_NAME"

STEP_END=$(date +%s)
STEP_TIME=$((STEP_END - STEP_START))
log_message "Step 1 completed in ${STEP_TIME}s"

# ==============================================================================
# STEP 2: READ TRIMMING
# ==============================================================================
echo ""
echo "🔹 Step 2: Read Trimming"
log_message "Starting Step 2: Read Trimming"
STEP_START=$(date +%s)

mkdir -p "$PROCESSED/trimmed_reads"

java -jar tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $THREADS \
    "$R1_FILE" "$R2_FILE" \
    "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_1_paired.fastq.gz" \
    "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_1_unpaired.fastq.gz" \
    "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_2_paired.fastq.gz" \
    "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || \
    handle_error "Trimmomatic failed for $SAMPLE_NAME"

# Log trimming statistics
TRIMMED_R1_SIZE=$(ls -lh "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_1_paired.fastq.gz" | awk '{print $5}')
TRIMMED_R2_SIZE=$(ls -lh "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_2_paired.fastq.gz" | awk '{print $5}')
log_message "Trimmed file sizes: R1=$TRIMMED_R1_SIZE, R2=$TRIMMED_R2_SIZE"

STEP_END=$(date +%s)
STEP_TIME=$((STEP_END - STEP_START))
log_message "Step 2 completed in ${STEP_TIME}s"

# ==============================================================================
# STEP 3: GENOME ALIGNMENT
# ==============================================================================
echo ""
echo "🔹 Step 3: Read Alignment"
log_message "Starting Step 3: Genome Alignment"
STEP_START=$(date +%s)

mkdir -p "$PROCESSED/alignment"

hisat2 -p $THREADS \
    -x reference/genome/hisat2_index/grch38/genome \
    -1 "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_1_paired.fastq.gz" \
    -2 "$PROCESSED/trimmed_reads/${SAMPLE_NAME}_2_paired.fastq.gz" \
    -S "$PROCESSED/alignment/${SAMPLE_NAME}.sam" \
    2> "$PROCESSED/alignment/${SAMPLE_NAME}_hisat2.log" || \
    handle_error "HISAT2 alignment failed for $SAMPLE_NAME"

# Convert to sorted BAM
echo "🔄 Converting to BAM and sorting..."
log_message "Converting SAM to sorted BAM"

samtools view -@ $THREADS -bS "$PROCESSED/alignment/${SAMPLE_NAME}.sam" | \
    samtools sort -@ $THREADS -o "$PROCESSED/alignment/${SAMPLE_NAME}.sorted.bam" || \
    handle_error "SAM to BAM conversion failed for $SAMPLE_NAME"

samtools index "$PROCESSED/alignment/${SAMPLE_NAME}.sorted.bam" || \
    handle_error "BAM indexing failed for $SAMPLE_NAME"

# Clean up SAM file to save space
rm "$PROCESSED/alignment/${SAMPLE_NAME}.sam"

# Log alignment statistics
BAM_SIZE=$(ls -lh "$PROCESSED/alignment/${SAMPLE_NAME}.sorted.bam" | awk '{print $5}')
log_message "BAM file size: $BAM_SIZE"

STEP_END=$(date +%s)
STEP_TIME=$((STEP_END - STEP_START))
log_message "Step 3 completed in ${STEP_TIME}s"

# ==============================================================================
# STEP 4: GENE QUANTIFICATION
# ==============================================================================
echo ""
echo "🔹 Step 4: Gene Quantification"
log_message "Starting Step 4: Gene Quantification"
STEP_START=$(date +%s)

mkdir -p "$PROCESSED/quantification"

stringtie "$PROCESSED/alignment/${SAMPLE_NAME}.sorted.bam" \
    -G reference/annotation/Homo_sapiens.GRCh38.110.gtf \
    -o "$PROCESSED/quantification/${SAMPLE_NAME}.gtf" \
    -A "$PROCESSED/quantification/${SAMPLE_NAME}_abundance.tab" \
    -e -B -p $THREADS || \
    handle_error "StringTie quantification failed for $SAMPLE_NAME"

# Log quantification results
ABUNDANCE_LINES=$(wc -l < "$PROCESSED/quantification/${SAMPLE_NAME}_abundance.tab")
log_message "Gene abundance table: $ABUNDANCE_LINES genes quantified"

STEP_END=$(date +%s)
STEP_TIME=$((STEP_END - STEP_START))
log_message "Step 4 completed in ${STEP_TIME}s"

# ==============================================================================
# PIPELINE COMPLETION
# ==============================================================================
TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - START_TIME))
HOURS=$((TOTAL_TIME / 3600))
MINUTES=$(((TOTAL_TIME % 3600) / 60))
SECONDS=$((TOTAL_TIME % 60))

echo ""
echo "✅ Pipeline completed successfully for $SAMPLE_NAME"
echo "⏰ Total processing time: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""
echo "📁 Results locations:"
echo "   Quality Control: $PROCESSED/quality_control/${SAMPLE_NAME}_*_fastqc.html"
echo "   BAM file: $PROCESSED/alignment/${SAMPLE_NAME}.sorted.bam"
echo "   Gene quantification: $PROCESSED/quantification/${SAMPLE_NAME}_abundance.tab"
echo "   Log file: $LOGFILE"
echo ""
echo "🎯 Status: READY FOR DOWNSTREAM ANALYSIS"
echo "   Note: Pipeline stops here before DEG analysis as requested"

log_message "Pipeline completed successfully for $SAMPLE_NAME"
log_message "Total processing time: ${HOURS}h ${MINUTES}m ${SECONDS}s"
log_message "Results ready for downstream analysis"

echo "╚══════════════════════════════════════════════════════════════════════════════╝"

# Close the if statement for single_sample mode check
fi
