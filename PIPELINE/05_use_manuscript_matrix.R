#!/usr/bin/env Rscript
# Use manuscript count matrix directly (Option 1 strategy)
# Verified: reproduces manuscript DEG values 100%
#
# Usage:
#   Rscript 05_use_manuscript_matrix.R [contrast] [out_dir]
#   contrast: "G-AA_vs_Control" (default), "U-AA_vs_Control", "G-AA_vs_U-AA"

suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble)})

args <- commandArgs(trailingOnly=TRUE)
CONTRAST <- if (length(args)>=1) args[1] else "G-AA_vs_Control"
OUT <- if (length(args)>=2) args[2] else "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/deseq2_replica"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

MATRIX <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"

cat("Loading manuscript count matrix...\n")
raw <- read_tsv(MATRIX, show_col_types=FALSE)
cat(sprintf("Loaded: %d genes × %d columns\n", nrow(raw), ncol(raw)))

# Manuscript 14 samples (eAT5.R order)
ms_samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
              "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18",
              "AA-HMH","AA-PJH")

# Extract count matrix
counts <- as.matrix(raw[, ms_samps])
rownames(counts) <- raw$EnsemblID
storage.mode(counts) <- "integer"

# Manuscript filter (eAT5.R Oct 20 2025): rowSums > 0
counts <- counts[rowSums(counts) > 0, ]
cat(sprintf("After rowSums>0: %d genes\n", nrow(counts)))

# Manuscript metadata
meta <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8)),
                                 levels=c("Control","U-AA","G-AA")))
rownames(meta) <- ms_samps

# Run DESeq2 (manuscript-style ~ group)
cat("Running DESeq2 (manuscript-style ~ group)...\n")
dds <- DESeqDataSetFromMatrix(counts, meta, ~ group)
dds <- DESeq(dds, parallel=FALSE)

# Helper: get contrast results + filter
get_contrast <- function(g1, g2, label) {
  res <- results(dds, contrast=c("group", g1, g2))
  df <- as.data.frame(res) %>% rownames_to_column("gene_id")

  # Add gene_type from GTF
  gtf_lines <- readLines(GTF); gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
  parsed <- regmatches(gene_lines, regexec(
    'gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
  gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4)
    data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
  gt <- gt %>% distinct(gene_id, .keep_all=TRUE)
  df$gene_name <- gt$gene_name[match(df$gene_id, gt$gene_id)]
  df$gene_type <- gt$gene_type[match(df$gene_id, gt$gene_id)]

  # Manuscript post-filter (eAT5.R line 154): padj != NA
  df <- df %>% filter(!is.na(padj))

  # Add significance label
  df$significance <- with(df, ifelse(padj<0.05 & log2FoldChange>1, "Upregulated",
                            ifelse(padj<0.05 & log2FoldChange < -1, "Downregulated", "Not significant")))

  out <- file.path(OUT, paste0("DE_", label, ".tsv"))
  write_tsv(df, out)
  cat(sprintf("Saved: %s\n", out))
  cat(sprintf("  Total: %d, Sig: %d (mRNA %d, lncRNA %d)\n",
              nrow(df), sum(df$significance != "Not significant"),
              sum(df$significance != "Not significant" & df$gene_type == "protein_coding"),
              sum(df$significance != "Not significant" & df$gene_type == "lncRNA")))
  df
}

# Run all 3 contrasts
cat("\n=== Running all 3 contrasts ===\n")
gvc <- get_contrast("G-AA", "Control", "G-AA_vs_Control")
uvc <- get_contrast("U-AA", "Control", "U-AA_vs_Control")
gvu <- get_contrast("G-AA", "U-AA", "G-AA_vs_U-AA")

saveRDS(dds, file.path(OUT, "dds_manuscript_replica.rds"))
cat(sprintf("\nAll outputs in: %s\n", OUT))
cat("Expected manuscript results: 2078 mRNA + 1167 lncRNA (G-AA vs Ctrl), 1315 + 992 (U-AA vs Ctrl), 4 + 4 (G-AA vs U-AA)\n")
