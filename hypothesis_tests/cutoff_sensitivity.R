#!/usr/bin/env Rscript
# Cutoff sensitivity analysis — multiple thresholds for manuscript vs our v3
suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble)})

MS_DIR <- "/Volumes/ExtremeSSD/ibmfs/03_original_analysis/deg_results"
OUR_COUNTS <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_manuscript_v3.txt"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"

# Our v3 DESeq2
raw <- read_tsv(OUR_COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(raw[, -(1:6)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"
colnames(counts) <- sub("_sorted$","",sub("^.*/","",sub(".bam$","",colnames(counts))))
samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
           "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH")
counts <- counts[, samps]; counts <- counts[rowSums(counts)>0,]
meta <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8)),
                                 levels=c("Control","U-AA","G-AA")))
rownames(meta) <- samps
dds <- DESeqDataSetFromMatrix(counts, meta, ~group); dds <- DESeq(dds, parallel=FALSE)

gtf_lines <- readLines(GTF); gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec('gene_id "([^"]+)".*?gene_type "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==3) data.frame(gene_id=x[2], gene_type=x[3]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE); gt$gene_id <- as.character(gt$gene_id)

get_our_deg <- function(g1, g2) {
  res <- results(dds, contrast=c("group", g1, g2))
  df <- as.data.frame(res) %>% rownames_to_column("gene_id")
  df$gene_type <- gt$gene_type[match(df$gene_id, gt$gene_id)]
  df %>% filter(!is.na(padj))
}

# Load manuscript
ms_GvC <- read_tsv(file.path(MS_DIR, "G-AA_vs_Control_results.txt"), show_col_types=FALSE)
ms_UvC <- read_tsv(file.path(MS_DIR, "U-AA_vs_Control_results.txt"), show_col_types=FALSE)
ms_GvU <- read_tsv(file.path(MS_DIR, "G-AA_vs_U-AA_results.txt"), show_col_types=FALSE)
our_GvC <- get_our_deg("G-AA","Control")
our_UvC <- get_our_deg("U-AA","Control")
our_GvU <- get_our_deg("G-AA","U-AA")

count_at_cutoff <- function(df, padj_thr, lfc_thr, biotype=NULL) {
  s <- !is.na(df$padj) & df$padj < padj_thr & abs(df$log2FoldChange) > lfc_thr
  if (!is.null(biotype)) s <- s & df$gene_type == biotype
  sum(s, na.rm=TRUE)
}

# Multi-cutoff sweep
cutoffs <- list(
  list(padj=0.05, lfc=1),
  list(padj=0.05, lfc=1.5),
  list(padj=0.05, lfc=2),
  list(padj=0.01, lfc=1),
  list(padj=0.01, lfc=1.5),
  list(padj=0.01, lfc=2),
  list(padj=0.001, lfc=1),
  list(padj=0.001, lfc=2)
)

cat("\n=== G-AA vs Control: mRNA DEG count at various cutoffs ===\n")
cat(sprintf("%-6s %-6s | %-10s %-10s %-10s %-10s\n", "padj", "|LFC|", "MS_mRNA", "Our_mRNA", "Ratio", "Δ"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (cf in cutoffs) {
  ms_n  <- count_at_cutoff(ms_GvC, cf$padj, cf$lfc, "protein_coding")
  our_n <- count_at_cutoff(our_GvC, cf$padj, cf$lfc, "protein_coding")
  ratio <- if (ms_n>0) round(our_n/ms_n,2) else NA
  delta <- our_n - ms_n
  cat(sprintf("<%-5g %-6g | %-10d %-10d %-10.2f %+d\n",
              cf$padj, cf$lfc, ms_n, our_n, ratio, delta))
}

cat("\n=== G-AA vs Control: lncRNA DEG count ===\n")
cat(sprintf("%-6s %-6s | %-10s %-10s %-10s %-10s\n", "padj", "|LFC|", "MS_lnc", "Our_lnc", "Ratio", "Δ"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (cf in cutoffs) {
  ms_n  <- count_at_cutoff(ms_GvC, cf$padj, cf$lfc, "lncRNA")
  our_n <- count_at_cutoff(our_GvC, cf$padj, cf$lfc, "lncRNA")
  ratio <- if (ms_n>0) round(our_n/ms_n,2) else NA
  delta <- our_n - ms_n
  cat(sprintf("<%-5g %-6g | %-10d %-10d %-10.2f %+d\n",
              cf$padj, cf$lfc, ms_n, our_n, ratio, delta))
}

cat("\n=== U-AA vs Control: mRNA DEG count ===\n")
cat(sprintf("%-6s %-6s | %-10s %-10s %-10s %-10s\n", "padj", "|LFC|", "MS_mRNA", "Our_mRNA", "Ratio", "Δ"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (cf in cutoffs) {
  ms_n  <- count_at_cutoff(ms_UvC, cf$padj, cf$lfc, "protein_coding")
  our_n <- count_at_cutoff(our_UvC, cf$padj, cf$lfc, "protein_coding")
  ratio <- if (ms_n>0) round(our_n/ms_n,2) else NA
  delta <- our_n - ms_n
  cat(sprintf("<%-5g %-6g | %-10d %-10d %-10.2f %+d\n",
              cf$padj, cf$lfc, ms_n, our_n, ratio, delta))
}

cat("\n=== U-AA vs Control: lncRNA DEG count ===\n")
cat(sprintf("%-6s %-6s | %-10s %-10s %-10s %-10s\n", "padj", "|LFC|", "MS_lnc", "Our_lnc", "Ratio", "Δ"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (cf in cutoffs) {
  ms_n  <- count_at_cutoff(ms_UvC, cf$padj, cf$lfc, "lncRNA")
  our_n <- count_at_cutoff(our_UvC, cf$padj, cf$lfc, "lncRNA")
  ratio <- if (ms_n>0) round(our_n/ms_n,2) else NA
  delta <- our_n - ms_n
  cat(sprintf("<%-5g %-6g | %-10d %-10d %-10.2f %+d\n",
              cf$padj, cf$lfc, ms_n, our_n, ratio, delta))
}

cat("\n=== G-AA vs U-AA total DEG ===\n")
cat(sprintf("%-6s %-6s | %-10s %-10s\n", "padj", "|LFC|", "MS_total", "Our_total"))
cat(paste(rep("-", 40), collapse=""), "\n")
for (cf in cutoffs) {
  ms_n  <- count_at_cutoff(ms_GvU, cf$padj, cf$lfc)
  our_n <- count_at_cutoff(our_GvU, cf$padj, cf$lfc)
  cat(sprintf("<%-5g %-6g | %-10d %-10d\n", cf$padj, cf$lfc, ms_n, our_n))
}

# Save TSV summary
results <- expand.grid(
  contrast = c("G-AA vs Ctrl","U-AA vs Ctrl","G-AA vs U-AA"),
  padj_thr = c(0.05, 0.01, 0.001),
  lfc_thr = c(1, 1.5, 2),
  biotype = c("protein_coding","lncRNA","all")
)
results$ms_n <- NA; results$our_n <- NA
for (i in 1:nrow(results)) {
  r <- results[i,]
  ms_df <- list("G-AA vs Ctrl"=ms_GvC, "U-AA vs Ctrl"=ms_UvC, "G-AA vs U-AA"=ms_GvU)[[r$contrast]]
  our_df<- list("G-AA vs Ctrl"=our_GvC, "U-AA vs Ctrl"=our_UvC, "G-AA vs U-AA"=our_GvU)[[r$contrast]]
  bt <- if (r$biotype=="all") NULL else as.character(r$biotype)
  results$ms_n[i]  <- count_at_cutoff(ms_df, r$padj_thr, r$lfc_thr, bt)
  results$our_n[i] <- count_at_cutoff(our_df, r$padj_thr, r$lfc_thr, bt)
}
results$ratio <- round(results$our_n / results$ms_n, 3)
results$delta <- results$our_n - results$ms_n
write_tsv(results, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/cutoff_sensitivity_table.tsv")
cat("\nSaved: docs/cutoff_sensitivity_table.tsv\n")
