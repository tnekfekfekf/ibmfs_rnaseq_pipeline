#!/usr/bin/env Rscript
# Compare DEG results: 2 controls (manuscript) vs 5 controls (with Child1/2/3 added)
# All 17 samples quantified with IDENTICAL pipeline (v3 options) → differences are biology, not pipeline.
suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble); library(sva)})

COUNTS <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_v3_17samples.txt"
GTF    <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"
OUT    <- "/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/control_comparison"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# Load counts
raw <- read_tsv(COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(raw[, -(1:2)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"
cat(sprintf("Loaded: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# Parse gene_type
gtf_lines <- readLines(GTF)
gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec('gene_id "([^"]+)".*?gene_type "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==3) data.frame(gene_id=x[2], gene_type=x[3]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE)
gt$gene_id <- as.character(gt$gene_id)

g_aa <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3")
ctrl_internal <- c("AA-PRO","AA-KEW")
ctrl_public   <- c("Child1","Child2","Child3")
u_aa  <- c("AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH")

# ===== Analysis A: 2 controls (manuscript replica) =====
cat("\n========== ANALYSIS A: 2 controls (manuscript replica, n=14) ==========\n")
samps_A <- c(g_aa, ctrl_internal, u_aa)
counts_A <- counts[, samps_A]; counts_A <- counts_A[rowSums(counts_A) > 0, ]
meta_A <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8)),
                                   levels=c("Control","U-AA","G-AA")))
rownames(meta_A) <- samps_A
dds_A <- DESeqDataSetFromMatrix(counts_A, meta_A, ~ group); dds_A <- DESeq(dds_A, parallel=FALSE)

# ===== Analysis B: 5 controls (2 internal + 3 public, simple) =====
cat("\n========== ANALYSIS B: 5 controls (n=17, simple ~group, NAIVE) ==========\n")
samps_B <- c(g_aa, ctrl_internal, ctrl_public, u_aa)
counts_B <- counts[, samps_B]; counts_B <- counts_B[rowSums(counts_B) > 0, ]
meta_B <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",5),rep("U-AA",8)),
                                   levels=c("Control","U-AA","G-AA")),
                      cohort=factor(c(rep("internal",4),rep("internal",2),rep("public",3),rep("internal",8)),
                                    levels=c("internal","public")))
rownames(meta_B) <- samps_B
dds_B <- DESeqDataSetFromMatrix(counts_B, meta_B, ~ group); dds_B <- DESeq(dds_B, parallel=FALSE)

# ===== Analysis C: 5 controls with batch-aware design (~ cohort + group) =====
cat("\n========== ANALYSIS C: 5 controls + batch-aware (~cohort + group) ==========\n")
dds_C <- DESeqDataSetFromMatrix(counts_B, meta_B, ~ cohort + group); dds_C <- DESeq(dds_C, parallel=FALSE)

# ===== Analysis D: 3 public controls only (Child1/2/3) =====
cat("\n========== ANALYSIS D: 3 Child controls only (n=15) ==========\n")
samps_D <- c(g_aa, ctrl_public, u_aa)
counts_D <- counts[, samps_D]; counts_D <- counts_D[rowSums(counts_D) > 0, ]
meta_D <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",3),rep("U-AA",8)),
                                   levels=c("Control","U-AA","G-AA")))
rownames(meta_D) <- samps_D
dds_D <- DESeqDataSetFromMatrix(counts_D, meta_D, ~ group); dds_D <- DESeq(dds_D, parallel=FALSE)

# Helper: count sig DEGs (filtered by gene_type)
sig_counts <- function(dds, label) {
  out <- list()
  for (cc in list(c("G-AA","Control"), c("U-AA","Control"), c("G-AA","U-AA"))) {
    res <- results(dds, contrast=c("group", cc[1], cc[2]))
    df <- as.data.frame(res) %>% rownames_to_column("gene_id")
    df$gene_type <- gt$gene_type[match(df$gene_id, gt$gene_id)]
    df <- df %>% filter(!is.na(padj))
    sig_total <- sum(df$padj<0.05 & abs(df$log2FoldChange)>1, na.rm=TRUE)
    sig_m <- sum(df$padj<0.05 & abs(df$log2FoldChange)>1 & df$gene_type=="protein_coding", na.rm=TRUE)
    sig_l <- sum(df$padj<0.05 & abs(df$log2FoldChange)>1 & df$gene_type=="lncRNA", na.rm=TRUE)
    out[[paste(cc, collapse="_vs_")]] <- df
    cat(sprintf("  %-22s tested=%d  total_sig=%d  mRNA=%d  lncRNA=%d\n",
                paste(cc[1],"vs",cc[2]), nrow(df), sig_total, sig_m, sig_l))
  }
  out
}

cat("\n=== Significance summary (padj<0.05, |log2FC|>1) ===\n")
cat("\n--- Analysis A: 2 internal controls (manuscript replica, n=14) ---\n")
resA <- sig_counts(dds_A, "A")
cat("\n--- Analysis B: 5 controls, naive ~group (n=17) ---\n")
resB <- sig_counts(dds_B, "B")
cat("\n--- Analysis C: 5 controls, batch-aware ~cohort+group (n=17) ---\n")
resC <- sig_counts(dds_C, "C")
cat("\n--- Analysis D: 3 public Child controls only (n=15) ---\n")
resD <- sig_counts(dds_D, "D")

# Save results
saveRDS(list(A=dds_A, B=dds_B, C=dds_C, D=dds_D), file.path(OUT, "dds_all.rds"))
saveRDS(list(A=resA, B=resB, C=resC, D=resD), file.path(OUT, "deg_tables.rds"))

# Detailed comparison: G-AA vs Control DEG overlap
cat("\n=== G-AA vs Control DEG overlap ===\n")
gvc <- list(
  A = resA[["G-AA_vs_Control"]],
  B = resB[["G-AA_vs_Control"]],
  C = resC[["G-AA_vs_Control"]],
  D = resD[["G-AA_vs_Control"]]
)
get_sig <- function(df) df$gene_id[!is.na(df$padj) & df$padj<0.05 & abs(df$log2FoldChange)>1]
sig_lists <- lapply(gvc, get_sig)
cat("Sig counts:\n")
for (n in names(sig_lists)) cat(sprintf("  %s: %d\n", n, length(sig_lists[[n]])))
cat("\nOverlaps (intersection):\n")
cat(sprintf("  A ∩ B = %d (n=%d sigs share between manuscript 2ctrl and naive 5ctrl)\n",
            length(intersect(sig_lists$A, sig_lists$B)), length(sig_lists$A)))
cat(sprintf("  A ∩ C = %d (n=%d, manuscript vs batch-aware 5ctrl)\n",
            length(intersect(sig_lists$A, sig_lists$C)), length(sig_lists$A)))
cat(sprintf("  A ∩ D = %d (n=%d, manuscript vs public-only 3ctrl)\n",
            length(intersect(sig_lists$A, sig_lists$D)), length(sig_lists$A)))
cat(sprintf("  B ∩ C = %d (n=%d, naive 5ctrl vs batch-aware)\n",
            length(intersect(sig_lists$B, sig_lists$C)), length(sig_lists$B)))

# Save summary table
summary_tab <- data.frame(
  Analysis = c("A: 2 internal (manuscript)", "B: 5 ctrl naive", "C: 5 ctrl batch-aware", "D: 3 public only"),
  n_samples = c(14, 17, 17, 15),
  GvC_mRNA  = sapply(list(resA$`G-AA_vs_Control`, resB$`G-AA_vs_Control`, resC$`G-AA_vs_Control`, resD$`G-AA_vs_Control`),
                     function(d) sum(!is.na(d$padj) & d$padj<0.05 & abs(d$log2FoldChange)>1 & d$gene_type=="protein_coding")),
  GvC_lncRNA= sapply(list(resA$`G-AA_vs_Control`, resB$`G-AA_vs_Control`, resC$`G-AA_vs_Control`, resD$`G-AA_vs_Control`),
                     function(d) sum(!is.na(d$padj) & d$padj<0.05 & abs(d$log2FoldChange)>1 & d$gene_type=="lncRNA")),
  UvC_mRNA  = sapply(list(resA$`U-AA_vs_Control`, resB$`U-AA_vs_Control`, resC$`U-AA_vs_Control`, resD$`U-AA_vs_Control`),
                     function(d) sum(!is.na(d$padj) & d$padj<0.05 & abs(d$log2FoldChange)>1 & d$gene_type=="protein_coding")),
  UvC_lncRNA= sapply(list(resA$`U-AA_vs_Control`, resB$`U-AA_vs_Control`, resC$`U-AA_vs_Control`, resD$`U-AA_vs_Control`),
                     function(d) sum(!is.na(d$padj) & d$padj<0.05 & abs(d$log2FoldChange)>1 & d$gene_type=="lncRNA")),
  GvU_total = sapply(list(resA$`G-AA_vs_U-AA`, resB$`G-AA_vs_U-AA`, resC$`G-AA_vs_U-AA`, resD$`G-AA_vs_U-AA`),
                     function(d) sum(!is.na(d$padj) & d$padj<0.05 & abs(d$log2FoldChange)>1))
)
write_tsv(summary_tab, file.path(OUT, "DEG_summary_4analyses.tsv"))
cat("\n=== Summary table ===\n")
print(summary_tab)
cat(sprintf("\nSaved: %s/DEG_summary_4analyses.tsv\n", OUT))
