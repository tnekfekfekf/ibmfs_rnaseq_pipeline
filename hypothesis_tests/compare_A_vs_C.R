#!/usr/bin/env Rscript
# Comprehensive comparison: Analysis A (2 ctrl, manuscript) vs Analysis C (5 ctrl, batch-aware)
# Goal: demonstrate that adding public controls (with proper design) does NOT change conclusions
suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble)})

DDS_RDS <- "/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/control_comparison/dds_all.rds"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"

dds_all <- readRDS(DDS_RDS)
dds_A <- dds_all$A  # 14 samples, 2 controls (manuscript replica)
dds_C <- dds_all$C  # 17 samples, 5 controls (~ cohort + group)

# Gene annotation
gtf_lines <- readLines(GTF); gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec(
  'gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4)
  data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE)
gt$gene_id <- as.character(gt$gene_id)

get_deg <- function(dds, g1, g2) {
  res <- results(dds, contrast=c("group", g1, g2))
  df <- as.data.frame(res) %>% rownames_to_column("gene_id")
  df$gene_name <- gt$gene_name[match(df$gene_id, gt$gene_id)]
  df$gene_type <- gt$gene_type[match(df$gene_id, gt$gene_id)]
  df %>% filter(!is.na(padj))
}

A_GvC <- get_deg(dds_A, "G-AA","Control")
A_UvC <- get_deg(dds_A, "U-AA","Control")
A_GvU <- get_deg(dds_A, "G-AA","U-AA")
C_GvC <- get_deg(dds_C, "G-AA","Control")
C_UvC <- get_deg(dds_C, "U-AA","Control")
C_GvU <- get_deg(dds_C, "G-AA","U-AA")

is_sig <- function(d, padj=0.05, lfc=1) !is.na(d$padj) & d$padj<padj & abs(d$log2FoldChange)>lfc

# =====================================================
cat("=========================================================\n")
cat("    ANALYSIS A vs ANALYSIS C — Comprehensive Comparison\n")
cat("    A: 14 samples, 2 internal controls (manuscript replica)\n")
cat("    C: 17 samples, 5 controls (~ cohort + group, batch-aware)\n")
cat("=========================================================\n\n")

# Table 1: Total DEG counts
cat("=== Table 1: DEG counts per contrast/biotype ===\n\n")
build_summary <- function(d, label, dataset) {
  list(
    contrast=label, dataset=dataset,
    total = sum(is_sig(d)),
    mRNA = sum(is_sig(d) & d$gene_type=="protein_coding"),
    lncRNA = sum(is_sig(d) & d$gene_type=="lncRNA"),
    up = sum(is_sig(d) & d$log2FoldChange>0),
    down = sum(is_sig(d) & d$log2FoldChange<0)
  )
}
results_list <- list(
  build_summary(A_GvC, "G-AA vs Ctrl", "A (14, 2 ctrl)"),
  build_summary(C_GvC, "G-AA vs Ctrl", "C (17, 5 ctrl)"),
  build_summary(A_UvC, "U-AA vs Ctrl", "A (14, 2 ctrl)"),
  build_summary(C_UvC, "U-AA vs Ctrl", "C (17, 5 ctrl)"),
  build_summary(A_GvU, "G-AA vs U-AA", "A (14, 2 ctrl)"),
  build_summary(C_GvU, "G-AA vs U-AA", "C (17, 5 ctrl)")
)
summary_tab <- do.call(rbind, lapply(results_list, function(x) data.frame(x)))
print(summary_tab, row.names=FALSE)

# Table 2: DEG overlap (gene-level)
cat("\n\n=== Table 2: Gene-level overlap between A and C ===\n\n")
overlap_table <- function(dA, dC, label) {
  sig_A <- dA$gene_id[is_sig(dA)]
  sig_C <- dC$gene_id[is_sig(dC)]
  overlap <- intersect(sig_A, sig_C)
  data.frame(
    contrast=label,
    A_sig=length(sig_A),
    C_sig=length(sig_C),
    overlap=length(overlap),
    A_minus_C=length(setdiff(sig_A, sig_C)),
    C_minus_A=length(setdiff(sig_C, sig_A)),
    pct_overlap_of_A=sprintf("%.1f%%", 100*length(overlap)/length(sig_A))
  )
}
overlap_tab <- rbind(
  overlap_table(A_GvC, C_GvC, "G-AA vs Ctrl"),
  overlap_table(A_UvC, C_UvC, "U-AA vs Ctrl"),
  overlap_table(A_GvU, C_GvU, "G-AA vs U-AA")
)
print(overlap_tab, row.names=FALSE)

# Table 3: log2FC correlation between A and C (continuous)
cat("\n\n=== Table 3: Continuous-value correlation (log2FC, padj) ===\n\n")
cor_table <- function(dA, dC, label) {
  m <- merge(dA[, c("gene_id","baseMean","log2FoldChange","padj")],
             dC[, c("gene_id","baseMean","log2FoldChange","padj")],
             by="gene_id", suffixes=c("_A","_C"))
  data.frame(
    contrast=label,
    common_genes=nrow(m),
    baseMean_r=round(cor(m$baseMean_A, m$baseMean_C, use="complete"),4),
    log2FC_r=round(cor(m$log2FoldChange_A, m$log2FoldChange_C, use="complete"),4),
    log2FC_diff_median=round(median(abs(m$log2FoldChange_A - m$log2FoldChange_C), na.rm=TRUE),4)
  )
}
cor_tab <- rbind(
  cor_table(A_GvC, C_GvC, "G-AA vs Ctrl"),
  cor_table(A_UvC, C_UvC, "U-AA vs Ctrl"),
  cor_table(A_GvU, C_GvU, "G-AA vs U-AA")
)
print(cor_tab, row.names=FALSE)

# Table 4: Manuscript-highlighted genes — same outcome in A and C?
cat("\n\n=== Table 4: Manuscript-highlighted genes — A vs C concordance ===\n\n")
key_genes <- c(
  "HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",  # 6 strong lncRNAs (manuscript Figure 2)
  "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1", # 5 moderate lncRNAs
  "TEN1-CDK3","SFT2D3","OR52K1","LRRC24",                   # G-AA vs U-AA 4 mRNAs
  "FANCA","FANCG","FANCD2","TERT"                            # FA-related
)

cat(sprintf("%-15s %-12s | %-22s | %-22s | %s\n",
            "Gene", "Type", "Analysis A (G-AA vs Ctrl)", "Analysis C (G-AA vs Ctrl)", "Match"))
cat(paste(rep("-", 110), collapse=""), "\n")
for (sym in key_genes) {
  gid <- gt$gene_id[gt$gene_name == sym][1]
  if (is.na(gid)) next
  gtype <- gt$gene_type[gt$gene_name == sym][1]
  ra <- A_GvC[A_GvC$gene_id == gid, ]
  rc <- C_GvC[C_GvC$gene_id == gid, ]
  if (nrow(ra) == 0 || nrow(rc) == 0) next

  a_sig <- isTRUE(ra$padj<0.05 && abs(ra$log2FoldChange)>1)
  c_sig <- isTRUE(rc$padj<0.05 && abs(rc$log2FoldChange)>1)
  match_str <- if (a_sig == c_sig) {
    if (a_sig) "✅ both sig" else "○ both n.s."
  } else "❌ DIFFER"

  cat(sprintf("%-15s %-12s | LFC=%6.2f padj=%8.2e | LFC=%6.2f padj=%8.2e | %s\n",
              sym, substr(gtype,1,12), ra$log2FoldChange, ra$padj,
              rc$log2FoldChange, rc$padj, match_str))
}

# Save table
write_tsv(summary_tab, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/A_vs_C_DEG_counts.tsv")
write_tsv(overlap_tab, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/A_vs_C_overlap.tsv")
write_tsv(cor_tab, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/A_vs_C_correlation.tsv")

cat("\n\nAll tables saved to docs/A_vs_C_*.tsv\n")
