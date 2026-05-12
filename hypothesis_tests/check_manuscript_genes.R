#!/usr/bin/env Rscript
# Check if manuscript-mentioned genes are detected in our analyses
# - Analysis A: 14 samples, 2 controls (manuscript replica)
# - Analysis C: 17 samples, 5 controls + batch-aware design
suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble)})

DDS_RDS <- "/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/control_comparison/dds_all.rds"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"

# Load DESeq2 objects
dds_all <- readRDS(DDS_RDS)
dds_A <- dds_all$A   # 14 samples, manuscript replica
dds_C <- dds_all$C   # 17 samples, batch-aware

# Build gene_name → gene_id (with version) lookup from GTF
cat("Building gene symbol → Ensembl lookup...\n")
gtf_lines <- readLines(GTF)
gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec(
  'gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4)
  data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE)
gt$gene_id <- as.character(gt$gene_id)
cat(sprintf("Loaded %d gene annotations\n", nrow(gt)))

# Manuscript-mentioned genes
manuscript_genes <- list(
  "G-AA vs U-AA (4 mRNAs)" = c("SFT2D3","OR52K1","LRRC24","TEN1-CDK3"),
  "Key 11 lncRNAs (avg CPM>=10, CPM>=1 all)" = c(
    "ATP1A1-AS1","LINC01036","HCG11","HCP5","SNHG32","PSMB8-AS1",
    "TAGAP-AS1","MALAT1","FAM30A","USP3-AS1","MIR22HG"),
  "RT-qPCR validated 3 (key biomarkers)" = c("ATP1A1-AS1","USP3-AS1","SNHG32"),
  "FA-related mutations (manuscript mentions)" = c("FANCG","FANCD2","FANCD1","FANCA","TERT")
)

# Helper: get gene_id from symbol
find_gene_id <- function(symbol) {
  matches <- gt$gene_id[gt$gene_name == symbol]
  if (length(matches) == 0) return(NA)
  matches[1]
}

# Helper: get DEG info for gene in specific contrast and DESeq2 object
get_gene_info <- function(dds, contrast_groups, gene_id) {
  if (is.na(gene_id)) return(list(baseMean=NA, log2FC=NA, padj=NA, sig=NA))
  res <- results(dds, contrast=c("group", contrast_groups[1], contrast_groups[2]))
  if (!gene_id %in% rownames(res)) return(list(baseMean=NA, log2FC=NA, padj=NA, sig=NA))
  r <- res[gene_id, ]
  list(
    baseMean = round(r$baseMean, 1),
    log2FC = round(r$log2FoldChange, 2),
    padj = signif(r$padj, 3),
    sig = !is.na(r$padj) && r$padj < 0.05 && abs(r$log2FoldChange) > 1
  )
}

# Build comprehensive table for each gene
cat("\n\n========== MANUSCRIPT GENES — CHECK IN OUR ANALYSES ==========\n\n")

for (group_label in names(manuscript_genes)) {
  cat(sprintf("\n=== %s ===\n", group_label))
  for (sym in manuscript_genes[[group_label]]) {
    gene_id <- find_gene_id(sym)
    if (is.na(gene_id)) {
      cat(sprintf("  %-15s — NOT FOUND in GTF\n", sym))
      next
    }
    gtype <- gt$gene_type[gt$gene_id == gene_id]
    cat(sprintf("  %-15s (%s, %s)\n", sym, gene_id, gtype))

    # G-AA vs Control
    info_A <- get_gene_info(dds_A, c("G-AA","Control"), gene_id)
    info_C <- get_gene_info(dds_C, c("G-AA","Control"), gene_id)
    cat(sprintf("    G-AA vs Ctrl  A(14):  log2FC=%6s  padj=%8s  %s\n",
                ifelse(is.na(info_A$log2FC), "NA", info_A$log2FC),
                ifelse(is.na(info_A$padj), "NA", info_A$padj),
                ifelse(isTRUE(info_A$sig), "✅sig", "❌nonsig")))
    cat(sprintf("    G-AA vs Ctrl  C(17):  log2FC=%6s  padj=%8s  %s\n",
                ifelse(is.na(info_C$log2FC), "NA", info_C$log2FC),
                ifelse(is.na(info_C$padj), "NA", info_C$padj),
                ifelse(isTRUE(info_C$sig), "✅sig", "❌nonsig")))

    # U-AA vs Control
    info_A2 <- get_gene_info(dds_A, c("U-AA","Control"), gene_id)
    info_C2 <- get_gene_info(dds_C, c("U-AA","Control"), gene_id)
    cat(sprintf("    U-AA vs Ctrl  A(14):  log2FC=%6s  padj=%8s  %s\n",
                ifelse(is.na(info_A2$log2FC), "NA", info_A2$log2FC),
                ifelse(is.na(info_A2$padj), "NA", info_A2$padj),
                ifelse(isTRUE(info_A2$sig), "✅sig", "❌nonsig")))
    cat(sprintf("    U-AA vs Ctrl  C(17):  log2FC=%6s  padj=%8s  %s\n",
                ifelse(is.na(info_C2$log2FC), "NA", info_C2$log2FC),
                ifelse(is.na(info_C2$padj), "NA", info_C2$padj),
                ifelse(isTRUE(info_C2$sig), "✅sig", "❌nonsig")))

    # G-AA vs U-AA
    info_A3 <- get_gene_info(dds_A, c("G-AA","U-AA"), gene_id)
    info_C3 <- get_gene_info(dds_C, c("G-AA","U-AA"), gene_id)
    cat(sprintf("    G-AA vs U-AA  A(14):  log2FC=%6s  padj=%8s  %s\n",
                ifelse(is.na(info_A3$log2FC), "NA", info_A3$log2FC),
                ifelse(is.na(info_A3$padj), "NA", info_A3$padj),
                ifelse(isTRUE(info_A3$sig), "✅sig", "❌nonsig")))
    cat(sprintf("    G-AA vs U-AA  C(17):  log2FC=%6s  padj=%8s  %s\n",
                ifelse(is.na(info_C3$log2FC), "NA", info_C3$log2FC),
                ifelse(is.na(info_C3$padj), "NA", info_C3$padj),
                ifelse(isTRUE(info_C3$sig), "✅sig", "❌nonsig")))
  }
}

# Build summary TSV
all_results <- data.frame()
for (group_label in names(manuscript_genes)) {
  for (sym in manuscript_genes[[group_label]]) {
    gene_id <- find_gene_id(sym)
    gtype <- if (!is.na(gene_id)) gt$gene_type[gt$gene_id == gene_id] else NA
    for (cc_label in c("G-AA_vs_Control", "U-AA_vs_Control", "G-AA_vs_U-AA")) {
      cc <- strsplit(cc_label, "_vs_")[[1]]
      a <- get_gene_info(dds_A, cc, gene_id)
      c <- get_gene_info(dds_C, cc, gene_id)
      all_results <- rbind(all_results, data.frame(
        gene_group=group_label, gene_symbol=sym, gene_id=gene_id, gene_type=gtype,
        contrast=cc_label,
        A_baseMean=a$baseMean, A_log2FC=a$log2FC, A_padj=a$padj, A_sig=a$sig,
        C_baseMean=c$baseMean, C_log2FC=c$log2FC, C_padj=c$padj, C_sig=c$sig
      ))
    }
  }
}
write_tsv(all_results, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/manuscript_genes_check.tsv")
cat(sprintf("\nSaved table: docs/manuscript_genes_check.tsv (%d rows)\n", nrow(all_results)))
