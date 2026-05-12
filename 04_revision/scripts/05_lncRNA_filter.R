#!/usr/bin/env Rscript
# Reviewer 1: re-do lncRNA filtering with relaxed/group-aware thresholds.
# Original (paper): mean CPM >= 10 AND CPM >= 1 across ALL 14 samples.
# This excludes biologically relevant low-expression lncRNAs.
#
# Relaxed strategies tested:
#   (R1) "any-group" filter:  CPM >= 1 in >= 50% of samples WITHIN ANY group
#   (R2) "lncRNA-aware":      mean CPM >= 1 in any group, CPM >= 0.5 in >= 50%
#                             of that group
#   (R3) untouched: only DESeq2's independent filtering (current best practice)
#
# Robust DE-lncRNA = called by all 3 strategies + DESeq2 padj<0.05 + |LFC|>1.

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(dplyr); library(readr)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
OUT  <- file.path(ROOT, "deseq2", "lncRNA"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dds  <- readRDS(file.path(ROOT, "deseq2/dds_full.rds"))
gm   <- readRDS(file.path(ROOT, "deseq2/gene_meta.rds"))
meta <- readRDS(file.path(ROOT, "deseq2/samples_meta.rds"))

# Restrict universe to lncRNAs (biotype) per Ensembl/GENCODE
lnc_ids <- gm$gene_id[gm$gene_biotype == "lncRNA"]
dds_l <- dds[rownames(dds) %in% lnc_ids, ]
message("[05] starting universe: ", nrow(dds_l), " lncRNA genes (post pre-filter)")

cpm_mat <- cpm(counts(dds_l), normalized.lib.sizes = FALSE)

# Strategy ORIGINAL (paper)
keep_orig <- (rowMeans(cpm_mat) >= 10) & (rowSums(cpm_mat >= 1) == ncol(cpm_mat))
# Strategy R1: >=50% samples >= 1 CPM within ANY group
g_combined <- meta$group_combined
keep_R1 <- sapply(rownames(cpm_mat), function(g) {
  any(sapply(levels(g_combined), function(grp) {
    idx <- which(g_combined == grp)
    sum(cpm_mat[g, idx] >= 1) / length(idx) >= 0.5
  }))
})
# Strategy R2: mean CPM >= 1 in any group AND >= 50% of that group >= 0.5
keep_R2 <- sapply(rownames(cpm_mat), function(g) {
  any(sapply(levels(g_combined), function(grp) {
    idx <- which(g_combined == grp)
    (mean(cpm_mat[g, idx]) >= 1) & (sum(cpm_mat[g, idx] >= 0.5) / length(idx) >= 0.5)
  }))
})
# Strategy R3: no extra filter — DESeq2 independent filtering only
keep_R3 <- rep(TRUE, nrow(cpm_mat))

message(sprintf("  ORIG: %d   R1: %d   R2: %d   R3: %d",
                sum(keep_orig), sum(keep_R1), sum(keep_R2), sum(keep_R3)))

run_de <- function(dds_in, keep_logical, label) {
  d <- dds_in[keep_logical, ]
  d <- DESeq(d, parallel = FALSE)
  out <- list()
  for (cc in list(c("group_combined","g_BMF","Control"),
                  c("group_combined","u_BMF","Control"),
                  c("group_combined","g_BMF","u_BMF"))) {
    r <- results(d, contrast = cc, alpha = 0.05)
    df <- as.data.frame(r) %>% tibble::rownames_to_column("gene_id") %>%
           left_join(gm, by = "gene_id") %>% arrange(padj)
    fname <- sprintf("lncRNA_DE_%s_%s_vs_%s.tsv", label, cc[2], cc[3])
    write_tsv(df, file.path(OUT, fname))
    sig <- sum(df$padj < 0.05 & abs(df$log2FoldChange) > 1, na.rm = TRUE)
    message("   ", label, " ", cc[2], "_vs_", cc[3], ": ", sig, " sig (FDR<0.05)")
    out[[paste(cc[2], cc[3], sep = "_vs_")]] <- df %>%
      filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) %>%
      pull(gene_id)
  }
  out
}

de_orig <- run_de(dds_l, keep_orig, "ORIG")
de_R1   <- run_de(dds_l, keep_R1,   "R1_anyGroup50")
de_R2   <- run_de(dds_l, keep_R2,   "R2_groupAware")
de_R3   <- run_de(dds_l, keep_R3,   "R3_noPreFilter")

# Robust intersection per contrast (called by all 4 strategies)
robust <- lapply(names(de_orig), function(ct) {
  Reduce(intersect, list(de_orig[[ct]], de_R1[[ct]], de_R2[[ct]], de_R3[[ct]]))
})
names(robust) <- names(de_orig)

stab_summary <- data.frame(
  contrast = names(de_orig),
  ORIG = sapply(de_orig, length),
  R1   = sapply(de_R1, length),
  R2   = sapply(de_R2, length),
  R3   = sapply(de_R3, length),
  ROBUST_intersect = sapply(robust, length)
)
print(stab_summary)
write.csv(stab_summary, file.path(OUT, "lncRNA_filter_strategy_summary.csv"), row.names = FALSE)
saveRDS(robust, file.path(OUT, "robust_lncRNA_DE_intersect.rds"))

message("[05] DONE. Output in ", OUT)
