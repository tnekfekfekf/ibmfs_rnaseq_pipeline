#!/usr/bin/env Rscript
# Hybrid analysis PCA: before and after ComBat-seq batch correction
# Combined matrix: 14 manuscript counts + 3 Child v3 counts
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap); library(RColorBrewer)
  library(sva); library(dplyr); library(readr); library(ggrepel); library(matrixStats)
})

MS_MATRIX <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt"
CHILD_V3_DIR <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3"
OUT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/hybrid_pca"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ===== Build hybrid count matrix =====
cat("Building hybrid matrix...\n")
ms <- read_tsv(MS_MATRIX, show_col_types=FALSE)
samps_14 <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
              "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18",
              "AA-HMH","AA-PJH")
ms14 <- as.matrix(ms[, samps_14]); rownames(ms14) <- ms$EnsemblID

# Load Child v3 counts
child_list <- list()
for (s in c("Child1","Child2","Child3")) {
  d <- read_tsv(file.path(CHILD_V3_DIR, paste0(s, ".counts.txt")), comment="#", show_col_types=FALSE)
  child_list[[s]] <- setNames(d[[7]], d$Geneid)
}
common <- Reduce(intersect, c(list(rownames(ms14)), lapply(child_list, names)))
ms14_c <- ms14[common, ]
child_mat <- do.call(cbind, lapply(child_list, function(x) x[common]))
colnames(child_mat) <- c("Child1","Child2","Child3")
counts <- cbind(ms14_c, child_mat)
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts) > 0, ]
cat(sprintf("Hybrid matrix: %d genes × %d samples\n", nrow(counts), ncol(counts)))

# Metadata
meta <- data.frame(
  sample_id = colnames(counts),
  group = factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8),rep("Control",3)),
                 levels=c("Control","U-AA","G-AA")),
  cohort = factor(c(rep("internal_MS",14), rep("public_v3",3)),
                  levels=c("internal_MS","public_v3")),
  control_type = factor(c(rep("G-AA",4),rep("internal_ctrl",2),rep("U-AA",8),rep("public_ctrl",3)),
                        levels=c("internal_ctrl","public_ctrl","U-AA","G-AA"))
)
rownames(meta) <- colnames(counts)

# ===== PCA before correction =====
cat("\n[1/2] VST + PCA — BEFORE batch correction...\n")
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ group)
vsd <- vst(dds, blind=TRUE)
mat_raw <- assay(vsd)

pca_plot <- function(mat, meta, ntop=2000, title="") {
  rv <- rowVars(mat); sel <- order(rv, decreasing=TRUE)[1:min(ntop, length(rv))]
  pca <- prcomp(t(mat[sel,]))
  pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=rownames(pca$x))
  d <- merge(d, meta, by.x="sample", by.y="sample_id")
  ggplot(d, aes(PC1, PC2, color=control_type, shape=cohort)) +
    geom_point(size=5, alpha=0.85) +
    geom_text_repel(aes(label=sample), size=3, max.overlaps=20) +
    scale_color_manual(values=c("internal_ctrl"="#2ECC71","public_ctrl"="#27AE60",
                                 "U-AA"="#F39C12","G-AA"="#E74C3C")) +
    scale_shape_manual(values=c("internal_MS"=16,"public_v3"=17)) +
    xlab(sprintf("PC1 (%.1f%%)", pv[1])) + ylab(sprintf("PC2 (%.1f%%)", pv[2])) +
    ggtitle(title) +
    theme_bw(12) + theme(plot.title=element_text(face="bold", size=14))
}

g1 <- pca_plot(mat_raw, meta, 2000, "PCA BEFORE batch correction\n(Hybrid: 14 manuscript + 3 Child)")
ggsave(file.path(OUT, "hybrid_PCA_before.pdf"), g1, width=10, height=7)
ggsave(file.path(OUT, "hybrid_PCA_before.png"), g1, width=10, height=7, dpi=150)
cat(sprintf("  Saved: %s/hybrid_PCA_before.{pdf,png}\n", OUT))

# ===== ComBat-seq correction =====
cat("\n[2/2] ComBat-seq → VST → PCA — AFTER correction...\n")
adj <- ComBat_seq(counts, batch=as.character(meta$cohort), group=as.character(meta$group))
saveRDS(adj, file.path(OUT, "hybrid_counts_combatseq.rds"))

dds_adj <- DESeqDataSetFromMatrix(adj, meta, design = ~ group)
vsd_adj <- vst(dds_adj, blind=TRUE)
mat_adj <- assay(vsd_adj)

g2 <- pca_plot(mat_adj, meta, 2000, "PCA AFTER ComBat-seq batch correction\n(Hybrid: 14 manuscript + 3 Child)")
ggsave(file.path(OUT, "hybrid_PCA_after.pdf"), g2, width=10, height=7)
ggsave(file.path(OUT, "hybrid_PCA_after.png"), g2, width=10, height=7, dpi=150)
cat(sprintf("  Saved: %s/hybrid_PCA_after.{pdf,png}\n", OUT))

# ===== ANOVA quantification =====
cat("\nANOVA — variance explained per PC:\n")
anova_pc <- function(mat, meta, label) {
  rv <- rowVars(mat); sel <- order(rv, decreasing=TRUE)[1:2000]
  pca <- prcomp(t(mat[sel,]))
  pv <- pca$sdev^2 / sum(pca$sdev^2)
  pca_df <- data.frame(pca$x[, 1:5])
  pca_df$cohort <- meta$cohort
  pca_df$group <- meta$group
  cat(sprintf("\n%s:\n", label))
  cat(sprintf("%-5s %8s %12s %12s\n", "PC", "var_pct", "p_cohort", "p_group"))
  for (pc in paste0("PC", 1:5)) {
    var_pct <- pv[as.numeric(sub("PC","",pc))]*100
    p_c <- summary(aov(as.formula(paste(pc, "~ cohort")), data=pca_df))[[1]][1, "Pr(>F)"]
    p_g <- summary(aov(as.formula(paste(pc, "~ group")), data=pca_df))[[1]][1, "Pr(>F)"]
    cat(sprintf("%-5s %7.1f%%  %12.3g  %12.3g\n", pc, var_pct, p_c, p_g))
  }
}
anova_pc(mat_raw, meta, "BEFORE correction")
anova_pc(mat_adj, meta, "AFTER correction")

cat(sprintf("\nDone. Output: %s\n", OUT))
list.files(OUT)
