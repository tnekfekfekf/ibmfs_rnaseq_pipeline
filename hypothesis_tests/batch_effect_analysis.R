#!/usr/bin/env Rscript
# Batch effect visualization & quantification
#   - PCA before/after ComBat-seq
#   - Sample-to-sample distance heatmap
#   - ANOVA: variance explained by cohort vs group
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap); library(RColorBrewer)
  library(sva); library(dplyr); library(readr); library(ggrepel); library(matrixStats)
})

COUNTS <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_v3_17samples.txt"
OUT    <- "/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/control_comparison/batch_effect"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

raw <- read_tsv(COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(raw[, -(1:2)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"

samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
           "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH",
           "Child1","Child2","Child3")
counts <- counts[, samps]
counts <- counts[rowSums(counts) > 0, ]

meta <- data.frame(
  sample_id = samps,
  group = factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8),rep("Control",3)),
                 levels=c("Control","U-AA","G-AA")),
  cohort = factor(c(rep("internal",14), rep("public",3)), levels=c("internal","public")),
  control_type = factor(c(rep("G-AA",4),rep("internal_ctrl",2),rep("U-AA",8),rep("public_ctrl",3)),
                        levels=c("internal_ctrl","public_ctrl","U-AA","G-AA"))
)
rownames(meta) <- samps

# --- 1. VST + PCA before correction ---
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
    geom_point(size=4, alpha=0.85) +
    geom_text_repel(aes(label=sample), size=2.7, max.overlaps=20) +
    scale_color_manual(values=c("internal_ctrl"="#2ECC71","public_ctrl"="#27AE60",
                                 "U-AA"="#F39C12","G-AA"="#E74C3C")) +
    scale_shape_manual(values=c("internal"=16,"public"=17)) +
    xlab(sprintf("PC1 (%.1f%%)", pv[1])) + ylab(sprintf("PC2 (%.1f%%)", pv[2])) +
    ggtitle(title) + theme_bw(11) + theme(plot.title=element_text(face="bold"))
}

cat("Generating PCA before correction...\n")
g1 <- pca_plot(mat_raw, meta, 2000, "PCA — Before batch correction (top 2000 var, VST)")
ggsave(file.path(OUT, "PCA_before_correction.pdf"), g1, width=8, height=6)
ggsave(file.path(OUT, "PCA_before_correction.png"), g1, width=8, height=6, dpi=150)

# --- 2. ComBat-seq batch correction ---
cat("Running ComBat-seq (cohort = batch, group preserved)...\n")
adj <- ComBat_seq(counts, batch=as.character(meta$cohort), group=as.character(meta$group))
saveRDS(adj, file.path(OUT, "counts_combatseq.rds"))

# --- 3. PCA after correction ---
dds_adj <- DESeqDataSetFromMatrix(adj, meta, design = ~ group)
vsd_adj <- vst(dds_adj, blind=TRUE)
mat_adj <- assay(vsd_adj)

cat("Generating PCA after correction...\n")
g2 <- pca_plot(mat_adj, meta, 2000, "PCA — After ComBat-seq correction (top 2000 var, VST)")
ggsave(file.path(OUT, "PCA_after_combatseq.pdf"), g2, width=8, height=6)
ggsave(file.path(OUT, "PCA_after_combatseq.png"), g2, width=8, height=6, dpi=150)

# --- 4. Sample-to-sample distance heatmap ---
cat("Generating sample distance heatmap...\n")
sd_mat <- as.matrix(dist(t(mat_raw)))
rownames(sd_mat) <- colnames(sd_mat) <- paste0(meta$sample_id, "_", meta$control_type)
ann <- data.frame(group=meta$control_type, cohort=meta$cohort, row.names=rownames(sd_mat))
ann_colors <- list(
  group=c(internal_ctrl="#2ECC71", public_ctrl="#27AE60", `U-AA`="#F39C12", `G-AA`="#E74C3C"),
  cohort=c(internal="#3498DB", public="#9B59B6"))
pdf(file.path(OUT, "sample_distance_heatmap.pdf"), width=10, height=9)
pheatmap(sd_mat, clustering_distance_rows=as.dist(sd_mat),
         clustering_distance_cols=as.dist(sd_mat),
         annotation_row=ann, annotation_col=ann, annotation_colors=ann_colors,
         col=colorRampPalette(rev(brewer.pal(9,"Blues")))(255),
         main="Sample-to-sample VST distance (before correction)")
dev.off()

# --- 5. ANOVA: PC variance explained by cohort vs group ---
cat("ANOVA: variance per PC explained by cohort vs group...\n")
rv <- rowVars(mat_raw); sel <- order(rv, decreasing=TRUE)[1:2000]
pca <- prcomp(t(mat_raw[sel,]))
pv <- pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(pca$x[, 1:5])
pca_df$cohort <- meta$cohort
pca_df$group <- meta$group
res <- list()
for (pc in paste0("PC", 1:5)) {
  m <- aov(get(pc) ~ cohort, data=pca_df)
  r_cohort <- summary(m)[[1]][1, "Pr(>F)"]
  m <- aov(get(pc) ~ group, data=pca_df)
  r_group <- summary(m)[[1]][1, "Pr(>F)"]
  res[[pc]] <- data.frame(PC=pc, var_pct=round(pv[as.numeric(sub("PC","",pc))]*100,1),
                           p_cohort=signif(r_cohort,3), p_group=signif(r_group,3))
}
anova_tab <- do.call(rbind, res)
cat("\n=== Per-PC ANOVA (BEFORE correction) ===\n")
print(anova_tab)
write_tsv(anova_tab, file.path(OUT, "PC_anova_before.tsv"))

# Same after correction
rv <- rowVars(mat_adj); sel <- order(rv, decreasing=TRUE)[1:2000]
pca <- prcomp(t(mat_adj[sel,]))
pv <- pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(pca$x[, 1:5])
pca_df$cohort <- meta$cohort
pca_df$group <- meta$group
res <- list()
for (pc in paste0("PC", 1:5)) {
  m <- aov(get(pc) ~ cohort, data=pca_df)
  r_cohort <- summary(m)[[1]][1, "Pr(>F)"]
  m <- aov(get(pc) ~ group, data=pca_df)
  r_group <- summary(m)[[1]][1, "Pr(>F)"]
  res[[pc]] <- data.frame(PC=pc, var_pct=round(pv[as.numeric(sub("PC","",pc))]*100,1),
                           p_cohort=signif(r_cohort,3), p_group=signif(r_group,3))
}
anova_tab_adj <- do.call(rbind, res)
cat("\n=== Per-PC ANOVA (AFTER ComBat-seq) ===\n")
print(anova_tab_adj)
write_tsv(anova_tab_adj, file.path(OUT, "PC_anova_after.tsv"))

cat(sprintf("\n[DONE] Outputs in %s\n", OUT))
list.files(OUT)
