#!/usr/bin/env Rscript
# Pre-correction QC: PCA, hierarchical clustering, sample-distance heatmap.
# Goal: visualize the batch effect (internal aspirate vs public MNC) before correction.
# Reviewer 1 also wants this hierarchical clustering of full transcriptome.

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap); library(RColorBrewer)
  library(sva); library(dplyr); library(ggrepel)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
OUT  <- file.path(ROOT, "deseq2"); FIG <- file.path(ROOT, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

dds <- readRDS(file.path(OUT, "dds_full.rds"))
meta <- readRDS(file.path(OUT, "samples_meta.rds"))
gene_meta <- readRDS(file.path(OUT, "gene_meta.rds"))

message("[02] VST-transforming ", ncol(dds), " samples")
vsd <- vst(dds, blind = TRUE)            # blind=TRUE: don't use design for variance stabilization (pure QC)
mat <- assay(vsd)
saveRDS(vsd, file.path(OUT, "vsd_blind.rds"))

# ---------- PCA on top variable genes (matches typical RNA-seq QC) ----------
pca_plot <- function(mat, meta, ntop = 500, color_var = "group", shape_var = "cohort", title = "") {
  rv <- rowVars(mat); sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(mat[sel, ]))
  pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  d <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample = rownames(pca$x))
  d <- merge(d, meta, by.x = "sample", by.y = "sample_id")
  ggplot(d, aes_string("PC1","PC2", color = color_var, shape = shape_var)) +
    geom_point(size = 3.5) + geom_text_repel(aes(label = sample), size = 3) +
    xlab(sprintf("PC1 (%.1f%%)", pv[1])) + ylab(sprintf("PC2 (%.1f%%)", pv[2])) +
    ggtitle(title) + theme_bw(11)
}

ggsave(file.path(FIG, "pca_top500_uncorrected.pdf"),
       pca_plot(mat, meta, 500, "group_combined", "cohort", "PCA, top 500 variable genes (VST, no batch correction)"),
       width = 8, height = 6)
ggsave(file.path(FIG, "pca_top2000_uncorrected.pdf"),
       pca_plot(mat, meta, 2000, "group_combined", "cohort", "PCA, top 2000 variable genes (VST, no batch correction)"),
       width = 8, height = 6)

# ---------- Sample-to-sample distance ----------
sd <- as.matrix(dist(t(mat)))
rownames(sd) <- colnames(sd) <- paste0(meta$sample_id, "_", meta$group_combined, "_", meta$cohort)
ann <- data.frame(group = meta$group_combined, cohort = meta$cohort, row.names = rownames(sd))
pdf(file.path(FIG, "sample_distance_heatmap.pdf"), width = 9, height = 8)
pheatmap(sd, clustering_distance_rows = as.dist(sd),
         clustering_distance_cols = as.dist(sd),
         annotation_row = ann, annotation_col = ann,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distance (VST, blind)")
dev.off()

# ---------- Full-transcriptome hierarchical clustering (Reviewer 1) ----------
# Use all expressed genes (already pre-filtered in 01) on VST scale
hc_col <- hclust(dist(t(mat)), method = "average")
hc_row <- hclust(dist(mat[order(rowVars(mat), decreasing = TRUE)[1:2000], ]), method = "average")
pdf(file.path(FIG, "full_transcriptome_clustering.pdf"), width = 9, height = 8)
pheatmap(mat[order(rowVars(mat), decreasing = TRUE)[1:2000], ],
         scale = "row",
         clustering_method = "average",
         show_rownames = FALSE,
         annotation_col = ann,
         main = "Hierarchical clustering, top 2000 variable genes (VST)")
dev.off()

# ---------- ComBat-seq batch correction (sensitivity) ----------
# Adjust counts for cohort, preserving group_combined biology
message("[02] ComBat-seq batch adjustment (cohort -> internal_aspirate vs public_MNC)")
counts_int <- counts(dds)
adj_counts <- ComBat_seq(counts_int,
                         batch = as.character(meta$cohort),
                         group = as.character(meta$group_combined))
saveRDS(adj_counts, file.path(OUT, "counts_combatseq.rds"))

# Re-PCA after ComBat-seq
dds_adj <- DESeqDataSetFromMatrix(adj_counts, colData = meta, design = ~ group_combined)
vsd_adj <- vst(dds_adj, blind = TRUE)
mat_adj <- assay(vsd_adj)
ggsave(file.path(FIG, "pca_top500_combatseq.pdf"),
       pca_plot(mat_adj, meta, 500, "group_combined", "cohort", "PCA, top 500 variable genes (after ComBat-seq)"),
       width = 8, height = 6)
ggsave(file.path(FIG, "pca_top2000_combatseq.pdf"),
       pca_plot(mat_adj, meta, 2000, "group_combined", "cohort", "PCA, top 2000 variable genes (after ComBat-seq)"),
       width = 8, height = 6)
saveRDS(vsd_adj, file.path(OUT, "vsd_combatseq.rds"))

# ---------- Diagnostic: how much variance does cohort explain on each PC? ----------
rv <- rowVars(mat); sel <- order(rv, decreasing = TRUE)[1:2000]
pca <- prcomp(t(mat[sel, ]))
pcs <- pca$x[, 1:5]
diag_anova <- sapply(1:5, function(i) {
  m <- summary(aov(pcs[, i] ~ meta$cohort))[[1]]
  c(F = m[1,"F value"], p = m[1,"Pr(>F)"])
})
colnames(diag_anova) <- paste0("PC", 1:5)
write.csv(t(diag_anova), file.path(OUT, "PC_vs_cohort_anova.csv"))
message("[02] Variance-by-cohort diagnostic:")
print(round(t(diag_anova), 4))

message("[02] DONE. Figures in ", FIG)
