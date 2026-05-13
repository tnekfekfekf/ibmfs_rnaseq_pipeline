#!/usr/bin/env Rscript
# Visualize batch correction effect on PCA for 17-sample analysis
suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(readr); library(dplyr); library(ggplot2); library(patchwork)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "pca_correction_comparison")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
cm  <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"
rownames(cm) <- raw$EnsemblID

s17 <- c("Child1","Child2","Child3","AA-PRO","AA-KEW",
         "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
         "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
         "AA-HMH","AA-PJH")
mat <- cm[, s17]
meta <- data.frame(
  sample=s17,
  site=factor(ifelse(s17 %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH"), "jinpyung",
               ifelse(s17 %in% c("Child1","Child2","Child3"), "public", "macrogen")),
              levels=c("macrogen","jinpyung","public")),
  cohort=factor(ifelse(s17 %in% c("Child1","Child2","Child3"), "public_MNC","internal_aspirate"),
                levels=c("internal_aspirate","public_MNC")),
  group=factor(ifelse(s17 %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
                ifelse(s17 %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF","u_BMF")),
              levels=c("Control","g_BMF","u_BMF")),
  row.names=s17
)

keep <- rowSums(mat) >= 10
dds <- DESeqDataSetFromMatrix(mat[keep,], meta, design= ~ 1)
vsd <- vst(dds, blind=TRUE)
mat_vst <- assay(vsd)

# 3 versions
mat_raw   <- mat_vst
mat_2lev  <- limma::removeBatchEffect(mat_vst, batch=meta$cohort,
                                       design=model.matrix(~ group, meta))
mat_3lev  <- limma::removeBatchEffect(mat_vst, batch=meta$site,
                                       design=model.matrix(~ group, meta))

make_pca_df <- function(mat, label) {
  p <- prcomp(t(mat))
  pct <- round(summary(p)$importance[2, 1:2]*100, 1)
  df <- data.frame(PC1=p$x[,1], PC2=p$x[,2], sample=rownames(p$x),
                   site=meta$site, group=meta$group, cohort=meta$cohort,
                   label=label, pct1=pct[1], pct2=pct[2])
  df
}

dfs <- bind_rows(
  make_pca_df(mat_raw,  "1) Raw VST (no correction)"),
  make_pca_df(mat_2lev, "2) After removing cohort (internal vs public, 2-level)"),
  make_pca_df(mat_3lev, "3) After removing site (macrogen / jinpyung / public, 3-level)")
)

# Make 3 plots side by side
plot_one <- function(d) {
  ttl <- d$label[1]
  pct1 <- d$pct1[1]; pct2 <- d$pct2[1]
  ggplot(d, aes(PC1, PC2, color=site, shape=group, label=sample)) +
    geom_point(size=4, alpha=0.85) +
    geom_text(vjust=-1, size=2.7, show.legend=FALSE) +
    scale_color_manual(values=c(macrogen="#1f77b4", jinpyung="#d62728", public="#2ca02c")) +
    labs(title=ttl,
         x=sprintf("PC1 (%.1f%%)", pct1),
         y=sprintf("PC2 (%.1f%%)", pct2)) +
    theme_bw(base_size=11) +
    theme(legend.position="bottom")
}

p1 <- plot_one(dfs %>% filter(label=="1) Raw VST (no correction)"))
p2 <- plot_one(dfs %>% filter(label=="2) After removing cohort (internal vs public, 2-level)"))
p3 <- plot_one(dfs %>% filter(label=="3) After removing site (macrogen / jinpyung / public, 3-level)"))

combined <- p1 + p2 + p3 + plot_layout(ncol=1)
ggsave(file.path(OUTDIR, "PCA_17sample_3panels.pdf"),
       combined, width=8, height=18, dpi=150)
ggsave(file.path(OUTDIR, "PCA_17sample_3panels.png"),
       combined, width=8, height=18, dpi=150)

# Quantify: site separation distance on PC1+PC2
for (lbl in unique(dfs$label)) {
  sub <- dfs %>% filter(label==lbl)
  centroids <- sub %>% group_by(site) %>% summarise(c1=mean(PC1), c2=mean(PC2), .groups="drop")
  d_pubvs_int <- sqrt((centroids$c1[centroids$site=="public"] - mean(centroids$c1[centroids$site!="public"]))^2 +
                       (centroids$c2[centroids$site=="public"] - mean(centroids$c2[centroids$site!="public"]))^2)
  d_macvjin   <- sqrt((centroids$c1[centroids$site=="macrogen"] - centroids$c1[centroids$site=="jinpyung"])^2 +
                       (centroids$c2[centroids$site=="macrogen"] - centroids$c2[centroids$site=="jinpyung"])^2)
  cat(sprintf("%s\n  centroid distance: public vs internal = %.2f, macrogen vs jinpyung = %.2f\n\n",
              lbl, d_pubvs_int, d_macvjin))
}

message("[DONE] Plot: ", file.path(OUTDIR, "PCA_17sample_3panels.pdf"))
