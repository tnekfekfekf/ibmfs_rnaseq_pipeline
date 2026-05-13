#!/usr/bin/env Rscript
# PCA visualization for the final winning design: 17-sample ~ W_1 + cohort + group
suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR); library(limma)
  library(readr); library(dplyr); library(ggplot2); library(patchwork)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "pca_final_visualization")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"; rownames(cm) <- raw$EnsemblID

s_use <- c("Child1","Child2","Child3","AA-PRO","AA-KEW",
           "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
           "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
           "AA-HMH","AA-PJH","AA-LES")
mat <- cm[, s_use]
meta <- data.frame(
  sample=s_use,
  site=factor(ifelse(s_use %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH","AA-LES"), "jinpyung",
               ifelse(s_use %in% c("Child1","Child2","Child3"), "public", "macrogen")),
              levels=c("macrogen","jinpyung","public")),
  cohort=factor(ifelse(s_use %in% c("Child1","Child2","Child3"), "public", "internal"),
                levels=c("internal","public")),
  group=factor(ifelse(s_use %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
                ifelse(s_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF","u_BMF")),
               levels=c("Control","g_BMF","u_BMF")),
  patient=c("Child1","Child2","Child3","PRO","KEW","DKC","FA","FA2","FA3",
            "LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES"),
  row.names=s_use
)

# RUVs with LES replicate
gmat <- matrix(-1, nrow=18, ncol=2); for (i in 1:18) gmat[i,1] <- i
les_idx <- which(meta$patient == "LES")
gmat[les_idx[1],2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]

keep <- rowSums(mat) >= 10
mat_f <- mat[keep, ]
set <- newSeqExpressionSet(mat_f, phenoData=meta)
set_uq <- betweenLaneNormalization(set, which="upper")
design0 <- model.matrix(~ group, data=pData(set_uq))
y <- DGEList(counts=counts(set_uq), group=meta$group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design0); y <- estimateGLMTagwiseDisp(y, design0)
fit <- glmFit(y, design0); lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set_uq))$table
empirical_neg <- rownames(top)[(nrow(top)-4999):nrow(top)]
set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=1, scIdx=gmat)
W <- pData(set_ruv)[, "W_1", drop=FALSE]

# Build DESeq2 for VST (just for normalization, not DE)
dds_norm <- DESeqDataSetFromMatrix(mat_f, meta, ~ 1)
vsd <- vst(dds_norm, blind=TRUE)
mat_vst <- assay(vsd)

# Three correction states
mat_raw     <- mat_vst
cohort_mat  <- model.matrix(~ cohort, meta)[, -1, drop=FALSE]
mat_W1only  <- limma::removeBatchEffect(mat_vst, covariates=as.matrix(W),
                                         design=model.matrix(~ group, meta))
mat_W1cohort <- limma::removeBatchEffect(mat_vst,
                                          covariates=cbind(as.matrix(W), cohort_mat),
                                          design=model.matrix(~ group, meta))

# PCA + plotting helper
make_pca_plot <- function(m, title_text) {
  p <- prcomp(t(m))
  pct <- round(summary(p)$importance[2, 1:2]*100, 1)
  df <- data.frame(PC1=p$x[,1], PC2=p$x[,2], sample=rownames(p$x),
                   site=meta$site, group=meta$group, cohort=meta$cohort)
  ggplot(df, aes(PC1, PC2, color=site, shape=group, label=sample)) +
    geom_point(size=4.5, alpha=0.85) +
    geom_text(vjust=-1, size=2.5, show.legend=FALSE) +
    scale_color_manual(values=c(macrogen="#1f77b4", jinpyung="#d62728", public="#2ca02c")) +
    scale_shape_manual(values=c(Control=16, g_BMF=17, u_BMF=15)) +
    labs(title=title_text,
         x=sprintf("PC1 (%.1f%%)", pct[1]),
         y=sprintf("PC2 (%.1f%%)", pct[2])) +
    theme_bw(base_size=11) +
    theme(legend.position="bottom",
          plot.title=element_text(size=10, face="bold"))
}

p1 <- make_pca_plot(mat_raw,
   "1) BEFORE correction (raw VST)\nthree sites separate strongly")
p2 <- make_pca_plot(mat_W1only,
   "2) AFTER W_1 only (RUV from LES replicate)\nmacrogen/jinpyung mixed, public still separate")
p3 <- make_pca_plot(mat_W1cohort,
   "3) AFTER W_1 + cohort (final winning design)\nall three sites mixed; biology (group) drives PC1/PC2")

combined <- p1 / p2 / p3
ggsave(file.path(OUTDIR, "PCA_final_3stages.pdf"), combined, width=8.5, height=20, dpi=150)
ggsave(file.path(OUTDIR, "PCA_final_3stages.png"), combined, width=8.5, height=20, dpi=150)

# Quantify
centroid_dist <- function(m, label) {
  p <- prcomp(t(m))
  df <- data.frame(PC1=p$x[,1], PC2=p$x[,2], site=meta$site, group=meta$group)
  ctr <- df %>% group_by(site) %>% summarise(c1=mean(PC1), c2=mean(PC2), .groups="drop")
  pairs <- combn(unique(as.character(ctr$site)), 2, simplify=FALSE)
  cat(sprintf("\n%s :\n", label))
  for (pair in pairs) {
    a <- ctr[ctr$site==pair[1], ]; b <- ctr[ctr$site==pair[2], ]
    cat(sprintf("  %s vs %s: %.2f\n", pair[1], pair[2],
                sqrt((a$c1-b$c1)^2 + (a$c2-b$c2)^2)))
  }
  ctr2 <- df %>% group_by(group) %>% summarise(c1=mean(PC1), c2=mean(PC2), .groups="drop")
  pairs2 <- combn(unique(as.character(ctr2$group)), 2, simplify=FALSE)
  for (pair in pairs2) {
    a <- ctr2[ctr2$group==pair[1], ]; b <- ctr2[ctr2$group==pair[2], ]
    cat(sprintf("  %s vs %s: %.2f\n", pair[1], pair[2],
                sqrt((a$c1-b$c1)^2 + (a$c2-b$c2)^2)))
  }
}
centroid_dist(mat_raw,      "Stage 1 (raw)")
centroid_dist(mat_W1only,   "Stage 2 (W_1 only)")
centroid_dist(mat_W1cohort, "Stage 3 (W_1 + cohort)")

message("\n[DONE] PDF: ", file.path(OUTDIR, "PCA_final_3stages.pdf"))
