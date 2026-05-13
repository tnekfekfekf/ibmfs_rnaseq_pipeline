#!/usr/bin/env Rscript
# 14-sample RUV analysis (matching manuscript cohort) using LES replicate
suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR)
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})
ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "ruv_14sample_analysis")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"; rownames(cm) <- raw$EnsemblID

samples_use <- c(
  "AA-PRO","AA-KEW","AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
  "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
  "AA-HMH","AA-PJH","AA-LES"
)
mat <- cm[, samples_use]
meta <- data.frame(
  sample=samples_use,
  site=factor(ifelse(samples_use %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH","AA-LES"),
                     "jinpyung","macrogen"), levels=c("macrogen","jinpyung")),
  group=factor(ifelse(samples_use %in% c("AA-PRO","AA-KEW"), "Control",
                ifelse(samples_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"),
                       "g_BMF","u_BMF")), levels=c("Control","g_BMF","u_BMF")),
  patient=c("PRO","KEW","DKC","FA","FA2","FA3","LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES"),
  row.names=samples_use
)
cat("=== Design (14-sample + LES replicate) ===\n"); print(table(meta$site, meta$group))

# RUVs replicate group matrix
gmat <- matrix(-1, nrow=nrow(meta), ncol=2)
for (i in seq_len(nrow(meta))) gmat[i, 1] <- i
les_idx <- which(meta$patient == "LES")
gmat[les_idx[1], 2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]

# Pre-filter
keep <- rowSums(mat) >= 10
mat_f <- mat[keep, ]
cat(sprintf("Genes after pre-filter: %d\n", nrow(mat_f)))

# RUVs
set <- newSeqExpressionSet(mat_f, phenoData=meta)
set_uq <- betweenLaneNormalization(set, which="upper")
design0 <- model.matrix(~group, data=pData(set_uq))
y <- DGEList(counts=counts(set_uq), group=meta$group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design0); y <- estimateGLMTagwiseDisp(y, design0)
fit <- glmFit(y, design0); lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set_uq))$table
empirical_neg <- rownames(top)[(nrow(top)-4999):nrow(top)]

set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=2, scIdx=gmat)
W <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop=FALSE]
cat(sprintf("\nRUVs %d W factors:\n", ncol(W))); print(round(W,3))

# Drop the AA-LES replicate from DE analysis (use it only for batch estimation)
samples_de <- setdiff(samples_use, "AA-LES")
mat_de <- mat_f[, samples_de]
meta_de <- meta[samples_de, ]
W_de <- W[samples_de, , drop=FALSE]
meta_dds <- cbind(meta_de, W_de)
form <- as.formula(paste("~", paste(c(colnames(W_de), "group"), collapse=" + ")))
cat("\nDESeq2 design:", deparse(form), "\n")

dds <- DESeqDataSetFromMatrix(countData=mat_de, colData=meta_dds, design=form)
dds <- DESeq(dds, quiet=TRUE)
saveRDS(dds, file.path(OUTDIR, "dds_ruv_14sample.rds"))

dump <- function(label, contrast) {
  r <- results(dds, contrast=contrast, alpha=0.05)
  d <- as.data.frame(r) %>% tibble::rownames_to_column("EnsemblID") %>%
    left_join(ann, by="EnsemblID") %>% arrange(padj)
  write_tsv(d, file.path(OUTDIR, sprintf("DE_%s.tsv", label)))
  cat(sprintf("\n  %s : n_DE = %d (up=%d, down=%d)\n", label,
              sum(d$padj<0.05 & abs(d$log2FoldChange)>1, na.rm=TRUE),
              sum(d$padj<0.05 & d$log2FoldChange> 1, na.rm=TRUE),
              sum(d$padj<0.05 & d$log2FoldChange< -1, na.rm=TRUE)))
  d
}
g_vs_c <- dump("gBMF_vs_Ctrl", c("group","g_BMF","Control"))
u_vs_c <- dump("uBMF_vs_Ctrl", c("group","u_BMF","Control"))
g_vs_u <- dump("gBMF_vs_uBMF", c("group","g_BMF","u_BMF"))

# 11 lncRNA audit
ms11 <- c("HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",
          "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1")
audit <- function(de, lbl) {
  cat(sprintf("\n--- %s ---\n", lbl))
  for (g in ms11) {
    r <- de %>% filter(GeneSymbol==g) %>% slice(1)
    if (nrow(r)==0) next
    sig <- !is.na(r$padj) && r$padj<0.05 && abs(r$log2FoldChange)>1
    cat(sprintf("  %-12s LFC=%6.2f padj=%.3e %s\n", g, r$log2FoldChange, r$padj,
                ifelse(sig,"SIG","ns")))
  }
}
audit(g_vs_c, "g-BMF vs Ctrl (14-sample RUV)")
audit(u_vs_c, "u-BMF vs Ctrl (14-sample RUV)")

# PCA before/after
vsd <- vst(dds, blind=TRUE)
mat_before <- assay(vsd)
W_aligned <- as.matrix(W_de[colnames(mat_before), , drop=FALSE])
mat_after <- limma::removeBatchEffect(mat_before, covariates=W_aligned)

plot_pca <- function(m, lbl, f) {
  pca <- prcomp(t(m))
  pct <- round(summary(pca)$importance[2,1:2]*100, 1)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=rownames(pca$x),
                   site=meta_de$site, group=meta_de$group)
  ggsave(f, ggplot(df, aes(PC1,PC2,color=site,shape=group,label=sample))+
           geom_point(size=4)+geom_text(vjust=-1,size=3)+
           labs(title=lbl, x=sprintf("PC1 (%.1f%%)",pct[1]), y=sprintf("PC2 (%.1f%%)",pct[2]))+
           theme_bw(), width=8, height=6, dpi=150)
}
plot_pca(mat_before, "Before RUV (VST blind)", file.path(OUTDIR,"PCA_before_RUV_14sample.pdf"))
plot_pca(mat_after,  "After RUV correction",   file.path(OUTDIR,"PCA_after_RUV_14sample.pdf"))

message("\n[DONE] Outputs in ", OUTDIR)
