#!/usr/bin/env Rscript
# 17-sample analysis with TWO covariates:
#   - W_1 from RUVs (site batch effect, learned via AA-RNA-1/AA-LES replicate)
#   - composition_PC1 (xCell-based cell composition: MNC vs whole BM)
# Goal: rescue RT-qPCR validated lncRNAs that were lost in plain 17-sample RUV.

suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR); library(limma)
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "ruv_celltype_17sample")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"; rownames(cm) <- raw$EnsemblID

samples_use <- c(
  "Child1","Child2","Child3","AA-PRO","AA-KEW",
  "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
  "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
  "AA-HMH","AA-PJH","AA-LES"
)
mat <- cm[, samples_use]
meta <- data.frame(
  sample=samples_use,
  site=factor(ifelse(samples_use %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH","AA-LES"), "jinpyung",
                     ifelse(samples_use %in% c("Child1","Child2","Child3"), "public", "macrogen")),
              levels=c("macrogen","jinpyung","public")),
  cohort=factor(ifelse(samples_use %in% c("Child1","Child2","Child3"), "public_MNC", "internal_aspirate"),
                levels=c("internal_aspirate","public_MNC")),
  group=factor(ifelse(samples_use %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
                ifelse(samples_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF", "u_BMF")),
              levels=c("Control","g_BMF","u_BMF")),
  patient=c("Child1","Child2","Child3","PRO","KEW","DKC","FA","FA2","FA3",
            "LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES"),
  row.names=samples_use
)

# ---- xCell scores (load and align) ----
xc <- read_tsv("/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/deconv/deconv_xCell.tsv",
               show_col_types=FALSE)
xc_mat <- as.matrix(xc[, -1]); rownames(xc_mat) <- xc$cell_type
# Use only the 17 samples we care about (AA-LES not in xCell, we'll handle via AA-RNA-1 duplicate)
xc_samples <- intersect(samples_use, colnames(xc_mat))
xc_mat_use <- xc_mat[, xc_samples]
cat("xCell samples available:", length(xc_samples), "/", length(samples_use), "\n")

# Filter to highly variable cell types (avoid noise)
cv <- apply(xc_mat_use, 1, function(x) sd(x)/(mean(x)+1e-6))
top_var <- head(rownames(xc_mat_use)[order(-cv)], 20)
xc_top <- xc_mat_use[top_var, ]
# Composition PCA across samples
pca <- prcomp(t(xc_top), scale.=TRUE)
comp_PC <- pca$x[, 1:2]
cat("\nVariance explained by xCell PCA:\n"); print(round(summary(pca)$importance[2, 1:3]*100, 1))
cat("\nPC1 PC2 of cell composition:\n"); print(round(comp_PC, 2))

# For AA-LES: use AA-RNA-1's composition score (same patient)
comp_full <- matrix(NA, nrow=length(samples_use), ncol=2,
                    dimnames=list(samples_use, c("comp_PC1","comp_PC2")))
for (s in samples_use) {
  if (s %in% rownames(comp_PC)) comp_full[s, ] <- comp_PC[s, ]
  else if (s == "AA-LES") comp_full[s, ] <- comp_PC["AA-RNA-1", ]
}
cat("\nFinal composition covariates per sample:\n"); print(round(comp_full, 2))

# ---- Pre-filter and RUVs (same as before) ----
keep <- rowSums(mat) >= 10
mat_f <- mat[keep, ]
gmat <- matrix(-1, nrow=18, ncol=2)
for (i in 1:18) gmat[i, 1] <- i
les_idx <- which(meta$patient == "LES")
gmat[les_idx[1], 2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]

set <- newSeqExpressionSet(mat_f, phenoData=meta)
set_uq <- betweenLaneNormalization(set, which="upper")
design0 <- model.matrix(~group, data=pData(set_uq))
y <- DGEList(counts=counts(set_uq), group=meta$group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design0); y <- estimateGLMTagwiseDisp(y, design0)
fit <- glmFit(y, design0); lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set_uq))$table
empirical_neg <- rownames(top)[(nrow(top)-4999):nrow(top)]
set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=1, scIdx=gmat)
W <- pData(set_ruv)[, "W_1", drop=FALSE]
cat("\nRUVs W_1:\n"); print(round(W, 3))

# ---- DESeq2 with ~ W_1 + comp_PC1 + comp_PC2 + group ----
# Drop AA-LES from DE to avoid double-counting patient
de_samples <- setdiff(samples_use, "AA-LES")
mat_de <- mat_f[, de_samples]
meta_de <- meta[de_samples, ]
W_de <- W[de_samples, , drop=FALSE]
comp_de <- comp_full[de_samples, , drop=FALSE]

meta_dds <- data.frame(meta_de, W_1=W_de[,1], comp_PC1=comp_de[,1], comp_PC2=comp_de[,2])
dds <- DESeqDataSetFromMatrix(countData=mat_de, colData=meta_dds,
                              design= ~ W_1 + comp_PC1 + comp_PC2 + group)
dds <- DESeq(dds, quiet=TRUE)
saveRDS(dds, file.path(OUTDIR, "dds_ruv_celltype.rds"))

cat("\n=== DESeq2 with RUV + cell composition ===\n")
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
audit(g_vs_c, "g-BMF vs Ctrl (17s + RUV + xCell)")
audit(u_vs_c, "u-BMF vs Ctrl (17s + RUV + xCell)")

# Visualization
vsd <- vst(dds, blind=TRUE)
mat_v <- assay(vsd)
covs <- as.matrix(meta_dds[, c("W_1","comp_PC1","comp_PC2")])
mat_after <- limma::removeBatchEffect(mat_v, covariates=covs)

plot_pca <- function(m, lbl, f) {
  pca <- prcomp(t(m))
  pct <- round(summary(pca)$importance[2,1:2]*100, 1)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=rownames(pca$x),
                   site=meta_de$site, group=meta_de$group, cohort=meta_de$cohort)
  ggsave(f, ggplot(df, aes(PC1,PC2,color=cohort,shape=group,label=sample))+
           geom_point(size=4)+geom_text(vjust=-1,size=3)+
           labs(title=lbl, x=sprintf("PC1 (%.1f%%)",pct[1]), y=sprintf("PC2 (%.1f%%)",pct[2]))+
           theme_bw(), width=8, height=6, dpi=150)
}
plot_pca(mat_v,      "Before W_1 + comp_PC correction (VST blind)",
         file.path(OUTDIR, "PCA_before.pdf"))
plot_pca(mat_after,  "After W_1 + comp_PC1 + comp_PC2 correction",
         file.path(OUTDIR, "PCA_after.pdf"))

message("\n[DONE] Outputs in ", OUTDIR)
