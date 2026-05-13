#!/usr/bin/env Rscript
# Pure RUV analysis: no site labels anywhere
# Only AA-RNA-1 <-> AA-LES technical replicate identity tells RUV about batch
suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR); library(limma)
  library(readr); library(dplyr); library(ggplot2)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "ruv_no_site_label")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"; rownames(cm) <- raw$EnsemblID

# 15 samples: 14 manuscript + AA-LES (technical replicate of AA-RNA-1)
s_use <- c("AA-PRO","AA-KEW",
           "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
           "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
           "AA-HMH","AA-PJH","AA-LES")
mat <- cm[, s_use]

# Meta has ONLY group and patient info, NO site label
meta <- data.frame(
  sample=s_use,
  group=factor(ifelse(s_use %in% c("AA-PRO","AA-KEW"),"Control",
               ifelse(s_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"),"g_BMF","u_BMF")),
               levels=c("Control","g_BMF","u_BMF")),
  patient=c("PRO","KEW","DKC","FA","FA2","FA3","LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES"),
  row.names=s_use
)
cat("=== Sample design (NO site labels) ===\n"); print(table(meta$group))
cat("Technical replicate pair (patient=LES): ", which(meta$patient=="LES"), "\n\n")

# RUVs replicate group matrix — only LES pair
gmat <- matrix(-1, nrow=15, ncol=2)
for (i in 1:15) gmat[i, 1] <- i
les_idx <- which(meta$patient == "LES")
gmat[les_idx[1], 2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]
cat("RUV replicate matrix (rows=replicate groups, cols=sample indices, -1=padding):\n")
print(gmat)

# Pre-filter
keep <- rowSums(mat) >= 10
mat_f <- mat[keep, ]
cat(sprintf("\nGenes after pre-filter: %d\n", nrow(mat_f)))

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

set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=1, scIdx=gmat)
W <- pData(set_ruv)[, "W_1", drop=FALSE]
cat("\nW_1 (RUV-estimated batch factor, learned only from LES replicate):\n")
print(round(W, 3))

# DESeq2 on 14 samples (drop AA-LES from DE to avoid double-counting)
s_de <- setdiff(s_use, "AA-LES")
mat_de <- mat_f[, s_de]
meta_de <- meta[s_de, ]
W_de <- W[s_de, , drop=FALSE]
meta_dds <- data.frame(meta_de, W_1=W_de[,1])

# Design = ~ W_1 + group ONLY. No site, no cohort, no extra covariates.
dds <- DESeqDataSetFromMatrix(mat_de, meta_dds, design = ~ W_1 + group)
dds <- DESeq(dds, quiet=TRUE)
saveRDS(dds, file.path(OUTDIR, "dds.rds"))

# Group-aware filter
cpm_raw <- cpm(mat_de)
group_aware <- function(cm, gr) {
  pass <- rep(FALSE, nrow(cm))
  for (g in unique(gr)) { idx <- which(gr==g)
    pass <- pass | (rowSums(cm[, idx, drop=FALSE] >= 1) >= ceiling(length(idx)*0.5)) }
  pass
}
fp <- rownames(cpm_raw)[group_aware(cpm_raw, as.character(meta_de$group))]
cat(sprintf("\nGenes passing group-aware filter (CPM>=1 in >=50%% of any group): %d\n", length(fp)))

dump <- function(label, contrast) {
  r <- results(dds, contrast=contrast, alpha=0.05)
  d <- as.data.frame(r) %>% tibble::rownames_to_column("EnsemblID") %>%
    left_join(ann, by="EnsemblID") %>% filter(EnsemblID %in% fp) %>% arrange(padj)
  write_tsv(d, file.path(OUTDIR, sprintf("DE_%s.tsv", label)))
  n_sig <- sum(d$padj<0.05 & abs(d$log2FoldChange)>1, na.rm=TRUE)
  n_lnc <- sum(d$padj<0.05 & abs(d$log2FoldChange)>1 & d$GeneType=="lncRNA", na.rm=TRUE)
  cat(sprintf("\n  %s : n_DE=%d (lncRNA=%d)\n", label, n_sig, n_lnc))
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
  rows <- list()
  for (g in ms11) {
    r <- de %>% filter(GeneSymbol==g) %>% slice(1)
    if (nrow(r)==0) next
    sig <- !is.na(r$padj) && r$padj<0.05 && abs(r$log2FoldChange)>1
    cat(sprintf("  %-12s LFC=%6.2f padj=%.2e %s\n", g, r$log2FoldChange, r$padj, ifelse(sig,"SIG","ns")))
    rows[[g]] <- data.frame(gene=g, LFC=r$log2FoldChange, padj=r$padj, sig=sig)
  }
  bind_rows(rows)
}
a1 <- audit(g_vs_c, "g-BMF vs Ctrl (NO site labels, only LES replicate)")
a2 <- audit(u_vs_c, "u-BMF vs Ctrl (NO site labels, only LES replicate)")
write_tsv(bind_rows(a1 %>% mutate(contrast="gBMF_vs_Ctrl"),
                    a2 %>% mutate(contrast="uBMF_vs_Ctrl")),
          file.path(OUTDIR, "ms11_audit_no_site.tsv"))

# PCA before/after
vsd <- vst(dds, blind=TRUE)
mat_vst <- assay(vsd)
mat_after <- limma::removeBatchEffect(mat_vst, covariates=as.matrix(W_de))

# Plot with KNOWN site labels for interpretation (not used in analysis)
known_site <- ifelse(s_de %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH"), "jinpyung", "macrogen")
plot_pca <- function(m, lbl, f) {
  pca <- prcomp(t(m))
  pct <- round(summary(pca)$importance[2,1:2]*100, 1)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=rownames(pca$x),
                   known_site=known_site, group=meta_de$group)
  ggsave(f, ggplot(df, aes(PC1,PC2,color=known_site,shape=group,label=sample))+
           geom_point(size=4)+geom_text(vjust=-1,size=3)+
           labs(title=lbl, x=sprintf("PC1 (%.1f%%)",pct[1]), y=sprintf("PC2 (%.1f%%)",pct[2]),
                color="actual site (not used)")+
           theme_bw(), width=8, height=6, dpi=150)
}
plot_pca(mat_vst,   "Before correction (raw VST)",                file.path(OUTDIR, "PCA_before.pdf"))
plot_pca(mat_after, "After removing W_1 (RUV from LES replicate)", file.path(OUTDIR, "PCA_after.pdf"))

message("\n[DONE] All outputs in ", OUTDIR)
