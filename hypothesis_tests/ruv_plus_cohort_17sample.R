#!/usr/bin/env Rscript
# 17-sample analysis with both RUV (for within-internal site batch via LES replicate)
# and cohort label (for internal vs public).
# Design: ~ W_1 + cohort + group
suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR); library(limma)
  library(readr); library(dplyr); library(ggplot2)
})
ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "ruv_plus_cohort_17sample")
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
  cohort=factor(ifelse(s_use %in% c("Child1","Child2","Child3"), "public", "internal"),
                levels=c("internal","public")),
  group=factor(ifelse(s_use %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
                ifelse(s_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF", "u_BMF")),
               levels=c("Control","g_BMF","u_BMF")),
  patient=c("Child1","Child2","Child3","PRO","KEW","DKC","FA","FA2","FA3",
            "LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES"),
  row.names=s_use
)
cat("Design:\n"); print(table(meta$cohort, meta$group))

# RUVs: LES replicate only
gmat <- matrix(-1, nrow=18, ncol=2)
for (i in 1:18) gmat[i, 1] <- i
les_idx <- which(meta$patient == "LES")
gmat[les_idx[1], 2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]

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
cat("\nW_1:\n"); print(round(W,3))

# Drop AA-LES from DE
s_de <- setdiff(s_use, "AA-LES")
mat_de <- mat_f[, s_de]; meta_de <- meta[s_de, ]; W_de <- W[s_de, , drop=FALSE]
meta_dds <- data.frame(meta_de, W_1=W_de[,1])

# Design: W_1 (within-internal site) + cohort (internal vs public) + group
dds <- DESeqDataSetFromMatrix(mat_de, meta_dds, design = ~ W_1 + cohort + group)
dds <- DESeq(dds, quiet=TRUE)

cpm_raw <- cpm(mat_de)
group_aware <- function(cm, gr) {
  pass <- rep(FALSE, nrow(cm))
  for (g in unique(gr)) { idx <- which(gr==g)
    pass <- pass | (rowSums(cm[, idx, drop=FALSE]>=1) >= ceiling(length(idx)*0.5)) }
  pass
}
fp <- rownames(cpm_raw)[group_aware(cpm_raw, as.character(meta_de$group))]
cat(sprintf("\nGenes passing group-aware filter: %d\n", length(fp)))

ms11 <- c("HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",
          "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1")

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

audit <- function(de, lbl) {
  cat(sprintf("\n--- %s ---\n", lbl))
  for (g in ms11) {
    r <- de %>% filter(GeneSymbol==g) %>% slice(1)
    if (nrow(r)==0) next
    sig <- !is.na(r$padj) && r$padj<0.05 && abs(r$log2FoldChange)>1
    cat(sprintf("  %-12s LFC=%6.2f padj=%.2e %s\n", g, r$log2FoldChange, r$padj, ifelse(sig,"SIG","ns")))
  }
}
audit(g_vs_c, "g-BMF vs Ctrl (17s, ~W_1+cohort+group)")
audit(u_vs_c, "u-BMF vs Ctrl (17s, ~W_1+cohort+group)")

message("\n[DONE] All outputs in ", OUTDIR)
