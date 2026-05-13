#!/usr/bin/env Rscript
# DE count comparison: manuscript vs our scenarios, ALL using rowSums(counts)>=10 filter
suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR)
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "de_count_comparison")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"; rownames(cm) <- raw$EnsemblID

run_scenario <- function(samples_use, design_formula, extra_covar=NULL, label) {
  mat <- cm[, samples_use]
  meta <- data.frame(
    sample=samples_use,
    cohort=factor(ifelse(samples_use %in% c("Child1","Child2","Child3"), "public", "internal"),
                  levels=c("internal","public")),
    group=factor(ifelse(samples_use %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
                  ifelse(samples_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF","u_BMF")),
                 levels=c("Control","g_BMF","u_BMF")),
    row.names=samples_use
  )

  # Filter: rowSums >= 10 (manuscript filter)
  keep <- rowSums(mat) >= 10
  mat_f <- mat[keep, ]

  meta_dds <- meta
  if (!is.null(extra_covar)) meta_dds <- cbind(meta_dds, extra_covar[samples_use, , drop=FALSE])

  dds <- DESeqDataSetFromMatrix(mat_f, meta_dds, design=design_formula)
  dds <- DESeq(dds, quiet=TRUE)

  get_counts <- function(contrast) {
    r <- results(dds, contrast=contrast, alpha=0.05)
    d <- as.data.frame(r) %>% tibble::rownames_to_column("EnsemblID") %>%
      left_join(ann, by="EnsemblID") %>%
      filter(padj < 0.05 & abs(log2FoldChange) > 1)
    list(
      total=nrow(d),
      pc   =sum(d$GeneType == "protein_coding", na.rm=TRUE),
      pc_up=sum(d$GeneType == "protein_coding" & d$log2FoldChange >  1, na.rm=TRUE),
      pc_dn=sum(d$GeneType == "protein_coding" & d$log2FoldChange < -1, na.rm=TRUE),
      lnc  =sum(d$GeneType == "lncRNA", na.rm=TRUE),
      lnc_up=sum(d$GeneType == "lncRNA" & d$log2FoldChange >  1, na.rm=TRUE),
      lnc_dn=sum(d$GeneType == "lncRNA" & d$log2FoldChange < -1, na.rm=TRUE)
    )
  }
  list(
    label=label,
    gBMF_vs_Ctrl=get_counts(c("group","g_BMF","Control")),
    uBMF_vs_Ctrl=get_counts(c("group","u_BMF","Control")),
    gBMF_vs_uBMF=get_counts(c("group","g_BMF","u_BMF"))
  )
}

# RUV W_1 for 14-sample
get_W1_14sample <- function() {
  s_use <- c("AA-PRO","AA-KEW","AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
             "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
             "AA-HMH","AA-PJH","AA-LES")
  mat <- cm[, s_use]
  patient <- c("PRO","KEW","DKC","FA","FA2","FA3","LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES")
  group <- factor(ifelse(s_use %in% c("AA-PRO","AA-KEW"),"Control",
                  ifelse(s_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"),"g_BMF","u_BMF")),
                  levels=c("Control","g_BMF","u_BMF"))
  gmat <- matrix(-1, nrow=15, ncol=2); for (i in 1:15) gmat[i,1] <- i
  les_idx <- which(patient=="LES")
  gmat[les_idx[1],2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]
  keep <- rowSums(mat) >= 10
  set <- newSeqExpressionSet(mat[keep,], phenoData=data.frame(group=group, row.names=s_use))
  set_uq <- betweenLaneNormalization(set, which="upper")
  design0 <- model.matrix(~ group)
  y <- DGEList(counts(set_uq), group=group)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design0); y <- estimateGLMTagwiseDisp(y, design0)
  fit <- glmFit(y, design0); lrt <- glmLRT(fit, coef=2)
  top <- topTags(lrt, n=nrow(set_uq))$table
  empirical_neg <- rownames(top)[(nrow(top)-4999):nrow(top)]
  set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=1, scIdx=gmat)
  W <- pData(set_ruv)[, "W_1", drop=FALSE]
  W
}

# RUV W_1 for 17-sample (including AA-LES)
get_W1_17sample <- function() {
  s_use <- c("Child1","Child2","Child3","AA-PRO","AA-KEW",
             "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
             "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
             "AA-HMH","AA-PJH","AA-LES")
  mat <- cm[, s_use]
  patient <- c("Child1","Child2","Child3","PRO","KEW","DKC","FA","FA2","FA3",
               "LES","RNA4","RNA5","RNA13","CSB","RNA16","HMH","PJH","LES")
  group <- factor(ifelse(s_use %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"),"Control",
                  ifelse(s_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"),"g_BMF","u_BMF")),
                  levels=c("Control","g_BMF","u_BMF"))
  gmat <- matrix(-1, nrow=18, ncol=2); for (i in 1:18) gmat[i,1] <- i
  les_idx <- which(patient=="LES")
  gmat[les_idx[1],2] <- les_idx[2]; gmat <- gmat[-les_idx[2], , drop=FALSE]
  keep <- rowSums(mat) >= 10
  set <- newSeqExpressionSet(mat[keep,], phenoData=data.frame(group=group, row.names=s_use))
  set_uq <- betweenLaneNormalization(set, which="upper")
  design0 <- model.matrix(~ group)
  y <- DGEList(counts(set_uq), group=group)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design0); y <- estimateGLMTagwiseDisp(y, design0)
  fit <- glmFit(y, design0); lrt <- glmLRT(fit, coef=2)
  top <- topTags(lrt, n=nrow(set_uq))$table
  empirical_neg <- rownames(top)[(nrow(top)-4999):nrow(top)]
  set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=1, scIdx=gmat)
  W <- pData(set_ruv)[, "W_1", drop=FALSE]
  W
}

W14 <- get_W1_14sample()
W17 <- get_W1_17sample()

# 4 scenarios
s14 <- c("AA-PRO","AA-KEW","AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
         "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
         "AA-HMH","AA-PJH")
s17 <- c("Child1","Child2","Child3", s14)

res1 <- run_scenario(s14, ~ group, NULL, "14-sample ~group (manuscript reproduction)")
res2 <- run_scenario(s14, ~ W_1 + group, W14[s14, , drop=FALSE], "14-sample + RUV (~W_1 + group)")
res3 <- run_scenario(s17, ~ cohort + group, NULL, "17-sample (~cohort + group)")
res4 <- run_scenario(s17, ~ W_1 + cohort + group, W17[s17, , drop=FALSE], "17-sample + RUV (~W_1 + cohort + group)")

# Manuscript numbers (from text and sTables)
ms_pub <- list(
  label="Manuscript published",
  gBMF_vs_Ctrl=list(total=NA, pc=2078, pc_up=1501, pc_dn=577, lnc=1169, lnc_up=NA, lnc_dn=NA),
  uBMF_vs_Ctrl=list(total=NA, pc=1315, pc_up=1201, pc_dn=114, lnc=994, lnc_up=NA, lnc_dn=NA),
  gBMF_vs_uBMF=list(total=NA, pc=4, pc_up=NA, pc_dn=NA, lnc=NA, lnc_up=NA, lnc_dn=NA)
)

# Build table
make_row <- function(res, contrast_name) {
  ct <- res[[contrast_name]]
  data.frame(
    scenario=res$label,
    contrast=contrast_name,
    n_DE_total=ct$total,
    n_PC=ct$pc, PC_up=ct$pc_up, PC_dn=ct$pc_dn,
    n_lncRNA=ct$lnc, lncRNA_up=ct$lnc_up, lncRNA_dn=ct$lnc_dn
  )
}

all_rows <- bind_rows(
  make_row(ms_pub, "gBMF_vs_Ctrl"), make_row(ms_pub, "uBMF_vs_Ctrl"), make_row(ms_pub, "gBMF_vs_uBMF"),
  make_row(res1, "gBMF_vs_Ctrl"),  make_row(res1, "uBMF_vs_Ctrl"),  make_row(res1, "gBMF_vs_uBMF"),
  make_row(res2, "gBMF_vs_Ctrl"),  make_row(res2, "uBMF_vs_Ctrl"),  make_row(res2, "gBMF_vs_uBMF"),
  make_row(res3, "gBMF_vs_Ctrl"),  make_row(res3, "uBMF_vs_Ctrl"),  make_row(res3, "gBMF_vs_uBMF"),
  make_row(res4, "gBMF_vs_Ctrl"),  make_row(res4, "uBMF_vs_Ctrl"),  make_row(res4, "gBMF_vs_uBMF")
)

cat("\n=========================================================\n")
cat("DE count comparison (filter: rowSums(counts)>=10, padj<0.05, |LFC|>1)\n")
cat("=========================================================\n")
print(all_rows, row.names=FALSE)
write_tsv(all_rows, file.path(OUTDIR, "DE_count_comparison_table.tsv"))

message("\n[DONE] Table saved to ", file.path(OUTDIR, "DE_count_comparison_table.tsv"))
