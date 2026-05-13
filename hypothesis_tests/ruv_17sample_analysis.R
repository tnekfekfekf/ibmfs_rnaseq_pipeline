#!/usr/bin/env Rscript
# =====================================================================
# RUV-based batch correction using cross-site technical replicates
# AA-RNA-1 (Macrogen) <-> AA-LES (Jinpyung) = same patient
# Then DESeq2 with ~ W + group on 17-sample raw counts
# =====================================================================

suppressPackageStartupMessages({
  library(RUVSeq); library(DESeq2); library(edgeR)
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "ruv_17sample_analysis")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"
rownames(cm) <- raw$EnsemblID

# Sample design — include the technical replicate AA-LES for RUVs
# Final analysis cohort = 17 samples (Child + 14 manuscript)
# Plus AA-LES included as a replicate of AA-RNA-1 to inform RUVs about batch
samples_use <- c(
  "Child1","Child2","Child3",                              # public controls (cohort=public)
  "AA-PRO","AA-KEW",                                       # internal controls (cohort=jinpyung)
  "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",      # g-BMF (Macrogen)
  "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13",
  "AA-RNA-18","AA-RNA-16",                                  # u-BMF (Macrogen)
  "AA-HMH","AA-PJH",                                        # u-BMF (Jinpyung)
  "AA-LES"                                                  # technical replicate of AA-RNA-1 at Jinpyung
)
mat <- cm[, samples_use]

meta <- data.frame(
  sample = samples_use,
  site = factor(
    ifelse(samples_use %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH","AA-LES"), "jinpyung",
           ifelse(samples_use %in% c("Child1","Child2","Child3"), "public", "macrogen")),
    levels = c("macrogen","jinpyung","public")
  ),
  group = factor(
    ifelse(samples_use %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
           ifelse(samples_use %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF", "u_BMF")),
    levels = c("Control","g_BMF","u_BMF")
  ),
  patient = c("Child1","Child2","Child3","PRO","KEW",
              "DKC","FA","FA2","FA3",
              "LES","RNA4","RNA5","RNA13","CSB","RNA16",
              "HMH","PJH",
              "LES"),  # AA-RNA-1 and AA-LES share patient ID "LES"
  row.names = samples_use
)
cat("=== Sample design ===\n"); print(table(meta$site, meta$group))
cat("\nPatient overlaps (>=2 samples = replicates):\n"); print(table(meta$patient))

# RUVs needs technical-replicate group matrix
# Each row = a "biological replicate group", entries are sample indices, -1 padding
patients_with_reps <- names(table(meta$patient))[table(meta$patient) > 1]
cat(sprintf("\nReplicate patient groups (n=%d): %s\n",
            length(patients_with_reps), paste(patients_with_reps, collapse=", ")))

# makeGroups: matrix where each row is one replicate group
make_groups <- function(meta) {
  groups <- list()
  for (p in unique(meta$patient)) {
    idx <- which(meta$patient == p)
    if (length(idx) >= 1) groups[[length(groups)+1]] <- idx
  }
  maxL <- max(sapply(groups, length))
  out <- matrix(-1, nrow=length(groups), ncol=maxL)
  for (i in seq_along(groups)) out[i, seq_along(groups[[i]])] <- groups[[i]]
  out
}
gmat <- make_groups(meta)
cat("\nGroup matrix for RUVs (-1 = padding):\n"); print(gmat)

# Pre-filter: standard rowSums>=10
keep <- rowSums(mat) >= 10
cat(sprintf("\nPre-filter rowSums>=10: %d / %d genes retained\n", sum(keep), length(keep)))
mat_f <- mat[keep, ]

# Step 1: Estimate empirical control genes (genes least likely to differ between groups)
# Use upper-quartile normalized data, then find low-variance/low-DE genes
set <- newSeqExpressionSet(mat_f, phenoData=meta)
set_uq <- betweenLaneNormalization(set, which="upper")

# Identify negative control genes: rank by p-value from naive edgeR test ~group, take bottom
design0 <- model.matrix(~group, data=pData(set_uq))
y <- DGEList(counts=counts(set_uq), group=meta$group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design0)
y <- estimateGLMTagwiseDisp(y, design0)
fit <- glmFit(y, design0)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set_uq))$table
empirical_neg <- rownames(top)[(nrow(top) - 4999):nrow(top)]   # bottom 5000
cat(sprintf("Empirical negative-control genes for RUVs: %d\n", length(empirical_neg)))

# Step 2: Run RUVs with technical replicates AND empirical negative controls
# k=2 requested; if only 1 replicate pair exists, RUVs may return fewer factors.
set_ruv <- RUVs(set_uq, cIdx=empirical_neg, k=2, scIdx=gmat)
W <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop=FALSE]
nW <- ncol(W)
cat(sprintf("\nRUVs returned %d W factor(s):\n", nW)); print(round(W, 3))

# Save corrected counts for visualization
norm_counts <- normCounts(set_ruv)
saveRDS(norm_counts, file.path(OUTDIR, "RUVs_normalized_counts.rds"))

# Step 3: DESeq2 with ~ W_1 (+W_2 if available) + group on RAW counts
meta_dds <- cbind(meta, W)
form_str <- paste("~", paste(c(colnames(W), "group"), collapse=" + "))
cat("DESeq2 design:", form_str, "\n")
dds <- DESeqDataSetFromMatrix(countData=mat_f, colData=meta_dds,
                              design= as.formula(form_str))
dds <- DESeq(dds, quiet=TRUE)
saveRDS(dds, file.path(OUTDIR, "dds_ruv.rds"))

cat("\n=== DESeq2 results ===\n")
dump <- function(label, contrast) {
  r <- results(dds, contrast=contrast, alpha=0.05)
  d <- as.data.frame(r) %>% tibble::rownames_to_column("EnsemblID") %>%
    left_join(ann, by="EnsemblID") %>% arrange(padj)
  write_tsv(d, file.path(OUTDIR, sprintf("DE_%s.tsv", label)))
  cat(sprintf("\n  %s : n_DE (padj<0.05, |LFC|>1) = %d (up=%d, down=%d)\n",
              label,
              sum(d$padj < 0.05 & abs(d$log2FoldChange) > 1, na.rm=TRUE),
              sum(d$padj < 0.05 & d$log2FoldChange >  1, na.rm=TRUE),
              sum(d$padj < 0.05 & d$log2FoldChange < -1, na.rm=TRUE)))
  d
}
g_vs_c <- dump("gBMF_vs_Ctrl", c("group","g_BMF","Control"))
u_vs_c <- dump("uBMF_vs_Ctrl", c("group","u_BMF","Control"))
g_vs_u <- dump("gBMF_vs_uBMF", c("group","g_BMF","u_BMF"))

# Step 4: Apply group-aware filter to RAW CPM
cpm_raw <- cpm(mat_f)
group_aware_pass <- function(cpm_mat, groups) {
  pass <- rep(FALSE, nrow(cpm_mat))
  for (g in unique(groups)) {
    idx <- which(groups == g)
    pass <- pass | (rowSums(cpm_mat[, idx, drop=FALSE] >= 1) >= ceiling(length(idx)*0.5))
  }
  pass
}
filter_pass <- group_aware_pass(cpm_raw, as.character(meta$group))
cat(sprintf("\nGroup-aware filter (CPM>=1 in >=50%% of any group): %d genes pass\n", sum(filter_pass)))

# Step 5: 11 manuscript lncRNA recovery
ms11 <- c("HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",
          "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1")
ms11_eid <- sapply(ms11, function(g) ann$EnsemblID[ann$GeneSymbol==g][1])
ms11_eid <- ms11_eid[!is.na(ms11_eid)]

de_lnc_audit <- function(de_df, lbl) {
  cat(sprintf("\n  --- %s : 11 manuscript lncRNA recovery ---\n", lbl))
  for (g in ms11) {
    row <- de_df %>% filter(GeneSymbol == g) %>% slice(1)
    if (nrow(row) == 0) { cat(sprintf("    %-12s : NOT IN MATRIX\n", g)); next }
    in_filter <- ms11_eid[g] %in% rownames(cpm_raw)[filter_pass]
    sig <- !is.na(row$padj) && row$padj < 0.05 && abs(row$log2FoldChange) > 1
    cat(sprintf("    %-12s : LFC=%6.2f padj=%.3e  filter=%s  DE=%s\n",
                g, row$log2FoldChange, row$padj,
                ifelse(in_filter, "PASS","FAIL"),
                ifelse(sig, "SIG", "ns")))
  }
}
de_lnc_audit(g_vs_c, "g-BMF vs Ctrl")
de_lnc_audit(u_vs_c, "u-BMF vs Ctrl")

# Save filter info
write_tsv(
  data.frame(EnsemblID=rownames(cpm_raw), pass_group_aware_filter=filter_pass),
  file.path(OUTDIR, "filter_pass_group_aware.tsv")
)

# Step 6: PCA visualization before/after correction
vsd_before <- vst(dds, blind=TRUE)   # before correction (but already DE-fit, so use blind)
mat_before <- assay(vsd_before)
W_aligned <- as.matrix(W[colnames(mat_before), , drop=FALSE])
mat_after  <- limma::removeBatchEffect(mat_before, covariates=W_aligned)

plot_pca <- function(mat, lbl, file) {
  pca <- prcomp(t(mat))
  pct <- round(summary(pca)$importance[2, 1:2]*100, 1)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=rownames(pca$x),
                   site=meta$site, group=meta$group)
  p <- ggplot(df, aes(PC1, PC2, color=site, shape=group, label=sample)) +
    geom_point(size=4) + geom_text(vjust=-1, size=3) +
    labs(title=lbl, x=sprintf("PC1 (%.1f%%)", pct[1]), y=sprintf("PC2 (%.1f%%)", pct[2])) +
    theme_bw()
  ggsave(file, p, width=8, height=6, dpi=150)
  invisible(pca)
}
plot_pca(mat_before, "Before RUV correction (VST, blind)",
         file.path(OUTDIR, "PCA_before_RUV.pdf"))
plot_pca(mat_after,  "After RUV correction (W1, W2 removed)",
         file.path(OUTDIR, "PCA_after_RUV.pdf"))

message("\n[DONE] Outputs in ", OUTDIR)
