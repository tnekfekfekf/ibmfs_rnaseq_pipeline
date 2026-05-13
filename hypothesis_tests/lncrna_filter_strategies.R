#!/usr/bin/env Rscript
# =====================================================================
# Reviewer 1, Issue 1: lncRNA filtering strategy sensitivity analysis
# ---------------------------------------------------------------------
# Tests 4 filter strategies on TWO cohort configurations:
#   (A) 14-sample manuscript matrix      | design = ~ group
#   (B) 17-sample hybrid (+3 Child v3)   | design = ~ cohort + group
#
# Strategies:
#   ORIG : mean CPM >= 10 AND CPM >= 1 in ALL samples  (manuscript)
#   R1   : CPM >= 1 in >= 50% samples of ANY group     (group-aware moderate)
#   R2   : mean CPM >= 1 in any group AND >=50% of that group >= 0.5
#   R3   : no pre-filter (DESeq2 IF only)              (most permissive)
#
# Outputs:
#   /Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/lncrna_filter_analysis/
#     - filter_summary_14sample.tsv
#     - filter_summary_17sample.tsv
#     - filter_summary_combined.tsv
#     - manuscript_11lncRNA_recovery.tsv
#     - de_lncRNA_{cohort}_{contrast}_{strategy}.tsv  (per cell)
#     - robust_lncRNAs_intersection.tsv
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(edgeR)   # for cpm()
})

# ---------- paths ----------
MATRIX  <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt"
OUTDIR  <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/lncrna_filter_analysis"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------- sample plan ----------
ctrl_int  <- c("AA-PRO","AA-KEW")
ctrl_pub  <- c("Child1","Child2","Child3")
g_BMF     <- c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3")
u_BMF     <- c("AA-RNA-5","AA-RNA-4","AA-RNA-13","AA-RNA-1","AA-RNA-18",
               "AA-RNA-16","AA-HMH","AA-PJH")

# ---------- 11 manuscript-highlighted DE-lncRNAs ----------
MS_LNCRNA_11 <- c("HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",
                  "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1")

# ---------- load matrix ----------
message("[01] Loading count matrix: ", MATRIX)
raw <- read_tsv(MATRIX, show_col_types = FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
mat <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
rownames(mat) <- raw$EnsemblID
storage.mode(mat) <- "integer"
message(sprintf("  matrix: %d genes x %d samples", nrow(mat), ncol(mat)))

# ---------- 4 filter strategies ----------
# ORIG  : the ACTUAL filter used in the manuscript DESeq2 script
#         (deseq2_analysis*.R: dds <- dds[rowSums(counts(dds)) >= 10, ])
# STRICT: the CPM-based filter as DESCRIBED in the manuscript methods text
#         (mean CPM >= 10 AND CPM >= 1 in all samples)  -- now in R1bis below
# R1/R2 : group-aware alternatives (biologically motivated)
# R3    : no pre-filter (DESeq2 IF only)
apply_filter <- function(counts, group_vec, strategy) {
  cpm_mat <- edgeR::cpm(counts)
  groups  <- unique(group_vec)
  if (strategy == "ORIG") {
    # Manuscript's actual DESeq2 pre-filter
    keep <- rowSums(counts) >= 10
  } else if (strategy == "STRICT") {
    # Manuscript's methods-text description (mean CPM>=10 AND CPM>=1 in all)
    keep <- (rowMeans(cpm_mat) >= 10) & (rowSums(cpm_mat >= 1) == ncol(cpm_mat))
  } else if (strategy == "R1") {
    pass_any <- rep(FALSE, nrow(cpm_mat))
    for (g in groups) {
      idx <- which(group_vec == g)
      pass_any <- pass_any | (rowSums(cpm_mat[, idx, drop=FALSE] >= 1) >= ceiling(length(idx) * 0.5))
    }
    keep <- pass_any
  } else if (strategy == "R2") {
    pass_any <- rep(FALSE, nrow(cpm_mat))
    for (g in groups) {
      idx <- which(group_vec == g)
      mean_ok <- rowMeans(cpm_mat[, idx, drop=FALSE]) >= 1
      half_ok <- rowSums(cpm_mat[, idx, drop=FALSE] >= 0.5) >= ceiling(length(idx) * 0.5)
      pass_any <- pass_any | (mean_ok & half_ok)
    }
    keep <- pass_any
  } else if (strategy == "R3") {
    keep <- rep(TRUE, nrow(cpm_mat))
  } else stop("unknown strategy: ", strategy)
  keep
}

# ---------- DESeq2 helper ----------
run_deseq <- function(counts, coldata, design_formula, contrasts) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  out <- list()
  for (nm in names(contrasts)) {
    cc <- contrasts[[nm]]
    res <- results(dds, contrast = cc, alpha = 0.05)
    out[[nm]] <- as.data.frame(res) %>%
      tibble::rownames_to_column("EnsemblID")
  }
  list(dds = dds, results = out)
}

# ---------- annotate & DE-lncRNA list ----------
de_lncRNA_table <- function(res_df, ann, padj_cut = 0.05, lfc_cut = 1) {
  res_df %>%
    left_join(ann, by = "EnsemblID") %>%
    filter(!is.na(padj), padj < padj_cut, abs(log2FoldChange) > lfc_cut) %>%
    filter(GeneType == "lncRNA") %>%
    arrange(padj)
}

# =====================================================================
# (A) 14-sample manuscript
# =====================================================================
message("\n[02] Building 14-sample manuscript subset")
samples14 <- c(ctrl_int, g_BMF, u_BMF)
mat14 <- mat[, samples14]
meta14 <- data.frame(
  sample = samples14,
  group  = factor(
    ifelse(samples14 %in% ctrl_int, "Control",
           ifelse(samples14 %in% g_BMF, "g_BMF", "u_BMF")),
    levels = c("Control","g_BMF","u_BMF")
  ),
  row.names = samples14
)
print(table(meta14$group))

contrasts_14 <- list(
  gBMF_vs_Ctrl = c("group","g_BMF","Control"),
  uBMF_vs_Ctrl = c("group","u_BMF","Control"),
  gBMF_vs_uBMF = c("group","g_BMF","u_BMF")
)

strategies <- c("ORIG","STRICT","R1","R2","R3")

summary14 <- list()
de_results_14 <- list()

for (s in strategies) {
  message(sprintf("  [14-sample] strategy = %s", s))
  keep <- apply_filter(mat14, as.character(meta14$group), s)
  message(sprintf("    genes kept: %d / %d", sum(keep), length(keep)))
  rr <- run_deseq(mat14[keep, ], meta14, ~ group, contrasts_14)
  de_results_14[[s]] <- rr$results
  for (nm in names(rr$results)) {
    de_lnc <- de_lncRNA_table(rr$results[[nm]], ann)
    n_up <- sum(de_lnc$log2FoldChange >  1)
    n_dn <- sum(de_lnc$log2FoldChange < -1)
    summary14[[length(summary14)+1]] <- tibble(
      cohort = "14-sample", strategy = s, contrast = nm,
      genes_kept = sum(keep),
      n_lncRNA_DE = nrow(de_lnc), n_up = n_up, n_down = n_dn
    )
    write_tsv(de_lnc, file.path(OUTDIR, sprintf("de_lncRNA_14sample_%s_%s.tsv", nm, s)))
  }
}
sum14 <- bind_rows(summary14)
write_tsv(sum14, file.path(OUTDIR, "filter_summary_14sample.tsv"))
print(sum14)

# =====================================================================
# (B) 17-sample hybrid (+3 Child) -- batch-aware
# =====================================================================
message("\n[03] Building 17-sample hybrid subset (+3 Child v3 controls)")
samples17 <- c(ctrl_int, ctrl_pub, g_BMF, u_BMF)
mat17 <- mat[, samples17]
meta17 <- data.frame(
  sample = samples17,
  group  = factor(
    ifelse(samples17 %in% c(ctrl_int, ctrl_pub), "Control",
           ifelse(samples17 %in% g_BMF, "g_BMF", "u_BMF")),
    levels = c("Control","g_BMF","u_BMF")
  ),
  cohort = factor(
    ifelse(samples17 %in% ctrl_pub, "public_MNC", "internal_aspirate"),
    levels = c("internal_aspirate","public_MNC")
  ),
  row.names = samples17
)
print(table(meta17$group, meta17$cohort))

contrasts_17 <- list(
  gBMF_vs_Ctrl = c("group","g_BMF","Control"),
  uBMF_vs_Ctrl = c("group","u_BMF","Control"),
  gBMF_vs_uBMF = c("group","g_BMF","u_BMF")
)

summary17 <- list()
de_results_17 <- list()

for (s in strategies) {
  message(sprintf("  [17-sample] strategy = %s", s))
  keep <- apply_filter(mat17, as.character(meta17$group), s)
  message(sprintf("    genes kept: %d / %d", sum(keep), length(keep)))
  rr <- run_deseq(mat17[keep, ], meta17, ~ cohort + group, contrasts_17)
  de_results_17[[s]] <- rr$results
  for (nm in names(rr$results)) {
    de_lnc <- de_lncRNA_table(rr$results[[nm]], ann)
    n_up <- sum(de_lnc$log2FoldChange >  1)
    n_dn <- sum(de_lnc$log2FoldChange < -1)
    summary17[[length(summary17)+1]] <- tibble(
      cohort = "17-sample", strategy = s, contrast = nm,
      genes_kept = sum(keep),
      n_lncRNA_DE = nrow(de_lnc), n_up = n_up, n_down = n_dn
    )
    write_tsv(de_lnc, file.path(OUTDIR, sprintf("de_lncRNA_17sample_%s_%s.tsv", nm, s)))
  }
}
sum17 <- bind_rows(summary17)
write_tsv(sum17, file.path(OUTDIR, "filter_summary_17sample.tsv"))
print(sum17)

# Combined summary table
sum_all <- bind_rows(sum14, sum17)
write_tsv(sum_all, file.path(OUTDIR, "filter_summary_combined.tsv"))

# =====================================================================
# (C) 11 manuscript lncRNA recovery across all conditions
# =====================================================================
message("\n[04] Manuscript 11-lncRNA recovery audit")
recovery <- list()
gather_recovery <- function(de_results, cohort_lbl) {
  for (s in strategies) {
    for (nm in names(de_results[[s]])) {
      tab <- de_lncRNA_table(de_results[[s]][[nm]], ann)
      for (gn in MS_LNCRNA_11) {
        row <- tab %>% filter(GeneSymbol == gn)
        if (nrow(row) >= 1) {
          recovery[[length(recovery)+1]] <<- tibble(
            cohort = cohort_lbl, strategy = s, contrast = nm,
            gene = gn, recovered = TRUE,
            log2FC = row$log2FoldChange[1], padj = row$padj[1]
          )
        } else {
          recovery[[length(recovery)+1]] <<- tibble(
            cohort = cohort_lbl, strategy = s, contrast = nm,
            gene = gn, recovered = FALSE,
            log2FC = NA_real_, padj = NA_real_
          )
        }
      }
    }
  }
}
gather_recovery(de_results_14, "14-sample")
gather_recovery(de_results_17, "17-sample")
recov_df <- bind_rows(recovery)
write_tsv(recov_df, file.path(OUTDIR, "manuscript_11lncRNA_recovery.tsv"))

# Wide pivot: gene x (cohort.strategy.contrast) -> recovered (TRUE/FALSE)
recov_wide <- recov_df %>%
  mutate(cell = paste(cohort, strategy, contrast, sep = "|")) %>%
  select(gene, cell, recovered) %>%
  pivot_wider(names_from = cell, values_from = recovered)
write_tsv(recov_wide, file.path(OUTDIR, "manuscript_11lncRNA_recovery_wide.tsv"))

# Per-gene recovery rate across the 24 cells (4 strategies x 3 contrasts x 2 cohorts)
recov_rate <- recov_df %>% group_by(gene) %>%
  summarise(n_total = n(), n_recovered = sum(recovered),
            recovery_rate = mean(recovered), .groups = "drop") %>%
  arrange(desc(recovery_rate))
write_tsv(recov_rate, file.path(OUTDIR, "manuscript_11lncRNA_recovery_rate.tsv"))
print(recov_rate)

# =====================================================================
# (D) Robust intersection set (gBMF_vs_Ctrl + uBMF_vs_Ctrl)
# =====================================================================
message("\n[05] Robust intersection across strategies")
get_de_set <- function(de_results, contrast_nm, strategy) {
  de_lncRNA_table(de_results[[strategy]][[contrast_nm]], ann)$GeneSymbol
}

robust_sets <- list()
for (cohort_lbl in c("14-sample","17-sample")) {
  de_obj <- if (cohort_lbl == "14-sample") de_results_14 else de_results_17
  for (nm in c("gBMF_vs_Ctrl","uBMF_vs_Ctrl","gBMF_vs_uBMF")) {
    sets <- lapply(strategies, function(s) get_de_set(de_obj, nm, s))
    names(sets) <- strategies
    all_genes <- unique(unlist(sets))
    member <- sapply(sets, function(x) all_genes %in% x)
    rownames(member) <- all_genes
    n_pass <- rowSums(member)
    robust_sets[[length(robust_sets)+1]] <- tibble(
      cohort = cohort_lbl, contrast = nm, gene = all_genes,
      pass_ORIG   = member[,"ORIG"],   pass_STRICT = member[,"STRICT"],
      pass_R1     = member[,"R1"],     pass_R2     = member[,"R2"],
      pass_R3     = member[,"R3"],
      n_strategies = n_pass
    )
  }
}
rob_df <- bind_rows(robust_sets)
write_tsv(rob_df, file.path(OUTDIR, "robust_lncRNAs_per_contrast.tsv"))

# Genes passing all 5 strategies in any contrast
rob_allK <- rob_df %>% filter(n_strategies == length(strategies)) %>% arrange(cohort, contrast, gene)
write_tsv(rob_allK, file.path(OUTDIR, "robust_lncRNAs_intersection_all_strategies.tsv"))
message(sprintf("  Genes passing all %d strategies (any contrast, any cohort): %d unique",
                length(strategies), length(unique(rob_allK$gene))))

# Ultra-robust: passes all strategies AND found in both cohort runs
ultra <- rob_allK %>% group_by(gene, contrast) %>%
  summarise(n_cohorts = n_distinct(cohort), .groups = "drop") %>%
  filter(n_cohorts == 2)
write_tsv(ultra, file.path(OUTDIR, "ultra_robust_lncRNAs.tsv"))
message(sprintf("  Ultra-robust (all strategies + both cohorts): %d gene-contrast pairs",
                nrow(ultra)))

# =====================================================================
# Done
# =====================================================================
message("\n[DONE] Outputs written to: ", OUTDIR)
message("Key files:")
message("  - filter_summary_combined.tsv       (Reviewer 1 Issue 1 main table)")
message("  - manuscript_11lncRNA_recovery_rate.tsv  (11-biomarker audit)")
message("  - robust_lncRNAs_intersection_all4.tsv")
message("  - ultra_robust_lncRNAs.tsv")
