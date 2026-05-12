#!/usr/bin/env Rscript
# Sensitivity analyses to address Reviewer 2 concerns:
#   (a) internal-only controls (n=2)         -> reproduces original analysis
#   (b) public-only controls (n=3)           -> independent replication
#   (c) combined (n=5) without Child1        -> drop age outlier (217mo)
#   (d) FA-only g-BMF vs Control             -> 3 FA, drop DKC
#   (e) DKC alone vs Control                 -> 1 DKC reported separately
# Each runs DESeq2 with its own appropriate design and dumps a summary.

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
OUT  <- file.path(ROOT, "deseq2", "sensitivity"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dds_full <- readRDS(file.path(ROOT, "deseq2/dds_full.rds"))
gm       <- readRDS(file.path(ROOT, "deseq2/gene_meta.rds"))
meta     <- readRDS(file.path(ROOT, "deseq2/samples_meta.rds"))

run_sub <- function(keep_ids, design_formula, contrasts, label) {
  message("\n[04] Subset: ", label, " (n=", length(keep_ids), ")")
  d <- dds_full[, keep_ids]
  colData(d) <- droplevels(colData(d))
  design(d) <- design_formula
  d <- DESeq(d, parallel = FALSE)
  for (cc in contrasts) {
    res <- results(d, contrast = cc[[1]], alpha = 0.05)
    fname <- sprintf("DE_%s_%s.tsv", label, cc[[2]])
    out <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id") %>%
           left_join(gm, by = "gene_id") %>% arrange(padj)
    write_tsv(out, file.path(OUT, fname))
    sig <- sum(out$padj < 0.05 & abs(out$log2FoldChange) > 1, na.rm = TRUE)
    message("   ", cc[[2]], ": ", sig, " DEGs (FDR<0.05, |LFC|>1)")
  }
}

# (a) Internal-only controls (n=2 IDA/HUS) --------------------
keep_a <- meta$sample_id[meta$cohort == "internal_aspirate"]
run_sub(keep_a,
        ~ group_combined,         # only one cohort, drop cohort term
        list(list(c("group_combined","g_BMF","Control"), "gBMF_vs_Control"),
             list(c("group_combined","u_BMF","Control"), "uBMF_vs_Control"),
             list(c("group_combined","g_BMF","u_BMF"),   "gBMF_vs_uBMF")),
        "INTERNAL_ONLY")

# (b) Public-only controls (n=3 PNBM) -----------------------
# In this subset: Controls are all public, BMF are all internal. cohort and
# group are perfectly confounded for the Control class, so cohort term must be
# dropped. The DE estimate then conflates cohort with group; this is reported
# as a disclosure for completeness of the sensitivity panel.
keep_b <- meta$sample_id[meta$group == "control_public" |
                         meta$group %in% c("u_BMF","g_BMF")]
run_sub(keep_b,
        ~ group_combined,        # cohort dropped (collinear with Control)
        list(list(c("group_combined","g_BMF","Control"), "gBMF_vs_Control"),
             list(c("group_combined","u_BMF","Control"), "uBMF_vs_Control"),
             list(c("group_combined","g_BMF","u_BMF"),   "gBMF_vs_uBMF")),
        "PUBLIC_ONLY_confounded")

# (c) Combined n=5, drop Child1 (age 217mo outlier) ---------
keep_c <- meta$sample_id[meta$sample_id != "Child1"]
run_sub(keep_c,
        ~ cohort + group_combined,
        list(list(c("group_combined","g_BMF","Control"), "gBMF_vs_Control"),
             list(c("group_combined","u_BMF","Control"), "uBMF_vs_Control"),
             list(c("group_combined","g_BMF","u_BMF"),   "gBMF_vs_uBMF")),
        "DROP_CHILD1")

# (d) FA-only g-BMF vs Control (drop DKC) -------------------
keep_d <- meta$sample_id[meta$sample_id != "AA-RNA-DKC"]
run_sub(keep_d,
        ~ cohort + subgroup,
        list(list(c("subgroup","FA",   "Control"), "FA_vs_Control"),
             list(c("subgroup","u_BMF","Control"), "uBMF_vs_Control"),
             list(c("subgroup","FA",   "u_BMF"),   "FA_vs_uBMF")),
        "FA_ONLY_gBMF")

# (e) DKC alone vs Control (1 DKC, all 5 controls) ---------
keep_e <- meta$sample_id[meta$subgroup %in% c("Control","DKC")]
run_sub(keep_e,
        ~ cohort + subgroup,
        list(list(c("subgroup","DKC","Control"), "DKC_vs_Control")),
        "DKC_ONLY")

message("\n[04] All sensitivity analyses written to ", OUT)
