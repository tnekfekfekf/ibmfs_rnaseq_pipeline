#!/usr/bin/env Rscript
# Final analysis: site-aware DESeq2 (~ site + group) on raw counts
# 14-sample (manuscript cohort): site = macrogen / jinpyung
# 17-sample (with Child): site = macrogen / jinpyung / public
# Group-aware filter on raw CPM

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma)
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS"
OUTDIR <- file.path(ROOT, "final_site_aware_analysis")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

raw <- read_tsv(file.path(ROOT, "manuscript_count_matrix_19samples.txt"), show_col_types=FALSE)
ann <- raw[, c("EnsemblID","GeneSymbol","GeneName","GeneType")]
cm  <- as.matrix(raw[, !(colnames(raw) %in% c("EnsemblID","GeneSymbol","GeneName","GeneType"))])
storage.mode(cm) <- "integer"; rownames(cm) <- raw$EnsemblID

ms11 <- c("HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",
          "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1")

site_of <- function(s) ifelse(s %in% c("AA-PRO","AA-KEW","AA-HMH","AA-PJH"), "jinpyung",
                       ifelse(s %in% c("Child1","Child2","Child3"), "public", "macrogen"))
group_of <- function(s) ifelse(s %in% c("Child1","Child2","Child3","AA-PRO","AA-KEW"), "Control",
                        ifelse(s %in% c("AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3"), "g_BMF", "u_BMF"))

# Group-aware filter (raw CPM)
group_aware_pass <- function(cpm_mat, groups) {
  pass <- rep(FALSE, nrow(cpm_mat))
  for (g in unique(groups)) {
    idx <- which(groups == g)
    pass <- pass | (rowSums(cpm_mat[, idx, drop=FALSE] >= 1) >= ceiling(length(idx)*0.5))
  }
  pass
}

run_one <- function(samples, label, design_formula) {
  cat("\n========================================================\n")
  cat("Running:", label, "\n")
  cat("Samples:", length(samples), " | Formula:", deparse(design_formula), "\n")

  mat <- cm[, samples]
  meta <- data.frame(
    sample=samples,
    site=factor(site_of(samples), levels=c("macrogen","jinpyung","public")),
    group=factor(group_of(samples), levels=c("Control","g_BMF","u_BMF")),
    row.names=samples
  )
  cat("Design table:\n"); print(table(meta$site, meta$group))

  keep <- rowSums(mat) >= 10
  mat_f <- mat[keep, ]
  cat(sprintf("Pre-filter rowSums>=10: %d genes\n", sum(keep)))

  dds <- DESeqDataSetFromMatrix(mat_f, meta, design=design_formula)
  dds <- DESeq(dds, quiet=TRUE)

  cpm_raw <- cpm(mat_f)
  filter_pass_ids <- rownames(cpm_raw)[group_aware_pass(cpm_raw, as.character(meta$group))]
  cat(sprintf("Group-aware filter (CPM>=1 in >=50%% of any group): %d genes pass\n",
              length(filter_pass_ids)))

  dump <- function(contrast_name, contrast_vec) {
    r <- results(dds, contrast=contrast_vec, alpha=0.05)
    d <- as.data.frame(r) %>% tibble::rownames_to_column("EnsemblID") %>%
      left_join(ann, by="EnsemblID")
    # Apply group-aware filter
    d <- d %>% filter(EnsemblID %in% filter_pass_ids) %>% arrange(padj)
    write_tsv(d, file.path(OUTDIR, sprintf("DE_%s_%s.tsv", label, contrast_name)))
    n_sig <- sum(d$padj<0.05 & abs(d$log2FoldChange)>1, na.rm=TRUE)
    n_up  <- sum(d$padj<0.05 & d$log2FoldChange>1, na.rm=TRUE)
    n_dn  <- sum(d$padj<0.05 & d$log2FoldChange< -1, na.rm=TRUE)
    n_lnc <- sum(d$padj<0.05 & abs(d$log2FoldChange)>1 & d$GeneType=="lncRNA", na.rm=TRUE)
    n_pc  <- sum(d$padj<0.05 & abs(d$log2FoldChange)>1 & d$GeneType=="protein_coding", na.rm=TRUE)
    cat(sprintf("  %s : n_DE=%d (up=%d, down=%d) | lncRNA=%d, protein_coding=%d\n",
                contrast_name, n_sig, n_up, n_dn, n_lnc, n_pc))
    d
  }

  g_vs_c <- dump("gBMF_vs_Ctrl", c("group","g_BMF","Control"))
  u_vs_c <- dump("uBMF_vs_Ctrl", c("group","u_BMF","Control"))
  g_vs_u <- dump("gBMF_vs_uBMF", c("group","g_BMF","u_BMF"))

  # 11 lncRNA audit
  cat("\n  --- 11 manuscript lncRNAs ---\n")
  audit_rows <- list()
  for (g in ms11) {
    r1 <- g_vs_c %>% filter(GeneSymbol==g) %>% slice(1)
    r2 <- u_vs_c %>% filter(GeneSymbol==g) %>% slice(1)
    if (nrow(r1)==0) next
    sig1 <- !is.na(r1$padj) && r1$padj<0.05 && abs(r1$log2FoldChange)>1
    sig2 <- !is.na(r2$padj) && r2$padj<0.05 && abs(r2$log2FoldChange)>1
    cat(sprintf("    %-12s | g-BMF LFC=%6.2f padj=%.2e %s | u-BMF LFC=%6.2f padj=%.2e %s\n",
                g, r1$log2FoldChange, r1$padj, ifelse(sig1,"SIG","ns"),
                r2$log2FoldChange, r2$padj, ifelse(sig2,"SIG","ns")))
    audit_rows[[g]] <- data.frame(
      gene=g, cohort=label,
      gBMF_LFC=r1$log2FoldChange, gBMF_padj=r1$padj, gBMF_sig=sig1,
      uBMF_LFC=r2$log2FoldChange, uBMF_padj=r2$padj, uBMF_sig=sig2
    )
  }
  audit_df <- bind_rows(audit_rows)
  write_tsv(audit_df, file.path(OUTDIR, sprintf("ms11_audit_%s.tsv", label)))

  list(dds=dds, g_vs_c=g_vs_c, u_vs_c=u_vs_c, g_vs_u=g_vs_u, audit=audit_df, meta=meta)
}

# -- Scenario A: 14-sample manuscript cohort, ~ site + group (2-level site) --
s14 <- c("AA-PRO","AA-KEW",
         "AA-RNA-DKC","AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3",
         "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-18","AA-RNA-16",
         "AA-HMH","AA-PJH")
res14 <- run_one(s14, "14sample_site_group", ~ site + group)

# -- Scenario B: 17-sample (+Child), ~ site + group (3-level site) --
s17 <- c("Child1","Child2","Child3", s14)
res17 <- run_one(s17, "17sample_site_group", ~ site + group)

# -- Scenario C: 14-sample, ~ group only (no batch covariate) for reference --
res14_nb <- run_one(s14, "14sample_group_only", ~ group)

# -- Combined summary table --
all_audit <- bind_rows(res14$audit, res17$audit, res14_nb$audit)
write_tsv(all_audit, file.path(OUTDIR, "ALL_ms11_audit_combined.tsv"))

# Pretty summary table
cat("\n\n=========================================================\n")
cat("SUMMARY: 11 manuscript lncRNA recovery across 3 scenarios\n")
cat("=========================================================\n")
summary_tbl <- all_audit %>%
  select(gene, cohort, gBMF_sig, uBMF_sig) %>%
  pivot_wider(names_from=cohort, values_from=c(gBMF_sig, uBMF_sig))
print(summary_tbl, n=Inf, width=Inf)
write_tsv(summary_tbl, file.path(OUTDIR, "SUMMARY_ms11_recovery_wide.tsv"))

# Counts per scenario
cnt <- all_audit %>% group_by(cohort) %>%
  summarise(g_BMF_sig=sum(gBMF_sig), u_BMF_sig=sum(uBMF_sig), n=n(), .groups="drop")
cat("\n=== Summary counts ===\n")
print(cnt)
write_tsv(cnt, file.path(OUTDIR, "SUMMARY_ms11_counts.tsv"))

# -- DE counts table (overall) --
de_counts <- list()
for (lbl in c("14sample_site_group","17sample_site_group","14sample_group_only")) {
  for (ct in c("gBMF_vs_Ctrl","uBMF_vs_Ctrl","gBMF_vs_uBMF")) {
    f <- file.path(OUTDIR, sprintf("DE_%s_%s.tsv", lbl, ct))
    if (!file.exists(f)) next
    d <- read_tsv(f, show_col_types=FALSE)
    sig <- d %>% filter(padj<0.05 & abs(log2FoldChange)>1)
    de_counts[[length(de_counts)+1]] <- data.frame(
      cohort=lbl, contrast=ct,
      n_DE=nrow(sig),
      n_lncRNA=sum(sig$GeneType=="lncRNA", na.rm=TRUE),
      n_pc=sum(sig$GeneType=="protein_coding", na.rm=TRUE)
    )
  }
}
de_counts_df <- bind_rows(de_counts)
write_tsv(de_counts_df, file.path(OUTDIR, "SUMMARY_DE_counts.tsv"))
cat("\n=== Overall DE counts (padj<0.05, |LFC|>1) ===\n")
print(de_counts_df)

message("\n[DONE] All outputs in ", OUTDIR)
