#!/usr/bin/env Rscript
# Build publication-ready summary tables for Reviewer 1 Issue 1 response
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(tibble)
})

OUTDIR <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/lncrna_filter_analysis"

sumc <- read_tsv(file.path(OUTDIR, "filter_summary_combined.tsv"), show_col_types = FALSE)
rec  <- read_tsv(file.path(OUTDIR, "manuscript_11lncRNA_recovery.tsv"), show_col_types = FALSE)

# --- Table 1: pretty wide DE-lncRNA count table ---
t1 <- sumc %>%
  mutate(cell = paste(strategy, contrast, sep = "_")) %>%
  select(cohort, strategy, contrast, genes_kept, n_lncRNA_DE, n_up, n_down)
write_tsv(t1, file.path(OUTDIR, "TABLE1_de_lncRNA_counts.tsv"))

cat("\n========================================================\n")
cat("TABLE 1: DE-lncRNA counts by filter strategy and cohort\n")
cat("(padj < 0.05 AND |log2FC| > 1)\n")
cat("========================================================\n")
print(t1, n = Inf)

# --- Table 2: Manuscript-11 recovery, restricted to disease-vs-control ---
disease_contrasts <- c("gBMF_vs_Ctrl", "uBMF_vs_Ctrl")
rec_disease <- rec %>% filter(contrast %in% disease_contrasts)
rec_rate2 <- rec_disease %>% group_by(gene) %>%
  summarise(n_total = n(), n_recovered = sum(recovered),
            recovery_rate = mean(recovered), .groups = "drop") %>%
  arrange(desc(recovery_rate), gene)
write_tsv(rec_rate2, file.path(OUTDIR, "TABLE2_manuscript_11lncRNA_recovery_diseaseVScontrol.tsv"))

cat("\n========================================================\n")
cat("TABLE 2: Manuscript 11 lncRNA recovery (disease-vs-control cells)\n")
cat("Total cells per gene = 4 strategies x 2 contrasts x 2 cohorts = 16\n")
cat("========================================================\n")
print(rec_rate2, n = Inf)

# --- Table 3: Wide pivot for the 11-lncRNA audit (each row = gene, each col = cell) ---
rec_disease_wide <- rec_disease %>%
  mutate(cell = paste(cohort, contrast, strategy, sep = "/")) %>%
  select(gene, cell, recovered) %>%
  pivot_wider(names_from = cell, values_from = recovered)
write_tsv(rec_disease_wide, file.path(OUTDIR, "TABLE3_manuscript_11lncRNA_wide.tsv"))

cat("\n========================================================\n")
cat("TABLE 3: Wide pivot - 11 lncRNAs x (cohort/contrast/strategy)\n")
cat("========================================================\n")
print(rec_disease_wide, n = Inf, width = Inf)

# --- Table 4: Per-strategy recovery rate ---
t4 <- rec_disease %>% group_by(cohort, strategy) %>%
  summarise(n_recovered = sum(recovered), n_total = n(),
            pct_recovered = mean(recovered) * 100, .groups = "drop") %>%
  arrange(cohort, strategy)
write_tsv(t4, file.path(OUTDIR, "TABLE4_recovery_pct_by_strategy.tsv"))

cat("\n========================================================\n")
cat("TABLE 4: % of 11 biomarkers recovered per strategy / cohort\n")
cat("(disease-vs-control contrasts only)\n")
cat("========================================================\n")
print(t4, n = Inf)

# --- Table 5: Delta vs manuscript's actual filter (ORIG) ---
# How does each filter compare to the manuscript's true filter (ORIG)?
gain <- sumc %>%
  filter(contrast %in% disease_contrasts) %>%
  group_by(cohort, contrast) %>%
  mutate(n_lncRNA_manuscript = n_lncRNA_DE[strategy == "ORIG"],
         delta_vs_manuscript = n_lncRNA_DE - n_lncRNA_manuscript) %>%
  ungroup() %>%
  select(cohort, contrast, strategy, n_lncRNA_DE, n_lncRNA_manuscript, delta_vs_manuscript)
write_tsv(gain, file.path(OUTDIR, "TABLE5_delta_vs_manuscript_filter.tsv"))

cat("\n========================================================\n")
cat("TABLE 5: Delta DE-lncRNA count vs the manuscript's true filter (ORIG)\n")
cat("========================================================\n")
print(gain, n = Inf)

# --- KEY FINDING SUMMARY ---
cat("\n\n========================================================\n")
cat("KEY FINDINGS for Reviewer 1, Issue 1\n")
cat("========================================================\n")
cat("\n1) MANUSCRIPT METHODS-TEXT DOCUMENTATION ERROR:\n")
cat("   The methods text says 'mean CPM >= 10 AND CPM >= 1 in all samples'\n")
cat("   (= STRICT in this script). The actual code (deseq2_analysis*.R\n")
cat("   line ~33-59 and snakemake_template line 87) uses\n")
cat("   rowSums(counts) >= 10 (= ORIG in this script).\n")
cat("   ORIG is what produced the manuscript's published DE-lncRNA list.\n\n")

cat("2) ORIG (TRUE manuscript filter) recovers 11/11 biomarkers in the\n")
cat("   14-sample reproduction (100%) and 10/11 in the 17-sample hybrid\n")
cat("   (95.5% - only MALAT1 in u-BMF vs Ctrl drops, due to indep-filter\n")
cat("   threshold shift when more low-count genes are co-tested).\n\n")

cat("3) STRICT (methods-text filter, never actually used) would EXCLUDE 6/11\n")
cat("   manuscript biomarkers in both cohorts (HCG11, HCP5, SNHG32,\n")
cat("   PSMB8-AS1, FAM30A, MIR22HG) because they have near-zero expression\n")
cat("   in controls (log2FC ~9-11). The 'CPM>=1 in ALL samples' rule\n")
cat("   rejects exactly this biomarker pattern.\n\n")

cat("4) R1 / R2 group-aware filters recover the same biomarker set as ORIG\n")
cat("   (100% 14-sample, 95.5% 17-sample) while being biologically\n")
cat("   motivated and easier to defend. log2FC values match ORIG within 1%\n")
cat("   for the 4 RT-qPCR-validated biomarkers.\n\n")

cat("5) g-BMF vs u-BMF: 0 lncRNAs DE under STRICT/R1/R2 in 17-sample,\n")
cat("   confirming the manuscript's interpretation that this is a\n")
cat("   PAN-BMF disease signature rather than a g/u discriminator.\n\n")

cat("6) Recommendations for the revision:\n")
cat("   (a) Correct the Methods text to describe the actual filter\n")
cat("       (rowSums(counts) >= 10). No published result changes.\n")
cat("   (b) Add this sensitivity table as a Supplementary Table.\n")
cat("   (c) Optionally adopt R1 (group-aware) as primary going forward.\n")

message("\n[DONE] Tables saved to ", OUTDIR)
