#!/usr/bin/env Rscript
# Bulk RNA-seq cellular deconvolution (Reviewer 2 #4 — composition confounder).
# Uses immunedeconv (xCell + MCP-counter + EPIC). xCell is signature-based,
# MCP-counter is markers-based, EPIC is constrained regression — three views.
#
# Input: VST-stabilized expression OR TPM-like; we use TPM derived from raw counts + gene length.
# Output: per-sample cell-type scores, plus statistical comparison vs group.

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  message("Installing immunedeconv from GitHub (Bioconductor removed it)...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(c("preprocessCore","limSolve","GSVA"), ask = FALSE, update = FALSE)
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("omnideconv/immunedeconv", upgrade = "never")
}
suppressPackageStartupMessages(library(immunedeconv))

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
OUT  <- file.path(ROOT, "deconv"); dir.create(OUT, showWarnings = FALSE)

dds <- readRDS(file.path(ROOT, "deseq2/dds_full.rds"))
gm  <- readRDS(file.path(ROOT, "deseq2/gene_meta.rds"))
meta <- readRDS(file.path(ROOT, "deseq2/samples_meta.rds"))

# Compute TPM from raw counts + gene length (mcols stored when GTF parsed not yet — use featureCounts output length)
fc_raw <- read_tsv(file.path(ROOT, "counts/featureCounts.cleaned.txt"), comment = "#", show_col_types = FALSE)
length_kb <- fc_raw$Length / 1000
names(length_kb) <- fc_raw$Geneid
counts <- counts(dds)
length_kb <- length_kb[rownames(counts)]
rate <- counts / length_kb
tpm <- t( t(rate) / colSums(rate, na.rm = TRUE) ) * 1e6
rownames(tpm) <- rownames(counts)

# immunedeconv expects HGNC symbols, not Ensembl IDs. Aggregate by gene_name.
sym_lookup <- setNames(gm$gene_name, gm$gene_id)
tpm_df <- as.data.frame(tpm) %>% tibble::rownames_to_column("gene_id") %>%
  mutate(gene = sym_lookup[gene_id]) %>%
  filter(!is.na(gene) & gene != "") %>%
  group_by(gene) %>% summarise(across(where(is.numeric), \(x) sum(x, na.rm=TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

# ---- xCell ----
message("[07] Running xCell on ", ncol(tpm_df), " samples")
res_xcell <- deconvolute(tpm_df, "xcell")
write_tsv(res_xcell, file.path(OUT, "deconv_xCell.tsv"))

# ---- MCP-counter ----
message("[07] Running MCP-counter")
res_mcp <- deconvolute(tpm_df, "mcp_counter")
write_tsv(res_mcp, file.path(OUT, "deconv_MCPcounter.tsv"))

# ---- EPIC (with HBM=FALSE for non-blood) — use default tumor signature ----
res_epic <- tryCatch({
  message("[07] Running EPIC")
  immunedeconv::deconvolute(tpm_df, "epic")
}, error = function(e) { message("EPIC failed: ", e$message); NULL })
if (!is.null(res_epic)) write_tsv(res_epic, file.path(OUT, "deconv_EPIC.tsv"))

# Diagnostic: stats per cell type, group_combined, pre/post batch correction
long <- res_xcell %>% pivot_longer(-cell_type, names_to = "sample_id", values_to = "score") %>%
  left_join(meta %>% select(sample_id, group_combined, cohort), by = "sample_id")

ggplot(long, aes(group_combined, score, fill = cohort)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 6) +
  theme_bw(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("xCell scores by group, split by cohort (internal vs public)")
ggsave(file.path(OUT, "deconv_xCell_by_group.pdf"), width = 18, height = 14)

# Test: does cohort drive cell composition more than group?
diff_tab <- long %>% group_by(cell_type) %>%
  summarise(
    p_cohort = tryCatch(kruskal.test(score ~ cohort)$p.value, error = \(e) NA),
    p_group  = tryCatch(kruskal.test(score ~ group_combined)$p.value, error = \(e) NA),
    .groups = "drop"
  ) %>%
  mutate(p_cohort_BH = p.adjust(p_cohort, "BH"),
         p_group_BH  = p.adjust(p_group,  "BH"))
write_tsv(diff_tab, file.path(OUT, "deconv_xCell_pvals.tsv"))
print(arrange(diff_tab, p_cohort_BH))

message("[07] DONE. Decomposition tables in ", OUT)
