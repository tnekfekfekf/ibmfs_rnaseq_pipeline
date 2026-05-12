#!/usr/bin/env Rscript
# Main DESeq2 contrasts using the batch-aware design (~ cohort + group_combined).
# Three contrasts: g_BMF vs Control, u_BMF vs Control, g_BMF vs u_BMF.
# Outputs full results CSVs (gene_id, gene_name, biotype, log2FC, padj, ...).

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr); library(ggplot2)
  library(EnhancedVolcano)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
OUT  <- file.path(ROOT, "deseq2"); FIG <- file.path(ROOT, "figures")
dir.create(file.path(OUT, "main"), showWarnings = FALSE)

dds  <- readRDS(file.path(OUT, "dds_full.rds"))
meta <- readRDS(file.path(OUT, "samples_meta.rds"))
gm   <- readRDS(file.path(OUT, "gene_meta.rds"))

# Confirm design used
print(design(dds))    # ~ cohort + group_combined
print(table(meta$group_combined, meta$cohort))

dds <- DESeq(dds, parallel = FALSE)
saveRDS(dds, file.path(OUT, "dds_fitted.rds"))

# Helper to dump a contrast with annotation
dump_contrast <- function(dds, contrast, label, ofile) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  cat("== ", label, " ==\n", sep = "")
  print(summary(res))
  d <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gm, by = "gene_id") %>%
    arrange(padj)
  write_tsv(d, ofile)
  d
}

# Three primary contrasts (treatment vs reference)
g_vs_c <- dump_contrast(dds, c("group_combined","g_BMF","Control"),
                        "g_BMF vs Control",
                        file.path(OUT, "main", "DE_gBMF_vs_Control.tsv"))
u_vs_c <- dump_contrast(dds, c("group_combined","u_BMF","Control"),
                        "u_BMF vs Control",
                        file.path(OUT, "main", "DE_uBMF_vs_Control.tsv"))
g_vs_u <- dump_contrast(dds, c("group_combined","g_BMF","u_BMF"),
                        "g_BMF vs u_BMF",
                        file.path(OUT, "main", "DE_gBMF_vs_uBMF.tsv"))

# Volcano plots (manuscript-ready)
volcano_one <- function(d, title, of) {
  pdf(of, width = 7, height = 6)
  p <- EnhancedVolcano(d, lab = d$gene_name, x = "log2FoldChange", y = "padj",
                       pCutoff = 0.05, FCcutoff = 1.0, title = title,
                       subtitle = sprintf("FDR<0.05, |log2FC|>1   (n_DE=%d)",
                                          sum(d$padj < 0.05 & abs(d$log2FoldChange) > 1, na.rm = TRUE)),
                       drawConnectors = TRUE, max.overlaps = 30)
  print(p); dev.off()
}
volcano_one(g_vs_c, "g-BMF vs Control",  file.path(FIG, "volcano_main_gBMF_vs_Control.pdf"))
volcano_one(u_vs_c, "u-BMF vs Control",  file.path(FIG, "volcano_main_uBMF_vs_Control.pdf"))
volcano_one(g_vs_u, "g-BMF vs u-BMF",    file.path(FIG, "volcano_main_gBMF_vs_uBMF.pdf"))

# Brief summary
sigtab <- function(d, fdr = 0.05, lfc = 1.0)
  c(up   = sum(d$padj < fdr & d$log2FoldChange >  lfc, na.rm = TRUE),
    down = sum(d$padj < fdr & d$log2FoldChange < -lfc, na.rm = TRUE),
    total= sum(d$padj < fdr & abs(d$log2FoldChange) > lfc, na.rm = TRUE))
sumtbl <- rbind(`g_BMF vs Control` = sigtab(g_vs_c),
                `u_BMF vs Control` = sigtab(u_vs_c),
                `g_BMF vs u_BMF`   = sigtab(g_vs_u))
print(sumtbl)
write.csv(sumtbl, file.path(OUT, "main", "DE_summary_FDR05_lfc1.csv"))

# Also: lncRNA-only and protein-coding-only summaries
for (bt in c("lncRNA","protein_coding")) {
  sub <- function(d) d[!is.na(d$gene_biotype) & d$gene_biotype == bt, ]
  cat("\n==", bt, "==\n")
  st <- rbind(`g_BMF vs Control` = sigtab(sub(g_vs_c)),
              `u_BMF vs Control` = sigtab(sub(u_vs_c)),
              `g_BMF vs u_BMF`   = sigtab(sub(g_vs_u)))
  print(st)
  write.csv(st, file.path(OUT, "main", paste0("DE_summary_", bt, "_FDR05_lfc1.csv")))
}

message("[03] DONE. Main contrasts written to ", file.path(OUT, "main"))
