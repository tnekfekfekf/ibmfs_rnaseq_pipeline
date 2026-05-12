#!/usr/bin/env Rscript
# GSEA at multiple FDR thresholds — Reviewer 1 asks how robust enrichments are
# under FDR<0.25 (exploratory) vs FDR<0.05 (stringent).
#
# Uses clusterProfiler/gseGO + gseKEGG with pre-ranked statistics from each contrast.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
})
need <- c("clusterProfiler","org.Hs.eg.db","msigdbr","enrichplot","DOSE")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  message("Installing: ", paste(miss, collapse=", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr); library(enrichplot)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
OUT  <- file.path(ROOT, "deseq2", "gsea"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Use a robust ranking metric: Wald stat (signed); fall back to log2FC * -log10(p)
rank_stat <- function(d) {
  s <- d$stat
  s[is.na(s)] <- (d$log2FoldChange[is.na(s)] * -log10(d$pvalue[is.na(s)] + 1e-300))
  s
}

run_gsea_one <- function(de_tsv, label) {
  cat("== GSEA: ", label, " ==\n", sep = "")
  d <- read_tsv(de_tsv, show_col_types = FALSE) %>% filter(!is.na(stat))
  d$ENTREZ <- mapIds(org.Hs.eg.db, keys = d$gene_name, column = "ENTREZID",
                     keytype = "SYMBOL", multiVals = "first")
  d <- d %>% filter(!is.na(ENTREZ)) %>% distinct(ENTREZ, .keep_all = TRUE) %>%
       mutate(stat = rank_stat(.))
  ranks <- setNames(d$stat, d$ENTREZ)
  ranks <- sort(ranks, decreasing = TRUE)

  results <- list()
  for (db in c("GO:BP","KEGG","H_HALLMARK")) {
    if (db == "GO:BP") {
      r <- gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, ont = "BP",
                 minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 0)
    } else if (db == "KEGG") {
      r <- tryCatch(gseKEGG(geneList = ranks, organism = "hsa",
                            minGSSize = 15, pvalueCutoff = 1, eps = 0),
                    error = function(e) { message("KEGG failed: ", e$message); NULL })
    } else {
      m <- msigdbr(species = "Homo sapiens", category = "H") %>%
           dplyr::select(gs_name, entrez_gene)
      r <- GSEA(geneList = ranks, TERM2GENE = m,
                minGSSize = 15, pvalueCutoff = 1, eps = 0)
    }
    if (is.null(r)) next
    df <- as.data.frame(r) %>% arrange(p.adjust)
    write_tsv(df, file.path(OUT, sprintf("GSEA_%s_%s.tsv", label, gsub("[: ]","_", db))))
    # Threshold-sensitivity counts
    for (thr in c(0.05, 0.10, 0.25)) {
      n <- sum(df$p.adjust < thr, na.rm = TRUE)
      cat(sprintf("   %s  FDR<%.2f: %d sig\n", db, thr, n))
    }
    results[[db]] <- df
  }
  results
}

main_dir <- file.path(ROOT, "deseq2", "main")
files <- list.files(main_dir, pattern = "^DE_.*\\.tsv$", full.names = TRUE)
for (f in files) {
  label <- sub("^DE_", "", sub("\\.tsv$", "", basename(f)))
  run_gsea_one(f, label)
}
message("[08] DONE. GSEA tables in ", OUT)
