#!/usr/bin/env Rscript
# Reproduce original manuscript DESeq2 analysis
# Based on: fix_gene_types_and_rerun_deseq.R (eAT5.R, Oct 20 2025 final version)
# Control = AA-PRO + AA-KEW
# Design = ~ group
# Separate mRNA / lncRNA analysis (filter by gene_type)
# Target: 2078 g-BMF mRNA DEGs, 1315 u-BMF mRNA DEGs, 4 g-BMF vs u-BMF mRNA

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr)
})

# Output
out_dir <- "/Users/jaeeunyoo/Desktop/star_workdir/deseq2_orig"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load featureCounts output
counts_file <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_orig_macrogen_v44.txt"
cat("Loading counts:", counts_file, "\n")
raw <- read_tsv(counts_file, comment = "#", show_col_types = FALSE)
gene_info <- raw[, 1:6]   # Geneid, Chr, Start, End, Strand, Length
counts <- as.matrix(raw[, -(1:6)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"

# Clean column names (remove path + suffix)
colnames(counts) <- sub("_sorted$", "", sub("^.*/", "", sub("\\.bam$", "", colnames(counts))))

cat(sprintf("Loaded: %d genes x %d samples\n", nrow(counts), ncol(counts)))
print(colnames(counts))

# Setup metadata (manuscript original — Oct 20 2025 final version)
g_aa  <- c("AA-RNA-FA", "AA-RNA-DKC", "AA-RNA-FA2", "AA-RNA-FA3")
ctrl  <- c("AA-PRO", "AA-KEW")
u_aa  <- c("AA-RNA-1", "AA-RNA-4", "AA-RNA-5", "AA-RNA-13", "AA-RNA-16", "AA-RNA-18", "AA-HMH", "AA-PJH")
all_samples <- c(g_aa, ctrl, u_aa)
metadata <- data.frame(
  sample_id = all_samples,
  group = c(rep("G-AA", length(g_aa)),
            rep("Control", length(ctrl)),
            rep("U-AA", length(u_aa))),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$sample_id

# Verify all samples present
available <- intersect(all_samples, colnames(counts))
cat("Available:", paste(available, collapse=", "), "\n")
counts <- counts[, available, drop=FALSE]
metadata <- metadata[available, , drop=FALSE]
metadata$group <- factor(metadata$group, levels = c("Control", "U-AA", "G-AA"))

# Remove all-zero genes (manuscript original — minimal filtering)
counts <- counts[rowSums(counts) > 0, ]
cat(sprintf("After zero-filter: %d genes\n", nrow(counts)))

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ group)
dds <- DESeq(dds, parallel = FALSE)

# Parse gene_type from GENCODE v44 GTF (for mRNA/lncRNA separation)
cat("Parsing gene_type from GTF...\n")
gtf_lines <- readLines("/Volumes/ExtremeSSD/download_ssd/gencode.v44.annotation.gtf", n = -1L)
gene_lines <- gtf_lines[grepl("\\tgene\\t", gtf_lines)]
gene_types <- data.frame(
  gene_id = character(0), gene_type = character(0), gene_name = character(0),
  stringsAsFactors = FALSE
)
parsed <- regmatches(gene_lines,
  regexec('gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
parsed_df <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4) data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
gene_types <- parsed_df %>% distinct(gene_id, .keep_all = TRUE)
cat(sprintf("Parsed %d genes from GTF\n", nrow(gene_types)))
print(table(gene_types$gene_type))

# Helper: get DEG table for a contrast + biotype filter
get_de <- function(grp1, grp2, biotype = "all") {
  res <- results(dds, contrast = c("group", grp1, grp2))
  df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  # Strip version for matching
  df$gene_id_clean <- sub("\\.[0-9]*$", "", df$gene_id)
  gt_clean <- gene_types %>% mutate(gene_id_clean = sub("\\.[0-9]*$", "", gene_id))
  df <- df %>% left_join(gt_clean %>% select(gene_id_clean, gene_type, gene_name), by = "gene_id_clean")
  df <- df %>% filter(!is.na(padj))
  if (biotype == "mrna") df <- df %>% filter(gene_type == "protein_coding")
  if (biotype == "lncrna") df <- df %>% filter(gene_type == "lncRNA")
  df
}

# Three contrasts (manuscript)
contrasts <- list(
  list("G-AA","Control"),
  list("U-AA","Control"),
  list("G-AA","U-AA")
)

results_summary <- data.frame()

for (c in contrasts) {
  g1 <- c[[1]]; g2 <- c[[2]]
  label <- paste0(g1, "_vs_", g2)
  cat("\n=== ", label, " ===\n")
  for (bt in c("all","mrna","lncrna")) {
    d <- get_de(g1, g2, bt)
    sig <- d %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
    up <- sum(sig$log2FoldChange > 0); down <- sum(sig$log2FoldChange < 0)
    cat(sprintf("  %s (%s): %d genes tested, %d sig (up: %d, down: %d)\n", label, bt, nrow(d), nrow(sig), up, down))
    results_summary <- rbind(results_summary, data.frame(
      contrast=label, biotype=bt, total=nrow(d), sig=nrow(sig), up=up, down=down
    ))
    # Save DEG table
    out <- file.path(out_dir, sprintf("DE_%s_%s.tsv", label, bt))
    write_tsv(d, out)
  }
}

cat("\n=== RESULTS SUMMARY ===\n")
print(results_summary)

cat("\n=== Manuscript target comparison ===\n")
cat("g-BMF vs Control mRNA:    target 2078 (sTable 2)\n")
cat("u-BMF vs Control mRNA:    target 1315 (sTable 3)\n")
cat("g-BMF vs u-BMF mRNA:      target    4 (sTable 4)\n")
cat("g-BMF vs Control lncRNA:  target 1167 (sTable 6)\n")
cat("u-BMF vs Control lncRNA:  target  992 (sTable 7)\n")

write_tsv(results_summary, file.path(out_dir, "results_summary.tsv"))
cat("\nDone.\n")
