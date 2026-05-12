#!/usr/bin/env Rscript
# Pipeline reproducibility validation
#   Compares our count matrix (as CPM) against manuscript CPM file
#   Per-sample Pearson + Spearman correlations
#   Independent of DESeq2 version — directly tests featureCounts pipeline
suppressPackageStartupMessages({library(readr); library(dplyr)})

args <- commandArgs(trailingOnly=TRUE)
OUR_COUNTS <- if (length(args)>=1) args[1] else "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_manuscript_v3.txt"
MS_CPM     <- if (length(args)>=2) args[2] else "/Users/jaeeunyoo/Downloads/all_samples_expression_table_CPM_with_entrez.txt"
OUT_LOG    <- if (length(args)>=3) args[3] else "/Volumes/ExtremeSSD/ibmfs/PIPELINE/validation_report.txt"

# Load our counts and compute CPM
our <- read_tsv(OUR_COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(our[, -(1:6)])
rownames(counts) <- our$Geneid
storage.mode(counts) <- "integer"
colnames(counts) <- sub("_sorted$","",sub("^.*/","",sub(".bam$","",colnames(counts))))
our_cpm <- t(t(counts) / colSums(counts)) * 1e6

# Load manuscript CPM
ms <- read_tsv(MS_CPM, show_col_types=FALSE)
ms$EnsemblID <- as.character(ms$EnsemblID)

# Common genes & samples
sample_cols <- intersect(colnames(ms), colnames(our_cpm))
common_ids <- intersect(rownames(our_cpm), ms$EnsemblID)
ms_sub <- ms[match(common_ids, ms$EnsemblID), sample_cols]
our_sub <- our_cpm[common_ids, sample_cols]

# Per-sample correlation
report <- c(
  sprintf("Pipeline Reproducibility Validation Report"),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Our counts: %s", OUR_COUNTS),
  sprintf("Manuscript CPM: %s", MS_CPM),
  sprintf("Common genes: %d, Common samples: %d", length(common_ids), length(sample_cols)),
  "",
  sprintf("%-15s %-10s %-10s %-15s %-12s", "Sample", "Pearson", "Spearman", "Median ratio", "Match"),
  paste(rep("-", 70), collapse="")
)
for (s in sample_cols) {
  v <- our_sub[, s]; m <- as.numeric(ms_sub[[s]])
  p <- cor(v, m, use="complete")
  sp <- cor(v, m, method="spearman", use="complete")
  both <- v > 1 & m > 1 & !is.na(v) & !is.na(m)
  r <- median(v[both]/m[both], na.rm=TRUE)
  match <- if (p > 0.999) "PERFECT" else if (p > 0.98) "Excellent" else if (p > 0.95) "Good" else "Partial"
  report <- c(report, sprintf("%-15s %-10.4f %-10.4f %-15.4f %s", s, p, sp, r, match))
}

# Print + save
cat(paste(report, collapse="\n"))
cat("\n")
writeLines(report, OUT_LOG)
cat(sprintf("\nReport saved: %s\n", OUT_LOG))
