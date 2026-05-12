#!/usr/bin/env Rscript
# Merge per-sample featureCounts outputs into a single count matrix
# Usage: Rscript 02_merge_counts.R [PER_SAMPLE_DIR] [OUT_FILE] [SAMPLE_NAMES...]
suppressPackageStartupMessages({library(dplyr); library(readr)})

args <- commandArgs(trailingOnly=TRUE)
DIR <- if (length(args) >= 1) args[1] else "/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3"
OUT <- if (length(args) >= 2) args[2] else "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_manuscript_v3.txt"

# Default: 14 manuscript samples in canonical order
default_samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
                    "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH")
samps <- if (length(args) >= 3) args[3:length(args)] else default_samps

cat(sprintf("Merging %d samples from %s\n", length(samps), DIR))
merged <- NULL
for (s in samps) {
  f <- file.path(DIR, paste0(s, ".counts.txt"))
  if (!file.exists(f)) { cat(sprintf("  MISSING: %s\n", f)); next }
  d <- read_tsv(f, comment="#", show_col_types=FALSE)
  if (is.null(merged)) { merged <- d[, 1:6]; merged[[s]] <- d[[7]] }
  else {
    stopifnot(all(merged$Geneid == d$Geneid))
    merged[[s]] <- d[[7]]
  }
}

write_tsv(merged, OUT)
cat(sprintf("\nMerged %d genes × %d samples\n", nrow(merged), length(samps)))
cat(sprintf("Library sizes (M): %s\n",
            paste(round(colSums(merged[, samps])/1e6, 1), collapse=", ")))
cat(sprintf("Saved: %s\n", OUT))
