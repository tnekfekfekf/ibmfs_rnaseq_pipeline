#!/usr/bin/env Rscript
# Merge 14 per-sample featureCounts files into combined matrix
suppressPackageStartupMessages({library(dplyr); library(readr)})

DIR <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample"
OUT <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_manuscript_v44_norRNA.txt"

samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
           "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH")

merged <- NULL
for (s in samps) {
  f <- file.path(DIR, paste0(s, ".counts.txt"))
  d <- read_tsv(f, comment="#", show_col_types=FALSE)
  if (is.null(merged)) {
    merged <- d[, 1:6]  # Geneid, Chr, Start, End, Strand, Length
    merged[[s]] <- d[[7]]
  } else {
    stopifnot(all(merged$Geneid == d$Geneid))
    merged[[s]] <- d[[7]]
  }
}
cat(sprintf("Merged: %d genes × %d samples\n", nrow(merged), length(samps)))
cat(sprintf("Lib sizes (M): %s\n", paste(round(colSums(merged[, samps])/1e6, 1), collapse=", ")))

# Write same format as featureCounts output
write_tsv(merged, OUT)
cat(sprintf("Saved: %s\n", OUT))

# Also write summary aggregate
sum_files <- list.files(DIR, pattern=".summary$", full.names=TRUE)
summary_all <- NULL
for (sf in sum_files) {
  s <- sub("\\.counts\\.txt\\.summary","", basename(sf))
  d <- read_tsv(sf, show_col_types=FALSE)
  colnames(d)[2] <- s
  if (is.null(summary_all)) summary_all <- d
  else summary_all <- merge(summary_all, d, by="Status")
}
write_tsv(summary_all, paste0(OUT, ".summary"))
cat("Summary saved\n")

# Print assignment rates
cat("\n=== Assignment rates (%) ===\n")
totals <- colSums(summary_all[, -1])
assigned <- as.numeric(summary_all[summary_all$Status=="Assigned", -1])
rates <- round(assigned / totals * 100, 1)
names(rates) <- colnames(summary_all)[-1]
print(rates[samps])
cat(sprintf("Mean: %.1f%%\n", mean(rates)))
