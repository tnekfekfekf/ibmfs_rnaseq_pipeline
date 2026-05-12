#!/usr/bin/env Rscript
# Compare v4 (-M --primary added) vs v3 (no -M) CPM correlation with manuscript
suppressPackageStartupMessages({library(readr); library(dplyr)})

V3_DIR <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3"
V4_DIR <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v4_Mprimary"
MS_CPM <- "/Users/jaeeunyoo/Downloads/all_samples_expression_table_CPM_with_entrez.txt"

# Load manuscript CPM
ms <- read_tsv(MS_CPM, show_col_types=FALSE)
ms$EnsemblID <- as.character(ms$EnsemblID)

# Load both v3 and v4 for the 10 patient samples
patient_samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3",
                   "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18")

load_counts <- function(dir, samps) {
  out <- NULL
  for (s in samps) {
    d <- read_tsv(file.path(dir, paste0(s, ".counts.txt")), comment="#", show_col_types=FALSE)
    if (is.null(out)) { out <- d[, "Geneid"]; }
    out[[s]] <- d[[7]]
  }
  m <- as.matrix(out[, -1])
  rownames(m) <- out$Geneid
  m
}

v3 <- load_counts(V3_DIR, patient_samps)
v4 <- load_counts(V4_DIR, patient_samps)

# Convert to CPM
v3_cpm <- t(t(v3) / colSums(v3)) * 1e6
v4_cpm <- t(t(v4) / colSums(v4)) * 1e6

# Match Ensembl IDs
common <- intersect(rownames(v3_cpm), ms$EnsemblID)
ms_sub <- ms[match(common, ms$EnsemblID), patient_samps]
v3_sub <- v3_cpm[common, ]
v4_sub <- v4_cpm[common, ]

cat("=== CPM correlation vs MANUSCRIPT — comparison ===\n")
cat(sprintf("%-15s %-15s %-15s %-15s %-15s %-15s\n",
            "Sample", "v3 Pearson", "v4 Pearson", "Δ Pearson",
            "v3 median ratio", "v4 median ratio"))
cat(paste(rep("-", 100), collapse=""), "\n")
results <- list()
for (s in patient_samps) {
  v3v <- v3_sub[, s]
  v4v <- v4_sub[, s]
  mv <- as.numeric(ms_sub[[s]])
  p3 <- cor(v3v, mv, use="complete")
  p4 <- cor(v4v, mv, use="complete")
  # Median ratios on non-zero pairs
  bp3 <- v3v > 1 & mv > 1; r3 <- median(v3v[bp3]/mv[bp3], na.rm=TRUE)
  bp4 <- v4v > 1 & mv > 1; r4 <- median(v4v[bp4]/mv[bp4], na.rm=TRUE)
  delta <- p4 - p3
  marker <- if (delta > 0.001) " ✓" else if (delta < -0.001) " ✗" else " ="
  cat(sprintf("%-15s %-15.4f %-15.4f %+.4f%s    %-15.4f %.4f\n",
              s, p3, p4, delta, marker, r3, r4))
  results[[s]] <- c(v3_pearson=p3, v4_pearson=p4, delta=delta, v3_ratio=r3, v4_ratio=r4)
}

cat("\n=== Summary ===\n")
all_p3 <- sapply(results, function(x) x["v3_pearson"])
all_p4 <- sapply(results, function(x) x["v4_pearson"])
cat(sprintf("Mean Pearson — v3: %.4f, v4: %.4f, Δ: %+.4f\n",
            mean(all_p3), mean(all_p4), mean(all_p4-all_p3)))
cat(sprintf("Median ratio — v3: %.4f, v4: %.4f\n",
            median(sapply(results, function(x) x["v3_ratio"])),
            median(sapply(results, function(x) x["v4_ratio"]))))

cat("\n=== Interpretation ===\n")
better <- sum(all_p4 > all_p3 + 0.001)
worse  <- sum(all_p4 < all_p3 - 0.001)
same   <- length(all_p3) - better - worse
cat(sprintf("v4 better (>+0.001 r): %d/10\n", better))
cat(sprintf("v4 worse  (<-0.001 r): %d/10\n", worse))
cat(sprintf("v4 same   (±0.001 r): %d/10\n", same))
