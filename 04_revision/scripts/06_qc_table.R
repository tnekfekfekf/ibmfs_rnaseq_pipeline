#!/usr/bin/env Rscript
# QC summary table — adapted for HISAT2 pipeline (parses HISAT2 summary + flagstat + featureCounts).

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })
ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
ALN  <- file.path(ROOT, "aligned"); CNT <- file.path(ROOT, "counts"); QC <- file.path(ROOT, "qc")
dir.create(QC, showWarnings = FALSE, recursive = TRUE)
meta <- read_tsv(file.path(ROOT, "metadata/samples.tsv"), show_col_types = FALSE)

parse_hisat2 <- function(s) {
  f <- file.path(ALN, paste0(s, ".hisat2_summary.txt"))
  if (!file.exists(f)) return(NULL)
  L <- readLines(f)
  num <- function(pat) {
    x <- L[grepl(pat, L)]
    if (length(x) == 0) return(NA_real_)
    as.numeric(sub("\\s+.*","", gsub("^\\s+","", x[1])))
  }
  pct <- function(pat) {
    x <- L[grepl(pat, L)]
    if (length(x) == 0) return(NA_real_)
    as.numeric(sub("%", "", regmatches(x[1], regexpr("[0-9.]+%", x[1]))))
  }
  data.frame(
    sample_id = s,
    input_pairs            = num("reads; of these"),
    conc_unique_pct        = pct("aligned concordantly exactly 1 time"),
    conc_multi_pct         = pct("aligned concordantly >1 times"),
    disc_pct               = pct("aligned discordantly 1 time"),
    overall_align_pct      = pct("overall alignment rate"),
    stringsAsFactors = FALSE
  )
}

parse_flagstat <- function(s) {
  f <- file.path(ALN, paste0(s, ".flagstat.txt"))
  if (!file.exists(f)) return(NULL)
  L <- readLines(f)
  num <- function(pat) as.numeric(sub("\\s+.*","", L[grepl(pat, L)][1]))
  total <- num("in total")
  primary <- num("primary$")     # samtools 1.22 has separate "primary" line
  if (is.na(primary)) primary <- total
  mapped <- num("primary mapped \\(")
  if (is.na(mapped)) mapped <- num("mapped \\(")
  dup <- num("duplicates")
  proper_paired_n <- num("properly paired")
  data.frame(
    sample_id = s,
    total_reads_bam = total,
    primary_mapped = mapped,
    duplicate_reads = dup,
    properly_paired = proper_paired_n,
    dup_pct = round(100 * dup / primary, 3),
    properly_paired_pct = round(100 * proper_paired_n / primary, 2),
    stringsAsFactors = FALSE
  )
}

# featureCounts summary
fc_sum_path <- file.path(CNT, "featureCounts.txt.summary")
fc <- if (file.exists(fc_sum_path)) {
  tab <- read_tsv(fc_sum_path, show_col_types = FALSE)
  colnames(tab)[-1] <- sub("\\.bam$", "", basename(colnames(tab)[-1]))
  tab <- as.data.frame(tab); rownames(tab) <- tab$Status; tab$Status <- NULL
  pick <- function(row) if (row %in% rownames(tab)) as.numeric(tab[row, ]) else rep(NA_real_, ncol(tab))
  data.frame(
    sample_id    = colnames(tab),
    fc_assigned  = pick("Assigned"),
    fc_unassigned_ambiguity = pick("Unassigned_Ambiguity"),
    fc_unassigned_nofeat    = pick("Unassigned_NoFeatures"),
    fc_unassigned_multi     = pick("Unassigned_MultiMapping"),
    fc_unassigned_chimera   = pick("Unassigned_Chimera"),
    stringsAsFactors = FALSE
  )
} else NULL

samples <- meta$sample_id
hs_df  <- bind_rows(lapply(samples, parse_hisat2))
flag_df <- bind_rows(lapply(samples, parse_flagstat))
combo <- meta %>%
  left_join(hs_df,  by = "sample_id") %>%
  left_join(flag_df, by = "sample_id")
if (!is.null(fc)) combo <- combo %>% left_join(fc, by = "sample_id")

# Derived metrics
combo <- combo %>% mutate(
  assigned_per_input_pair = round(fc_assigned / input_pairs, 3),
  assignment_rate_pct     = round(100 * fc_assigned / (fc_assigned + fc_unassigned_ambiguity + fc_unassigned_nofeat + fc_unassigned_multi + fc_unassigned_chimera), 2)
)

write_tsv(combo, file.path(QC, "qc_summary.tsv"))
write.csv(combo, file.path(QC, "qc_summary.csv"), row.names = FALSE)

# Print
cat("\n=== Per-sample QC summary ===\n")
print(combo %>% select(sample_id, group, cohort, input_pairs, overall_align_pct,
                       fc_assigned, assignment_rate_pct, dup_pct, properly_paired_pct),
      n = nrow(combo))

message("\n[06] QC summary written to ", file.path(QC, "qc_summary.tsv"))
