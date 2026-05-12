#!/usr/bin/env Rscript
# Load featureCounts matrix + samples.tsv, attach gene biotype, build DESeq2 dataset.
# Output: revision_analysis/deseq2/dds_full.rds, gene_meta.rds

suppressPackageStartupMessages({
  library(DESeq2)
  library(rtracklayer)
  library(tibble)
  library(dplyr)
  library(readr)
})

ROOT <- "/Volumes/ExtremeSSD/ibmfs/revision_analysis"
COUNTS <- file.path(ROOT, "counts/featureCounts.cleaned.txt")
META   <- file.path(ROOT, "metadata/samples.tsv")
GTF    <- "/Volumes/ExtremeSSD/ibmfs/ibmfs_fastq_raw_data/reference/gencode/gencode.v45.annotation.gtf"
OUT    <- file.path(ROOT, "deseq2"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

message("[01] Loading counts: ", COUNTS)
raw <- read_tsv(COUNTS, comment = "#", show_col_types = FALSE)
gene_info <- raw[, c("Geneid","Chr","Start","End","Strand","Length")]
counts    <- as.matrix(raw[, -(1:6)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"
message("  counts: ", nrow(counts), " genes x ", ncol(counts), " samples")
message("  samples: ", paste(colnames(counts), collapse=", "))

message("[01] Loading metadata: ", META)
meta <- read_tsv(META, show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$sample_id
stopifnot(all(colnames(counts) %in% meta$sample_id))
meta <- meta[colnames(counts), ]
meta$group  <- factor(meta$group,  levels = c("control_internal","control_public","u_BMF","g_BMF"))
meta$cohort <- factor(meta$cohort, levels = c("internal_aspirate","public_MNC"))
meta$source <- factor(meta$source)

# Combined "control" factor for primary contrasts (treats both control sets as one level)
meta$group_combined <- factor(
  ifelse(meta$group %in% c("control_internal","control_public"), "Control",
         as.character(meta$group)),
  levels = c("Control","u_BMF","g_BMF")
)

# Subgroup labels for FA-only / DKC sensitivity
meta$subgroup <- factor(
  case_when(
    meta$sample_id %in% c("AA-RNA-FA","AA-RNA-FA2","AA-RNA-FA3") ~ "FA",
    meta$sample_id == "AA-RNA-DKC" ~ "DKC",
    meta$group == "u_BMF" ~ "u_BMF",
    TRUE ~ "Control"
  ),
  levels = c("Control","u_BMF","FA","DKC")
)
print(table(meta$group, meta$cohort))
print(table(meta$subgroup, meta$cohort))

message("[01] Parsing GTF for gene biotype: ", GTF)
gtf <- import(GTF)
g <- gtf[gtf$type == "gene"]
gene_meta <- data.frame(
  gene_id        = g$gene_id,
  gene_name      = g$gene_name,
  gene_biotype   = g$gene_type,
  chr            = as.character(seqnames(g)),
  stringsAsFactors = FALSE
)
gene_meta <- gene_meta[!duplicated(gene_meta$gene_id), ]
rownames(gene_meta) <- gene_meta$gene_id
gene_meta <- gene_meta[rownames(counts), ]   # align order; missing -> NA
message("  biotype counts:")
print(table(gene_meta$gene_biotype, useNA = "ifany"))

message("[01] Building DESeqDataSet (design = ~ cohort + group_combined)")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ cohort + group_combined
)
mcols(dds) <- DataFrame(mcols(dds), gene_meta)

# Pre-filter: keep genes with >=10 reads in at least 3 samples (reasonable for 17-sample design)
keep <- rowSums(counts(dds) >= 10) >= 3
message(sprintf("  pre-filter: keeping %d / %d genes", sum(keep), length(keep)))
dds <- dds[keep, ]

saveRDS(dds, file.path(OUT, "dds_full.rds"))
saveRDS(gene_meta, file.path(OUT, "gene_meta.rds"))
saveRDS(meta, file.path(OUT, "samples_meta.rds"))
write_tsv(as.data.frame(meta), file.path(OUT, "samples_meta.tsv"))

message("[01] DONE. Saved to ", OUT)
