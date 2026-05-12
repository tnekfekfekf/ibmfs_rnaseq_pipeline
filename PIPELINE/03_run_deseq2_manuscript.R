#!/usr/bin/env Rscript
# Manuscript-style DESeq2 analysis (eAT5.R replica)
#   Design: ~ group
#   Controls: AA-PRO + AA-KEW (original 2)
#   Filter: rowSums > 0
#   Post-filter: padj != NA (DESeq2 independent filtering)
#   Output: per-contrast DEG tables (all + mRNA + lncRNA subset)
#
# Usage: Rscript 03_run_deseq2_manuscript.R [COUNTS_FILE] [GTF] [OUT_DIR]

suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble)})

args <- commandArgs(trailingOnly=TRUE)
COUNTS <- if (length(args)>=1) args[1] else "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_manuscript_v3.txt"
GTF    <- if (length(args)>=2) args[2] else "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"
OUT    <- if (length(args)>=3) args[3] else "/Users/jaeeunyoo/Desktop/star_workdir/deseq2_manuscript_style"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# Load counts
raw <- read_tsv(COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(raw[, -(1:6)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"
colnames(counts) <- sub("_sorted$","",sub("^.*/","",sub(".bam$","",colnames(counts))))

# Parse gene_type from GTF
cat("Parsing gene_type from GTF...\n")
gtf_lines <- readLines(GTF)
gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec('gene_id "([^"]+)".*?gene_type "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==3) data.frame(gene_id=x[2], gene_type=x[3]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE)
gt$gene_id <- as.character(gt$gene_id)

# Manuscript sample assignment (eAT5.R Oct 20 2025)
g_aa_samples    <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3")
control_samples <- c("AA-PRO","AA-KEW")
u_aa_samples    <- c("AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH")
samps <- c(g_aa_samples, control_samples, u_aa_samples)

stopifnot(all(samps %in% colnames(counts)))
counts <- counts[, samps]
meta <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8)),
                                 levels=c("Control","U-AA","G-AA")))
rownames(meta) <- samps

# Manuscript filter & DESeq2 (eAT5.R style)
counts <- counts[rowSums(counts) > 0, ]
cat(sprintf("After rowSums>0: %d genes × %d samples\n", nrow(counts), ncol(counts)))
dds <- DESeqDataSetFromMatrix(counts, meta, ~ group)
dds <- DESeq(dds, parallel=FALSE)
gene_types <- gt$gene_type[match(rownames(dds), gt$gene_id)]

# Helper: write DEG table per contrast + biotype
write_deg <- function(cc, label, biotype="all") {
  res <- results(dds, contrast=c("group", cc[1], cc[2]))
  df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% mutate(gene_type=gene_types)
  df <- df %>% filter(!is.na(padj))
  if (biotype == "mrna")   df <- df %>% filter(gene_type=="protein_coding")
  if (biotype == "lncrna") df <- df %>% filter(gene_type=="lncRNA")
  df$significance <- with(df, ifelse(padj<0.05 & log2FoldChange>1, "Upregulated",
                            ifelse(padj<0.05 & log2FoldChange < -1, "Downregulated", "Not significant")))
  out <- file.path(OUT, sprintf("DE_%s_%s.tsv", label, biotype))
  write_tsv(df, out)
  list(file=out, sig=sum(df$significance != "Not significant"),
       up=sum(df$significance=="Upregulated"),
       down=sum(df$significance=="Downregulated"),
       total=nrow(df))
}

cat("\n=== Results summary (|log2FC|>1, padj<0.05) ===\n")
cat(sprintf("%-22s %-8s %5s %5s %5s %6s\n", "Contrast", "Biotype", "Total", "Up", "Down", "Sig"))
for (cc in list(c("G-AA","Control"), c("U-AA","Control"), c("G-AA","U-AA"))) {
  label <- paste(cc[1], "vs", cc[2], sep="_")
  for (bt in c("all","mrna","lncrna")) {
    r <- write_deg(cc, label, bt)
    cat(sprintf("%-22s %-8s %5d %5d %5d %6d\n", label, bt, r$total, r$up, r$down, r$sig))
  }
}

cat(sprintf("\nAll DEG tables saved in: %s\n", OUT))
saveRDS(dds, file.path(OUT, "dds.rds"))
