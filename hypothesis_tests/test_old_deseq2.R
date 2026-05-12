#!/usr/bin/env Rscript
# Test reproducibility with DESeq2 1.48.2 (Bioc 3.21, Apr 2025 release)
# Same v3 count matrix, same eAT5.R-style analysis
# Check if ATP1A1-AS1, USP3-AS1, etc. become sig

# Force old library
LIB <- "/Users/jaeeunyoo/Desktop/star_workdir/R_libs_old"
.libPaths(c(LIB, .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2, lib.loc=LIB)
  library(dplyr); library(readr); library(tibble)
})
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
cat("Library:", normalizePath(find.package("DESeq2")), "\n\n")

COUNTS <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_manuscript_v3.txt"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"

# Load counts
raw <- read_tsv(COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(raw[, -(1:6)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"
colnames(counts) <- sub("_sorted$","",sub("^.*/","",sub(".bam$","",colnames(counts))))

# Manuscript samples
samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3","AA-PRO","AA-KEW",
           "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18","AA-HMH","AA-PJH")
counts <- counts[, samps]; counts <- counts[rowSums(counts) > 0, ]
meta <- data.frame(group=factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8)),
                                 levels=c("Control","U-AA","G-AA")))
rownames(meta) <- samps

cat(sprintf("After rowSums>0: %d genes × %d samples\n", nrow(counts), ncol(counts)))
cat("Running DESeq2 1.48.2...\n")
dds <- DESeqDataSetFromMatrix(counts, meta, ~ group)
dds <- DESeq(dds, parallel=FALSE)
cat("Done.\n\n")

# Parse gene_type
gtf_lines <- readLines(GTF); gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec(
  'gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4)
  data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE)
gt$gene_id <- as.character(gt$gene_id)

# Check specific genes
key_genes <- list(
  G_AA_vs_UAA_4 = c("SFT2D3","OR52K1","LRRC24","TEN1-CDK3"),
  Key_11_lncRNAs = c("ATP1A1-AS1","LINC01036","HCG11","HCP5","SNHG32","PSMB8-AS1",
                     "TAGAP-AS1","MALAT1","FAM30A","USP3-AS1","MIR22HG"),
  FA_genes = c("FANCG","FANCD2","FANCA","TERT")
)

# Manuscript DEG for comparison
ms_GvC <- read_tsv("/Volumes/ExtremeSSD/ibmfs/03_original_analysis/deg_results/G-AA_vs_Control_results.txt",
                   show_col_types=FALSE)

# Get our v1.50 results for the same genes (run prior comparison)
# We'll just compare old DESeq2 with manuscript directly
get_res <- function(g1, g2) {
  res <- results(dds, contrast=c("group", g1, g2))
  df <- as.data.frame(res) %>% rownames_to_column("gene_id")
  df$gene_name <- gt$gene_name[match(df$gene_id, gt$gene_id)]
  df$gene_type <- gt$gene_type[match(df$gene_id, gt$gene_id)]
  df %>% filter(!is.na(padj))
}

cat("========== COMPARISON: Manuscript vs DESeq2 1.48 (our v3 counts) ==========\n")
gvc_148 <- get_res("G-AA","Control")
uvc_148 <- get_res("U-AA","Control")
gvu_148 <- get_res("G-AA","U-AA")

# Compare key genes
cat("\n=== Critical genes (manuscript values vs DESeq2 1.48) ===\n")
cat(sprintf("%-15s | %-22s | %-22s | Status\n", "Gene", "Manuscript (G-AA)", "Our v3 + DESeq2 1.48"))
cat(paste(rep("-", 90), collapse=""), "\n")

for (cat_name in names(key_genes)) {
  cat(sprintf("\n[%s]\n", cat_name))
  for (sym in key_genes[[cat_name]]) {
    gene_id <- gt$gene_id[gt$gene_name == sym][1]
    if (is.na(gene_id)) { cat(sprintf("  %-15s NOT FOUND\n", sym)); next }
    ms_row <- ms_GvC[ms_GvC$GeneSymbol == sym, ][1, ]
    our_row <- gvc_148[gvc_148$gene_id == gene_id, ][1, ]
    if (nrow(ms_row) == 0 || is.na(our_row$padj)) next

    ms_sig <- !is.na(ms_row$padj) && ms_row$padj < 0.05 && abs(ms_row$log2FoldChange) > 1
    our_sig <- !is.na(our_row$padj) && our_row$padj < 0.05 && abs(our_row$log2FoldChange) > 1
    status <- ifelse(ms_sig == our_sig, "✅ MATCH", "❌ DIFFER")
    if (ms_sig && our_sig && abs(ms_row$log2FoldChange - our_row$log2FoldChange) < 0.5) {
      status <- "✅ PERFECT"
    }

    cat(sprintf("  %-15s | LFC=%6.2f padj=%8.2e | LFC=%6.2f padj=%8.2e | %s\n",
                sym, ms_row$log2FoldChange, ms_row$padj,
                our_row$log2FoldChange, our_row$padj, status))
  }
}

# Overall summary
ms_sig_all <- ms_GvC$gene_id[!is.na(ms_GvC$padj) & ms_GvC$padj<0.05 & abs(ms_GvC$log2FoldChange)>1]
our_sig_all <- gvc_148$gene_id[gvc_148$padj<0.05 & abs(gvc_148$log2FoldChange)>1]
overlap <- intersect(ms_sig_all, our_sig_all)
cat(sprintf("\n=== G-AA vs Control overall ===\n"))
cat(sprintf("Manuscript sig: %d\n", length(ms_sig_all)))
cat(sprintf("Our v3 + DESeq2 1.48 sig: %d\n", length(our_sig_all)))
cat(sprintf("Overlap: %d (recall %.1f%%, precision %.1f%%)\n",
            length(overlap), 100*length(overlap)/length(ms_sig_all),
            100*length(overlap)/length(our_sig_all)))

# By gene_type
sig_148 <- gvc_148[gvc_148$padj<0.05 & abs(gvc_148$log2FoldChange)>1, ]
mrna_n <- sum(sig_148$gene_type == "protein_coding", na.rm=TRUE)
lnc_n  <- sum(sig_148$gene_type == "lncRNA", na.rm=TRUE)
cat(sprintf("Our v3 + DESeq2 1.48 by biotype: mRNA=%d, lncRNA=%d\n", mrna_n, lnc_n))
cat(sprintf("Manuscript: mRNA=2078, lncRNA=1167\n"))

# Save DEG tables
write_tsv(gvc_148, "/Users/jaeeunyoo/Desktop/star_workdir/DE_GvC_DESeq2_148.tsv")
write_tsv(uvc_148, "/Users/jaeeunyoo/Desktop/star_workdir/DE_UvC_DESeq2_148.tsv")
write_tsv(gvu_148, "/Users/jaeeunyoo/Desktop/star_workdir/DE_GvU_DESeq2_148.tsv")
cat("\nSaved DEG tables to /Users/jaeeunyoo/Desktop/star_workdir/\n")
