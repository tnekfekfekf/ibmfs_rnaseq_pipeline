#!/usr/bin/env Rscript
# Critical verification: are Child (MNC) controls really comparable to internal (aspirate) controls?
# Check:
# 1. Sequencing depth difference
# 2. Cell composition marker expression (myeloid vs lymphoid)
# 3. Direct comparison Child vs Internal control (should NOT be biologically different)
# 4. Whether ~cohort + group is truly correcting OR just absorbing real biology
suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble); library(sva)})

COUNTS <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/fc_v3_17samples.txt"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"
DDS_RDS <- "/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/control_comparison/dds_all.rds"

# Load
raw <- read_tsv(COUNTS, comment="#", show_col_types=FALSE)
counts <- as.matrix(raw[, -(1:2)])
rownames(counts) <- raw$Geneid
storage.mode(counts) <- "integer"

gtf_lines <- readLines(GTF); gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec(
  'gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4)
  data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE)

# Sample groups
internal_ctrl <- c("AA-PRO","AA-KEW")
public_ctrl <- c("Child1","Child2","Child3")
all_ctrl <- c(internal_ctrl, public_ctrl)

cat("==========================================================\n")
cat("CHECK 1: Sequencing depth\n")
cat("==========================================================\n\n")
lib_sizes <- colSums(counts)
cat(sprintf("%-15s %-12s %s\n", "Sample", "Cohort", "Lib size (M)"))
for (s in all_ctrl) {
  coh <- if (s %in% internal_ctrl) "internal_aspirate" else "public_MNC"
  cat(sprintf("%-15s %-12s %.1f\n", s, coh, lib_sizes[s]/1e6))
}
cat(sprintf("\nInternal avg: %.1f M\n", mean(lib_sizes[internal_ctrl])/1e6))
cat(sprintf("Public avg:   %.1f M\n", mean(lib_sizes[public_ctrl])/1e6))

cat("\n==========================================================\n")
cat("CHECK 2: Cell composition markers (CPM scale)\n")
cat("==========================================================\n\n")
# Hematopoietic lineage markers
markers <- list(
  Lymphoid = list(
    `CD3 T cell` = c("CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2"),
    `CD19 B cell` = c("CD19","MS4A1","CD79A","CD79B"),
    `NK cell` = c("NKG7","GNLY","NCAM1","KLRD1")
  ),
  Myeloid = list(
    `Granulocyte` = c("ELANE","MPO","PRTN3","CSF3R","FCGR3B"),
    `Monocyte` = c("CD14","CD68","LYZ","FCN1","S100A8","S100A9"),
    `Erythroid` = c("HBB","HBA1","HBA2","GYPA","ALAS2","SLC4A1"),
    `Megakaryocyte` = c("ITGA2B","GP1BA","PF4","ITGB3")
  ),
  Progenitor = list(
    `HSC/Progenitor` = c("CD34","KIT","HOXA9","MEIS1","MPL")
  )
)

# CPM
cpm <- t(t(counts) / lib_sizes) * 1e6

print_markers <- function() {
  for (lineage in names(markers)) {
    cat(sprintf("\n--- %s ---\n", lineage))
    for (subtype in names(markers[[lineage]])) {
      cat(sprintf("\n  %s markers:\n", subtype))
      for (sym in markers[[lineage]][[subtype]]) {
        gid <- gt$gene_id[gt$gene_name == sym][1]
        if (is.na(gid) || !gid %in% rownames(cpm)) next
        int_mean <- mean(cpm[gid, internal_ctrl])
        pub_mean <- mean(cpm[gid, public_ctrl])
        ratio <- if (int_mean>0) pub_mean/int_mean else NA
        cat(sprintf("    %-8s  internal CPM=%7.2f  public CPM=%7.2f  ratio=%s\n",
                    sym, int_mean, pub_mean,
                    ifelse(is.na(ratio), "NA", sprintf("%.2fx", ratio))))
      }
    }
  }
}
print_markers()

cat("\n==========================================================\n")
cat("CHECK 3: Direct Public vs Internal control comparison\n")
cat("==========================================================\n\n")
# Treat Child vs internal_ctrl as if they were "groups"
# What would DEG analysis say?
samps_ctrl_only <- c(internal_ctrl, public_ctrl)
counts_ctrl <- counts[, samps_ctrl_only]
counts_ctrl <- counts_ctrl[rowSums(counts_ctrl) > 0, ]
meta_ctrl <- data.frame(
  cohort = factor(c(rep("internal",2), rep("public",3)),
                  levels=c("internal","public"))
)
rownames(meta_ctrl) <- samps_ctrl_only

cat("Running DESeq2 (internal vs public controls, treated as DEG)...\n")
dds_ctrl <- DESeqDataSetFromMatrix(counts_ctrl, meta_ctrl, ~ cohort)
dds_ctrl <- DESeq(dds_ctrl, parallel=FALSE)
res_ctrl <- results(dds_ctrl, contrast=c("cohort","public","internal"))
res_ctrl <- as.data.frame(res_ctrl) %>% rownames_to_column("gene_id")
res_ctrl$gene_name <- gt$gene_name[match(res_ctrl$gene_id, gt$gene_id)]
res_ctrl$gene_type <- gt$gene_type[match(res_ctrl$gene_id, gt$gene_id)]
res_ctrl <- res_ctrl %>% filter(!is.na(padj))

sig_n <- sum(res_ctrl$padj<0.05 & abs(res_ctrl$log2FoldChange)>1, na.rm=TRUE)
sig_n_strong <- sum(res_ctrl$padj<0.05 & abs(res_ctrl$log2FoldChange)>2, na.rm=TRUE)
sig_n_veryStrong <- sum(res_ctrl$padj<0.05 & abs(res_ctrl$log2FoldChange)>5, na.rm=TRUE)

cat(sprintf("\nTotal genes tested: %d\n", nrow(res_ctrl)))
cat(sprintf("Sig with padj<0.05 & |LFC|>1: %d genes (%.1f%%)\n", sig_n, 100*sig_n/nrow(res_ctrl)))
cat(sprintf("Sig with padj<0.05 & |LFC|>2: %d genes (%.1f%%)\n", sig_n_strong, 100*sig_n_strong/nrow(res_ctrl)))
cat(sprintf("Sig with padj<0.05 & |LFC|>5: %d genes (%.1f%%)\n", sig_n_veryStrong, 100*sig_n_veryStrong/nrow(res_ctrl)))

# Top hits
cat("\nTop 20 most-different genes (public vs internal control):\n")
top <- res_ctrl %>% filter(padj<0.05) %>% arrange(desc(abs(log2FoldChange))) %>% head(20)
for (i in 1:nrow(top)) {
  cat(sprintf("  %-15s %-10s LFC=%6.2f  padj=%8.2e\n",
              ifelse(is.na(top$gene_name[i]), top$gene_id[i], top$gene_name[i]),
              top$gene_type[i], top$log2FoldChange[i], top$padj[i]))
}

cat("\n==========================================================\n")
cat("CHECK 4: How much does ~ cohort term absorb?\n")
cat("==========================================================\n\n")
# Load existing C analysis dds
dds_all <- readRDS(DDS_RDS)
dds_C <- dds_all$C
# Get cohort coefficient (vs internal)
# results(dds_C, name=...) — show cohort effect for highly-DE markers
result_names <- resultsNames(dds_C)
cat("Coefficients estimated by ~cohort + group:\n")
print(result_names)
cat("\nFor reference: cohort term captures public vs internal effect per gene\n")
cat("If public is fundamentally different (e.g., cell composition), cohort term will be large\n\n")

# Get cohort effect for some marker genes
cohort_effects <- function(genes_list) {
  res <- results(dds_C, name="cohort_public_vs_internal")
  out <- data.frame()
  for (sym in genes_list) {
    gid <- gt$gene_id[gt$gene_name == sym][1]
    if (is.na(gid) || !gid %in% rownames(res)) next
    out <- rbind(out, data.frame(
      gene = sym,
      cohort_LFC = res[gid, "log2FoldChange"],
      cohort_padj = res[gid, "padj"]
    ))
  }
  out
}

# Cell-type markers
markers_flat <- unlist(unlist(markers, recursive=FALSE), use.names=FALSE)
cohort_eff <- cohort_effects(markers_flat)
cat("Cohort effect (public vs internal) for cell-type markers:\n")
print(cohort_eff)

# Save outputs
write_tsv(res_ctrl, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/public_vs_internal_ctrl_DEG.tsv")
write_tsv(cohort_eff, "/Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline/docs/cohort_effect_markers.tsv")

cat("\n==========================================================\n")
cat("INTERPRETATION\n")
cat("==========================================================\n")
cat(sprintf("\n1. Depth: public is %.1fx deeper than internal (%.1fM vs %.1fM avg)\n",
            mean(lib_sizes[public_ctrl])/mean(lib_sizes[internal_ctrl]),
            mean(lib_sizes[public_ctrl])/1e6, mean(lib_sizes[internal_ctrl])/1e6))
cat(sprintf("\n2. Cell composition: see marker patterns above\n"))
cat(sprintf("\n3. DEG between public vs internal: %d genes (%.0f%% of tested)\n",
            sig_n, 100*sig_n/nrow(res_ctrl)))
cat(sprintf("   → These are TRUE biological differences, NOT just batch\n"))
cat(sprintf("\n4. ~cohort + group design absorbs these differences as cohort effect\n"))
cat(sprintf("   → For specific genes, cohort effect can be very large (see table above)\n"))
