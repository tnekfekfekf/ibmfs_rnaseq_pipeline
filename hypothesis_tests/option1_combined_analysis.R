#!/usr/bin/env Rscript
# Option 1 analysis: 14 manuscript samples (from authentic count matrix)
# + 3 Child samples (from our v3 quantification)
# Run 4 designs and compare results.
suppressPackageStartupMessages({library(DESeq2); library(dplyr); library(readr); library(tibble)})

MS_MATRIX <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt"
CHILD_V3_DIR <- "/Users/jaeeunyoo/Desktop/star_workdir/counts/per_sample_v3"
GTF <- "/Users/jaeeunyoo/Desktop/star_workdir/gencode.v44.annotation.no_rRNA.gtf"
OUT <- "/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/option1_analyses"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ===== Step 1: Build combined count matrix =====
cat("[1/4] Building combined matrix: 14 manuscript + 3 Child(v3)...\n")

# Load manuscript matrix
ms <- read_tsv(MS_MATRIX, show_col_types=FALSE)
manuscript_samps <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3",
                      "AA-PRO","AA-KEW",
                      "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18",
                      "AA-HMH","AA-PJH")
ms_counts <- as.matrix(ms[, manuscript_samps])
rownames(ms_counts) <- ms$EnsemblID
storage.mode(ms_counts) <- "integer"
cat(sprintf("  Manuscript counts: %d genes × %d samples\n", nrow(ms_counts), ncol(ms_counts)))

# Load our v3 Child counts (separately quantified)
child_samps <- c("Child1","Child2","Child3")
child_list <- list()
for (s in child_samps) {
  f <- file.path(CHILD_V3_DIR, paste0(s, ".counts.txt"))
  d <- read_tsv(f, comment="#", show_col_types=FALSE)
  child_list[[s]] <- setNames(d[[7]], d$Geneid)
}
common_genes <- Reduce(intersect, c(list(rownames(ms_counts)), lapply(child_list, names)))
ms_counts <- ms_counts[common_genes, ]
child_counts <- do.call(cbind, lapply(child_list, function(x) x[common_genes]))
colnames(child_counts) <- child_samps
cat(sprintf("  Child v3 counts: %d genes × %d samples\n", nrow(child_counts), ncol(child_counts)))

# Combine
counts_all <- cbind(ms_counts, child_counts)
storage.mode(counts_all) <- "integer"
cat(sprintf("  COMBINED: %d genes × %d samples\n", nrow(counts_all), ncol(counts_all)))

# Lib size compare
lib_sizes <- colSums(counts_all)/1e6
cat("  Library sizes (M):\n")
for (s in colnames(counts_all)) {
  cohort <- if (s %in% child_samps) "child_v3" else "manuscript"
  cat(sprintf("    %-12s %s  %.1fM\n", s, cohort, lib_sizes[s]))
}

# Save combined matrix
write_tsv(as.data.frame(cbind(EnsemblID=rownames(counts_all), counts_all)),
          file.path(OUT, "combined_14ms_3child_v3.txt"))

# ===== Step 2: Gene annotation =====
cat("\n[2/4] Loading gene annotation...\n")
gtf_lines <- readLines(GTF); gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
parsed <- regmatches(gene_lines, regexec(
  'gene_id "([^"]+)".*?gene_type "([^"]+)".*?gene_name "([^"]+)"', gene_lines))
gt <- do.call(rbind, lapply(parsed, function(x) if(length(x)==4)
  data.frame(gene_id=x[2], gene_type=x[3], gene_name=x[4]) else NULL))
gt <- gt %>% distinct(gene_id, .keep_all=TRUE); gt$gene_id <- as.character(gt$gene_id)

# ===== Step 3: 4 analyses =====
cat("\n[3/4] Running 4 DESeq2 analyses...\n")

run_deseq <- function(samp_subset, meta_df, design_formula, label) {
  cat(sprintf("  -- %s --\n", label))
  cts <- counts_all[, samp_subset]
  cts <- cts[rowSums(cts) > 0, ]
  dds <- DESeqDataSetFromMatrix(cts, meta_df, design=design_formula)
  dds <- DESeq(dds, parallel=FALSE)
  dds
}

# A: 14 manuscript samples only (authentic), ~ group — this IS manuscript replica
meta_A <- data.frame(
  group = factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8)),
                 levels=c("Control","U-AA","G-AA")),
  row.names = manuscript_samps
)
dds_A <- run_deseq(manuscript_samps, meta_A, ~ group, "A: 14 manuscript (authentic)")

# B: 14 ms + 3 child, naive ~ group
samps_BC <- c(manuscript_samps, child_samps)
meta_BC <- data.frame(
  group = factor(c(rep("G-AA",4),rep("Control",2),rep("U-AA",8),rep("Control",3)),
                 levels=c("Control","U-AA","G-AA")),
  cohort = factor(c(rep("manuscript",14), rep("child_v3",3)),
                  levels=c("manuscript","child_v3")),
  row.names = samps_BC
)
dds_B <- run_deseq(samps_BC, meta_BC, ~ group, "B: 14ms + 3child naive (n=17)")

# C: 14 ms + 3 child, batch-aware ~ cohort + group
dds_C <- run_deseq(samps_BC, meta_BC, ~ cohort + group, "C: 14ms + 3child batch-aware (n=17)")

# D: 12 manuscript patients + 3 child (no internal controls)
samps_D <- c("AA-RNA-FA","AA-RNA-DKC","AA-RNA-FA2","AA-RNA-FA3",
             "AA-RNA-1","AA-RNA-4","AA-RNA-5","AA-RNA-13","AA-RNA-16","AA-RNA-18",
             "AA-HMH","AA-PJH",
             "Child1","Child2","Child3")
meta_D <- data.frame(
  group = factor(c(rep("G-AA",4),rep("U-AA",8),rep("Control",3)),
                 levels=c("Control","U-AA","G-AA")),
  row.names = samps_D
)
dds_D <- run_deseq(samps_D, meta_D, ~ group, "D: 3 Child only as controls (n=15)")

# ===== Step 4: Extract DEG + compare =====
cat("\n[4/4] Computing DEG metrics + writing tables...\n")

get_res <- function(dds, g1, g2) {
  res <- results(dds, contrast=c("group", g1, g2))
  df <- as.data.frame(res) %>% rownames_to_column("gene_id")
  df$gene_name <- gt$gene_name[match(df$gene_id, gt$gene_id)]
  df$gene_type <- gt$gene_type[match(df$gene_id, gt$gene_id)]
  df %>% filter(!is.na(padj))
}

is_sig <- function(d, padj=0.05, lfc=1) !is.na(d$padj) & d$padj<padj & abs(d$log2FoldChange)>lfc

# Get all contrasts for all analyses
res_list <- list()
for (lab in c("A","B","C","D")) {
  dds <- get(paste0("dds_", lab))
  for (cc in list(c("G-AA","Control"), c("U-AA","Control"), c("G-AA","U-AA"))) {
    res_list[[paste(lab, cc[1], cc[2], sep="_")]] <- get_res(dds, cc[1], cc[2])
  }
}

# === TABLE 1: DEG counts (padj<0.05, |LFC|>1) ===
cat("\n=== TABLE 1: DEG counts per analysis & contrast ===\n")
tab1 <- data.frame()
for (lab in c("A","B","C","D")) {
  for (cc_str in c("G-AA_Control","U-AA_Control","G-AA_U-AA")) {
    k <- paste(lab, cc_str, sep="_")
    d <- res_list[[k]]
    sig <- is_sig(d)
    tab1 <- rbind(tab1, data.frame(
      Analysis=lab, Contrast=gsub("_"," vs ", cc_str),
      n_tested=nrow(d),
      n_sig=sum(sig),
      n_mRNA=sum(sig & d$gene_type=="protein_coding"),
      n_lncRNA=sum(sig & d$gene_type=="lncRNA"),
      n_up=sum(sig & d$log2FoldChange>0),
      n_down=sum(sig & d$log2FoldChange<0)
    ))
  }
}
print(tab1, row.names=FALSE)
write_tsv(tab1, file.path(OUT, "table1_deg_counts.tsv"))

# === TABLE 2: Gene-level overlap with Analysis A ===
cat("\n=== TABLE 2: Gene-level overlap with Analysis A (reference) ===\n")
tab2 <- data.frame()
for (cc_str in c("G-AA_Control","U-AA_Control","G-AA_U-AA")) {
  A_sig <- res_list[[paste("A", cc_str, sep="_")]]$gene_id[is_sig(res_list[[paste("A", cc_str, sep="_")]])]
  for (lab in c("B","C","D")) {
    other_sig <- res_list[[paste(lab, cc_str, sep="_")]]$gene_id[is_sig(res_list[[paste(lab, cc_str, sep="_")]])]
    overlap <- intersect(A_sig, other_sig)
    tab2 <- rbind(tab2, data.frame(
      Contrast=gsub("_"," vs ", cc_str),
      Compared=paste("A vs", lab),
      A_sig=length(A_sig), other_sig=length(other_sig), Overlap=length(overlap),
      Pct_A_overlap=ifelse(length(A_sig)>0, sprintf("%.1f%%", 100*length(overlap)/length(A_sig)), "NA")
    ))
  }
}
print(tab2, row.names=FALSE)
write_tsv(tab2, file.path(OUT, "table2_overlap_with_A.tsv"))

# === TABLE 3: Key manuscript genes — sig status across analyses (G-AA vs Control) ===
cat("\n=== TABLE 3: Key manuscript-highlighted genes (G-AA vs Ctrl) ===\n")
key_genes <- c("HCG11","HCP5","SNHG32","PSMB8-AS1","FAM30A","MIR22HG",
               "ATP1A1-AS1","USP3-AS1","TAGAP-AS1","LINC01036","MALAT1",
               "TEN1-CDK3","FANCA","FANCG","TERT")
tab3 <- data.frame(Gene=key_genes)
for (lab in c("A","B","C","D")) {
  d <- res_list[[paste(lab, "G-AA_Control", sep="_")]]
  out_lfc <- c()
  out_padj <- c()
  out_sig <- c()
  for (sym in key_genes) {
    gid <- gt$gene_id[gt$gene_name == sym][1]
    if (is.na(gid) || !gid %in% d$gene_id) {
      out_lfc <- c(out_lfc, NA); out_padj <- c(out_padj, NA); out_sig <- c(out_sig, NA)
    } else {
      r <- d[d$gene_id == gid, ]
      out_lfc <- c(out_lfc, round(r$log2FoldChange, 2))
      out_padj <- c(out_padj, signif(r$padj, 3))
      out_sig <- c(out_sig, ifelse(is_sig(r), "✓", "·"))
    }
  }
  tab3[[paste0(lab,"_LFC")]] <- out_lfc
  tab3[[paste0(lab,"_padj")]] <- out_padj
  tab3[[paste0(lab,"_sig")]] <- out_sig
}
print(tab3, row.names=FALSE)
write_tsv(tab3, file.path(OUT, "table3_key_genes_GvC.tsv"))

cat(sprintf("\nAll outputs saved to: %s\n", OUT))
saveRDS(list(A=dds_A, B=dds_B, C=dds_C, D=dds_D, res=res_list),
        file.path(OUT, "all_dds_objects.rds"))
