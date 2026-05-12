#!/usr/bin/env Rscript
# ============================================================
# Enhanced DESeq2 + Pathway Analysis for Aplastic Anemia RNA-seq
# Samples: AA-HMH (healthy_control), AA-RNA-1 (aplastic_anemia), AA-RNA-FA (fanconi_anemia)
# Author: Analysis pipeline for IBMFS study
# Date: 2026-03-09
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) args[[1]] else "/Users/jaeeunyoo/Desktop/researches/ibmfs_fastq_raw_data"

counts_path <- file.path(root_dir, "counts", "featureCounts.cleaned.txt")
meta_path   <- file.path(root_dir, "metadata", "samples.tsv")
out_dir     <- file.path(root_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_path)) stop("Missing counts file: ", counts_path)
if (!file.exists(meta_path))   stop("Missing metadata file: ", meta_path)

# ── 1. Load libraries ─────────────────────────────────────────
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

cat(">>> Loading count matrix...\n")
counts_tbl <- read.table(counts_path, header = TRUE, sep = "\t",
                         quote = "", comment.char = "", check.names = FALSE)
gene_id    <- counts_tbl[[1]]
count_mat  <- as.matrix(counts_tbl[, 7:ncol(counts_tbl)])
rownames(count_mat) <- gene_id
storage.mode(count_mat) <- "integer"

# ── 2. Load metadata ──────────────────────────────────────────
cat(">>> Loading metadata...\n")
meta <- read.table(meta_path, header = TRUE, sep = "\t",
                   quote = "", comment.char = "", stringsAsFactors = FALSE)

# Align samples
colnames(count_mat) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(count_mat))
meta <- meta[match(colnames(count_mat), meta$sample), , drop = FALSE]
if (any(is.na(meta$sample))) stop("Metadata samples do not match count columns. Check names.")

# Set reference level to healthy_control
meta$condition <- factor(meta$condition, levels = c("healthy_control", "aplastic_anemia", "fanconi_anemia"))

cat(sprintf(">>> %d genes x %d samples\n", nrow(count_mat), ncol(count_mat)))
cat("    Conditions:", paste(meta$condition, collapse = ", "), "\n")

# ── 3. Build DESeqDataSet ─────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                               colData   = meta,
                               design    = ~ condition)
# Pre-filter: keep genes with >= 10 total counts
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat(sprintf(">>> After pre-filtering: %d genes retained\n", nrow(dds)))

# ── 4. VST normalization ──────────────────────────────────────
cat(">>> Computing VST...\n")
vsd <- vst(dds, blind = TRUE)
write.csv(as.data.frame(assay(vsd)),
          file.path(out_dir, "vst_matrix.csv"))

# ── 5. PCA plot ───────────────────────────────────────────────
cat(">>> Generating PCA plot...\n")
pca_data    <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.5, fontface = "bold", max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("healthy_control" = "#2ecc71",
                                "aplastic_anemia"  = "#e74c3c",
                                "fanconi_anemia"   = "#3498db")) +
  ggtitle("PCA — VST-normalized counts") +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(out_dir, "pca_vst.png"), p_pca, width = 7, height = 6, dpi = 150)

# ── 6. Sample-to-sample distance heatmap ─────────────────────
cat(">>> Generating sample distance heatmap...\n")
sampleDists <- dist(t(assay(vsd)))
sampleDistMat <- as.matrix(sampleDists)
colnames(sampleDistMat) <- rownames(sampleDistMat) <- meta$condition

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png(file.path(out_dir, "sample_distance_heatmap.png"), width = 700, height = 600, res = 100)
pheatmap(sampleDistMat, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Euclidean Distances (VST)")
dev.off()

# ── 7. DESeq2 differential expression ────────────────────────
cat(">>> Running DESeq2...\n")
dds <- DESeq(dds)

# Helper: save DE results with annotation
save_deseq_results <- function(dds, contrast, label, out_dir) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  write.csv(res_df, file.path(out_dir, paste0("deseq2_", label, ".csv")), row.names = FALSE)

  # Volcano plot
  res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                       ifelse(res_df$log2FoldChange > 1, "Up", "Down"), "NS")
  top_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ] |>
    arrange(padj) |> head(20)

  p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pmax(padj, 1e-300)), color = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(data = top_genes,
                    aes(label = gene_id), size = 2.5, max.overlaps = 15,
                    color = "black", fontface = "italic") +
    scale_color_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "grey60")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    xlab("log2 Fold Change") + ylab("-log10(adjusted p-value)") +
    ggtitle(paste0("Volcano: ", gsub("_", " ", label))) +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(out_dir, paste0("volcano_", label, ".png")), p_vol,
         width = 7, height = 6, dpi = 150)

  # MA plot
  png(file.path(out_dir, paste0("MA_", label, ".png")), width = 700, height = 550, res = 100)
  plotMA(res, main = paste0("MA Plot: ", gsub("_", " ", label)), alpha = 0.05)
  dev.off()

  invisible(res_df)
}

# Contrast 1: aplastic_anemia vs healthy_control
cat("  Contrast 1: aplastic_anemia vs healthy_control\n")
res1 <- save_deseq_results(dds,
  contrast = c("condition", "aplastic_anemia", "healthy_control"),
  label    = "AA_vs_HC",
  out_dir  = out_dir)

# Contrast 2: fanconi_anemia vs healthy_control
cat("  Contrast 2: fanconi_anemia vs healthy_control\n")
res2 <- save_deseq_results(dds,
  contrast = c("condition", "fanconi_anemia", "healthy_control"),
  label    = "FA_vs_HC",
  out_dir  = out_dir)

# Contrast 3: fanconi_anemia vs aplastic_anemia
cat("  Contrast 3: fanconi_anemia vs aplastic_anemia\n")
res3 <- save_deseq_results(dds,
  contrast = c("condition", "fanconi_anemia", "aplastic_anemia"),
  label    = "FA_vs_AA",
  out_dir  = out_dir)

# ── 8. Top DEG heatmap ────────────────────────────────────────
cat(">>> Generating top DEG heatmap...\n")
sig_genes <- unique(c(
  res1[!is.na(res1$padj) & res1$padj < 0.05 & abs(res1$log2FoldChange) > 1, "gene_id"][1:50],
  res2[!is.na(res2$padj) & res2$padj < 0.05 & abs(res2$log2FoldChange) > 1, "gene_id"][1:50]
))
sig_genes <- sig_genes[!is.na(sig_genes)]

if (length(sig_genes) > 5) {
  vst_mat <- assay(vsd)[sig_genes, , drop = FALSE]
  vst_mat <- vst_mat - rowMeans(vst_mat)

  anno_col <- data.frame(Condition = meta$condition)
  rownames(anno_col) <- colnames(vst_mat)
  anno_colors <- list(
    Condition = c(healthy_control = "#2ecc71",
                  aplastic_anemia  = "#e74c3c",
                  fanconi_anemia   = "#3498db"))

  png(file.path(out_dir, "top_DEG_heatmap.png"), width = 900, height = 1000, res = 100)
  pheatmap(vst_mat,
           annotation_col  = anno_col,
           annotation_colors = anno_colors,
           show_rownames   = TRUE,
           show_colnames   = TRUE,
           fontsize_row    = 7,
           scale           = "row",
           color           = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
           main            = "Top Differentially Expressed Genes (row-scaled VST)")
  dev.off()
}

# ── 9. Gene Ontology enrichment (requires clusterProfiler + org.Hs.eg.db) ────
cat(">>> Attempting Gene Ontology enrichment...\n")
go_enrichment_attempted <- FALSE

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db",    quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  go_enrichment_attempted <- TRUE

  run_go <- function(gene_list, label, out_dir) {
    # Convert Ensembl IDs to Entrez (strip version suffix if present)
    ensembl_ids <- sub("\\..*$", "", gene_list)
    eg <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db, drop = TRUE)
    if (nrow(eg) == 0) {
      cat("  [WARN] No Entrez IDs mapped for", label, "\n"); return(invisible(NULL))
    }
    # GO Biological Process
    ego <- enrichGO(gene          = eg$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE)
    if (!is.null(ego) && nrow(ego) > 0) {
      write.csv(as.data.frame(ego), file.path(out_dir, paste0("GO_BP_", label, ".csv")), row.names = FALSE)
      p_go <- dotplot(ego, showCategory = 20, title = paste0("GO BP — ", gsub("_", " ", label)))
      ggsave(file.path(out_dir, paste0("GO_BP_dotplot_", label, ".png")), p_go,
             width = 10, height = 8, dpi = 150)
    }
    # KEGG
    ekegg <- enrichKEGG(gene         = eg$ENTREZID,
                        organism     = "hsa",
                        pvalueCutoff = 0.05)
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      write.csv(as.data.frame(ekegg), file.path(out_dir, paste0("KEGG_", label, ".csv")), row.names = FALSE)
      p_kegg <- dotplot(ekegg, showCategory = 20, title = paste0("KEGG — ", gsub("_", " ", label)))
      ggsave(file.path(out_dir, paste0("KEGG_dotplot_", label, ".png")), p_kegg,
             width = 10, height = 8, dpi = 150)
    }
  }

  # Upregulated genes in AA vs HC
  up_AA_HC <- res1[!is.na(res1$padj) & res1$padj < 0.05 & res1$log2FoldChange > 1, "gene_id"]
  if (length(up_AA_HC) > 5) run_go(up_AA_HC, "UP_AA_vs_HC", out_dir)

  # Upregulated genes in FA vs HC
  up_FA_HC <- res2[!is.na(res2$padj) & res2$padj < 0.05 & res2$log2FoldChange > 1, "gene_id"]
  if (length(up_FA_HC) > 5) run_go(up_FA_HC, "UP_FA_vs_HC", out_dir)

  cat("  GO/KEGG enrichment done.\n")
} else {
  cat("  [INFO] clusterProfiler or org.Hs.eg.db not installed. Skipping GO/KEGG.\n")
  cat("  Install with: BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))\n")
}

# ── 10. GSEA ranked analysis (clusterProfiler) ───────────────
cat(">>> Attempting GSEA...\n")
if (go_enrichment_attempted) {
  run_gsea <- function(res_df, label, out_dir) {
    res_df <- res_df[!is.na(res_df$stat), ]
    gene_list <- res_df$stat
    names(gene_list) <- sub("\\..*$", "", res_df$gene_id)   # strip Ensembl version
    eg_map <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db, drop = TRUE)
    gene_list2 <- gene_list[names(gene_list) %in% eg_map$ENSEMBL]
    names(gene_list2) <- eg_map$ENTREZID[match(names(gene_list2), eg_map$ENSEMBL)]
    # Remove duplicate Entrez IDs — keep highest absolute stat value
    gene_list2 <- gene_list2[!duplicated(names(gene_list2))]
    gene_list2 <- sort(gene_list2, decreasing = TRUE)

    gsea_res <- gseKEGG(geneList     = gene_list2,
                        organism     = "hsa",
                        minGSSize    = 10,
                        pvalueCutoff = 0.2,
                        verbose      = FALSE)
    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      write.csv(as.data.frame(gsea_res),
                file.path(out_dir, paste0("GSEA_KEGG_", label, ".csv")), row.names = FALSE)
      p_gsea <- dotplot(gsea_res, showCategory = 15,
                        title = paste0("GSEA KEGG — ", gsub("_", " ", label)), split = ".sign") +
        facet_grid(. ~ .sign)
      ggsave(file.path(out_dir, paste0("GSEA_KEGG_dotplot_", label, ".png")), p_gsea,
             width = 12, height = 8, dpi = 150)
    }
  }
  run_gsea(res1, "AA_vs_HC", out_dir)
  run_gsea(res2, "FA_vs_HC", out_dir)
  cat("  GSEA done.\n")
}

# ── 11. Hematopoiesis & Bone Marrow Failure gene signature ───
cat(">>> Checking key IBMFS / bone marrow failure signature genes...\n")

ibmfs_genes <- list(
  Fanconi_Anemia_Core  = c("FANCA", "FANCB", "FANCC", "FANCD1", "FANCD2",
                            "FANCE", "FANCF", "FANCG", "FANCI", "FANCJ",
                            "FANCL", "FANCM", "FANCN", "FANCP", "FANCQ",
                            "FANCR", "FANCS", "FANCT", "FANCU", "FANCV", "FANCW"),
  DNA_Damage_Response  = c("TP53", "ATM", "ATR", "BRCA1", "BRCA2", "CHEK1", "CHEK2",
                            "RAD51", "H2AFX", "CDKN1A", "MDM2"),
  Hematopoiesis_TF     = c("GATA1", "GATA2", "TAL1", "KLF1", "NFE2",
                            "RUNX1", "CEBPA", "SPI1", "FLI1", "ETV6"),
  Cytokine_Signaling   = c("IFNG", "TNF", "IL2", "IL6", "IL10", "TGFB1",
                            "TGFB2", "TGFB3", "CXCL10", "CCL5"),
  HSC_Markers          = c("CD34", "THY1", "KIT", "FLT3", "PROM1", "PTPRC",
                            "HOXA9", "MEIS1", "MPL"),
  Apoptosis            = c("BCL2", "BAX", "BCL2L1", "CASP3", "CASP9",
                            "CYCS", "APAF1", "FAS", "FASLG")
)

# Extract VST values for signature genes (match by gene symbol – may need conversion)
# NOTE: featureCounts output uses Ensembl gene IDs; if symbols are used, matching depends on GTF
# This section will work if rownames are gene symbols or if you convert them first.
sig_check <- lapply(names(ibmfs_genes), function(group) {
  genes <- ibmfs_genes[[group]]
  found <- intersect(genes, rownames(assay(vsd)))
  data.frame(Group = group,
             Gene  = genes,
             Found_in_data = genes %in% rownames(assay(vsd)))
})
sig_check_df <- do.call(rbind, sig_check)
write.csv(sig_check_df, file.path(out_dir, "ibmfs_signature_genes_check.csv"), row.names = FALSE)

found_genes <- sig_check_df[sig_check_df$Found_in_data, "Gene"]
if (length(found_genes) >= 3) {
  vst_sig <- assay(vsd)[found_genes, , drop = FALSE]
  vst_sig <- vst_sig - rowMeans(vst_sig)

  anno_col <- data.frame(Condition = meta$condition)
  rownames(anno_col) <- colnames(vst_sig)

  png(file.path(out_dir, "ibmfs_signature_heatmap.png"), width = 900, height = 800, res = 100)
  pheatmap(vst_sig,
           annotation_col  = anno_col,
           show_rownames   = TRUE,
           fontsize_row    = 8,
           scale           = "row",
           color           = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
           main            = "IBMFS / Hematopoiesis Signature Genes (row-scaled VST)")
  dev.off()
}

# ── 12. Summary report ────────────────────────────────────────
cat(">>> Writing summary report...\n")
n_sig_AA_HC <- sum(!is.na(res1$padj) & res1$padj < 0.05, na.rm = TRUE)
n_sig_FA_HC <- sum(!is.na(res2$padj) & res2$padj < 0.05, na.rm = TRUE)
n_sig_FA_AA <- sum(!is.na(res3$padj) & res3$padj < 0.05, na.rm = TRUE)

summary_lines <- c(
  "=== DESeq2 Analysis Summary ===",
  sprintf("Samples analyzed: %s", paste(meta$sample, collapse = ", ")),
  sprintf("Conditions: %s", paste(levels(meta$condition), collapse = " | ")),
  sprintf("Genes tested (post-filter): %d", nrow(dds)),
  "",
  sprintf("DEGs (padj < 0.05): AA vs Healthy Control: %d", n_sig_AA_HC),
  sprintf("DEGs (padj < 0.05): FA vs Healthy Control: %d", n_sig_FA_HC),
  sprintf("DEGs (padj < 0.05): FA vs AA:               %d", n_sig_FA_AA),
  "",
  "Output files:",
  "  vst_matrix.csv               - VST-normalized expression matrix",
  "  pca_vst.png                  - PCA plot",
  "  sample_distance_heatmap.png  - Sample similarity heatmap",
  "  volcano_AA_vs_HC.png         - Volcano plot: AA vs HC",
  "  volcano_FA_vs_HC.png         - Volcano plot: FA vs HC",
  "  volcano_FA_vs_AA.png         - Volcano plot: FA vs AA",
  "  deseq2_AA_vs_HC.csv          - Full DE results: AA vs HC",
  "  deseq2_FA_vs_HC.csv          - Full DE results: FA vs HC",
  "  deseq2_FA_vs_AA.csv          - Full DE results: FA vs AA",
  "  top_DEG_heatmap.png          - Heatmap of top significant DEGs",
  "  ibmfs_signature_heatmap.png  - IBMFS key gene heatmap",
  "  ibmfs_signature_genes_check.csv - Signature gene presence check",
  if (go_enrichment_attempted) "  GO_BP_*.csv / KEGG_*.csv    - Enrichment results" else
    "  [GO/KEGG: install clusterProfiler + org.Hs.eg.db to enable]"
)

writeLines(summary_lines, file.path(out_dir, "analysis_summary.txt"))
cat(">>> Done! Results saved in:", out_dir, "\n")
