# R script for RNA-seq Differential Expression Analysis (AA vs HD)

# Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 1. Load Count Matrix to see structure and clean colnames
countData <- read.table("counts_mini.txt", header = TRUE, row.names = 1, skip = 1)

# Clean colnames: remove path and extension
# Expected format: "AA.CSB_sorted.bam", "HD_10_sorted.bam" etc. (R replaces '-' with '.')
clean_names <- colnames(countData)
clean_names <- gsub(".*AA", "AA", clean_names) # Fix path prefix if any
clean_names <- gsub(".*HD", "HD", clean_names) # Fix path prefix if any
clean_names <- gsub("_sorted\\.bam", "", clean_names)
clean_names <- gsub("\\.sorted\\.bam", "", clean_names)
colnames(countData) <- clean_names

# Filter: We only want column 6 onwards (featureCounts has 6 info columns: Chr, Start, End, Strand, Length)
# Wait, featureCounts output with metadata?
# Yes, first 6 cols are metadata. So count columns start from 7.
countData_counts <- countData[, 6:ncol(countData)]
# Rename again just to be safe
colnames(countData_counts) <- clean_names[6:length(clean_names)]

# 2. Create Metadata
sampleNames <- colnames(countData_counts)
condition <- ifelse(grepl("^AA", sampleNames), "Patient",
  ifelse(grepl("^HD", sampleNames), "Control", "Other")
)

# Check if any "Other" (should not be if we only input AA and HD)
colData <- data.frame(row.names = sampleNames, Condition = factor(condition))

# Relevel to set Control as reference
colData$Condition <- relevel(colData$Condition, ref = "Control")

# 3. Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = countData_counts,
  colData = colData,
  design = ~Condition
)

# Filter low counts (row sum < 10)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 4. Run DESeq
dds <- DESeq(dds)

# 5. Results
res <- results(dds, contrast = c("Condition", "Patient", "Control"))
summary(res)

# Save Results
write.csv(as.data.frame(res), file = "results_AA_vs_HD.csv")

# Filter significant results (p-adj < 0.05)
resSig <- subset(res, padj < 0.05)
write.csv(as.data.frame(resSig), file = "results_AA_vs_HD_significant.csv")

# 6. Plots (Enhanced)

# --- PCA Plot ---
# Use Variance Stabilizing Transformation (VST)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Custom Colors
group_colors <- c("Patient" = "#E41A1C", "Control" = "#377EB8")

png("pca_plot.png", width = 1200, height = 1000, res = 150)
ggplot(pcaData, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 5, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  coord_fixed() +
  ggtitle("Principal Component Analysis")
dev.off()

# --- Volcano Plot (ggplot2) ---
res_df <- as.data.frame(res)
# Add significance column
res_df$Significance <- "NS"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
res_df$Significance <- factor(res_df$Significance, levels = c("Up", "Down", "NS"))

png("volcano_plot.png", width = 1200, height = 1000, res = 150)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot (Patient vs Control)",
    x = "Log2 Fold Change",
    y = "-Log10 P-value"
  ) +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
dev.off()

# --- Heatmap (Enhanced) ---
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), min(50, nrow(assay(vsd))))
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(dds)[, c("Condition")])
colnames(df) <- "Condition"
rownames(df) <- colnames(mat)

# Annotation Key Coloring
ann_colors <- list(Condition = group_colors)

png("heatmap_top50.png", width = 1200, height = 1200, res = 150)
pheatmap(mat,
  annotation_col = df,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  show_rownames = FALSE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 10,
  main = "Top 50 Variable Genes",
  treeheight_row = 30,
  treeheight_col = 30
)
dev.off()
