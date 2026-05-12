args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) args[[1]] else "/Users/jaeeunyoo/Desktop/researches/ibmfs_fastq_raw_data"

counts_path <- file.path(root_dir, "counts", "featureCounts.cleaned.txt")
meta_path <- file.path(root_dir, "metadata", "samples.tsv")
out_dir <- file.path(root_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_path)) stop("Missing counts file: ", counts_path)
if (!file.exists(meta_path)) stop("Missing metadata file: ", meta_path)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

counts_tbl <- read.table(counts_path, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
gene_id <- counts_tbl[[1]]
count_mat <- as.matrix(counts_tbl[, 7:ncol(counts_tbl)])
rownames(count_mat) <- gene_id
storage.mode(count_mat) <- "integer"

meta <- read.table(meta_path, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
req_cols <- c("sample", "condition")
if (!all(req_cols %in% colnames(meta))) stop("metadata/samples.tsv must have columns: sample, condition")

colnames(count_mat) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(count_mat))
meta <- meta[match(colnames(count_mat), meta$sample), , drop = FALSE]
if (any(is.na(meta$sample))) stop("Metadata samples do not match count columns")

meta$condition <- factor(meta$condition)
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

vsd <- vst(dds, blind = TRUE)
write.csv(as.data.frame(assay(vsd)), file.path(out_dir, "vst_matrix.csv"))

pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()
ggsave(filename = file.path(out_dir, "pca_vst.png"), plot = p_pca, width = 6, height = 5, dpi = 150)

if (nlevels(meta$condition) < 2) {
  writeLines("Only one condition present; skipping differential expression.", con = file.path(out_dir, "DE_SKIPPED.txt"))
  quit(save = "no", status = 0)
}

dds <- DESeq(dds)

conds <- levels(meta$condition)
if (length(conds) != 2) {
  writeLines("More than 2 conditions present; edit script to choose contrasts.", con = file.path(out_dir, "DE_SKIPPED.txt"))
  quit(save = "no", status = 0)
}

res <- results(dds, contrast = c("condition", conds[[2]], conds[[1]]))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)

res_df$neglog10padj <- -log10(pmax(res_df$padj, 1e-300))
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(alpha = 0.4, size = 1) +
  theme_bw()
ggsave(filename = file.path(out_dir, "volcano.png"), plot = p_volcano, width = 6, height = 5, dpi = 150)

