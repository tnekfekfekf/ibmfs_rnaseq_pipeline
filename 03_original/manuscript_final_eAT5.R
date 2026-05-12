#!/usr/bin/env Rscript

# Fix gene types and re-run DESeq2 analysis
# Extract proper gene types from Gencode GTF

library(DESeq2)
library(ggplot2)
library(dplyr)

# Set working directory
setwd("/Users/jaeeunyoo/Desktop/researches/rna_analysis_desktop")

# Create output directories
dir.create("processed_data/deg_latest_grouping_fixed", showWarnings = FALSE, recursive = TRUE)
dir.create("processed_data/deg_latest_grouping_fixed/lncrna", showWarnings = FALSE, recursive = TRUE)
dir.create("processed_data/deg_latest_grouping_fixed/mrna", showWarnings = FALSE, recursive = TRUE)

# Function to extract gene types from GTF
extract_gene_types <- function(gtf_file) {
  cat("📊 Extracting gene types from GTF file...\n")
  
  # Read GTF file
  gtf_lines <- readLines(gtf_file)
  gene_lines <- gtf_lines[grepl('\tgene\t', gtf_lines)]
  
  gene_info <- data.frame(
    gene_id = character(),
    gene_type = character(),
    stringsAsFactors = FALSE
  )
  
  for (line in gene_lines) {
    # Extract gene_id
    gene_id_match <- regexpr('gene_id "([^"]*)"', line)
    if (gene_id_match > 0) {
      gene_id <- regmatches(line, gene_id_match)[1]
      gene_id <- gsub('gene_id "([^"]*)"', '\\1', gene_id)
      
      # Extract gene_type
      gene_type_match <- regexpr('gene_type "([^"]*)"', line)
      if (gene_type_match > 0) {
        gene_type <- regmatches(line, gene_type_match)[1]
        gene_type <- gsub('gene_type "([^"]*)"', '\\1', gene_type)
        
        # Remove version number from gene_id
        gene_id_clean <- gsub("\\.[0-9]*$", "", gene_id)
        
        gene_info <- rbind(gene_info, data.frame(
          gene_id = gene_id_clean,
          gene_type = gene_type,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Remove duplicates
  gene_info <- gene_info[!duplicated(gene_info$gene_id), ]
  
  cat("Extracted gene types for", nrow(gene_info), "genes\n")
  cat("Gene type distribution:\n")
  print(table(gene_info$gene_type))
  
  return(gene_info)
}

# Extract gene types
gene_types <- extract_gene_types("reference_data/annotations/gencode.v44.annotation.no_rRNA.gtf")

# Define sample groupings
cat("📊 Setting up sample groupings...\n")

g_aa_samples <- c("AA-RNA-FA", "AA-RNA-DKC", "AA-RNA-FA2", "AA-RNA-FA3")  # 4 G-AA
control_samples <- c("AA-PRO", "AA-KEW")  # 2 Controls
u_aa_samples <- c("AA-RNA-1", "AA-RNA-4", "AA-RNA-5", "AA-RNA-13", "AA-RNA-16", "AA-RNA-18", "AA-HMH", "AA-PJH")  # 8 U-AA

all_samples <- c(g_aa_samples, control_samples, u_aa_samples)

# Create metadata
metadata <- data.frame(
  sample_id = all_samples,
  group = c(rep("G-AA", length(g_aa_samples)),
            rep("Control", length(control_samples)),
            rep("U-AA", length(u_aa_samples))),
  stringsAsFactors = FALSE
)

# Read count data
cat("📊 Reading count data...\n")

# Read the count file properly
header_row <- readLines("processed_data/deg_analysis_final_results/expression_tables/all_samples_expression_table_counts_final_correct.txt", n = 1)
col_names <- strsplit(header_row, "\t")[[1]]

count_data <- read.table("processed_data/deg_analysis_final_results/expression_tables/all_samples_expression_table_counts_final_correct.txt", 
                        header = FALSE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, skip = 2)
colnames(count_data) <- col_names[-1]

# Filter to only include our samples
available_samples <- intersect(all_samples, colnames(count_data))
cat("Available samples:", paste(available_samples, collapse = ", "), "\n")

# Filter metadata to only include available samples
metadata <- metadata[metadata$sample_id %in% available_samples, ]

# Extract count matrix
count_matrix <- count_data[, available_samples, drop = FALSE]
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

# Remove genes with zero counts across all samples
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

cat("Count matrix dimensions:", dim(count_matrix), "\n")

# Create DESeq2 object
cat("📊 Creating DESeq2 object...\n")

rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(count_matrix), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ group
)

# Run DESeq2
cat("📊 Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Function to perform DEG analysis with proper gene type filtering
perform_deg_analysis <- function(group1, group2, gene_type = "all") {
  cat("📊 Analyzing", group1, "vs", group2, "(", gene_type, ")\n")
  
  results <- results(dds, contrast = c("group", group1, group2))
  
  # Convert to data frame
  results_df <- as.data.frame(results)
  results_df$gene_id <- rownames(results_df)
  
  # Remove version numbers from gene IDs for matching
  results_df$gene_id_clean <- gsub("\\.[0-9]*$", "", results_df$gene_id)
  
  # Add gene type information
  results_df <- merge(results_df, gene_types, by.x = "gene_id_clean", by.y = "gene_id", all.x = TRUE)
  
  # Add significance categories
  results_df$significance <- "Not significant"
  results_df$significance[results_df$padj < 0.05 & results_df$log2FoldChange > 1] <- "Upregulated"
  results_df$significance[results_df$padj < 0.05 & results_df$log2FoldChange < -1] <- "Downregulated"
  
  # Filter out NA values
  results_df <- results_df[!is.na(results_df$padj), ]
  
  # Filter by gene type if specified
  if (gene_type == "mrna") {
    results_df <- results_df[results_df$gene_type == "protein_coding", ]
  } else if (gene_type == "lncrna") {
    results_df <- results_df[results_df$gene_type == "lncRNA", ]
  }
  
  # Clean up columns
  results_df <- results_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", 
                               "stat", "pvalue", "padj", "significance", "gene_type")]
  
  return(results_df)
}

# Perform all comparisons
comparisons <- list(
  c("G-AA", "Control"),
  c("U-AA", "Control"),
  c("G-AA", "U-AA")
)

gene_types_list <- c("all", "mrna", "lncrna")

# Run analyses
for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  comparison_name <- paste0(group1, "_vs_", group2)
  
  for (gene_type in gene_types_list) {
    cat("📊 Running analysis:", comparison_name, "(", gene_type, ")\n")
    
    results_df <- perform_deg_analysis(group1, group2, gene_type)
    
    # Save results
    if (gene_type == "all") {
      output_file <- paste0("processed_data/deg_latest_grouping_fixed/", comparison_name, "_results.txt")
    } else {
      output_file <- paste0("processed_data/deg_latest_grouping_fixed/", gene_type, "/", comparison_name, "_", gene_type, "_results.txt")
    }
    
    write.table(results_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("✅ Saved:", output_file, "\n")
    cat("   Total genes:", nrow(results_df), "\n")
    cat("   Significant genes:", sum(results_df$significance != "Not significant"), "\n")
    cat("   Upregulated:", sum(results_df$significance == "Upregulated"), "\n")
    cat("   Downregulated:", sum(results_df$significance == "Downregulated"), "\n")
    cat("\n")
  }
}

cat("🎉 DESeq2 analysis with proper gene types completed successfully!\n")
cat("📁 Results saved in: processed_data/deg_latest_grouping_fixed/\n")
