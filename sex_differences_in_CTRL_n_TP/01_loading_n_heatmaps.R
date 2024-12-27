# Load required libraries
library(pheatmap)

# Create output directory if it doesn't exist
output_dir <- "output_DEG_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
}

# Create subdirectories for different output types
plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")

for (dir in c(plots_dir, tables_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    cat("Created directory:", dir, "\n")
  }
}

# Function to load headerless DESeq2 results
load_deseq_results <- function(file_path) {
  # Define column names for DESeq2 results
  deseq2_colnames <- c("GeneID", "baseMean", "log2FoldChange", 
                       "lfcSE", "stat", "pvalue", "padj")
  
  # Load the data with defined column names
  data <- read.table(file_path, 
                     sep = "\t", 
                     header = FALSE,
                     col.names = deseq2_colnames,
                     stringsAsFactors = FALSE)
  
  # Print summary of loaded data
  cat("\nLoaded", nrow(data), "genes from", file_path, "\n")
  cat("First few rows:\n")
  print(head(data))
  
  return(data)
}

# Load DESeq2 results for both conditions
ctrl_DESeq2_results <- load_deseq_results("data/DESeq2_results_CTRLs.tabular")
tp_DESeq2_results <- load_deseq_results("data/DESeq2_results_TPs.tabular")

# Load additional data files (if needed)
asd_gene_list <- read.csv("data/SFARI-Gene_genes_08-19-2024release_09-19-2024export.csv")

# Remove rows with NA in log2FoldChange
ctrl_DESeq2_results <- ctrl_DESeq2_results[!is.na(ctrl_DESeq2_results$log2FoldChange), ]
tp_DESeq2_results <- tp_DESeq2_results[!is.na(tp_DESeq2_results$log2FoldChange), ]

# Set adjusted p-value threshold
pvalue_cutoff <- 0.05

# Filter DEGs for CTRLs using adjusted p-value
deg_ctrl <- ctrl_DESeq2_results[ctrl_DESeq2_results$padj < pvalue_cutoff, ]
deg_ctrl <- deg_ctrl[!is.na(deg_ctrl$GeneID), ]
deg_GeneID_ctrl <- deg_ctrl$GeneID

# Filter DEGs for TPs using adjusted p-value
deg_tp <- tp_DESeq2_results[tp_DESeq2_results$padj < pvalue_cutoff, ]
deg_tp <- deg_tp[!is.na(deg_tp$GeneID), ]
deg_GeneID_tp <- deg_tp$GeneID

# Print summary statistics
cat("\nDifferentially Expressed Genes Summary:\n")
cat("CTRL Sex Differences DEGs:", length(deg_GeneID_ctrl), "\n")
cat("TP Sex Differences DEGs:", length(deg_GeneID_tp), "\n")

# Write DEG lists to files
write.table(deg_GeneID_ctrl, 
            file = file.path(tables_dir, "deg_GeneID_ctrl.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

write.table(deg_GeneID_tp, 
            file = file.path(tables_dir, "deg_GeneID_tp.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

# Save complete DEG results
write.table(deg_ctrl,
            file = file.path(tables_dir, "complete_DEG_results_ctrl.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(deg_tp,
            file = file.path(tables_dir, "complete_DEG_results_tp.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Function to create heatmap with simplified column names
create_heatmap <- function(deg_genes, normalized_counts_file, condition, output_png) {
  if(length(deg_genes) < 2) {
    cat("Not enough DEGs to create heatmap for", condition, "\n")
    return(NULL)
  }
  
  # Load normalized counts
  normalized_counts <- try(
    read.table(normalized_counts_file, header = TRUE, row.names = 1, sep = "\t")
  )
  
  if(inherits(normalized_counts, "try-error")) {
    cat("Error loading normalized counts for", condition, "\n")
    return(NULL)
  }
  
  # Subset normalized counts to only include DEGs
  deg_counts <- normalized_counts[deg_genes, ]
  
  # Create simplified column names
  old_names <- colnames(deg_counts)
  new_names <- character(length(old_names))
  
  # Counter for females and males in each condition
  f_ctrl <- 1
  f_tp <- 1
  m_ctrl <- 1
  m_tp <- 1
  
  for(i in seq_along(old_names)) {
    if(grepl("CTRL.*F", old_names[i])) {
      new_names[i] <- sprintf("CTRL F%d", f_ctrl)
      f_ctrl <- f_ctrl + 1
    } else if(grepl("TP.*F", old_names[i])) {
      new_names[i] <- sprintf("TP F%d", f_tp)
      f_tp <- f_tp + 1
    } else if(grepl("CTRL.*M", old_names[i])) {
      new_names[i] <- sprintf("CTRL M%d", m_ctrl)
      m_ctrl <- m_ctrl + 1
    } else if(grepl("TP.*M", old_names[i])) {
      new_names[i] <- sprintf("TP M%d", m_tp)
      m_tp <- m_tp + 1
    } else {
      new_names[i] <- old_names[i]
    }
  }
  
  colnames(deg_counts) <- new_names
  
  # Scale the data
  scaled_deg_counts <- t(scale(t(deg_counts)))
  
  # Create heatmap and save to PNG
  png(output_png, width = 1500, height = 2400, res = 300)
  pheatmap(scaled_deg_counts,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize_col = 14,
           color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  
  cat("Heatmap saved to:", output_png, "\n")
}

# Create heatmaps for both conditions
cat("\nCreating heatmaps...\n")

# CTRL heatmap
create_heatmap(deg_GeneID_ctrl, 
               "data/normalized_counts_DESeq2_CTRLs.tabular",
               "CTRL",
               file.path(plots_dir, "heatmap_ctrl_DEGs.png"))

# TP heatmap
create_heatmap(deg_GeneID_tp,
               "data/normalized_counts_DESeq2_TPs.tabular",
               "TP",
               file.path(plots_dir, "heatmap_tp_DEGs.png"))

# Save analysis summary
summary_file <- file.path(output_dir, "analysis_summary.txt")
sink(summary_file)
cat("DEG Analysis Summary\n")
cat("===================\n\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("P-value cutoff:", pvalue_cutoff, "\n\n")
cat("Results:\n")
cat("- CTRL Sex Differences DEGs:", length(deg_GeneID_ctrl), "\n")
cat("- TP Sex Differences DEGs:", length(deg_GeneID_tp), "\n")
sink()

cat("\nAnalysis complete! Results saved in:", output_dir, "\n")
