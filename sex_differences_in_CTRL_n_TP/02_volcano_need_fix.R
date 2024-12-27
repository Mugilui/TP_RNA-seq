# Load required libraries
library(ggplot2)
library(ggrepel)

# Create output directory if it doesn't exist
output_dir <- "output_DEG_analysis"
plots_dir <- file.path(output_dir, "plots")

for (dir in c(output_dir, plots_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
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

# Function to create volcano plot data using UNFILTERED data
prepare_volcano_data <- function(deseq_data, pvalue_cutoff) {
  # Create a copy of the data
  plot_data <- deseq_data
  
  # Handle NA values in padj - set them to 1 (not significant)
  plot_data$padj[is.na(plot_data$padj)] <- 1
  
  # Add Expression column for coloring
  plot_data$Expression <- "Not Significant"
  
  # Label up and down regulated genes
  plot_data$Expression[plot_data$log2FoldChange > 0 & plot_data$padj < pvalue_cutoff] <- "Upregulated"
  plot_data$Expression[plot_data$log2FoldChange < 0 & plot_data$padj < pvalue_cutoff] <- "Downregulated"
  
  # Calculate -log10(p-value), handle NA values in pvalue
  plot_data$pvalue[is.na(plot_data$pvalue)] <- 1
  plot_data$log10_pvalue <- -log10(plot_data$pvalue)
  
  # Remove infinite values
  plot_data <- subset(plot_data, is.finite(log10_pvalue) & is.finite(log2FoldChange))
  
  # Get top 20 genes by p-value (excluding NA p-values)
  top_genes <- head(plot_data[plot_data$padj < pvalue_cutoff, ][
    order(-plot_data[plot_data$padj < pvalue_cutoff, ]$log10_pvalue), ], 20)
  
  return(list(plot_data = plot_data, top_genes = top_genes))
}

# Function to create volcano plot
create_volcano_plot <- function(plot_data, top_genes, title) {
  # Calculate axis limits with some padding
  max_x <- max(abs(plot_data$log2FoldChange), na.rm = TRUE) * 1.1
  max_y <- max(plot_data$log10_pvalue, na.rm = TRUE) * 1.1
  
  volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = log10_pvalue, color = Expression)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Upregulated" = "tomato", 
                                  "Downregulated" = "cornflowerblue", 
                                  "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 P-value") +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    coord_cartesian(xlim = c(-max_x, max_x), ylim = c(0, max_y)) +  # Set axis limits
    geom_text_repel(data = top_genes, 
                    aes(label = GeneID),
                    size = 5,  # Reduced text size for better PDF rendering
                    color = "black",
                    box.padding = unit(1.2, "lines"), # 0.7
                    point.padding = unit(1.0, "lines"), # 0.5
                    segment.size = 0.15,
                    min.segment.length = 0.1,
                    max.overlaps = Inf)
  
  return(volcano_plot)
}

# Main execution
main <- function() {
  # Set adjusted p-value threshold
  pvalue_cutoff <- 0.05
  
  # Load DESeq2 results
  ctrl_DESeq2_results <- load_deseq_results("data/DESeq2_results_CTRLs.tabular")
  tp_DESeq2_results <- load_deseq_results("data/DESeq2_results_TPs.tabular")
  
  # Create volcano plots
  ctrl_volcano_data <- prepare_volcano_data(ctrl_DESeq2_results, pvalue_cutoff)
  ctrl_volcano <- create_volcano_plot(ctrl_volcano_data$plot_data, 
                                      ctrl_volcano_data$top_genes,
                                      "CTRL sex-specific DEGs volcano plot")
  
  tp_volcano_data <- prepare_volcano_data(tp_DESeq2_results, pvalue_cutoff)
  tp_volcano <- create_volcano_plot(tp_volcano_data$plot_data,
                                    tp_volcano_data$top_genes,
                                    "TP sex-specific DEGs volcano plot")
  
  # Save plots as PDFs using standard PDF device
  ggsave(file.path(plots_dir, "volcano_ctrl.pdf"), 
         ctrl_volcano, 
         width = 6, 
         height = 5)
  
  ggsave(file.path(plots_dir, "volcano_tp.pdf"), 
         tp_volcano, 
         width = 6, 
         height = 5)
  
  # Print summary of points in each category
  cat("\nVolcano Plot Categories:\n")
  cat("CTRL genes:\n")
  print(table(ctrl_volcano_data$plot_data$Expression))
  cat("\nTP genes:\n")
  print(table(tp_volcano_data$plot_data$Expression))
}

# Run the main function
main()