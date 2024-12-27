# Load the necessary libraries for volcano plots
library(ggplot2)
library(ggrepel)

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
  
  # Get top 20 genes by p-value (excluding NA p-values)
  top_genes <- head(plot_data[!is.na(plot_data$pvalue) & 
                                plot_data$padj < pvalue_cutoff, ][
                                  order(-plot_data[!is.na(plot_data$pvalue) & 
                                                     plot_data$padj < pvalue_cutoff, ]$log10_pvalue), ], 20)
  
  return(list(plot_data = plot_data, top_genes = top_genes))
}

# Function to create volcano plot
create_volcano_plot <- function(plot_data, top_genes, title) {
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
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    geom_text_repel(data = top_genes, 
                    aes(label = GeneID),
                    size = 5,
                    color = "black",
                    box.padding = unit(0.7, "lines"),
                    point.padding = unit(0.5, "lines"),
                    segment.size = 0.15,
                    min.segment.length = 0.1,
                    max.overlaps = Inf)
  
  return(volcano_plot)
}

# Create volcano plots using the ORIGINAL unfiltered data
female_volcano_data <- prepare_volcano_data(female_DESeq2_results, pvalue_cutoff)
female_volcano <- create_volcano_plot(female_volcano_data$plot_data, 
                                      female_volcano_data$top_genes,
                                      "Female DEGs Volcano Plot")

male_volcano_data <- prepare_volcano_data(male_DESeq2_results, pvalue_cutoff)
male_volcano <- create_volcano_plot(male_volcano_data$plot_data,
                                    male_volcano_data$top_genes,
                                    "Male DEGs Volcano Plot")

# Save plots
ggsave(file.path(plots_dir, "volcano_female.png"), 
       female_volcano, 
       width = 6, 
       height = 5, 
       dpi = 300)

ggsave(file.path(plots_dir, "volcano_male.png"), 
       male_volcano, 
       width = 6, 
       height = 5, 
       dpi = 300)

# Print summary of points in each category
cat("\nVolcano Plot Categories:\n")
cat("Female genes:\n")
print(table(female_volcano_data$plot_data$Expression))
cat("\nMale genes:\n")
print(table(male_volcano_data$plot_data$Expression))