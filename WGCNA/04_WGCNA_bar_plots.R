# Load required libraries
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Source common functions
source("scripts/WGCNA_all/00_WGCNA_common.R")

# Load saved workspace
load("wgcna_output/all/complete_analysis.RData")

# Check data consistency
check_data_consistency()

# Create directory for eigengene plots
dir.create("wgcna_output/all/eigengene_plots", showWarnings = FALSE)

# Get standardized sample ordering
sample_order <- create_sample_order()
sample_info <- sample_order$info

# Get standardized module colors
module_color_map <- create_module_color_map()

# Function to create eigengene plot with fixed labels
create_eigengene_plot <- function(module_name, eigengene_values, sample_info) {
  # Reorder eigengene values to match desired sample order
  eigengene_values <- eigengene_values[match(sample_info$SampleID, rownames(all_datExpr))]
  
  # Create labels in the correct order
  plot_data <- data.frame(
    SampleID = sample_info$SampleID,
    Group = sample_info$Group,
    Sample_Number = sample_info$Sample_Number,
    Eigengene = eigengene_values
  )
  
  # Force the order of the x-axis labels
  plot_data$x_labels <- factor(
    paste(plot_data$Group, plot_data$Sample_Number),
    levels = c(
      paste("CTRL_F", 1:4),
      paste("TP_F", 1:4),
      paste("CTRL_M", 1:4),
      paste("TP_M", 1:3)
    )
  )
  
  # Add position column for spacing
  plot_data$position <- as.numeric(factor(plot_data$Group, 
                                          levels = c("CTRL_F", "TP_F", "CTRL_M", "TP_M")))
  plot_data$x_position <- seq_along(plot_data$x_labels) + 
    (plot_data$position - 1) * 0.5  # Add space between groups
  
  p <- ggplot(plot_data, aes(x = x_position, y = Eigengene, fill = Group)) +
    geom_bar(stat = "identity", width = 0.8) +  # Made bars slightly wider
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c(
      "CTRL_F" = "pink",
      "TP_F" = "pink4",
      "CTRL_M" = "skyblue",
      "TP_M" = "skyblue4"
    )) +
    scale_x_continuous(breaks = plot_data$x_position) +  # Custom x-axis positions
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.title.x = element_blank(), # Remove x-axis title
      plot.title = element_text(size = 20, hjust = 0.5),  # Larger title
      axis.text.y = element_text(size = 14),  # Larger y-axis text
      axis.title.y = element_text(size = 16),  # Larger y-axis title
      legend.position = "none",  # Remove legend
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(
      title = paste0(module_name, " module"),  # Modified title format
      y = "Eigengene Value"
    )
  
  return(p)
}

# Get module colors and create plots
module_colors <- unique(all_moduleColors)
all_plots <- list()

# Create and save individual plots
for(module in module_colors) {
  ME_name <- paste0("ME", module)
  if(ME_name %in% colnames(MEs)) {
    cat("Creating plot for module:", module, "\n")
    eigengene_values <- MEs[, ME_name]
    p <- create_eigengene_plot(module, eigengene_values, sample_info)
    
    # Save individual plot as PNG with high resolution
    filename <- paste0("wgcna_output/all/eigengene_plots/Module_", module, "_eigengene.png")
    ggsave(filename, plot = p, width = 7, height = 3, dpi = 300, device = "png")
    
    all_plots[[module]] <- p
  }
}

# Create summary data frame with reordered eigengene values
eigengene_summary <- data.frame(
  SampleID = sample_info$SampleID,
  Group = sample_info$Group,
  Sample_Number = sample_info$Sample_Number
)

for(module in module_colors) {
  ME_name <- paste0("ME", module)
  if(ME_name %in% colnames(MEs)) {
    eigengene_values <- MEs[, ME_name]
    eigengene_summary[[paste0("Module_", module)]] <- 
      eigengene_values[match(sample_info$SampleID, rownames(all_datExpr))]
  }
}

# Calculate summary statistics
summary_stats <- eigengene_summary %>%
  gather(Module, Value, starts_with("Module_")) %>%
  group_by(Module, Group) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    N = n(),
    SE = SD/sqrt(N),
    .groups = 'drop'
  )

# Save summary statistics
write.csv(summary_stats, 
          "wgcna_output/all/eigengene_plots/module_eigengene_summary_stats.csv", 
          row.names = FALSE)

# Create a combined PNG with all plots
png("wgcna_output/all/eigengene_plots/all_module_eigengenes.png", 
    width = 15, height = 8, units = "in", res = 300)
for(i in seq(1, length(all_plots), 2)) {
  plots_to_combine <- all_plots[i:min(i+1, length(all_plots))]
  do.call(grid.arrange, c(plots_to_combine, ncol = 1))
}
dev.off()

# Create heatmap with reordered data
eigengene_matrix <- as.matrix(eigengene_summary[, startsWith(colnames(eigengene_summary), 
                                                             "Module_")])
rownames(eigengene_matrix) <- paste(eigengene_summary$Group, eigengene_summary$Sample_Number)

png("wgcna_output/all/eigengene_plots/eigengene_heatmap.png", width = 12, height = 8,
    units = "in", res = 300)
heatmap(eigengene_matrix,
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margin = c(10, 10),
        main = "Module Eigengenes Across Samples",
        xlab = "Modules",
        ylab = "Samples")
dev.off()

# Create correlation heatmap between modules
module_correlations <- cor(eigengene_matrix)
png("wgcna_output/all/eigengene_plots/module_correlation_heatmap.png", 
    width = 12, height = 10, units = "in", res = 300)
heatmap(module_correlations,
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margin = c(10, 10),
        main = "Module-Module Correlations",
        xlab = "Modules",
        ylab = "Modules")
dev.off()

# Calculate average eigengene values by group
group_averages <- eigengene_summary %>%
  group_by(Group) %>%
  summarise(across(starts_with("Module_"), mean)) %>%
  column_to_rownames("Group")

# Create group comparison heatmap
png("wgcna_output/all/eigengene_plots/group_comparison_heatmap.png", 
    width = 12, height = 6, units = "in", res = 300)
heatmap(as.matrix(group_averages),
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margin = c(10, 10),
        main = "Average Module Expression by Group",
        xlab = "Modules",
        ylab = "Groups")
dev.off()

# Print summary information
cat("\nEigengene Analysis Summary:\n")
cat("Number of modules analyzed:", length(module_colors), "\n")
cat("Number of samples per group:\n")
print(table(sample_info$Group))

# Save complete analysis results
save(eigengene_summary, summary_stats, module_correlations, group_averages,
     file = "wgcna_output/all/eigengene_plots/eigengene_analysis_results.RData")

# Create summary file
sink("wgcna_output/all/eigengene_plots/analysis_summary.txt")
cat("Eigengene Analysis Summary\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nNumber of modules:", length(module_colors), "\n")
cat("\nSample distribution:\n")
print(table(sample_info$Group))
cat("\nTop module correlations:\n")
print(head(sort(module_correlations[upper.tri(module_correlations)], 
                decreasing = TRUE), 10))
sink()

print("Analysis complete. Results saved in wgcna_output/all/eigengene_plots/")