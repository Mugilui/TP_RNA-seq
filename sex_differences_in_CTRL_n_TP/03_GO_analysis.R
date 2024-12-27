# Load required libraries
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(stringr)

# Create output directories
dir.create("GO_DEGs/CTRL", recursive = TRUE, showWarnings = FALSE)
dir.create("GO_DEGs/TP", recursive = TRUE, showWarnings = FALSE)

# Function to separate up and downregulated genes
separate_up_down_genes <- function(deg_data, log2FC_threshold = 0) {
  up_genes <- deg_data$GeneID[deg_data$log2FoldChange > log2FC_threshold]
  down_genes <- deg_data$GeneID[deg_data$log2FoldChange < -log2FC_threshold]
  return(list(all = deg_data$GeneID, up = up_genes, down = down_genes))
}

# Separate up and down-regulated genes
ctrl_genes <- separate_up_down_genes(deg_ctrl)
tp_genes <- separate_up_down_genes(deg_tp)

# Function to perform GO enrichment
perform_go_enrichment <- function(gene_list, ont = "BP") {
  if(length(gene_list) == 0) {
    warning("No genes provided for enrichment analysis")
    return(NULL)
  }
  
  ego <- try(enrichGO(gene = gene_list,
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = ont,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05))
  
  if(inherits(ego, "try-error")) {
    warning("Error in GO enrichment analysis")
    return(NULL)
  }
  
  return(ego)
}

# Perform GO enrichment for both conditions
ctrl_all_ego <- perform_go_enrichment(ctrl_genes$all)
ctrl_up_ego <- perform_go_enrichment(ctrl_genes$up)
ctrl_down_ego <- perform_go_enrichment(ctrl_genes$down)

tp_all_ego <- perform_go_enrichment(tp_genes$all)
tp_up_ego <- perform_go_enrichment(tp_genes$up)
tp_down_ego <- perform_go_enrichment(tp_genes$down)

# Dotplot function with significance filtering
custom_dotplot <- function(ego_up, ego_down, title, base_size = 15, dot_size_range = c(4, 12), width = 40, p_cutoff = 0.05) {
  if(is.null(ego_up) || is.null(ego_down) || 
     nrow(ego_up@result) == 0 || nrow(ego_down@result) == 0) {
    warning("Insufficient enrichment results for plotting")
    return(NULL)
  }
  
  # Add regulation information
  ego_up@result$Regulation <- "UP"
  ego_down@result$Regulation <- "DOWN"
  
  # Combine and get significant terms
  combined_result <- rbind(ego_up@result, ego_down@result) %>%
    # Filter for significant terms
    filter(p.adjust < p_cutoff) %>%
    # Remove any redundant terms
    distinct(Description, Regulation, .keep_all = TRUE)
  
  # Print number of significant terms for each regulation
  cat("\nNumber of significant GO terms (p.adjust <", p_cutoff, "):\n")
  print(table(combined_result$Regulation))
  
  # Continue with visualization for top 5 terms
  combined_result <- combined_result %>%
    # Group by regulation direction
    group_by(Regulation) %>%
    # Sort by adjusted p-value within each group
    arrange(p.adjust, .by_group = TRUE) %>%
    # Take top 5 most significant terms for each regulation
    slice_head(n = 10) %>%
    ungroup()
  
  # Check fot significant terms
  if(nrow(combined_result) == 0) {
    warning("No significant GO terms found with p.adjust < ", p_cutoff)
    return(NULL)
  }
  
  # Wrap the GO term descriptions
  combined_result$Description <- str_wrap(combined_result$Description, width = width)
  
  # Create factor levels based on p.adjust to order the terms
  combined_result$Description <- factor(combined_result$Description,
                                        levels = rev(unique(combined_result$Description)))
  
  # Create the plot
  ggplot(combined_result, 
         aes(x = Regulation, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    scale_size(range = dot_size_range, name = "Gene Count") +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted\np-value") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = base_size),
      axis.text.y = element_text(size = base_size + 3, lineheight = 0.8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 4, face = "bold")
    ) +
    labs(title = title,
         y = "GO Terms",
         size = "Gene Count",
         color = "Adjusted p-value")
}

# Create and save dotplots for both conditions
if(!is.null(ctrl_up_ego) && !is.null(ctrl_down_ego)) {
  ctrl_dotplot <- custom_dotplot(ctrl_up_ego, ctrl_down_ego, 
                                 "GO enrichment of sex-specific DEGs in CTRL", 
                                 base_size = 20, 
                                 dot_size_range = c(10, 30),
                                 width = 30,
                                 p_cutoff = 0.05)
  
  if(!is.null(ctrl_dotplot)) {
    ggsave("GO_DEGs/CTRL/ctrl_custom_dotplot.png", 
           ctrl_dotplot, 
           width = 13,
           height = 13,
           dpi = 300)
  }
}

if(!is.null(tp_up_ego) && !is.null(tp_down_ego)) {
  tp_dotplot <- custom_dotplot(tp_up_ego, tp_down_ego, 
                               "GO enrichment of sex-specific DEGs in TP", 
                               base_size = 20, 
                               dot_size_range = c(10, 30),
                               width = 30,
                               p_cutoff = 0.05)
  
  if(!is.null(tp_dotplot)) {
    ggsave("GO_DEGs/TP/tp_custom_dotplot.png", 
           tp_dotplot, 
           width = 13,
           height = 13,
           dpi = 300)
  }
}

# Function to save GO results as CSV
save_go_results <- function(ego, filename) {
  if(!is.null(ego) && nrow(ego@result) > 0) {
    write.csv(ego@result, file = filename, row.names = FALSE)
  } else {
    warning(paste("No results to save for", filename))
  }
}

# Save GO results for CTRL condition
save_go_results(ctrl_all_ego, "GO_DEGs/CTRL/ctrl_all_DEGs_GO.csv")
save_go_results(ctrl_up_ego, "GO_DEGs/CTRL/ctrl_up_regulated_GO.csv")
save_go_results(ctrl_down_ego, "GO_DEGs/CTRL/ctrl_down_regulated_GO.csv")

# Save GO results for TP condition
save_go_results(tp_all_ego, "GO_DEGs/TP/tp_all_DEGs_GO.csv")
save_go_results(tp_up_ego, "GO_DEGs/TP/tp_up_regulated_GO.csv")
save_go_results(tp_down_ego, "GO_DEGs/TP/tp_down_regulated_GO.csv")

# Print summary of results
cat("\nGO enrichment analysis completed.\n")
cat("CTRL condition:\n")
cat("Number of all DEGs:", length(ctrl_genes$all), "\n")
cat("Number of up-regulated genes:", length(ctrl_genes$up), "\n")
cat("Number of down-regulated genes:", length(ctrl_genes$down), "\n")

cat("\nTP condition:\n")
cat("Number of all DEGs:", length(tp_genes$all), "\n")
cat("Number of up-regulated genes:", length(tp_genes$up), "\n")
cat("Number of down-regulated genes:", length(tp_genes$down), "\n")
cat("\nResults saved in GO_DEGs/CTRL and GO_DEGs/TP directories.\n")