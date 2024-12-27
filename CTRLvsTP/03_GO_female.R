# Load required libraries
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(stringr)

# Create output directory if it doesn't exist
dir.create("GO_DEGs/female", recursive = TRUE, showWarnings = FALSE)

# Function to separate up and down-regulated genes
separate_up_down_genes <- function(deg_data, log2FC_threshold = 0) {
  # Modified to use log2FoldChange column name from DESeq2 results
  up_genes <- deg_data$GeneID[deg_data$log2FoldChange > log2FC_threshold]
  down_genes <- deg_data$GeneID[deg_data$log2FoldChange < -log2FC_threshold]
  return(list(all = deg_data$GeneID, up = up_genes, down = down_genes))
}

# Separate up and down-regulated genes using the filtered DEG results
female_genes <- separate_up_down_genes(deg_female)

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

# Perform GO enrichment for all, up and downregulated genes
female_all_ego <- perform_go_enrichment(female_genes$all)
female_up_ego <- perform_go_enrichment(female_genes$up)
female_down_ego <- perform_go_enrichment(female_genes$down)

custom_dotplot <- function(ego_up, ego_down, title, base_size = 15, dot_size_range = c(4, 12), width = 40) {
  # Check if both enrichment results exist and have results
  if(is.null(ego_up) || is.null(ego_down) || 
     nrow(ego_up@result) == 0 || nrow(ego_down@result) == 0) {
    warning("Insufficient enrichment results for plotting")
    return(NULL)
  }
  
  # Combine results
  ego_up@result$Regulation <- "UP"
  ego_down@result$Regulation <- "DOWN"
  combined_result <- rbind(ego_up@result, ego_down@result)
  
  # Select top 5 terms for each regulation based on adjusted p-value
  top_terms <- combined_result %>%
    group_by(Regulation) %>%
    top_n(5, wt = -p.adjust) %>%
    ungroup()
  
  # Wrap the text in the Description column
  top_terms$Description <- str_wrap(top_terms$Description, width = width)
  
  # Create the plot
  ggplot(top_terms, aes(x = Regulation, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    scale_size(range = dot_size_range) +
    scale_color_gradient(low = "red", high = "blue") +
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

# Create and save dotplot
if(!is.null(female_up_ego) && !is.null(female_down_ego)) {
  female_dotplot <- custom_dotplot(female_up_ego, female_down_ego, 
                                   "Female DEGs GO enrichment", 
                                   base_size = 20, 
                                   dot_size_range = c(10, 30),
                                   width = 30)
  
  if(!is.null(female_dotplot)) {
    ggsave("GO_DEGs/female/female_custom_dotplot.png", 
           female_dotplot, 
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

# Save GO results
save_go_results(female_all_ego, "GO_DEGs/female/female_all_DEGs_GO.csv")
save_go_results(female_up_ego, "GO_DEGs/female/female_up_regulated_GO.csv")
save_go_results(female_down_ego, "GO_DEGs/female/female_down_regulated_GO.csv")

# Print summary of results
cat("\nGO enrichment analysis for female DEGs completed.\n")
cat("Number of all DEGs:", length(female_genes$all), "\n")
cat("Number of up-regulated genes:", length(female_genes$up), "\n")
cat("Number of down-regulated genes:", length(female_genes$down), "\n")
cat("Results saved in GO_DEGs/female directory.\n")
