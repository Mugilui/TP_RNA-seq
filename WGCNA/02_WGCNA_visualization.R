# Load required libraries
library(igraph)
library(WGCNA)

# Source common functions
source("scripts/WGCNA_all/00_WGCNA_common.R")

# Load saved workspace
load("wgcna_output/all/complete_analysis.RData")

# Check data consistency
check_data_consistency()

# Get standardized module colors
module_color_map <- create_module_color_map()
unique_modules <- names(module_color_map)

# Create directory for saving the plots
output_dir <- "wgcna_output/all/modules_visualizations"
dir.create(output_dir, recursive = TRUE)

# Load the TOM (assuming it was saved as a single file)
TOM_file <- "all_TOM-block.1.RData"

if (file.exists(TOM_file)) {
  load(TOM_file)
  module_TOM <- as.matrix(TOM)
  print(paste("Dimensions of module_TOM:", dim(module_TOM)[1], "x", dim(module_TOM)[2]))
  
  # Loop over each module
  for (module in unique_modules) {
    print(paste("Processing module:", module))
    
    # Get genes for this module
    module_genes <- which(all_moduleColors == module)
    print(paste("Number of genes in module:", length(module_genes)))
    
    # Ensure module_genes indices are within the TOM dimensions
    valid_genes <- module_genes[module_genes <= nrow(module_TOM)]
    
    if (length(valid_genes) > 1) {
      # Subset the TOM for the module genes
      module_TOM_subset <- module_TOM[valid_genes, valid_genes]
      print(paste("Dimensions of module_TOM_subset:", dim(module_TOM_subset)[1], "x", dim(module_TOM_subset)[2]))
      
      # Calculate the connectivity (degree) within the module based on TOM
      gene_connectivity <- rowSums(module_TOM_subset)
      
      # Ensure that we don't attempt to select more genes than available
      n_genes_to_select <- min(10, length(gene_connectivity))
      
      # Rank genes by connectivity and select the top N hub genes (at most 10)
      top_hub_genes <- order(gene_connectivity, decreasing = TRUE)[1:n_genes_to_select]
      
      # Extract gene names for the top hub genes
      gene_names <- colnames(all_datExpr)[valid_genes[top_hub_genes]]
      
      # Extract the TOM for the top hub genes
      top_genes_TOM <- module_TOM_subset[top_hub_genes, top_hub_genes]
      diag(top_genes_TOM) <- 0
      
      # Create the network graph
      network <- graph_from_adjacency_matrix(top_genes_TOM, 
                                             mode = "undirected", 
                                             weighted = TRUE, 
                                             diag = FALSE)
      network <- simplify(network)
      
      # Set vertex color to match module color
      vertex_colors <- rep(module, length(gene_names))
      
      # Create the file name for saving the plot
      plot_filename <- paste0(output_dir, "/", module, "_module_network.png")
      
      # Save the plot as a PNG with high resolution and transparency
      png(plot_filename, width = 2100, height = 2100, res = 300, bg = "transparent")
      
      # Set transparent background for the plot
      par(bg = NA)
      
      # Plot the network with larger text
      plot(network, 
           vertex.size = 10,           # Slightly larger vertices
           vertex.label = gene_names,
           vertex.label.cex = 1.7,     # Increased text size for gene names
           vertex.label.font = 2,
           vertex.label.color = "black",
           vertex.color = vertex_colors,
           vertex.frame.color = "black",
           vertex.label.dist = -1.5,   # Adjusted label distance for better visibility
           vertex.label.degree = pi/2,
           vertex.label.family = "sans",
           main = "")
      
      # Close the PNG device
      dev.off()
      
    } else {
      print(paste("Warning: Not enough valid genes for module", module))
      print("Skipping network creation for this module.")
    }
  }
} else {
  warning("TOM file not found!")
}

# Create a more detailed traits data frame
all_traits <- data.frame(
  SampleID = rownames(all_datExpr),
  Group = factor(c(rep("CTRL_F", 4), rep("CTRL_M", 4), rep("TP_F", 4), rep("TP_M", 3)),
                 levels = c("CTRL_F", "CTRL_M", "TP_F", "TP_M")),
  Treatment = factor(c(rep("CTRL", 8), rep("TP", 7))),
  Sex = factor(c(rep("F", 4), rep("M", 4), rep("F", 4), rep("M", 3)))
)

# Create binary variables for all comparisons
all_traits$Treatment_Numeric <- ifelse(all_traits$Treatment == "CTRL", 0, 1)
all_traits$Sex_Numeric <- ifelse(all_traits$Sex == "F", 0, 1)
all_traits$CTRL_F_vs_Others <- ifelse(all_traits$Group == "CTRL_F", 1, 0)
all_traits$CTRL_M_vs_Others <- ifelse(all_traits$Group == "CTRL_M", 1, 0)
all_traits$TP_F_vs_Others <- ifelse(all_traits$Group == "TP_F", 1, 0)
all_traits$TP_M_vs_Others <- ifelse(all_traits$Group == "TP_M", 1, 0)

# Correlate all module eigengenes with traits
trait_columns <- c("Treatment_Numeric", "Sex_Numeric", 
                   "CTRL_F_vs_Others", "CTRL_M_vs_Others", 
                   "TP_F_vs_Others", "TP_M_vs_Others")

all_moduleTraitCor <- cor(MEs, all_traits[, trait_columns], 
                          use = "pairwise.complete.obs")

# Calculate p-values
all_nSamples <- nrow(all_datExpr)
all_moduleTraitPvalue <- corPvalueStudent(all_moduleTraitCor, all_nSamples)

# Create text matrix for the heatmap
all_textMatrix <- paste(signif(all_moduleTraitCor, 2), "\n(",
                        signif(all_moduleTraitPvalue, 1), ")", sep = "")

# Create better labels for the heatmap
trait_labels <- c("Treatment (CTRL vs TP)", 
                  "Sex (F vs M)",
                  "CTRL Female vs Others",
                  "CTRL Male vs Others",
                  "TP Female vs Others",
                  "TP Male vs Others")

# Create the heatmap visualization as PNG
png(file.path(output_dir, "Module_Trait_Relationships.png"), 
    width = 3000, height = 3600, res = 300)
labeledHeatmap(Matrix = all_moduleTraitCor,
               xLabels = trait_labels,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = all_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,  # Increased text size for correlation values
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

# Save results
save(all_moduleTraitCor, all_moduleTraitPvalue, all_traits, 
     file = "wgcna_output/all/trait_correlations.RData")

print("Network visualization and trait correlation analysis complete.")