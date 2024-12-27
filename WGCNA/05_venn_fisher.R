# Disable scientific notation
options(scipen = 999)

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(WGCNA)
  library(homologene)
  library(VennDiagram)
  library(grid)
  library(gridExtra)
})

# Create output directory structure
output_dir <- "wgcna_output/all/SFARI_analysis"
venn_dir <- file.path(output_dir, "venn_diagrams")
dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

# Create a log file
log_file <- file.path(output_dir, "analysis_log.txt")
sink(log_file)
cat("WGCNA-SFARI Analysis Log\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Load WGCNA results
load("wgcna_output/all/complete_analysis.RData")

# Get module assignments
module_df <- data.frame(
  gene = colnames(all_datExpr),
  module = all_moduleColors
)

# Convert mouse genes to human genes using homologene
human_homologs <- homologene(module_df$gene, inTax = 10090, outTax = 9606)
names(human_homologs) <- c("mouse_gene", "human_gene", "mouse_id", "human_id")

# Clean up homolog mapping
gene_map <- human_homologs %>%
  dplyr::select(mouse_gene, human_gene) %>%
  dplyr::distinct()

# Get background gene set (only successfully mapped genes)
background_genes <- unique(gene_map$human_gene[!is.na(gene_map$human_gene)])
n_background <- length(background_genes)

# 3. Read SFARI genes
sfari_genes <- read.csv("data/SFARI-Gene_genes_08-19-2024release_09-19-2024export.csv")
sfari_symbols <- unique(sfari_genes$gene.symbol)

# Print detailed mapping statistics
cat("Gene mapping statistics:\n")
cat("Starting mouse genes from WGCNA:", length(unique(module_df$gene)), "\n")
cat("Successfully mapped to human:", n_background, "\n")
cat("Total SFARI genes:", length(sfari_symbols), "\n")
cat("SFARI genes found in our background set:", length(intersect(background_genes, sfari_symbols)), "\n\n")

cat("This means our Fisher's test will use:\n")
cat("- Total background genes (universe):", n_background, "\n")
cat("- These are only the genes that were successfully mapped from mouse to human\n\n")

# Function to perform Fisher's exact test for enrichment
perform_enrichment_test <- function(module_genes, sfari_genes, universe_genes) {
  # Create contingency table
  in_module <- universe_genes %in% module_genes
  is_sfari <- universe_genes %in% sfari_genes
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(table(in_module, is_sfari))
  
  # Get overlapping genes
  overlapping_genes <- intersect(module_genes, sfari_genes)
  
  return(list(
    p_value = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    n_overlap = length(overlapping_genes),
    overlapping_genes = overlapping_genes
  ))
}

# Find overlaps and create Venn diagrams for each module
results <- list()

# Set up colors for Venn diagrams
module_colors <- setNames(
  colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))(length(unique(module_df$module))),
  unique(module_df$module)
)

for (module in unique(module_df$module)) {
  cat("\nProcessing module:", module, "\n")
  
  # Get module genes
  module_mouse_genes <- module_df$gene[module_df$module == module]
  
  # Convert to human genes
  module_human_genes <- gene_map$human_gene[gene_map$mouse_gene %in% module_mouse_genes]
  module_human_genes <- module_human_genes[!is.na(module_human_genes)]
  
  # Calculate statistics
  stats <- perform_enrichment_test(module_human_genes, sfari_symbols, background_genes)
  
  # Print detailed information for this module
  cat("Module size (mapped to human):", length(module_human_genes), "\n")
  cat("SFARI genes in module:", stats$n_overlap, "\n")
  cat("P-value:", format(stats$p_value, scientific = FALSE, digits = 4), "\n")
  cat("Odds ratio:", round(stats$odds_ratio, 2), "\n\n")
  
  # Store results
  results[[module]] <- data.frame(
    module = module,
    total_genes = length(module_human_genes),
    n_asd_genes = stats$n_overlap,
    percent_overlap = round(stats$n_overlap / length(module_human_genes) * 100, 2),
    p_value = stats$p_value,
    odds_ratio = stats$odds_ratio,
    asd_genes = paste(stats$overlapping_genes, collapse = ", ")
  )
  
  # Create Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(module_human_genes),
    area2 = length(sfari_symbols),
    cross.area = stats$n_overlap,
    category = c(paste("Module", module), "SFARI Genes"),
    fill = c(module_colors[module], "#999999"),
    alpha = 0.5,
    scaled = TRUE,
    cat.pos = c(0, 0),
    cat.dist = c(0.025, 0.025),
    euler.d = TRUE,
    sep.dist = 0.03
  )
  
  # Save Venn diagram
  pdf(file = file.path(venn_dir, paste0("module_", module, "_venn.pdf")),
      width = 8, height = 8)
  grid.draw(venn.plot)
  dev.off()
  
  # Save gene lists for each module
  module_gene_file <- file.path(output_dir, paste0("module_", module, "_genes.txt"))
  writeLines(c(
    "Module genes (human symbols):",
    module_human_genes,
    "\nOverlapping SFARI genes:",
    stats$overlapping_genes
  ), module_gene_file)
}

# Combine all results
results_df <- do.call(rbind, results)

# Sort by p-value
results_df <- results_df[order(results_df$p_value), ]

# Save results
write.csv(results_df, 
          file.path(output_dir, "module_asd_overlaps.csv"), 
          row.names = FALSE)

# Create a summary table of all overlaps
overlap_summary <- results_df %>%
  dplyr::select(module, total_genes, n_asd_genes, percent_overlap, p_value, odds_ratio) %>%
  dplyr::mutate(
    significant = ifelse(p_value < 0.05, "Yes", "No"),
    p_value = format(p_value, scientific = FALSE, digits = 4),
    odds_ratio = round(odds_ratio, 2)
  ) %>%
  dplyr::arrange(p_value)

# Save summary table as PDF
pdf(file.path(output_dir, "overlap_summary_table.pdf"), 
    width = 12, height = length(unique(module_df$module))/2)
grid.table(overlap_summary)
dev.off()

# Close the log file
sink()

cat("\nAnalysis complete! All results saved in:", output_dir, "\n")