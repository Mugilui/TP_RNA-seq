# Load required libraries
library(homologene)
library(VennDiagram)
library(grid)

# Create output directories
asd_dir <- "ASD_analysis"
dir.create(asd_dir, showWarnings = FALSE)
female_dir <- file.path(asd_dir, "female_ASD_analysis")
male_dir <- file.path(asd_dir, "male_ASD_analysis")
dir.create(female_dir, showWarnings = FALSE)
dir.create(male_dir, showWarnings = FALSE)

# Function to perform Fisher's exact test for overlap significance
perform_overlap_test <- function(deg_genes, asd_genes, background_size) {
  # Input validation
  if(background_size < length(deg_genes) || background_size < length(asd_genes)) {
    stop("Background size must be greater than both gene sets")
  }
  
  # Create contingency table
  overlap <- length(intersect(deg_genes, asd_genes))
  deg_only <- length(deg_genes) - overlap
  asd_only <- length(asd_genes) - overlap
  neither <- background_size - length(deg_genes) - length(asd_genes) + overlap
  
  # Create contingency matrix
  contingency_table <- matrix(c(overlap, deg_only, asd_only, neither), 
                              nrow = 2,
                              dimnames = list(c("In_ASD", "Not_in_ASD"),
                                              c("In_DEG", "Not_in_DEG")))
  
  # Perform Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  
  return(list(
    overlap = overlap,
    odds_ratio = fisher_test$estimate,
    p_value = fisher_test$p.value,
    confidence_interval = fisher_test$conf.int,
    background_size = background_size,
    deg_size = length(deg_genes),
    asd_size = length(asd_genes)
  ))
}

# Female Analysis
# Get all tested genes from female dataset
female_tested_genes <- unique(female_DESeq2_results$GeneID)

# Translate all female tested genes to human orthologs
female_all_human_orthologs <- homologene::homologene(female_tested_genes, 
                                                     inTax = 10090, 
                                                     outTax = 9606)
print("Female analysis:")
print(paste("Number of female genes with human orthologs:", nrow(female_all_human_orthologs)))

# Female background set
female_background_genes <- nrow(female_all_human_orthologs)

# Get human orthologs for female DEGs
female_deg_orthologs <- female_all_human_orthologs[female_all_human_orthologs$`10090` %in% deg_female$GeneID, ]
print(paste("Number of female DEGs with human orthologs:", nrow(female_deg_orthologs)))

# Merge female DEG information with orthologs
female_merged_data <- merge(deg_female, female_deg_orthologs,
                            by.x = "GeneID", by.y = "10090")

# Filter for ASD-related genes in females
female_asd_related_degs <- female_merged_data[female_merged_data$`9606` %in% asd_gene_list$gene.symbol, ]

# Get top 10 female ASD-related DEGs
female_top_10_asd_degs <- female_asd_related_degs[order(female_asd_related_degs$padj), ]
female_top_10_asd_degs <- female_top_10_asd_degs[1:min(10, nrow(female_top_10_asd_degs)), ]

# Save female results
write.csv(female_top_10_asd_degs, 
          file = file.path(female_dir, "top_10_asd_degs.csv"), 
          row.names = FALSE)

# Remove duplicates from female data
female_merged_data_no_duplicates <- female_merged_data[!duplicated(female_merged_data$GeneID), ]

# Perform statistical test for females
female_stats <- perform_overlap_test(female_merged_data_no_duplicates$`9606`, 
                                     asd_gene_list$gene.symbol,
                                     female_background_genes)

# Create female summary
female_summary <- data.frame(
  Overlap_Size = female_stats$overlap,
  DEG_Size = female_stats$deg_size,
  ASD_Size = female_stats$asd_size,
  Background_Size = female_stats$background_size,
  Odds_Ratio = female_stats$odds_ratio,
  P_Value = female_stats$p_value,
  CI_Lower = female_stats$confidence_interval[1],
  CI_Upper = female_stats$confidence_interval[2]
)

# Save female statistics
write.csv(female_summary, 
          file = file.path(female_dir, "statistics.csv"), 
          row.names = FALSE)

# Print female summary statistics
cat("\nFemale Summary Statistics:\n")
cat("Total background genes (tested AND have human orthologs):", female_background_genes, "\n")
cat("DEGs with human orthologs:", nrow(female_merged_data_no_duplicates), "\n")
cat("ASD-related DEGs:", nrow(female_asd_related_degs), "\n")
cat("Overlap significance: p =", format(female_stats$p_value, digits = 3),
    ", OR =", format(female_stats$odds_ratio, digits = 3), "\n")

# Create female Venn diagram
pdf(file.path(female_dir, "venn_diagram.pdf"))
venn.plot_female <- venn.diagram(
  x = list(
    "DEGs" = female_merged_data_no_duplicates$`9606`,
    "ASD Genes" = asd_gene_list$gene.symbol
  ),
  category.names = c("DEGs", "ASD Genes"),
  fill = c("pink", "yellow"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  filename = NULL
)
grid.draw(venn.plot_female)
grid.text(paste("p =", format(female_stats$p_value, digits = 3),
                "\nOR =", format(female_stats$odds_ratio, digits = 3)),
          x = 0.5, y = 0.1, gp = gpar(fontsize = 12))
dev.off()

# Save female intersection genes
female_asd_intersect <- female_merged_data_no_duplicates$`9606`[female_merged_data_no_duplicates$`9606` %in% asd_gene_list$gene.symbol]
write.csv(data.frame(ASD_Genes = female_asd_intersect),
          file = file.path(female_dir, "asd_intersecting_genes.csv"),
          row.names = FALSE)


# Male Analysis
# Get all tested genes from male dataset
male_tested_genes <- unique(male_DESeq2_results$GeneID)

# Translate all male tested genes to human orthologs
male_all_human_orthologs <- homologene::homologene(male_tested_genes, 
                                                   inTax = 10090, 
                                                   outTax = 9606)
print("\nMale analysis:")
print(paste("Number of male genes with human orthologs:", nrow(male_all_human_orthologs)))

# Male background set
male_background_genes <- nrow(male_all_human_orthologs)

# Get human orthologs for male DEGs
male_deg_orthologs <- male_all_human_orthologs[male_all_human_orthologs$`10090` %in% deg_male$GeneID, ]
print(paste("Number of male DEGs with human orthologs:", nrow(male_deg_orthologs)))

# Merge male DEG information with orthologs
male_merged_data <- merge(deg_male, male_deg_orthologs,
                          by.x = "GeneID", by.y = "10090")

# Filter for ASD-related genes in males
male_asd_related_degs <- male_merged_data[male_merged_data$`9606` %in% asd_gene_list$gene.symbol, ]

# Get top 10 male ASD-related DEGs
male_top_10_asd_degs <- male_asd_related_degs[order(male_asd_related_degs$padj), ]
male_top_10_asd_degs <- male_top_10_asd_degs[1:min(10, nrow(male_top_10_asd_degs)), ]

# Save male results
write.csv(male_top_10_asd_degs, 
          file = file.path(male_dir, "top_10_asd_degs.csv"), 
          row.names = FALSE)

# Remove duplicates from male data
male_merged_data_no_duplicates <- male_merged_data[!duplicated(male_merged_data$GeneID), ]

# Perform statistical test for males
male_stats <- perform_overlap_test(male_merged_data_no_duplicates$`9606`, 
                                   asd_gene_list$gene.symbol,
                                   male_background_genes)

# Create male summary
male_summary <- data.frame(
  Overlap_Size = male_stats$overlap,
  DEG_Size = male_stats$deg_size,
  ASD_Size = male_stats$asd_size,
  Background_Size = male_stats$background_size,
  Odds_Ratio = male_stats$odds_ratio,
  P_Value = male_stats$p_value,
  CI_Lower = male_stats$confidence_interval[1],
  CI_Upper = male_stats$confidence_interval[2]
)

# Save male statistics
write.csv(male_summary, 
          file = file.path(male_dir, "statistics.csv"), 
          row.names = FALSE)

# Print male summary statistics
cat("\nMale Summary Statistics:\n")
cat("Total background genes (tested AND have human orthologs):", male_background_genes, "\n")
cat("DEGs with human orthologs:", nrow(male_merged_data_no_duplicates), "\n")
cat("ASD-related DEGs:", nrow(male_asd_related_degs), "\n")
cat("Overlap significance: p =", format(male_stats$p_value, digits = 3),
    ", OR =", format(male_stats$odds_ratio, digits = 3), "\n")

# Create male Venn diagram
pdf(file.path(male_dir, "venn_diagram.pdf"))
venn.plot_male <- venn.diagram(
  x = list(
    "DEGs" = male_merged_data_no_duplicates$`9606`,
    "ASD Genes" = asd_gene_list$gene.symbol
  ),
  category.names = c("DEGs", "ASD Genes"),
  fill = c("lightblue", "yellow"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  filename = NULL
)
grid.draw(venn.plot_male)
grid.text(paste("p =", format(male_stats$p_value, digits = 3),
                "\nOR =", format(male_stats$odds_ratio, digits = 3)),
          x = 0.5, y = 0.1, gp = gpar(fontsize = 12))
dev.off()

# Save male intersection genes
male_asd_intersect <- male_merged_data_no_duplicates$`9606`[male_merged_data_no_duplicates$`9606` %in% asd_gene_list$gene.symbol]
write.csv(data.frame(ASD_Genes = male_asd_intersect),
          file = file.path(male_dir, "asd_intersecting_genes.csv"),
          row.names = FALSE)
