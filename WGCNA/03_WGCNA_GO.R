# Load required libraries
library(WGCNA)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)

# Source common functions
source("scripts/WGCNA_all/00_WGCNA_common.R")

# Load saved workspace
load("wgcna_output/all/complete_analysis.RData")

# Check data consistency
check_data_consistency()

# Create output directory
dir.create("wgcna_output/all/enrichment", recursive = TRUE, showWarnings = FALSE)

# Get standardized sample ordering
sample_order <- create_sample_order()
sample_info <- sample_order$info

# Create sample information with new order
sample_info <- data.frame(
  SampleID = sample_order$order,
  Treatment = factor(c(
    rep("CTRL", 8),  # Controls
    rep("TP", 7)     # Treatment
  ), levels = c("CTRL", "TP")),
  Sample_Number = c(1:8, 1:7)
)

# Reorder the expression data and traits
all_datExpr <- all_datExpr[sample_order$order,]
treatment <- factor(ifelse(grepl("^CTRL", sample_order$order), "Control", "Treatment"))

# Calculate module eigengenes
MEs <- moduleEigengenes(all_datExpr, colors = all_moduleColors)$eigengenes

# Create a mapping between numeric labels and color names
module_color_map <- create_module_color_map()

# Calculate correlations with treatment
ME_correlation_treatment <- cor(MEs, as.numeric(treatment == "Treatment"), use = "p")
ME_pvalue_treatment <- corPvalueStudent(ME_correlation_treatment, length(treatment))

# Create a data frame for plotting
plot_data <- data.frame(
  Module = module_color_map[gsub("ME", "", rownames(ME_correlation_treatment))],
  Correlation_Treatment = ME_correlation_treatment[,1],
  PValue_Treatment = ME_pvalue_treatment[,1]
)

# Create the treatment effect plot
p_treatment <- ggplot(plot_data, 
                      aes(x = reorder(Module, Correlation_Treatment), 
                          y = Correlation_Treatment)) +
  geom_bar(stat = "identity", fill = "grey50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module Correlations",
       x = "Modules", 
       y = "Correlation with Treatment")

# Save the plot
ggsave("wgcna_output/all/enrichment/module_correlation_plot.png", 
       p_treatment, width = 12, height = 8)

# Function to perform enrichment analysis for a single module
perform_enrichment <- function(module_genes, module_color) {
  cat("\nPerforming enrichment for module:", module_color, "\n")
  cat("Number of genes in module:", length(module_genes), "\n")
  
  if (length(module_genes) < 5) {
    message("Module ", module_color, " has fewer than 5 genes. Skipping enrichment analysis.")
    return(NULL)
  }
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(module_genes, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)
  
  if (nrow(entrez_ids) == 0) {
    message("No valid Entrez IDs found for module: ", module_color)
    return(NULL)
  }
  
  # Log the number of successfully converted genes
  message(paste0("Module ", module_color, ": ", 
                 nrow(entrez_ids), " out of ", 
                 length(module_genes), 
                 " genes successfully converted to Entrez IDs"))
  
  # Perform GO enrichment analysis
  go_enrichment <- tryCatch({
    enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    message("Error in GO enrichment for module ", module_color, ": ", e$message)
    return(NULL)
  })
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- tryCatch({
    enrichKEGG(
      gene = entrez_ids$ENTREZID,
      organism = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    message("Error in KEGG enrichment for module ", module_color, ": ", e$message)
    return(NULL)
  })
  
  return(list(go = go_enrichment, kegg = kegg_enrichment))
}

# Perform enrichment analysis for each module
enrichment_results <- list()
for (module_num in unique(all_moduleColors)) {
  module_color <- module_color_map[as.character(module_num)]
  module_genes <- colnames(all_datExpr)[all_moduleColors == module_num]
  
  result <- tryCatch({
    perform_enrichment(module_genes, module_color)
  }, error = function(e) {
    message("Error in enrichment analysis for module ", module_color, ": ", e$message)
    NULL
  })
  
  if (!is.null(result)) {
    enrichment_results[[module_color]] <- result
  }
}

# Function to create dotplots using ggplot2
create_dotplot <- function(enrichment_data, title, output_file, wrap_width = 50) {
  if (!is.null(enrichment_data) && nrow(enrichment_data@result) > 0) {
    # Check if there are significant terms after p-value adjustment
    sig_terms <- sum(enrichment_data@result$p.adjust < 0.05)
    
    if (sig_terms > 0) {
      tryCatch({
        # Convert enrichment results to data frame and select top terms
        plot_data <- as.data.frame(enrichment_data@result) %>%
          arrange(p.adjust) %>%
          head(min(10, sig_terms)) %>%
          mutate(
            Description = stringr::str_wrap(Description, width = wrap_width),
            Description = factor(Description, levels = rev(Description)),
            GeneRatio = sapply(GeneRatio, function(x) {
              nums <- as.numeric(strsplit(x, "/")[[1]])
              nums[1] / nums[2]
            })
          )
        
        # Create ggplot
        p <- ggplot(plot_data, 
                    aes(x = GeneRatio, 
                        y = Description,
                        size = Count,
                        color = p.adjust)) +
          geom_point() +
          scale_color_gradient(
            low = "red",
            high = "blue",
            trans = "log10",
            name = "Adjusted\np-value"
          ) +
          scale_size_continuous(
            name = "Gene Count",
            range = c(4, 12)
          ) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 18),
            axis.text.x = element_text(size = 15),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
          #  panel.grid.major = element_line(color = "grey90"),
          #  panel.grid.minor = element_line(color = "grey95")
          #  plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
          ) +
          labs(
            title = title,
            x = "Gene Ratio",
            y = ""
          )
        
        # Save plot
        ggsave(output_file, p, width = 10, height = 10, dpi = 300)
        message("Successfully created ggplot dotplot: ", output_file)
        
      }, error = function(e) {
        message("Error in creating dotplot: ", e$message)
      })
    } else {
      message("No significant terms (p.adjust < 0.05) found for: ", title)
    }
  } else {
    message("No enrichment data available for: ", title)
  }
}

# Create dotplots for each module's enrichment results
for (module_color in names(enrichment_results)) {
  if (!is.null(enrichment_results[[module_color]])) {
    go_results <- enrichment_results[[module_color]]$go
    kegg_results <- enrichment_results[[module_color]]$kegg
    
    create_dotplot(
      go_results, 
      paste("GO Enrichment -", module_color, "Module"),
      paste0("wgcna_output/all/enrichment/", module_color, "_GO_dotplot.pdf"),
      wrap_width = 30
    )
    
    create_dotplot(
      kegg_results, 
      paste("KEGG Enrichment -", module_color, "Module"),
      paste0("wgcna_output/all/enrichment/", module_color, "_KEGG_dotplot.pdf"),
      wrap_width = 30
    )
  }
}

# Save results
save(enrichment_results, plot_data, 
     file = "wgcna_output/all/enrichment/enrichment_results.RData")

# Create summary file
sink("wgcna_output/all/enrichment/enrichment_summary.txt")
cat("Enrichment Analysis Summary\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

for (module_color in names(enrichment_results)) {
  cat("\nModule:", module_color, "\n")
  
  if (!is.null(enrichment_results[[module_color]]$go)) {
    cat("GO terms found:", nrow(enrichment_results[[module_color]]$go@result), "\n")
  } else {
    cat("No significant GO terms found\n")
  }
  
  if (!is.null(enrichment_results[[module_color]]$kegg)) {
    cat("KEGG pathways found:", nrow(enrichment_results[[module_color]]$kegg@result), "\n")
  } else {
    cat("No significant KEGG pathways found\n")
  }
}
sink()

print("Enrichment analysis complete. Check wgcna_output/all/enrichment/ for results.")