# Load required libraries
library(DESeq2)
library(VennDiagram)
library(ggplot2)
library(ggrepel)

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
  
  return(data)
}

# Load DESeq2 results for both conditions
ctrl_DESeq2_results <- load_deseq_results("data/DESeq2_results_CTRLs.tabular")
tp_DESeq2_results <- load_deseq_results("data/DESeq2_results_TPs.tabular")

# Function to format DEG tables
format_deg_table <- function(data) {
  return(data.frame(
    GeneID = data$GeneID,
    log2FoldChange = round(data$log2FoldChange, 3),
    pvalue = format(data$pvalue, scientific = TRUE, digits = 3),
    padj = format(data$padj, scientific = TRUE, digits = 3)
  ))
}

# Function to get significant genes
get_sig_genes <- function(data) {
  # Remove NA values and filter by significance
  return(data[!is.na(data$padj) & !is.na(data$log2FoldChange) & data$padj < 0.05, ])
}

# Get significant genes for both conditions
ctrl_sig <- get_sig_genes(ctrl_DESeq2_results)
tp_sig <- get_sig_genes(tp_DESeq2_results)

# Find genes that are only significant in TP condition
new_degs <- tp_sig[!(tp_sig$GeneID %in% ctrl_sig$GeneID), ]
new_degs <- new_degs[order(new_degs$padj), ]

# Save new DEGs
write.table(format_deg_table(new_degs),
            file = file.path(tables_dir, "1_new_DEGs_in_TP.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Find genes that lost significance in TP condition
lost_degs <- ctrl_sig[!(ctrl_sig$GeneID %in% tp_sig$GeneID), ]
lost_degs <- lost_degs[order(lost_degs$padj), ]

# Save genes that lost significance
write.table(format_deg_table(lost_degs),
            file = file.path(tables_dir, "2_lost_DEGs_in_TP.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Create summary file
sink(file.path(tables_dir, "deg_analysis_summary.txt"))
cat("Analysis of Differential Expression Changes between CTRL and TP\n")
cat("========================================================\n\n")

cat("1. DEG Changes:\n")
cat("   Total DEGs in CTRL condition:", nrow(ctrl_sig), "\n")
cat("   Total DEGs in TP condition:", nrow(tp_sig), "\n\n")

cat("2. New and Lost DEGs:\n")
cat("   New DEGs in TP:", nrow(new_degs), "\n")
cat("   Lost DEGs in TP:", nrow(lost_degs), "\n")

sink()

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(
    CTRL = ctrl_sig$GeneID,
    TP = tp_sig$GeneID
  ),
  filename = file.path(plots_dir, "venn_diagram_DEGs.png"),
  category.names = c("CTRL", "TP"),
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("cornflowerblue", "tomato"),
  fill = c(alpha("cornflowerblue", 0.3), alpha("tomato", 0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# Create visualization of top changed DEGs
# Prepare data for new DEGs
top_new_degs <- head(new_degs[order(abs(new_degs$log2FoldChange), decreasing = TRUE), ], 10)
top_new_degs$change_type <- "New in TP"

# Prepare data for lost DEGs
top_lost_degs <- head(lost_degs[order(abs(lost_degs$log2FoldChange), decreasing = TRUE), ], 10)
top_lost_degs$change_type <- "Lost in TP"

# Combine the data
plot_data <- rbind(top_new_degs, top_lost_degs)
plot_data$GeneID <- factor(plot_data$GeneID, 
                           levels = plot_data$GeneID[order(plot_data$log2FoldChange)])

# Create the plot
changed_degs_plot <- ggplot(plot_data, 
                            aes(x = log2FoldChange, y = GeneID, fill = change_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Lost in TP" = "cornflowerblue", "New in TP" = "tomato")) +
  theme_minimal() +
  labs(
    title = "Top 10 Changed DEGs in TP Treatment",
    x = "log2 Fold Change",
    y = "Gene ID",
    fill = NULL  # Remove legend title
  ) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave(file.path(plots_dir, "top_changed_DEGs.png"), 
       changed_degs_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Print final summary to console
cat("\nAnalysis complete! Created the following files:\n")
cat("\nIn", tables_dir, ":\n")
cat("1. New DEGs in TP condition (1_new_DEGs_in_TP.txt)\n")
cat("2. Lost DEGs in TP condition (2_lost_DEGs_in_TP.txt)\n")
cat("3. Summary file (deg_analysis_summary.txt)\n")
cat("\nIn", plots_dir, ":\n")
cat("1. Venn diagram (venn_diagram_DEGs.png)\n")
cat("2. Volcano plots (volcano_plots.png)\n")