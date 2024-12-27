# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(WGCNA)
  library(cluster)
  library(ggplot2)
  library(ggpubr)
})

# Source common functions
source("scripts/WGCNA_all/00_WGCNA_common.R")

# Enable multi-threading for better performance
disableWGCNAThreads()

# Create output directory structure
dir.create("wgcna_output/all", recursive = TRUE, showWarnings = FALSE)

# Load normalized counts data
all_data <- read.table("data/normalized_counts_DESeq2_all.tabular", 
                       header = TRUE, row.names = 1, sep = "\t")

# Create sample metadata
sample_info <- data.frame(
  SampleID = colnames(all_data),
  Group = factor(c(rep("CTRL_F", 4), rep("CTRL_M", 4), 
                   rep("TP_F", 4), rep("TP_M", 3)),
                 levels = c("CTRL_F", "CTRL_M", "TP_F", "TP_M")),
  Treatment = factor(c(rep("CTRL", 8), rep("TP", 7))),
  Sex = factor(c(rep("F", 4), rep("M", 4), rep("F", 4), rep("M", 3)))
)

# Check if samples match between metadata and expression data
if (!all(sample_info$SampleID == colnames(all_data))) {
  stop("Sample IDs in metadata don't match expression data")
}

# Filter low expression genes
mean_counts <- rowMeans(all_data)
all_data_filt <- all_data[mean_counts > 50, ]
print(paste("Genes after filtering:", nrow(all_data_filt)))

# Calculate coefficient of variation
cv <- apply(all_data_filt, 1, function(x) sd(x)/mean(x))
hist(cv, breaks = 50, main = "Distribution of Coefficient of Variation")

# Log2 transformation
all_datExpr <- log2(all_data_filt + 1)

# Transpose for WGCNA (samples as rows)
all_datExpr <- t(all_datExpr)

# Sample clustering to detect outliers
sampleTree <- hclust(dist(all_datExpr), method = "average")

# Plot sample clustering
pdf("wgcna_output/all/sample_clustering.pdf", width = 12, height = 8)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample Clustering", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Add colored bars for groups
group_colors <- c(CTRL_F = "skyblue", CTRL_M = "blue", 
                  TP_F = "pink", TP_M = "red")
plotColorUnderTree(sampleTree, colors = group_colors[sample_info$Group],
                   groupLabels = "Treatment-Sex Groups")
dev.off()

# Check for missing values
if(sum(is.na(all_datExpr)) > 0) {
  stop("Error: Missing values detected in the expression data")
}

# Soft thresholding power analysis
powers <- c(1:20)
sft <- pickSoftThreshold(all_datExpr, 
                         powerVector = powers,
                         verbose = 5,
                         networkType = "unsigned")

# Plot power analysis results
pdf("wgcna_output/all/power_analysis_plot.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
# Scale independence plot
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     main = "Scale independence",
     type = "n")
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, 
     col = "red")
abline(h = 0.80, col = "red")

# Mean connectivity plot
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels = powers, 
     col = "red")
dev.off()

# Select power
power <- min(sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > 0.80])
if(is.infinite(power)) {
  warning("No power reached the 0.85 threshold. Choosing the power with the highest R^2.")
  power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
}
print(paste("Selected power:", power))

# Network construction and module detection
all_net <- blockwiseModules(all_datExpr, 
                            power = power,
                            TOMType = "unsigned", 
                            minModuleSize = 30,
                            reassignThreshold = 0, 
                            mergeCutHeight = 0.5,
                            numericLabels = TRUE, 
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "all_TOM",
                            verbose = 3,
                            maxBlockSize = 8000)

# Save network
save(all_net, file = "wgcna_output/all/all_net.RData")

# Convert labels to colors
all_moduleColors <- labels2colors(all_net$colors)

# Get the gene ordering from the dendrogram
all_geneTreeOrder <- all_net$dendrograms[[1]]$order
all_moduleColorsOrdered <- all_moduleColors[all_geneTreeOrder]

# Plot the dendrogram with just module colors
pdf("wgcna_output/all/module_dendrogram.pdf", width = 15, height = 3)
plotDendroAndColors(all_net$dendrograms[[1]],
                    colors = all_moduleColorsOrdered,
                    groupLabels = "Modules",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

# Calculate Module Eigengenes
MEs <- moduleEigengenes(all_datExpr, colors = all_moduleColors)$eigengenes

# Plot eigenmodule relationships
pdf("wgcna_output/all/eigengene_dendrogram.pdf", width = 10, height = 6)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", 
                      marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

# Save the module color assignments
write.table(data.frame(Gene = colnames(all_datExpr),
                       Module = all_moduleColors),
            file = "wgcna_output/all/module_assignments.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

print("Dendrogram plots have been created successfully.")

# Save workspace
save.image(file = "wgcna_output/all/complete_analysis.RData")

# Create summary file
sink("wgcna_output/all/analysis_summary.txt")
cat("WGCNA Analysis Summary\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Number of samples:", nrow(all_datExpr), "\n")
cat("Number of genes:", ncol(all_datExpr), "\n")
cat("Selected power:", power, "\n")
cat("\nSample groups:\n")
print(table(sample_info$Group))
cat("\nModule sizes:\n")
print(table(all_moduleColors))
sink()

print("Analysis complete. Check the output directory for results.")