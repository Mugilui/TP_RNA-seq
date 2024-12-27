# Common functions and utilities for WGCNA analysis
# Save this as 00_WGCNA_common.R

# Function to create standardized sample ordering
create_sample_order <- function() {
  # Create the desired order (females first, then males)
  desired_sample_order <- c(
    grep("CTRL_F\\d+", rownames(all_datExpr), value = TRUE),
    grep("TP_F\\d+", rownames(all_datExpr), value = TRUE),
    grep("CTRL_M\\d+", rownames(all_datExpr), value = TRUE),
    grep("TP_M\\d+", rownames(all_datExpr), value = TRUE)
  )
  
  # Create sample information with new order
  sample_info <- data.frame(
    SampleID = desired_sample_order,
    Group = factor(c(
      rep("CTRL_F", 4), # Female controls
      rep("TP_F", 4),   # Female treatment
      rep("CTRL_M", 4), # Male controls
      rep("TP_M", 3)    # Male treatment
    ), levels = c("CTRL_F", "TP_F", "CTRL_M", "TP_M")),
    Sample_Number = c(1:4, 1:4, 1:4, 1:3)
  )
  
  return(list(
    order = desired_sample_order,
    info = sample_info
  ))
}

# Function to check data consistency
check_data_consistency <- function() {
  required_objects <- c("all_datExpr", "all_moduleColors", "MEs")
  missing_objects <- required_objects[!sapply(required_objects, exists)]
  
  if (length(missing_objects) > 0) {
    stop("Missing required objects: ", paste(missing_objects, collapse = ", "), 
         ". Please ensure complete_analysis.RData is loaded.")
  }
  
  if (!file.exists("wgcna_output/all/complete_analysis.RData")) {
    stop("Analysis data file not found. Please run 01_WGCNA.R first.")
  }
}

# Function to create standardized module color mapping
create_module_color_map <- function() {
  module_nums <- unique(all_moduleColors)
  module_colors <- labels2colors(module_nums)
  names(module_colors) <- module_nums
  return(module_colors)
}