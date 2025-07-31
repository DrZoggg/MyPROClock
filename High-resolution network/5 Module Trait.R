




# Define traits for correlation analysis
selected_traits <- c('plaque progression')

# Perform correlation analysis between modules and traits
seurat_obj <- ComputeModuleTraitCorrelation(
  seurat_obj,
  traits = selected_traits,
  grouping_var = 'cell_type'
)


#--------------------------
# Custom Function Definitions
#--------------------------

ComputeModuleTraitCorrelation <- function (seurat_obj, traits, grouping_var = NULL, module_features = "harmonized_eigengenes",
                                           correlation_method = "pearson", subset_var = NULL, subset_values = NULL,
                                           analysis_name = NULL, ...) {
  
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis
  
  # Extract module eigengenes based on the specified features
  if (module_features == "harmonized_eigengenes") {
    eigengenes <- seurat_obj@misc$harmonized_eigengenes
  } else if (module_features == "original_eigengenes") {
    eigengenes <- seurat_obj@misc$original_eigengenes
  } else if (module_features == "scores") {
    eigengenes <- seurat_obj@misc$module_gene_scores
  } else {
    stop("Invalid module_features selection. Options: harmonized_eigengenes, original_eigengenes, scores")
  }
  
  # Subset Seurat object if required
  if (!is.null(subset_var)) {
    cells_to_keep <- seurat_obj@meta.data[[subset_var]] %in% subset_values
    eigengenes <- eigengenes[cells_to_keep, ]
    seurat_obj <- seurat_obj[, cells_to_keep]
  }
  
  # Validate trait availability
  missing_traits <- setdiff(traits, colnames(seurat_obj@meta.data))
  if (length(missing_traits) > 0) {
    stop(paste("Traits not found in metadata:", paste(missing_traits, collapse = ", ")))
  }
  
  # Check data types of traits
  valid_classes <- c("numeric", "factor", "integer")
  trait_classes <- sapply(traits, function(trait) class(seurat_obj@meta.data[[trait]]))
  
  if (!all(trait_classes %in% valid_classes)) {
    invalid_traits <- traits[!(trait_classes %in% valid_classes)]
    stop(paste("Invalid trait types:", paste(invalid_traits, collapse = ", "), ". Must be numeric, factor, or integer."))
  }
  
  # Convert factor traits to numeric with warning
  factor_traits <- traits[trait_classes == "factor"]
  if (length(factor_traits) > 0) {
    for (ft in factor_traits) {
      warning(sprintf("Converting factor trait '%s' levels (%s) to numeric.", ft,
                      paste(levels(seurat_obj@meta.data[[ft]]), collapse = ", ")))
      seurat_obj@meta.data[[ft]] <- as.numeric(seurat_obj@meta.data[[ft]])
    }
  }
  
  # Retrieve module names excluding 'grey'
  module_assignments <- seurat_obj@misc$module_assignments
  modules <- unique(module_assignments$module)
  modules <- modules[modules != "grey"]
  
  # Prepare trait and module eigengene data
  trait_data <- seurat_obj@meta.data[, traits, drop = FALSE]
  if (length(traits) == 1) {
    trait_data <- data.frame(trait_data)
    colnames(trait_data) <- traits
  }
  
  # Initialize lists for storing results
  correlation_results <- list()
  pvalue_results <- list()
  adjusted_pvalue_results <- list()
  
  # Compute correlation and p-values
  correlation_matrix <- Hmisc::rcorr(as.matrix(trait_data), as.matrix(eigengenes), type = correlation_method)
  current_cor <- correlation_matrix$r[traits, modules, drop = FALSE]
  current_pvalues <- correlation_matrix$P[traits, modules, drop = FALSE]
  
  # Adjust p-values (FDR correction)
  pvalues_long <- reshape2::melt(current_pvalues)
  if (length(traits) == 1) {
    pvalues_long$Var1 <- traits
    pvalues_long$Var2 <- factor(modules, levels = modules)
    pvalues_long <- pvalues_long[, c("Var1", "Var2", "value")]
  }
  pvalues_long$fdr <- p.adjust(pvalues_long$value, method = "fdr")
  adjusted_pvalues <- reshape2::dcast(pvalues_long, Var1 ~ Var2, value.var = "fdr")
  rownames(adjusted_pvalues) <- adjusted_pvalues$Var1
  adjusted_pvalues <- adjusted_pvalues[, -1]
  
  # Store results
  correlation_results[["all_cells"]] <- current_cor
  pvalue_results[["all_cells"]] <- current_pvalues
  adjusted_pvalue_results[["all_cells"]] <- adjusted_pvalues
  
  # Group-wise correlations if grouping_var provided
  if (!is.null(grouping_var)) {
    trait_data$group <- seurat_obj@meta.data[[grouping_var]]
    eigengenes$group <- seurat_obj@meta.data[[grouping_var]]
    
    group_levels <- unique(seurat_obj@meta.data[[grouping_var]])
    
    for (group in group_levels) {
      trait_subset <- trait_data[trait_data$group == group, traits, drop = FALSE]
      eigengene_subset <- eigengenes[eigengenes$group == group, modules, drop = FALSE]
      
      group_cor_matrix <- Hmisc::rcorr(as.matrix(trait_subset), as.matrix(eigengene_subset), type = correlation_method)
      group_cor <- group_cor_matrix$r[traits, modules, drop = FALSE]
      group_pvalues <- group_cor_matrix$P[traits, modules, drop = FALSE]
      
      group_pvalues_long <- reshape2::melt(group_pvalues)
      if (length(traits) == 1) {
        group_pvalues_long$Var1 <- traits
        group_pvalues_long$Var2 <- factor(modules, levels = modules)
        group_pvalues_long <- group_pvalues_long[, c("Var1", "Var2", "value")]
      }
      group_pvalues_long$fdr <- p.adjust(group_pvalues_long$value, method = "fdr")
      group_adjusted_pvalues <- reshape2::dcast(group_pvalues_long, Var1 ~ Var2, value.var = "fdr")
      rownames(group_adjusted_pvalues) <- group_adjusted_pvalues$Var1
      group_adjusted_pvalues <- group_adjusted_pvalues[, -1]
      
      correlation_results[[group]] <- group_cor
      pvalue_results[[group]] <- group_pvalues
      adjusted_pvalue_results[[group]] <- as.matrix(group_adjusted_pvalues)
    }
  }
  
  module_trait_correlation <- list(
    correlation = correlation_results,
    p_values = pvalue_results,
    adjusted_p_values = adjusted_pvalue_results
  )
  
  # Save correlation results to Seurat object
  seurat_obj@misc[[analysis_name]]$module_trait_correlation <- module_trait_correlation
  
  return(seurat_obj)
}









