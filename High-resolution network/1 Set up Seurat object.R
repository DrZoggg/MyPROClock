# High-resolution network analysis using custom Seurat functions
# Step 1: Initialize Seurat object for custom network analysis

# Ensure the Seurat object is up-to-date
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)

# Prepare Seurat object for custom network analysis
seurat_obj <- InitializeNetworkAnalysis(
  seurat_obj,
  gene_method = "fraction", # Strategy for gene selection
  fraction_threshold = 0.05, # Minimum fraction of cells expressing genes
  analysis_name = "custom_network" # Identifier for the custom analysis
)


#--------------------------
# Custom Function Definitions
#--------------------------

InitializeNetworkAnalysis <- function(seurat_obj, analysis_name, selected_features = NULL, metacell_ref = NULL, ...) {
  # Set active analysis session
  seurat_obj <- ActivateNetworkAnalysis(seurat_obj, analysis_name)

  # Select genes for analysis based on user input or default strategy
  if (is.null(selected_features)) {
    seurat_obj <- FilterNetworkGenes(seurat_obj, analysis_name = analysis_name, ...)
  } else {
    seurat_obj <- FilterNetworkGenes(
      seurat_obj,
      analysis_name = analysis_name,
      gene_method = "custom",
      feature_list = selected_features
    )
  }

  # Warn user if too few genes are selected
  genes_selected <- RetrieveAnalysisGenes(seurat_obj, analysis_name)
  if (length(genes_selected) < 500) {
    warning(sprintf("%d genes selected. Consider adjusting your parameters to include more genes.", length(genes_selected)))
  }

  # Incorporate metacell information, if provided
  if (!is.null(metacell_ref)) {
    if (!(metacell_ref %in% names(seurat_obj@misc))) {
      stop("Specified metacell_ref not found in seurat_obj@misc.")
    }
    seurat_obj <- AssignMetacellReference(seurat_obj, metacell_ref, analysis_name)
  }

  return(seurat_obj)
}


ActivateNetworkAnalysis <- function(seurat_obj, analysis_name) {
  seurat_obj@misc$active_analysis <- analysis_name
  if (!(analysis_name %in% names(seurat_obj@misc))) {
    seurat_obj@misc[[analysis_name]] <- list()
  }
  return(seurat_obj)
}


FilterNetworkGenes <- function(seurat_obj, gene_method = "fraction", fraction_threshold = 0.05, group_variable = NULL, feature_list = NULL, assay_type = NULL, analysis_name = NULL) {
  # Validate gene selection strategy
  valid_methods <- c("variable", "fraction", "all", "custom")
  if (!(gene_method %in% valid_methods)) {
    stop(sprintf('Invalid gene_method option: "%s". Choose from: %s', gene_method, paste(valid_methods, collapse = ", ")))
  }

  # Determine active assay
  assay_type <- assay_type %||% DefaultAssay(seurat_obj)

  # Extract expression matrix
  expr_matrix <- if (CheckSeurat5()) {
    SeuratObject::LayerData(seurat_obj, layer = "counts", assay = assay_type)
  } else {
    Seurat::GetAssayData(seurat_obj, slot = "counts", assay = assay_type)
  }

  # Gene selection logic
  if (gene_method == "fraction") {
    # Binarize counts to reduce memory usage
    expr_binary <- expr_matrix
    expr_binary[expr_binary > 0] <- 1

    if (!is.null(group_variable)) {
      groups <- unique(seurat_obj@meta.data[[group_variable]])
      feature_list <- unique(unlist(lapply(groups, function(grp) {
        cur_expr <- expr_binary[, seurat_obj@meta.data[[group_variable]] == grp]
        rownames(cur_expr)[Matrix::rowSums(cur_expr) >= fraction_threshold * ncol(cur_expr)]
      })))
    } else {
      feature_list <- rownames(expr_binary)[Matrix::rowSums(expr_binary) >= fraction_threshold * ncol(seurat_obj)]
    }
  } else if (gene_method == "variable") {
    feature_list <- VariableFeatures(seurat_obj)
  } else if (gene_method == "all") {
    feature_list <- rownames(seurat_obj)
  } else if (gene_method == "custom") {
    feature_list <- unique(feature_list)
    missing_features <- setdiff(feature_list, rownames(seurat_obj))
    if (length(missing_features) > 0) {
      stop(sprintf("The following genes are missing from seurat_obj: %s", paste(missing_features, collapse = ", ")))
    }
  }

  if (length(feature_list) == 0) stop("No genes selected. Adjust selection criteria.")
  if (length(feature_list) <= 100) warning(sprintf("Only %d genes selected; consider alternative selection strategies.", length(feature_list)))

  seurat_obj <- AssignAnalysisGenes(seurat_obj, feature_list, analysis_name)
  return(seurat_obj)
}


RetrieveAnalysisGenes <- function(seurat_obj, analysis_name = NULL) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis
  seurat_obj@misc[[analysis_name]]$analysis_genes
}


AssignMetacellReference <- function(seurat_obj, metacell_obj, analysis_name = NULL) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis
  seurat_obj@misc[[analysis_name]]$analysis_metacell_obj <- metacell_obj
  return(seurat_obj)
}


