# Set expression data for specific group analysis
seurat_obj <- PrepareExpressionMatrix(
  seurat_obj,
  target_group = "SC", # specific group to analyze
  grouping_variable = "cell_type", # metadata column for grouping
  assay_type = "RNA", # assay to use
  data_layer = "data" # use normalized expression data
)

# Evaluate optimal soft-thresholding powers
seurat_obj <- EvaluateNetworkPowers(
  seurat_obj,
  network_type = "signed" # options: "unsigned", "signed", "signed hybrid"
)

# Visualize soft-thresholding power analysis results
power_plots <- GeneratePowerPlots(seurat_obj)

# Assemble plots into a cohesive layout
CombinePlots(power_plots, columns = 2)

# Retrieve soft-threshold power evaluation summary
power_results <- ExtractPowerResults(seurat_obj)

# Construct gene co-expression network
seurat_obj <- GenerateCoexpressionNetwork(
  seurat_obj,
  matrix_identifier = "SC" # identifier for the generated overlap matrix
)

# Visualize the dendrogram of gene clustering
VisualizeNetworkDendrogram(seurat_obj, plot_title = "SC Gene Co-expression Network Dendrogram")

# Obtain Topological Overlap Matrix (TOM)
TOM_matrix <- ExtractOverlapMatrix(seurat_obj)


#--------------------------
# Custom Function Definitions
#--------------------------

PrepareExpressionMatrix <- function(seurat_obj, target_group, use_metacells = TRUE, grouping_variable = NULL,
                                    additional_grouping = NULL, additional_group_name = NULL, return_seurat = TRUE,
                                    analysis_name = NULL, assay_type = NULL, data_layer = "data", ...) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis

  genes_selected <- RetrieveAnalysisGenes(seurat_obj, analysis_name)
  metacell_obj <- RetrieveMetacellReference(seurat_obj, analysis_name)

  analysis_obj <- if (use_metacells && !is.null(metacell_obj)) metacell_obj else seurat_obj

  assay_type <- assay_type %||% DefaultAssay(analysis_obj)
  if (!(assay_type %in% names(analysis_obj@assays))) stop("Specified assay not found in the Seurat object.")

  metadata_subset <- analysis_obj@meta.data
  if (!is.null(grouping_variable)) {
    if (!(grouping_variable %in% colnames(metadata_subset))) stop("grouping_variable not found in metadata.")
    metadata_subset <- subset(metadata_subset, metadata_subset[[grouping_variable]] %in% target_group)
  }

  selected_cells <- rownames(metadata_subset)
  expr_matrix <- t(as.data.frame(Seurat::GetAssayData(analysis_obj, assay = assay_type, slot = data_layer)[genes_selected, selected_cells]))

  valid_genes <- genes_selected[WGCNA::goodGenes(expr_matrix, ...)]
  expr_matrix <- expr_matrix[, valid_genes]

  if (return_seurat) {
    seurat_obj <- AssignAnalysisGenes(seurat_obj, valid_genes, analysis_name)
    seurat_obj@misc[[analysis_name]]$expression_matrix <- expr_matrix
    return(seurat_obj)
  }
  return(expr_matrix)
}


EvaluateNetworkPowers <- function(seurat_obj, power_sequence = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
                                  use_metacells = TRUE, network_type = "signed", correlation_method = "bicor",
                                  recalculate_matrix = FALSE, grouping_variable = NULL, target_group = NULL, ...) {
  if (!"expression_matrix" %in% names(seurat_obj@misc[[seurat_obj@misc$active_analysis]]) || recalculate_matrix) {
    seurat_obj <- PrepareExpressionMatrix(seurat_obj,
      use_metacells = use_metacells,
      grouping_variable = grouping_variable, target_group = target_group
    )
  }

  expr_matrix <- seurat_obj@misc[[seurat_obj@misc$active_analysis]]$expression_matrix
  power_analysis <- WGCNA::pickSoftThreshold(expr_matrix,
    powerVector = power_sequence, verbose = 100,
    networkType = network_type, corFnc = correlation_method, ...
  )[[2]]

  seurat_obj@misc[[seurat_obj@misc$active_analysis]]$power_analysis <- power_analysis
  return(seurat_obj)
}


GeneratePowerPlots <- function(seurat_obj, chosen_power = NULL) {
  power_data <- seurat_obj@misc[[seurat_obj@misc$active_analysis]]$power_analysis

  chosen_power <- chosen_power %||% min(power_data$Power[power_data$SFT.R.sq >= 0.8])

  plot_fit <- power_data %>%
    ggplot(aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = chosen_power, linetype = "dashed") +
    labs(y = "Scale-Free Topology Fit", x = "Soft-Threshold Power")

  plot_connectivity <- power_data %>%
    ggplot(aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = chosen_power, linetype = "dashed") +
    labs(y = "Mean Connectivity", x = "Soft-Threshold Power")

  return(list(plot_fit, plot_connectivity))
}


CombinePlots <- function(plot_list, columns = 2) {
  patchwork::wrap_plots(plot_list, ncol = columns)
}


ExtractPowerResults <- function(seurat_obj) {
  seurat_obj@misc[[seurat_obj@misc$active_analysis]]$power_analysis
}


GenerateCoexpressionNetwork <- function(seurat_obj, soft_power = NULL, matrix_identifier = NULL,
                                        analysis_name = NULL, overwrite_matrix = FALSE, ...) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis

  soft_power <- soft_power %||% min(seurat_obj@misc[[analysis_name]]$power_analysis$Power[
    seurat_obj@misc[[analysis_name]]$power_analysis$SFT.R.sq >= 0.8
  ])

  expr_matrix <- seurat_obj@misc[[analysis_name]]$expression_matrix

  network <- WGCNA::blockwiseModules(expr_matrix, power = soft_power, TOMType = "signed", ...)

  seurat_obj@misc[[analysis_name]]$network <- network
  return(seurat_obj)
}


VisualizeNetworkDendrogram <- function(seurat_obj, plot_title = "Gene Dendrogram") {
  network_data <- seurat_obj@misc[[seurat_obj@misc$active_analysis]]$network
  WGCNA::plotDendroAndColors(network_data$dendrograms[[1]], network_data$colors, main = plot_title)
}


ExtractOverlapMatrix <- function(seurat_obj) {
  network_data <- seurat_obj@misc[[seurat_obj@misc$active_analysis]]$network
  TOM <- WGCNA::TOMsimilarityFromExpr(network_data$data, power = network_data$powerEstimate)
  return(TOM)
}
