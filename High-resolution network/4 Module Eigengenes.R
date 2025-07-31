# Section 4: Module Eigengene and Connectivity Analysis

# First, perform data scaling to prevent downstream analysis errors
# seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# Compute module eigengenes across the dataset
dataset_obj <- ComputeAllModuleEigengenes(
  seurat_obj,
  harmonize_var = "Sample"
)

# Extract harmonized module eigengenes
harmonized_eigengenes <- GetModuleEigengenes(dataset_obj)

# Extract original (non-harmonized) module eigengenes
original_eigengenes <- GetModuleEigengenes(dataset_obj, harmonized = FALSE)

# Compute gene connectivity based on module eigengenes (kME)
dataset_obj <- CalculateGeneConnectivity(
  dataset_obj,
  group_var = 'cell_type',
  selected_group = 'SC'
)

# Score cells based on the expression of the top 25 genes within each module using UCell
library(UCell)
dataset_obj <- ScoreTopModuleGenes(
  dataset_obj,
  top_gene_count = 25,
  scoring_algorithm = 'UCell'
)

# Plot FeaturePlots for each module eigengene
module_feature_plots <- GenerateModuleFeaturePlots(
  dataset_obj,
  eigengene_type = 'harmonized_eigengenes',
  order_cells = TRUE
)

# Combine FeaturePlots into a cohesive figure
CombinePlots(module_feature_plots, columns = 6)


#--------------------------
# Custom Function Definitions
#--------------------------

ComputeAllModuleEigengenes <- function(seurat_obj, harmonize_var = NULL, ...) {
  modules <- GetModuleAssignments(seurat_obj)
  unique_modules <- unique(modules$module)
  
  eigengene_results <- list()
  harmonized_results <- list()
  
  for(mod in unique_modules){
    seurat_obj <- CalculateModuleEigengene(
      seurat_obj,
      module = mod,
      module_data = modules,
      harmonize_var = harmonize_var,
      ...
    )
    
    eigengene_results[[mod]] <- seurat_obj@reductions$ModuleEigengene@cell.embeddings[,1]
    
    if(!is.null(harmonize_var)){
      harmonized_results[[mod]] <- seurat_obj@reductions$HarmonizedEigengene@cell.embeddings[,1]
    }
  }
  
  seurat_obj <- SaveEigengenes(seurat_obj, eigengene_results, harmonized=FALSE)
  
  if(!is.null(harmonize_var)){
    seurat_obj <- SaveEigengenes(seurat_obj, harmonized_results, harmonized=TRUE)
  }
  
  return(seurat_obj)
}

GetModuleEigengenes <- function(seurat_obj, harmonized = TRUE) {
  if(harmonized){
    return(seurat_obj@misc$harmonized_eigengenes)
  } else {
    return(seurat_obj@misc$original_eigengenes)
  }
}

CalculateGeneConnectivity <- function(seurat_obj, group_var, selected_group, ...) {
  cells_to_use <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_var]] == selected_group, ])
  expr_data <- GetAssayData(seurat_obj, slot='data')[, cells_to_use]
  eigengene_data <- seurat_obj@misc$harmonized_eigengenes[cells_to_use, ]
  
  connectivity_matrix <- cor(as.matrix(t(expr_data)), as.matrix(eigengene_data))
  rownames(connectivity_matrix) <- rownames(expr_data)
  
  seurat_obj@misc$gene_connectivity <- connectivity_matrix
  
  return(seurat_obj)
}

ScoreTopModuleGenes <- function(seurat_obj, top_gene_count=25, scoring_algorithm='UCell', ...) {
  modules <- GetModuleAssignments(seurat_obj)
  module_list <- unique(modules$module)
  
  gene_sets <- lapply(module_list, function(mod){
    module_genes <- modules[modules$module == mod, ]
    top_genes <- module_genes[order(-module_genes$kME), 'gene_name'][1:top_gene_count]
    return(top_genes)
  })
  names(gene_sets) <- module_list
  
  if(scoring_algorithm == 'UCell'){
    scored_data <- UCell::AddModuleScore_UCell(seurat_obj, features=gene_sets, ...)@meta.data
  } else if(scoring_algorithm == 'Seurat'){
    scored_data <- Seurat::AddModuleScore(seurat_obj, features=gene_sets, ...)@meta.data
  }
  
  seurat_obj@misc$module_gene_scores <- scored_data
  return(seurat_obj)
}

GenerateModuleFeaturePlots <- function(seurat_obj, eigengene_type='harmonized_eigengenes', order_cells=TRUE, ...) {
  eigengene_data <- if(eigengene_type == 'harmonized_eigengenes') seurat_obj@misc$harmonized_eigengenes else seurat_obj@misc$original_eigengenes
  
  plots <- lapply(colnames(eigengene_data), function(feature){
    FeaturePlot(seurat_obj, features=feature, order=order_cells, ...) + ggtitle(feature)
  })
  
  return(plots)
}

CombinePlots <- function(plot_list, columns=6){
  patchwork::wrap_plots(plot_list, ncol=columns)
}

GetModuleAssignments <- function(seurat_obj){
  seurat_obj@misc$module_assignments
}

SaveEigengenes <- function(seurat_obj, eigengene_list, harmonized=FALSE){
  eigengene_matrix <- do.call(cbind, eigengene_list)
  if(harmonized){
    seurat_obj@misc$harmonized_eigengenes <- eigengene_matrix
  } else {
    seurat_obj@misc$original_eigengenes <- eigengene_matrix
  }
  return(seurat_obj)
}

CalculateModuleEigengene <- function(seurat_obj, module, module_data, harmonize_var=NULL, ...){
  genes_in_module <- module_data$gene_name[module_data$module == module]
  seurat_obj <- AddModuleScore(seurat_obj, features=list(genes_in_module), name=module, ...)
  
  if(!is.null(harmonize_var)){
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars=harmonize_var, reduction.save='HarmonizedEigengene')
  }
  
  return(seurat_obj)
}
