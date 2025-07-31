



# Construct metacells within each group
seurat_obj <- BuildGroupMetacells(
  seurat_obj = seurat_obj,
  grouping_vars = c("cell_type", "Sample"), # columns in metadata used for grouping
  dimension_reduction = "harmony", # dimensionality reduction method
  neighbors_k = 25, # number of nearest neighbors
  max_overlap = 10, # max shared cells between two metacells
  identity_var = "cell_type" # metadata column for cell identities
)

# Normalize metacell expression matrix
seurat_obj <- NormalizeGroupMetacells(seurat_obj)


#--------------------------
# Custom Function Definitions
#--------------------------

BuildGroupMetacells <- function(seurat_obj, grouping_vars = c("seurat_clusters"), identity_var = "seurat_clusters",
                                neighbors_k = 25, dimension_reduction = "pca", assay_type = NULL, selected_cells = NULL,
                                data_slot = "counts", aggregation_mode = "average", min_cells_group = 100, max_overlap = 15,
                                target_metacell_number = 1000, iteration_limit = 5000, verbose_output = FALSE,
                                analysis_name = NULL) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis

  if (any(grepl("#", grouping_vars))) {
    stop("Grouping variable names must not contain '#'.")
  }
  if (!(identity_var %in% grouping_vars)) {
    stop("identity_var must be included in grouping_vars")
  }
  if (!(aggregation_mode %in% c("sum", "average"))) {
    stop("aggregation_mode must be either 'sum' or 'average'.")
  }
  if (!(dimension_reduction %in% names(seurat_obj@reductions))) {
    stop(sprintf("Reduction method '%s' not found. Available methods: %s", dimension_reduction, paste(names(seurat_obj@reductions), collapse = ", ")))
  }
  assay_type <- assay_type %||% DefaultAssay(seurat_obj)
  if (!(data_slot %in% c("counts", "data", "scale.data"))) {
    stop("data_slot must be one of 'counts', 'data', or 'scale.data'.")
  }
  if (!is.null(selected_cells)) {
    original_obj <- seurat_obj
    seurat_obj <- seurat_obj[, selected_cells]
  }

  seurat_meta <- seurat_obj@meta.data[, grouping_vars, drop = FALSE]
  seurat_meta[] <- lapply(seurat_meta, as.character)
  seurat_obj$metacell_group <- apply(seurat_meta, 1, paste, collapse = "#")

  group_counts <- table(seurat_obj$metacell_group)
  valid_groups <- names(group_counts[group_counts >= min_cells_group])
  if (length(valid_groups) == 0) stop("No groups have sufficient cells (min_cells_group requirement unmet).")
  seurat_list <- lapply(valid_groups, function(grp) seurat_obj[, seurat_obj$metacell_group == grp])

  metacell_objs <- mapply(CreateMetacells,
    seurat_obj = seurat_list,
    group_name = valid_groups,
    MoreArgs = list(
      k = neighbors_k,
      reduction = dimension_reduction,
      assay = assay_type,
      slot = data_slot,
      return_metacell = TRUE,
      mode = aggregation_mode,
      max_shared = max_overlap,
      max_iter = iteration_limit,
      target_metacells = target_metacell_number,
      verbose = verbose_output,
      analysis_name = analysis_name
    ), SIMPLIFY = FALSE
  )

  metacell_objs <- metacell_objs[!sapply(metacell_objs, is.null)]

  run_stats <- do.call(rbind, lapply(metacell_objs, function(obj) obj@misc$run_stats))
  run_stats[grouping_vars] <- do.call(rbind, strsplit(run_stats$name, "#"))

  combined_metacell_obj <- Reduce(function(x, y) merge(x, y), metacell_objs)

  Idents(combined_metacell_obj) <- combined_metacell_obj@meta.data[[identity_var]]

  if (!is.null(selected_cells)) {
    seurat_obj <- original_obj
  }

  seurat_obj <- AssignMetacellReference(seurat_obj, combined_metacell_obj, analysis_name)
  seurat_obj <- AssignAnalysisParameters(seurat_obj, params = list(
    metacell_neighbors = neighbors_k,
    metacell_reduction = dimension_reduction,
    metacell_slot = data_slot,
    metacell_assay = assay_type,
    metacell_run_stats = run_stats
  ), analysis_name)
  return(seurat_obj)
}


AssignMetacellReference <- function(seurat_obj, metacell_obj, analysis_name = NULL) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis
  seurat_obj@misc[[analysis_name]]$metacell_obj <- metacell_obj
  return(seurat_obj)
}


AssignAnalysisParameters <- function(seurat_obj, params, analysis_name = NULL) {
  analysis_name <- analysis_name %||% seurat_obj@misc$active_analysis
  if (is.null(seurat_obj@misc[[analysis_name]]$analysis_params)) {
    seurat_obj@misc[[analysis_name]]$analysis_params <- params
  } else {
    seurat_obj@misc[[analysis_name]]$analysis_params <- modifyList(seurat_obj@misc[[analysis_name]]$analysis_params, params)
  }
  return(seurat_obj)
}
