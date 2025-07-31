# Example demonstrating how to run Module UMAP using 50 hub genes per module:
set.seed(12345)
seurat_obj9 <- GenerateModuleUMAP(
  seurat_obj,
  hub_gene_count = 50, # number of hub genes per module
  neighbors = 40, # UMAP neighbors parameter
  min_distance = 0.3 # minimum distance parameter for UMAP
)

# Plotting with ggplot
umap_results <- RetrieveModuleUMAP(seurat_obj)
table(umap_results$module)
table(umap_results$hub, umap_results$module)

# Initial basic ggplot
library(ggplot2)
ggplot(umap_results, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    color = umap_results$color,
    size = umap_results$kME * 2
  ) +
  umap_theme()

# Prepare data frame for plotting with annotations
plot_df <- umap_results

# Calculate cluster centroids
centroids <- data.frame()
for (current_module in unique(plot_df$module)) {
  module_data <- plot_df[plot_df$module == current_module, ]
  centroid <- data.frame(
    module = current_module,
    UMAP1 = mean(module_data$UMAP1),
    UMAP2 = mean(module_data$UMAP2)
  )
  centroids <- rbind(centroids, centroid)
}

# Retrieve top 3 hub genes per module for annotation
hub_genes <- RetrieveHubGenes(seurat_obj, top_n = 3)
annotation_genes <- hub_genes$gene_name
plot_df$label <- ifelse(plot_df$gene %in% annotation_genes, plot_df$gene, "")
annotated_plot_df <- subset(plot_df, label != "")

# Enhanced plot with rasterized points and annotations
library(ggrepel)
library(ggrastr)
p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = module)) +
  rasterise(
    geom_point(
      inherit.aes = FALSE,
      data = plot_df,
      aes(x = UMAP1, y = UMAP2),
      color = plot_df$color,
      size = plot_df$kME * 4
    ),
    dpi = 500, dpi_scale = 0.5
  ) +
  geom_point(
    inherit.aes = FALSE,
    data = annotated_plot_df,
    shape = 21, color = "black",
    fill = annotated_plot_df$color,
    size = annotated_plot_df$kME * 2,
    aes(x = UMAP1, y = UMAP2)
  ) +
  geom_text_repel(
    data = centroids,
    aes(label = module),
    color = "black", max.overlaps = Inf, size = 3, fontface = "bold"
  ) +
  geom_text_repel(
    aes(label = label),
    max.overlaps = Inf,
    color = "black", fontface = "italic", size = 2
  ) +
  umap_theme() +
  NoLegend() +
  coord_equal() +
  theme(plot.margin = margin(0, 0, 0, 0))

#--------------------------
# Custom Function Definitions
#--------------------------

GenerateModuleUMAP <- function(seurat_obj, feature_set = "TOM", hub_gene_count = 10, exclude_grey = TRUE,
                               analysis_label = NULL, neighbors = 25, metric_type = "cosine", spread_factor = 1,
                               min_distance = 0.4, supervised = FALSE, ...) {
  analysis_label <- analysis_label %||% seurat_obj@misc$active_analysis

  modules <- RetrieveModuleAssignments(seurat_obj, analysis_label)
  module_names <- unique(modules$module)

  if (!all(paste0("kME_", module_names) %in% colnames(modules))) {
    stop("kME values missing; ensure ModuleEigengenes and ModuleConnectivity have been run.")
  }

  if (exclude_grey) {
    module_names <- module_names[module_names != "grey"]
    modules <- subset(modules, module != "grey")
  }

  hub_genes <- lapply(module_names, function(mod) {
    mod_genes <- subset(modules, module == mod)
    mod_genes[order(-mod_genes[[paste0("kME_", mod)]]), ][1:hub_gene_count, "gene_name"]
  })

  names(hub_genes) <- module_names
  all_selected_genes <- modules$gene_name[modules$module %in% module_names]

  TOM_matrix <- GetTOM(seurat_obj, analysis_label)
  feature_matrix <- TOM_matrix[all_selected_genes, unlist(hub_genes)]

  if (supervised) {
    umap_result <- uwot::umap(
      X = feature_matrix,
      min_dist = min_distance,
      n_neighbors = neighbors,
      metric = metric_type,
      spread = spread_factor,
      y = modules$module,
      ...
    )
  } else {
    umap_result <- uwot::umap(
      X = feature_matrix,
      min_dist = min_distance,
      n_neighbors = neighbors,
      metric = metric_type,
      spread = spread_factor,
      ...
    )
  }

  plot_df <- as.data.frame(umap_result)
  colnames(plot_df) <- c("UMAP1", "UMAP2")
  plot_df$gene <- rownames(feature_matrix)
  plot_df$module <- modules$module[match(plot_df$gene, modules$gene_name)]
  plot_df$color <- modules$color[match(plot_df$gene, modules$gene_name)]
  plot_df$hub <- ifelse(plot_df$gene %in% unlist(hub_genes), "hub", "other")

  kME_df <- do.call(rbind, lapply(module_names, function(mod) {
    mod_subset <- subset(modules, module == mod)
    kME_data <- data.frame(
      gene_name = mod_subset$gene_name,
      kME = scale01(mod_subset[[paste0("kME_", mod)]])
    )
    return(kME_data)
  }))

  plot_df$kME <- kME_df$kME[match(plot_df$gene, kME_df$gene_name)]
  seurat_obj@misc[[analysis_label]]$module_umap <- plot_df

  return(seurat_obj)
}

RetrieveHubGenes <- function(seurat_obj, top_n = 10, modules = NULL, analysis_label = NULL) {
  analysis_label <- analysis_label %||% seurat_obj@misc$active_analysis
  module_data <- subset(RetrieveModuleAssignments(seurat_obj, analysis_label), module != "grey")

  if (is.null(modules)) {
    modules <- unique(module_data$module)
  }

  hub_df <- do.call(rbind, lapply(modules, function(mod) {
    mod_subset <- subset(module_data, module == mod)
    mod_subset <- mod_subset[order(-mod_subset[[paste0("kME_", mod)]]), ][1:top_n, ]
    mod_subset
  }))

  rownames(hub_df) <- NULL
  return(hub_df)
}
