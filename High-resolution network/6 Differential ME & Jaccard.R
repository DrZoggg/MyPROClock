# Identify Differential Module Eigengenes (DMEs)
DMEs <- IdentifyDifferentialModules(
  seurat_obj,
  group_1_barcodes = group1,
  group_2_barcodes = group2,
  statistical_test = "wilcox",
  analysis_label = "SC"
)

head(DMEs)

# Generate volcano plot for Differential Module Eigengenes
PlotDifferentialModulesVolcano(
  seurat_obj,
  DMEs,
  analysis_label = "SC"
)


#--------------------------
# Custom Function Definitions
#--------------------------

IdentifyDifferentialModules <- function(seurat_obj, group_1_barcodes, group_2_barcodes, harmonized = TRUE,
                                        analysis_label = NULL, include_missing = FALSE, statistical_test = "wilcox",
                                        positive_only = FALSE, logfc_threshold = 0, min_pct = 0, verbose = FALSE,
                                        pseudocount = 0, ...) {
  analysis_label <- analysis_label %||% seurat_obj@misc$active_analysis

  # Barcode validation
  if (!all(group_1_barcodes %in% colnames(seurat_obj))) {
    stop("Some barcodes in group_1_barcodes are not present in seurat_obj.")
  }
  if (!all(group_2_barcodes %in% colnames(seurat_obj))) {
    stop("Some barcodes in group_2_barcodes are not present in seurat_obj.")
  }
  if (length(intersect(group_1_barcodes, group_2_barcodes)) > 0) {
    stop("Overlapping barcodes detected between group_1_barcodes and group_2_barcodes.")
  }

  # Retrieve module eigengenes
  eigengenes <- GetModuleEigengenes(seurat_obj, harmonized, analysis_label)
  eigengenes <- eigengenes[, colnames(eigengenes) != "grey"]

  eigengenes[eigengenes < 0] <- 0
  eigengenes_transposed <- t(eigengenes)

  ME_assay <- Seurat::CreateAssayObject(eigengenes_transposed)
  differential_modules <- FindMarkers(ME_assay,
    cells.1 = group_1_barcodes, cells.2 = group_2_barcodes,
    slot = "counts", test.use = statistical_test, only.pos = positive_only,
    logfc.threshold = logfc_threshold, min.pct = min_pct,
    verbose = verbose, pseudocount.use = pseudocount,
    ...
  )

  differential_modules$module <- rownames(differential_modules)

  if (include_missing) {
    missing_modules <- setdiff(rownames(eigengenes), differential_modules$module)
    for (mod in missing_modules) {
      differential_modules[mod, ] <- NA
      differential_modules[mod, "module"] <- mod
    }
  }

  return(differential_modules)
}

PlotDifferentialModulesVolcano <- function(seurat_obj, differential_modules, show_labels = TRUE, module_point_size = 4,
                                           label_text_size = 4, display_cutoff = TRUE, analysis_label = NULL) {
  differential_modules <- na.omit(differential_modules)
  adjusted_pval_min <- min(differential_modules$p_val_adj[differential_modules$p_val_adj != 0])

  differential_modules$p_val_adj[differential_modules$p_val_adj == 0] <- adjusted_pval_min
  finite_fc <- differential_modules$avg_log2FC[is.finite(differential_modules$avg_log2FC)]
  max_fc <- max(abs(finite_fc))

  differential_modules$avg_log2FC[is.infinite(differential_modules$avg_log2FC) & differential_modules$avg_log2FC < 0] <- -max_fc
  differential_modules$avg_log2FC[is.infinite(differential_modules$avg_log2FC) & differential_modules$avg_log2FC > 0] <- max_fc

  module_colors <- seurat_obj@misc$module_assignments %>%
    dplyr::filter(module != "grey") %>%
    dplyr::select(module, color) %>%
    distinct()

  colors <- setNames(module_colors$color, module_colors$module)

  differential_modules$label <- ifelse(differential_modules$p_val_adj < 0.05, differential_modules$module, "")

  plot <- ggplot(differential_modules, aes(x = avg_log2FC, y = -log10(p_val_adj), fill = module, color = module))

  if (display_cutoff) {
    plot <- plot + geom_vline(xintercept = 0, linetype = "dashed", color = "grey75", alpha = 0.8) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)), fill = "grey75", alpha = 0.8)
  }

  plot <- plot + geom_point(size = module_point_size, shape = 21, color = "black")

  if (show_labels) {
    plot <- plot + ggrepel::geom_text_repel(aes(label = label),
      color = "black",
      min.segment.length = 0, max.overlaps = Inf, size = label_text_size
    )
  }

  plot <- plot + scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    xlim(c(-max_fc - 0.1, max_fc + 0.1)) +
    labs(x = expression("Average log"[2] ~ "(Fold Change)"), y = expression("-log"[10] ~ "(Adj. P-value)")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(), axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5), legend.position = "none"
    )

  return(plot)
}
