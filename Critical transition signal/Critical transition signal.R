# ====================
# DATA PREPROCESSING
# ====================

# Factorize group variable with specified levels
cli$group <- factor(
  cli$group,
  levels = c("CCS", "UA", "NSTEMI", "STEMI", "PCI")
)

# Split sample IDs by group
sample_groups <- split(cli[, "ID"], f = cli$group)

# Subset expression matrix for each group (only MF genes)
group_expression <- lapply(sample_groups, function(samples) {
  df_sub <- df[, samples] %>% as.data.frame()
  df_sub[rownames(df_sub) %in% MF_genes, ]
})

# ====================
# NETWORK CONSTRUCTION
# ====================

#' Construct Correlation-Based Networks
#'
#' @param expression_list List of expression matrices per group
#' @param fdr_threshold False discovery rate cutoff for edges
#' @return List of igraph objects representing gene networks
construct_correlation_networks <- function(expression_list, fdr_threshold = 0.05) {
  # Validate input
  if (!all(sapply(expression_list, is.matrix) |
    !all(sapply(expression_list, is.data.frame)))) {
    stop("Input must be a list of matrices or data frames")
  }

  # Initialize result list
  network_list <- list()

  # Calculate correlation matrices
  correlation_matrices <- lapply(expression_list, function(expr) {
    corr.test(t(expr), adjust = "fdr", ci = FALSE)$r
  })

  # Calculate p-value matrices
  pvalue_matrices <- lapply(expression_list, function(expr) {
    corr.test(t(expr), adjust = "fdr", ci = FALSE)$p
  })

  # Process each group
  for (group_name in names(expression_list)) {
    cor_mat <- correlation_matrices[[group_name]]
    pval_mat <- pvalue_matrices[[group_name]]

    # Handle gene names with special characters
    rownames(cor_mat) <- gsub(".", "+", rownames(cor_mat), fixed = TRUE)
    rownames(pval_mat) <- gsub(".", "+", rownames(pval_mat), fixed = TRUE)

    # Prepare for edge extraction
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
    significant_edges <- list()

    # Identify significant correlations
    for (i in 1:nrow(cor_mat)) {
      sig_cols <- which(pval_mat[i, ] < fdr_threshold)
      significant_edges[[i]] <- cor_mat[i, sig_cols]
      names(significant_edges[[i]]) <- names(sig_cols)
    }
    names(significant_edges) <- rownames(cor_mat)

    # Create edge dataframe
    edges_df <- stack(do.call(c, significant_edges))
    edges_df <- edges_df[!is.na(edges_df$values), ]

    # Parse node pairs
    node_pairs <- str_split_fixed(edges_df$ind, "\\.", 2)
    edges_df$node1 <- node_pairs[, 1]
    edges_df$node2 <- node_pairs[, 2]

    # Revert gene name formatting
    edges_df$node1 <- gsub("+", ".", edges_df$node1, fixed = TRUE)
    edges_df$node2 <- gsub("+", ".", edges_df$node2, fixed = TRUE)

    # Final edge data
    edges_df <- edges_df[, c("node1", "node2", "values")]
    colnames(edges_df)[3] <- "correlation"
    edges_df$weight <- abs(edges_df$correlation)

    # Create nodes
    nodes_df <- data.frame(id = unique(c(edges_df$node1, edges_df$node2)))

    message(paste0(group_name, ": ", nrow(nodes_df), " nodes"))

    # Construct igraph object
    network_list[[group_name]] <- graph_from_data_frame(
      d = edges_df,
      vertices = nodes_df,
      directed = FALSE
    )
  }

  return(network_list)
}

# Build networks
gene_networks <- construct_correlation_networks(group_expression, fdr = 0.01)

# ====================
# NETWORK CLUSTERING
# ====================

#' Cluster Network Using Specified Method
#'
#' @param network_input List of networks or correlation matrices
#' @param method Clustering method ("rw", "hcm", "km", "pam", "natural")
#' @param cluster_cutoff Cluster size cutoff or number of clusters
#' @return List of cluster assignments
cluster_networks <- function(network_input,
                             method = c("rw", "hcm", "km", "pam", "natural"),
                             cluster_cutoff = NULL) {
  method <- match.arg(method)

  # Validate input based on method
  if (method == "rw" || method == "natural") {
    if (!all(sapply(network_input, inherits, "igraph"))) {
      stop("This method requires igraph objects")
    }
  } else {
    if (!all(sapply(network_input, is.matrix) |
      !all(sapply(network_input, is.data.frame)))) {
      stop("This method requires matrices or data frames")
    }
  }

  # Method-specific clustering
  if (method == "rw") {
    if (is.null(cluster_cutoff)) cluster_cutoff <- 4
    if (!is.integer(cluster_cutoff)) {
      warning("Using integer cutoff for random walk clustering")
      cluster_cutoff <- as.integer(cluster_cutoff)
    }
    cluster_results <- detect_communities(network_input, cluster_cutoff)
  } else if (method == "hcm") {
    if (is.null(cluster_cutoff)) {
      stop("Hierarchical clustering requires cluster cutoff")
    }
    correlation_matrices <- lapply(network_input, function(x) {
      corr.test(t(x), adjust = "fdr", ci = FALSE)$r
    })
    cluster_results <- perform_hierarchical_clustering(correlation_matrices, cluster_cutoff)
  } else if (method %in% c("km", "pam")) {
    if (is.null(cluster_cutoff)) {
      stop("This method requires cluster cutoff")
    }
    correlation_matrices <- lapply(network_input, function(x) {
      corr.test(t(x), adjust = "fdr", ci = FALSE)$r
    })
    cluster_results <- perform_partition_clustering(correlation_matrices, method, cluster_cutoff)
  } else if (method == "natural") {
    cluster_results <- lapply(network_input, function(g) {
      components(g)$membership
    })
  }

  return(cluster_results)
}

# Helper function for random walk clustering
detect_communities <- function(network_list, steps) {
  lapply(network_list, function(g) {
    cluster_walktrap(g, steps = steps)$membership
  })
}

# Helper function for hierarchical clustering
perform_hierarchical_clustering <- function(cor_matrices, k) {
  clusters <- lapply(cor_matrices, function(mat) {
    hc <- hclust(dist(mat), method = "complete")
    cutree(hc, k = k)
  })
  # Visualization
  par(mfrow = c(1, length(clusters)))
  lapply(clusters, function(hc) plot(hc))
  return(clusters)
}

# Helper function for partition-based clustering
perform_partition_clustering <- function(cor_matrices, method, k) {
  lapply(cor_matrices, function(mat) {
    if (method == "km") {
      kmeans(mat, centers = k)$cluster
    } else {
      pam(mat, k = k, metric = "euclidean")$clustering
    }
  })
}

# Cluster networks using random walk
network_clusters <- cluster_networks(
  gene_networks,
  method = "rw"
)

# ====================
# MODULE CRITICAL INDEX CALCULATION
# ====================

#' Calculate Module Critical Index (MCI)
#'
#' @param cluster_assignments List of cluster assignments
#' @param expression_data List of expression matrices
#' @param adjust_size Whether to adjust for module size
#' @param correlation_method Correlation calculation method ("cor", "BioTIP")
#' @param precomputed_corr Precomputed correlation matrix (optional)
#' @return List containing MCI and related statistics
calculate_module_critical_index <- function(cluster_assignments,
                                            expression_data,
                                            adjust_size = FALSE,
                                            correlation_method = c("cor", "BioTIP"),
                                            precomputed_corr = NULL) {
  correlation_method <- match.arg(correlation_method)
  results <- list(
    members = list(),
    MCI = list(),
    sd = list(),
    PCC = list(),
    PCCo = list()
  )

  # Validate input names
  if (is.null(names(cluster_assignments)) || is.null(names(expression_data))) {
    warning("Missing names in input lists")
  }

  # Process each group
  for (group_name in names(cluster_assignments)) {
    clusters <- cluster_assignments[[group_name]]
    expr_matrix <- expression_data[[group_name]]

    # Skip invalid groups
    if (all(is.na(clusters))) {
      results$MCI[[group_name]] <- NA
      results$sd[[group_name]] <- NA
      results$PCC[[group_name]] <- NA
      results$PCCo[[group_name]] <- NA
      next
    }

    # Initialize storage
    module_MCI <- numeric(max(clusters))
    module_sd <- list()
    module_pcc <- numeric(max(clusters))
    module_pcco <- numeric(max(clusters))

    # Process each cluster
    for (cluster_id in 1:max(clusters)) {
      # Extract module genes
      module_genes <- names(clusters[clusters == cluster_id])
      module_expr <- expr_matrix[rownames(expr_matrix) %in% module_genes, ]

      # Extract non-module genes
      non_module_expr <- expr_matrix[!rownames(expr_matrix) %in% module_genes, ]

      # Calculate correlations
      if (correlation_method == "cor") {
        # Intra-module correlation
        if (nrow(module_expr) > 1) {
          module_cor <- abs(cor(t(module_expr)))
          diag(module_cor) <- NA
          module_pcc[cluster_id] <- mean(module_cor, na.rm = TRUE)
        } else {
          module_pcc[cluster_id] <- NA
        }

        # Extra-module correlation
        if (nrow(module_expr) > 0 && nrow(non_module_expr) > 0) {
          cross_cor <- abs(cor(t(non_module_expr), t(module_expr)))
          module_pcco[cluster_id] <- mean(cross_cor, na.rm = TRUE)
        } else {
          module_pcco[cluster_id] <- NA
        }
      } else if (correlation_method == "BioTIP") {
        # Use precomputed correlations if available
        if (!is.null(precomputed_corr)) {
          # Intra-module correlation
          if (length(module_genes) > 1) {
            mod_cor <- precomputed_corr[module_genes, module_genes]
            module_pcc[cluster_id] <- mean(abs(mod_cor[upper.tri(mod_cor)]))
          }

          # Extra-module correlation
          if (length(module_genes) > 0 && nrow(non_module_expr) > 0) {
            cross_cor <- precomputed_corr[rownames(non_module_expr), module_genes]
            module_pcco[cluster_id] <- mean(abs(cross_cor))
          }
        } else {
          # Calculate correlations on the fly
          module_pcc[cluster_id] <- mean(
            abs(cor.shrink(module_expr, MARGIN = 1, target = 0)),
            na.rm = TRUE
          )
          module_pcco[cluster_id] <- mean(
            abs(cor.shrink(non_module_expr, module_expr, MARGIN = 1, target = 0)),
            na.rm = TRUE
          )
        }
      }

      # Calculate gene standard deviations
      gene_sds <- apply(module_expr, 1, sd, na.rm = TRUE)
      module_sd[[cluster_id]] <- gene_sds

      # Calculate MCI with optional size adjustment
      if (adjust_size) {
        module_MCI[cluster_id] <- mean(gene_sds) *
          (module_pcc[cluster_id] / module_pcco[cluster_id]) *
          sqrt(nrow(module_expr))
      } else {
        module_MCI[cluster_id] <- mean(gene_sds) *
          (module_pcc[cluster_id] / module_pcco[cluster_id])
      }
    }

    # Store results
    results$members[[group_name]] <- clusters
    results$MCI[[group_name]] <- module_MCI
    results$sd[[group_name]] <- module_sd
    results$PCC[[group_name]] <- module_pcc
    results$PCCo[[group_name]] <- module_pcco
  }

  return(results)
}

# Calculate MCI
mci_results <- calculate_module_critical_index(
  network_clusters,
  group_expression,
  correlation_method = "BioTIP"
)

# Plot MCI values
plot_module_critical_index(mci_results$MCI, ylim = c(0, 1))

# ====================
# TOP MODULE IDENTIFICATION
# ====================

#' Identify Top Critical Modules
#'
#' @param cluster_assignments List of cluster assignments
#' @param mci_values List of MCI values
#' @param min_size Minimum module size
#' @param top_n Number of top modules to return
#' @return List containing top module indices and members
identify_top_modules <- function(cluster_assignments, mci_values, min_size = 1, top_n = 1) {
  # Validate input
  if (top_n < 1) stop("top_n must be >= 1")
  if (min_size < 1) stop("min_size must be >= 1")

  results <- list(
    indices = list(),
    members = list()
  )

  # Process each group
  for (group_name in names(cluster_assignments)) {
    clusters <- cluster_assignments[[group_name]]
    mci <- mci_values[[group_name]]

    # Filter modules by size
    module_sizes <- table(clusters)
    valid_modules <- names(module_sizes)[module_sizes > min_size]

    if (length(valid_modules) == 0) {
      results$indices[[group_name]] <- NA
      results$members[[group_name]] <- NA
      next
    }

    # Select top modules
    top_indices <- order(mci[as.numeric(valid_modules)], decreasing = TRUE)[1:min(top_n, length(valid_modules))]
    top_module_ids <- as.numeric(valid_modules)[top_indices]

    # Extract module members
    module_members <- lapply(top_module_ids, function(mod_id) {
      names(clusters)[clusters == mod_id]
    })

    # Store results
    results$indices[[group_name]] <- top_module_ids
    results$members[[group_name]] <- module_members
  }

  return(results)
}

# Identify top modules
top_modules <- identify_top_modules(
  mci_results$members,
  mci_results$MCI,
  min_size = 10,
  top_n = 3
)

# Extract module statistics
extract_module_statistics <- function(stat_list, module_indices) {
  result <- list()
  for (group_name in names(module_indices)) {
    group_stats <- stat_list[[group_name]]
    group_indices <- module_indices[[group_name]]
    result[[group_name]] <- sapply(group_indices, function(idx) {
      if (is.list(group_stats)) {
        mean(group_stats[[idx]])
      } else {
        group_stats[idx]
      }
    })
  }
  return(result)
}

# Extract various statistics for top modules
max_mci <- extract_module_statistics(mci_results$MCI, top_modules$indices)
max_sd <- extract_module_statistics(mci_results$sd, top_modules$indices)
max_pcc <- extract_module_statistics(mci_results$PCC, top_modules$indices)
max_pcco <- extract_module_statistics(mci_results$PCCo, top_modules$indices)

# Combine into DNB score matrix
dnb_scores <- rbind(
  MCI = unlist(max_mci),
  SD = unlist(max_sd),
  PCC = unlist(max_pcc),
  PCCo = unlist(max_pcco)
) %>% as.data.frame()

# ====================
# CRITICAL TRANSITION STATE
# ====================

#' Identify Critical Transition State (CTS) Genes
#'
#' @param top_mci Named vector of top MCI values
#' @param top_modules List of top module members
#' @return List of CTS genes per state
identify_cts_genes <- function(top_mci, top_modules) {
  # Validate input
  if (is.null(names(top_mci))) stop("top_mci must be named")
  if (!all(names(top_mci) %in% names(top_modules))) {
    stop("Mismatched names between top_mci and top_modules")
  }

  cts_genes <- lapply(names(top_mci), function(state) {
    top_modules[[state]][[which.max(top_mci[[state]])]]
  })
  names(cts_genes) <- names(top_mci)

  return(cts_genes)
}

# Identify CTS genes
cts_genes <- identify_cts_genes(top_mci, top_modules$members)

# ====================
# SIMULATION AND VISUALIZATION
# ====================

#' Simulate MCI Distribution
#'
#' @param module_size Size of modules to simulate
#' @param sample_groups List of sample groupings
#' @param expression_matrix Full expression matrix
#' @param permutations Number of permutations
#' @param correlation_method Correlation calculation method
#' @return Matrix of simulated MCI values
simulate_mci_distribution <- function(module_size, sample_groups, expression_matrix,
                                      permutations = 100,
                                      correlation_method = c("cor", "BioTIP")) {
  correlation_method <- match.arg(correlation_method)
  precomputed_corr <- NULL

  # Precompute correlations if using BioTIP
  if (correlation_method == "BioTIP") {
    precomputed_corr <- cor.shrink(expression_matrix, MARGIN = 1, target = 0)
  }

  # Initialize result matrix
  simulated_mci <- matrix(
    nrow = length(sample_groups),
    ncol = permutations,
    dimnames = list(names(sample_groups), NULL)
  )

  # Perform permutations
  pb <- txtProgressBar(min = 0, max = permutations, style = 3)
  for (i in 1:permutations) {
    # Create random clusters
    random_clusters <- lapply(sample_groups, function(samples) {
      genes <- sample(rownames(expression_matrix), module_size)
      setNames(rep(1, length(genes)), genes)
    })

    # Calculate MCI for random modules
    random_mci <- calculate_module_critical_index(
      random_clusters,
      lapply(sample_groups, function(s) expression_matrix[, s]),
      correlation_method = correlation_method,
      precomputed_corr = precomputed_corr
    )$MCI

    # Store max MCI per group
    for (group_name in names(random_mci)) {
      if (all(is.na(random_mci[[group_name]]))) next
      simulated_mci[group_name, i] <- max(random_mci[[group_name]], na.rm = TRUE)
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(simulated_mci)
}

# Run simulation
set.seed(12345)
simulated_mci <- simulate_mci_distribution(
  module_size = 8,
  sample_groups = sample_groups,
  expression_matrix = df,
  permutations = 50,
  correlation_method = "BioTIP"
)

#' Visualize MCI Simulation Results
#'
#' @param observed_mci Named vector of observed MCI values
#' @param simulated_mci Matrix of simulated MCI values
#' @param state_order Order of states for plotting
#' @param y_axis_limits Custom y-axis limits
#' @param main_title Plot title
#' @param reference_state State for reference lines
plot_mci_simulation <- function(observed_mci, simulated_mci,
                                state_order = NULL, y_axis_limits = NULL,
                                main_title = "MCI Simulation Results",
                                reference_state = NULL) {
  # Set state order
  if (!is.null(state_order)) {
    if (!all(state_order %in% rownames(simulated_mci))) {
      stop("Invalid states in state_order")
    }
    simulated_mci <- simulated_mci[state_order, ]
  }

  # Prepare plot data
  plot_data <- as.data.frame(t(simulated_mci))

  # Set y-axis limits
  if (is.null(y_axis_limits)) {
    y_axis_limits <- range(c(observed_mci, simulated_mci), na.rm = TRUE)
  }

  # Create boxplot
  boxplot(plot_data,
    col = "lightgray", ylab = "Module Critical Index (MCI)",
    main = main_title, ylim = y_axis_limits, outline = FALSE
  )

  # Add observed values
  for (state in names(observed_mci)) {
    state_idx <- which(colnames(plot_data) == state)
    points(state_idx, observed_mci[state], col = "red", pch = 19, cex = 1.5)
  }

  # Add reference lines
  if (!is.null(reference_state)) {
    ref_values <- simulated_mci[reference_state, ]
    abline(h = min(ref_values), lty = 2, col = "gray")
    abline(h = max(ref_values), lty = 2, col = "gray")
    abline(h = mean(ref_values) + 2 * sd(ref_values), lty = 3, col = "darkgray")
  }

  # Add significance indicators
  for (state in names(observed_mci)) {
    state_idx <- which(colnames(plot_data) == state)
    p_value <- mean(simulated_mci[state, ] > observed_mci[state])
    if (p_value < 0.05) {
      text(state_idx, max(simulated_mci[state, ]) * 0.95, "*", cex = 2)
    }
  }
}

# Visualize results
plot_mci_simulation(
  observed_mci = unlist(max_mci),
  simulated_mci = simulated_mci,
  state_order = c("CCS", "UA", "NSTEMI", "STEMI", "PCI"),
  main_title = "MCI Simulation Results"
)
