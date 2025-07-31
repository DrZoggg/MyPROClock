




#Advanced randomization function for computing randomized Clock scores
compute_randomized_clock <- function(data_input,
                                     gene_bins,
                                     significant_genes,
                                     permutations = 1000) {
  # Calculate the frequency of significant genes in each bin
  significant_bins <- as.matrix(table(gene_bins[significant_genes]))
  bin_levels <- rownames(significant_bins)
  
  # Initialize matrix for randomized background generation
  randomized_background <- matrix(FALSE,
                                  nrow = length(gene_bins),
                                  ncol = permutations)
  
  # Perform bin-specific random sampling
  for (bin_idx in seq_len(nrow(significant_bins))) {
    genes_per_bin <- significant_bins[bin_idx]
    if (genes_per_bin > 0) {
      bin_indices <- which(gene_bins %in% bin_levels[bin_idx])
      for (perm in seq_len(permutations)) {
        sampled_indices <- sample(bin_indices, genes_per_bin)
        randomized_background[sampled_indices, perm] <- TRUE
      }
    }
  }
  
  # Compute randomized scores by averaging selected genes
  randomized_scores <- apply(randomized_background, 2, function(selection) {
    selected_expression <- data_input$scaled_expression[selection, , drop = FALSE]
    colMeans(selected_expression)
  })
  
  # Average across permutations to obtain final randomized scores
  averaged_randomized_scores <- rowMeans(randomized_scores)
  
  return(averaged_randomized_scores)
}






# Comprehensive calculation function for the MyPROClock metric
calculate_MyPROClock <- function(expression_data,
                                 circadian_genes = NULL,
                                 analysis_mode = c("scRNAseq", "bulk_RNAseq"),
                                 binning_method = "frequency",
                                 permutations = 1000,
                                 bins = 50,
                                 reproducible_seed = TRUE) {
  
  analysis_mode <- match.arg(analysis_mode)
  options(warn = -1)
  
  # Validate input parameters
  if (!analysis_mode %in% c("scRNAseq", "bulk_RNAseq")) {
    stop("Invalid analysis_mode provided. Must be either 'scRNAseq' or 'bulk_RNAseq'.")
  }
  
  # Prepare input data structure
  data_input <- list(
    expression_matrix = as.matrix(expression_data),
    gene_names = rownames(expression_data)
  )
  
  # Data preprocessing and scaling based on analysis type
  if (analysis_mode == "scRNAseq") {
    gene_means <- rowMeans(data_input$expression_matrix)
    data_input$scaled_expression <- sweep(data_input$expression_matrix,
                                          MARGIN = 1,
                                          STATS = gene_means,
                                          FUN = "-")
    
    adjusted_expression <- 10 * (2 ^ data_input$expression_matrix - 1)
    data_input$gene_distribution <- log2(rowMeans(adjusted_expression, na.rm = TRUE) + 1)
    
  } else if (analysis_mode == "bulk_RNAseq") {
    data_input$gene_distribution <- rowMeans(data_input$expression_matrix)
    data_input$scaled_expression <- sweep(data_input$expression_matrix,
                                          MARGIN = 1,
                                          STATS = data_input$gene_distribution,
                                          FUN = "-")
  }
  
  # Discretize gene expression distribution into specified bins
  discretized_bins <- arules::discretize(data_input$gene_distribution,
                                         method = binning_method,
                                         breaks = bins)
  
  numeric_bins <- plyr::mapvalues(discretized_bins,
                                  from = levels(discretized_bins),
                                  to = seq_along(levels(discretized_bins)))
  
  data_input$gene_bins <- as.numeric(as.character(numeric_bins))
  
  # Determine which genes are circadian-related
  circadian_flags <- data_input$gene_names %in% circadian_genes
  
  if (sum(circadian_flags) <= 5) {
    stop("Insufficient circadian genes identified to reliably compute MyPROClock.")
  }
  
  # Set reproducible seed for consistent randomization
  if (reproducible_seed) {
    set.seed(123)
  }
  
  # Generate randomized background and calculate MyPROClock
  randomized_scores <- compute_randomized_clock(
    data_input = data_input,
    gene_bins = data_input$gene_bins,
    significant_genes = circadian_flags,
    permutations = permutations
  )
  
  raw_circadian_scores <- colMeans(data_input$scaled_expression[circadian_flags, , drop = FALSE])
  
  data_input$MyPROClock <- randomized_scores - raw_circadian_scores
  
  return(data_input$MyPROClock)
}



