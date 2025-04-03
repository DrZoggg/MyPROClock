

get_random_Clock <- function(input.list, 
                             genes.dist.bins, 
                             b.sign, 
                             num.rounds = 1000
                             ) {
  sign.bins <- as.matrix(table(genes.dist.bins[b.sign]))
  q <- rownames(sign.bins)
  bg <- matrix(data = F, nrow = length(genes.dist.bins), ncol = num.rounds)
  for (i in 1:nrow(sign.bins)) {
    num.genes <- sign.bins[i]
    if (num.genes > 0) {
      idx <- which(is.element(genes.dist.bins, q[i]))
      for (j in 1:num.rounds) {
        idxj <- sample(idx, num.genes)
        bg[idxj, j] <- T
      }
    }
  }
  r.scores <- apply(bg, 2, function(x) colMeans(input.list$expr.scaled[x, ]))
  r.scores <- rowMeans(r.scores)
  return(r.scores)
}



cal_MyPROClock <- function(expr, 
                            circadians = NULL, 
                            study.type = NULL,
                            discretization.method = NULL,
                            num.rounds = 1000, 
                            n.bins = 50, 
                            seed = TRUE) {
  options(warn = -1)
  if (is.null(study.type)) {
    stop("please set your study types as 'scRNAseq' or 'bulk_RNAseq'")
  }

  # create input data
  input.list <- list(expr = as.matrix(expr), genes = rownames(expr))

  # set study types
  # scRNAseq
  if (study.type == "scRNAseq") {
    input.list$genes.mean <- rowMeans(input.list$expr)
    input.list$expr.scaled <- sweep(input.list$expr, 1, input.list$genes.mean) # scaled expr by mean(genes)

    dist <- 10 * ((2^input.list$expr) - 1)
    input.list$genes.dist <- log2(rowMeans(dist, na.rm = T) + 1)
  }

  # bulk_RNAseq
  if (study.type == "bulk_RNAseq") {
    input.list$genes.dist <- rowMeans(input.list$expr)
    input.list$expr.scaled <- sweep(input.list$expr, 1, input.list$genes.dist) # scaled expr by mean(genes)
  }

  # set bins, n.bins = 50
  if (is.null(discretization.method)) {
    input.list$genes.dist.bins <- arules::discretize(input.list$genes.dist,
      method = "frequency",
      breaks = n.bins
    ) %>%
      plyr::mapvalues(levels(.), 1:length(levels(.))) %>%
      as.numeric() %>%
      as.matrix()
  }

  if (!is.null(discretization.method)) {
    input.list$genes.dist.bins <- arules::discretize(input.list$genes.dist,
      method = discretization.method,
      breaks = n.bins
    ) %>%
      plyr::mapvalues(levels(.), 1:length(levels(.))) %>%
      as.numeric() %>%
      as.matrix()
  }

  # calculate MyPROClock
  b.sign <- is.element(input.list$genes, circadians)

  if (sum(b.sign) > 5) {
    if (seed) {
      set.seed(123)
      r.scores <- get_random_Clock(input.list, input.list$genes.dist.bins, b.sign, num.rounds = num.rounds)
    } else {
      r.scores <- get_random_Clock(input.list, input.list$genes.dist.bins, b.sign, num.rounds = num.rounds)
    }

    raw.scores <- colMeans(input.list$expr.scaled[b.sign, ])
    input.list$MyPROClock <- r.scores - raw.scores

    
  } else {
    stop("Non-enough-overlapping genes to calculate MyPROClock")
  }
  return(input.list$MyPROClock)
}
