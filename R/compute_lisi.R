# FILE: R/compute_lisi.R

#' Compute LISI in Parallel
#'
#' Calculates the Local Inverse Simpson's Index (LISI) for multiple metadata
#' labels in parallel.
#'
#' @param X A numeric matrix of embeddings (e.g., PCs) for which to compute LISI.
#' @param meta_data A data frame with cell metadata, where rows correspond to rows of X.
#' @param label_colnames A character vector of column names in `meta_data` for which to compute LISI.
#' @param perplexity The perplexity parameter used to determine the neighborhood size.
#' @param nn_eps A parameter for the nearest neighbor search.
#' @param n_thread The number of parallel threads to use.
#' @return A data frame with LISI scores for each cell (rows) and each label (columns).
#'
#' @importFrom RANN nn2
#' @importFrom parallel mclapply
#' @importFrom lisi compute_simpson_index
#'
#' @export
compute_lisi_parallel <- function(
    X, meta_data, label_colnames, perplexity = 30, nn_eps = 0, n_thread = 1
) {

  N <- nrow(meta_data)
  dknn <- nn2(X, k = perplexity * 3, eps = nn_eps)
  lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
  lisi_df <- Reduce(cbind, mclapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message('Cannot compute LISI on missing values')      
      return(rep(NA, N))
    } else {
      ## don't count yourself in your neighborhood
      dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
      dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
      labels <- as.integer(factor(labels)) - 1
      n_batches <- length(unique(labels))
      simpson <- compute_simpson_index(
        t(dknn$nn.dists),
        t(dknn$nn.idx) - 1, 
        labels,
        n_batches,
        perplexity
      )
      return(1 / simpson)
    }
  }, mc.cores = n_thread))
  lisi_df <- as.data.frame(lisi_df)  
  colnames(lisi_df) <- label_colnames
  row.names(lisi_df) <- row.names(meta_data)
  return(lisi_df)
}

