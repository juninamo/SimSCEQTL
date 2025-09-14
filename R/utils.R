# FILE: R/utils.R

#' Normalize data using Seurat's log-transform method
#'
#' A standalone implementation of Seurat's LogNormalize method.
#' @param A A sparse matrix of counts.
#' @param scaling_factor A scaling factor, defaults to 10000.
#' @param do_ftt Logical, if TRUE, applies a Freeman-Tukey stabilizing transformation.
#' @return A log-normalized sparse matrix.
#' @importFrom Matrix colSums
#' @export
NormalizeDataSeurat <- function(A, scaling_factor = 1e4, do_ftt = FALSE) {
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  A@x <- scaling_factor * A@x
  if (do_ftt) {
    A@x <- sqrt(A@x) + sqrt(1 + A@x)
  } else {
    A@x <- log(1 + A@x)
  }
  return(A)
}

#' Cosine Normalize a Matrix
#'
#' Normalizes vectors in a matrix to unit length (cosine normalization).
#' @param x A numeric matrix.
#' @param dim The margin to apply normalization over (1 for rows, 2 for columns).
#' @return A cosine-normalized matrix.
#' @export
cosine_normalize = function(x, dim){ apply(x, MARGIN = dim, function(a) a / sqrt(sum(a^2))) }

#' Unregister a doParallel Backend
#'
#' A helper function to clean up and unregister a parallel backend created by
#' the `foreach` package.
#' @return Invisible
#' @importFrom foreach .foreachGlobals
#' @export
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

