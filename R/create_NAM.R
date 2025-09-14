# FILE: R/create_NAM.R

#' Create a Neighborhood Abundance Matrix (NAM) and Decompose it
#'
#' This function implements the core logic of the CNA method. It builds a
#' neighborhood graph, diffuses sample information across it to create a NAM,
#' performs quality control, residualizes covariates, and finally decomposes
#' the NAM using SVD to produce NAM PCs.
#'
#' @param metadata A data frame with cell-level metadata.
#' @param pcs A numeric matrix of principal components or other embeddings (cells x PCs).
#' @param samplem_key The column name in `metadata` that identifies the sample/individual.
#' @param graph_use The name of the graph in the Seurat object to use.
#' @param batches A character vector of batch variable names in `metadata`.
#' @param covs A character vector of other covariate names in `metadata`.
#' @param nsteps The number of diffusion steps. If NULL, it's determined automatically.
#' @param verbose Logical, if TRUE, print progress messages.
#' @param assay The assay name to use when creating the Seurat object.
#' @param key A prefix for the resulting DimReduc object.
#'
#' @return A list containing the SVD of the NAM: `sampleXpc` (sample loadings),
#'   `svs` (squared singular values), `nbhdXpc` (neighborhood loadings), and
#'   `varexp` (variance explained).
#'
#' @import Seurat
#' @importFrom Matrix colSums Diagonal t
#' @importFrom tibble remove_rownames rownames_to_column
#' @importFrom dplyr select one_of
#' @importFrom glue glue
#' @importFrom moments kurtosis
#' @importFrom stats model.matrix as.formula sd svd
#' @importFrom purrr imap map map_dbl reduce
#' @importFrom RSpectra svds
#'
#' @export
create_NAM = function(metadata,
                      pcs,
                      samplem_key = NULL,
                      graph_use = 'RNA_snn',
                      batches = NULL,
                      covs = NULL ,
                      nsteps = NULL,
                      verbose=TRUE ,
                      assay=NULL,
                      key='NAMPC_'
){
  
  meta = metadata
  rownames(meta) = 1:nrow(meta)
  
  m <- as(t(pcs), "dgTMatrix") # by default, Matrix() returns dgCMatrix
  colnames(m) = 1:ncol(m)
  
  obj <- Seurat::CreateSeuratObject(
    counts = m, ## Subset expression matrix to cells in metadata
    meta.data = meta, 
    assay = 'RNA', 
    names.field = 1
  )
  
  harmony_embeddings_all = pcs
  rownames(harmony_embeddings_all) = 1:nrow(harmony_embeddings_all)
  
  obj@reductions$harmony = Seurat::CreateDimReducObject(
    embeddings = harmony_embeddings_all,
    stdev = as.numeric(apply(harmony_embeddings_all, 2, stats::sd)),
    assay = "RNA",
    key = Seurat::Key("HARMONY", quiet = TRUE)
  )
  
  obj <- obj %>%
    Seurat::FindNeighbors(verbose = TRUE, reduction = 'harmony', dims = 1:20, k.param = 30, nn.eps = 0)
  
  seurat_object = obj
  
  ## (1) format data 
  covs_keep <- NULL
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs)) covs_keep <- c(covs_keep, covs)
  
  if (length(names(seurat_object@graphs)) == 0) {
    stop('Must precompute graph in Seurat with FindNeighbors()')
  }
  if (is.null(graph_use)) {
    graph_use <- names(seurat_object@graphs)[[1]]
    message('Graph not specified. Using graph {graph_use}')
  } else {
    if (!graph_use %in% names(seurat_object@graphs)) {
      stop('Graph {graph_use} not in seurat object')
    }
  }
  covs_keep <- c(covs_keep, samplem_key)
  samplem_df <- tibble::remove_rownames(unique(dplyr::select(seurat_object@meta.data, one_of(covs_keep))))
  obs_df <- tibble::rownames_to_column(seurat_object@meta.data, 'CellID')
  if (nrow(samplem_df) == nrow(obs_df)) {
    stop(
      'Sample-level metadata is same length as cell-level metadata.       
             Please check that samplem_vars are sample-level covariates.'
    )
  }
  
  rcna_data <- list(
    samplem = samplem_df,
    obs = obs_df, 
    connectivities = seurat_object@graphs[[graph_use]],
    samplem_key = samplem_key,
    obs_key = 'CellID',
    N = nrow(samplem_df)
  )
  data = rcna_data
  suffix=''
  
  # formatting and error checking
  ## For association, batches needs to be a numeric vector
  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(as.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }
  max_frac_pcs=0.15
  res <- list()
  covs_mat <- as.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))
  
  f <- as.formula(as.character(glue('~0+{data$samplem_key}')))    
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]    
  s <- s[, data$samplem[[data$samplem_key]]] ## Necessary? 
  # s = Matrix indicating which sample the cells that represent a row are from
  ## row: cell
  ## column: individual
  ## 1=the cell from that sample
  ### prior knowledge of connectivity among cells
  
  prevmedkurt <- Inf
  maxnsteps=15L
  
  # NOTE: this function ignores distances and only uses unweighted connectivities
  # This function defines a process to diffuse (or propagate) information across a graph. 
  # This type of diffusion step is often used in tasks like clustering or dimensionality reduction for network data.
  diffuse_step <- function(data, s) {
    # 1. Get the connectivity matrix 'a' from `data$connectivities`. This matrix represents the graph's connectivity, 
    #    where rows and columns are vertices (cells), and values are the connection strengths (weights).
    a <- data$connectivities  
    # 2. Calculate the degree of each vertex (number of connected edges). The degree is a measure of a vertex's 'importance'. 
    #    Here, 1 is added to prevent division by zero for isolated nodes.
    degrees <- Matrix::colSums(a) + 1 
    # 3. Normalize the input vector 's' by the degree of each vertex.
    s_norm <- s / degrees  
    # 4. Multiply the connectivity matrix 'a' by the normalized vector 's_norm' and add the original 's_norm' back. 
    #    This updates each vertex's value with the sum of its neighbors' normalized values, simulating information 'diffusion'.
    res <- (a %*% s_norm) + s_norm  
    # 5. Finally, return the result as a dense matrix.
    return(as.matrix(res)) 
  }
  
  # This loop iterates through diffusion steps, evaluates the result, and stops when a condition is met.
  # It provides a mechanism to stop diffusion when it has reached an appropriate stage, checked by monitoring the median kurtosis.
  for (i in seq_len(maxnsteps)) { 
    s <- diffuse_step(data, s) 
    # Calculate the median kurtosis. Kurtosis measures the "tailedness" of the distribution.
    medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis)) 
    if (is.null(nsteps)) { 
      prevmedkurt <- medkurt
      # Stop if the change in kurtosis is small and at least 3 steps have been performed.
      if (prevmedkurt - medkurt < 3 & i > 3) { 
        message(glue::glue('stopping after {i} steps'))
        break 
      }            
    } else if (i == nsteps) {
      break
    }
  }  
  
  snorm <- t(prop.table(s, 2)) # normalization
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]
  
  NAM=snorm
  
  N <- nrow(NAM)
  ## NOTE: added NULL check 
  if (is.null(batches_vec) | length(unique(batches_vec)) == 1) {
    message('only one unique batch supplied to qc')
    keep <- rep(TRUE, ncol(NAM))
  } else {
    message('filtering based on batches kurtosis')  
  }
  
  # This function calculates the kurtosis of values for each batch.
  # In essence, it computes the kurtosis of the mean values of rows within each batch in the NAM matrix. 
  # This can be used to analyze batch effects.
  .batch_kurtosis <- function(NAM, batches_vec) {
    purrr::imap(split(seq_len(length(batches_vec)), batches_vec), function(i, b) {
      if (length(i)>1) {
        Matrix::colMeans(NAM[i, ])
      } else if (length(i)==1) {
        Matrix::colMeans(t(NAM[i, ]))
      }
    }) %>% 
      dplyr::bind_cols() %>% 
      apply(1, moments::kurtosis)
  }
  
  # This section filters data based on the kurtosis calculated for each batch.
  kurtoses <- .batch_kurtosis(NAM, batches_vec) 
  # Set a threshold for filtering. It's the greater of 6 or twice the median kurtosis.
  # This dynamic threshold helps ensure robustness against outliers and noise.
  threshold <- max(6, 2*median(kurtoses)) 
  message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}')) 
  # Get indices of batches with kurtosis below the threshold.
  keep <- which(kurtoses < threshold) 
  
  .res_qc_nam = list(NAM = NAM, keep = keep)
  
  res <- list()
  res[[paste0('NAM.T', suffix)]] <- t(.res_qc_nam[[1]])
  res[[paste0('keptcells', suffix)]] <- .res_qc_nam[[2]]
  res[[paste0('_batches', suffix)]] <- batches_vec
  
  # "Residualize" the NAM matrix, meaning to remove the effects of covariates.
  if (verbose) message('Residualize NAM')
  N <- nrow(NAM)
  # Center the NAM matrix (subtract the mean) but do not scale.
  NAM_ <- scale(NAM, center = TRUE, scale = FALSE) 
  ncols_C <- 0
  if (!is.null(covs_mat)) { 
    covs_mat <- scale(covs_mat)
    ncols_C <- ncols_C + ncol(covs_mat)
  }
  
  # This section creates a "residualizing matrix" M. Multiplying the data by M is equivalent to 
  # performing a linear regression against the covariates and taking the residuals.
  if (is.null(covs_mat)) {
    M <- Matrix::Diagonal(n = N) 
  } else {
    # M = I - X(X'X)^-1 X' , where X is the covariate matrix.
    M <- Matrix::Diagonal(n = N) - covs_mat %*% solve(t(covs_mat) %*% covs_mat, t(covs_mat))  
  }
  
  # Apply the residualizing matrix to NAM.
  NAM_ <- M %*% NAM_  
  
  .res_resid_nam = list(
    NAM_=scale(NAM_, center=FALSE, scale=TRUE), 
    M=M, 
    r=ncols_C
  )
  res[[paste0('_M', suffix)]] <- .res_resid_nam$M
  res[[paste0('_r', suffix)]] <- .res_resid_nam$r
  
  if (verbose) message('Decompose NAM')
  npcs <- max(10, round(max_frac_pcs * nrow(data$samplem)))
  
  npcs <- min(npcs, nrow(data$samplem) - 1) ## make sure you don't compute all SVs
  if (missing(npcs) | npcs > .5 * min(dim(NAM_))) {
    svd_res <- svd(scale(NAM_, center=FALSE, scale=TRUE))
  } else {
    svd_res <- RSpectra::svds(scale(NAM_, center=FALSE, scale=TRUE), k = npcs)
  }
  
  # SVD decomposition: A = U D V^T
  # d: singular values
  # u, v: orthogonal matrices
  
  U_df <- svd_res$u[, seq_len(npcs)]
  colnames(U_df) <- paste0('PC', seq_len(npcs))
  rownames(U_df) <- rownames(NAM_)
  V_df <- svd_res$v[, seq_len(npcs)]
  colnames(V_df) <- paste0('PC', seq_len(npcs))
  rownames(V_df) <- colnames(NAM_)
  .res_svd_nam <- list(U=U_df, svs=svd_res$d^2, V=V_df)
  
  res[[paste0('NAM_sampleXpc', suffix)]] <- .res_svd_nam$U
  res[[paste0('NAM_svs', suffix)]] <- .res_svd_nam$svs
  res[[paste0('NAM_varexp', suffix)]] <- .res_svd_nam$svs / nrow(.res_svd_nam$U) / nrow(.res_svd_nam$V)
  res[[paste0('NAM_nbhdXpc', suffix)]] <- .res_svd_nam$V
  
  nam_res = res
  
  NAMsvd=list(
    nam_res$NAM_sampleXpc,
    nam_res$NAM_svs,
    nam_res$NAM_nbhdXpc,
    nam_res$NAM_varexp
  )
  names(NAMsvd) = c("sampleXpc","svs","nbhdXpc","varexp")
  return(NAMsvd)
}