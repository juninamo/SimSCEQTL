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

#' Conditional Permutation
#'
#' Generates permutations of a vector, conditioned on a grouping factor.
#' The permutation is performed within each group defined by `B`.
#' @param B A vector defining the groups.
#' @param Y The numeric vector to be permuted.
#' @param num The number of permutations to generate.
#' @return A matrix with `num` columns, each representing a conditional permutation of `Y`.
#' @importFrom purrr map reduce
#' @importFrom dplyr bind_rows arrange
#' @importFrom Matrix cbind2
#' @export
conditional_permutation <- function(B, Y, num) {
  purrr::map(seq_len(num), function(i) {
    split(seq_len(length(Y)), B) %>% purrr::map(function(idx) {
      data.frame(idx, val=sample(Y[idx]))
    }) %>% dplyr::bind_rows() %>% 
      dplyr::arrange(idx) %>% 
      with(val)   
  }) %>% 
    purrr::reduce(Matrix::cbind2)
}

#' Calculate Empirical False Discovery Rates
#'
#' Computes empirical FDRs by comparing observed statistics to a set of null statistics
#' at given thresholds.
#' @param z A vector of observed statistics.
#' @param znull A matrix of null statistics from permutations.
#' @param thresholds A numeric vector of thresholds for calculating FDR.
#' @return A vector of empirical FDRs.
#' @importFrom Matrix colMeans
#' @export
empirical_fdrs <- function(z, znull, thresholds) {   
  n <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:n, ])
  ranks <- t(tail_counts(thresholds, z)[1:n, ])
  
  fdp <- sweep(tails, 2, ranks, '/')
  fdr <- Matrix::colMeans(fdp)
  
  return(fdr)
}

#' Tail Counts for FDR Calculation
#'
#' A helper function to count the number of null statistics exceeding a set of thresholds.
#' This is used in the empirical FDR calculation.
#' @param z A vector of thresholds (derived from observed statistics).
#' @param znull A matrix of null statistics.
#' @return A matrix where each column corresponds to a column of `znull` and rows
#' represent the count of null values exceeding the squared thresholds in `z`.
#' @export
tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))      
  })
}

#' Neighborhood Abundance Map (NAM) Association Test
#'
#' Performs an association test between a sample-level variable and cell state abundance,
#' as captured by a Neighborhood Abundance Map (NAM). The method involves diffusing
#' sample information across a cell-cell graph, performing dimensionality reduction on the
#' resulting NAM, and then using a permutation-based approach to test for significant
#' associations between the principal components of the NAM and the variable of interest.
#'
#' @param seurat_object A Seurat object containing the single-cell data and pre-computed graph.
#' @param metadata Optional. A data frame with cell metadata, used if `seurat_object` is NULL.
#' @param pcs Optional. A matrix of pre-computed principal components, used if `seurat_object` is NULL.
#' @param test_var The name of the numeric variable in sample-level metadata to test for association.
#' @param samplem_key The column name in the metadata that identifies unique samples or individuals.
#' @param graph_use The name of the graph to use from the Seurat object (e.g., 'RNA_snn').
#' @param batches Optional. The column name(s) for batch variables, used for conditional permutation.
#' @param covs Optional. The column name(s) for other covariates to regress out from the NAM.
#' @param nsteps Optional. The number of diffusion steps. If NULL, it is determined automatically.
#' @param verbose Logical. If TRUE, print progress messages.
#' @param assay The assay to store the results in.
#' @param key The key for the new DimReduc object.
#' @param maxnsteps The maximum number of diffusion steps allowed when `nsteps` is NULL.
#' @param max_frac_pcs The maximum fraction of samples to use as the number of NAM PCs.
#' @param ks A numeric vector specifying the number of PCs to test for the association model. If NULL, this is determined automatically.
#' @param Nnull The number of permutations to generate for the null distribution.
#' @param force_permute_all Logical. If TRUE, ignore batches and permute across all samples.
#' @param local_test Logical. If TRUE, computes neighborhood-level FDRs to identify specific significant neighborhoods.
#' @param seed An integer for setting the random seed for reproducibility.
#' @param return_nam Logical. If TRUE, stores the NAM embeddings and loadings in the result.
#' @return A Seurat object with a new `DimReduc` object (default key 'cna') containing the NAM results.
#' The `misc` slot of this object holds the detailed association test statistics. New metadata columns with
#' neighborhood correlations (ncorrs) are also added.
#' @importFrom glue glue
#' @importFrom moments kurtosis
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject FindNeighbors Key
#' @importFrom tibble remove_rownames rownames_to_column
#' @importFrom dplyr select one_of
#' @importFrom stats as.formula model.matrix median pf sd
#' @importFrom RSpectra svds
#' @export
association_nam = function(seurat_object=NULL,
                           metadata=NULL,
                           pcs=NULL,
                           test_var=NULL,
                           samplem_key = NULL,
                           graph_use = 'RNA_snn',
                           batches = NULL,
                           covs = NULL ,
                           nsteps = NULL,
                           verbose=TRUE ,
                           assay=NULL,
                           key='NAMPC_',
                           maxnsteps=15L,
                           max_frac_pcs=0.15, 
                           ks=NULL,
                           Nnull=1000,
                           force_permute_all=FALSE,
                           local_test=TRUE,
                           seed=1234,
                           return_nam=TRUE
){
  if (!is.null(seurat_object) && is.null(metadata) && is.null(pcs)) {
    message(paste0("will use Seurat object following analysis..."))
  } else if (is.null(seurat_object) && !is.null(metadata) && !is.null(pcs)){
    meta = metadata
    rownames(meta) = 1:nrow(meta)
    
    m <- as(t(pcs), "dgTMatrix") 
    colnames(m) = 1:ncol(m)
    
    obj <- Seurat::CreateSeuratObject(
      counts = m, 
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
  } else if ((is.null(seurat_object) && is.null(metadata) && !is.null(pcs))|
             (is.null(seurat_object) && !is.null(metadata) && is.null(pcs))){
    stop('Must provide both metadata and precomuputed PCs')
  } 
  
  covs_keep <- test_var
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
  
  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(data.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }
  
  res <- list()
  if(!is.null(covs)){
    covs_mat <- data.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))
  } else {
    covs_mat <- NULL
  }
  
  f <- as.formula(as.character(glue('~0+{data$samplem_key}')))   
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]    
  s <- s[, data$samplem[[data$samplem_key]]] 
  
  prevmedkurt <- Inf
  
  diffuse_step <- function(data, s) {
    a <- data$connectivities
    degrees <- Matrix::colSums(a) + 1 
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(as.matrix(res)) 
  }
  for (i in seq_len(maxnsteps)) { 
    s <- diffuse_step(data, s) 
    medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis)) 
    if (is.null(nsteps)) { 
      prevmedkurt <- medkurt
      if (prevmedkurt - medkurt < 3 & i > 3) { 
        message(glue::glue('stopping after {i} steps'))
        break 
      }         
    } else if (i == nsteps) {
      break
    }
  }  
  snorm <- t(prop.table(s, 2)) 
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]
  
  NAM=snorm
  
  N <- nrow(NAM)
  if (is.null(batches_vec) | length(unique(batches_vec)) == 1) {
    message('only one unique batch supplied to qc')
    keep <- rep(TRUE, ncol(NAM))
  } else {
    message('filtering based on batches kurtosis')  
  }
  
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
  
  kurtoses <- .batch_kurtosis(NAM, batches_vec) 
  threshold <- max(6, 2*median(kurtoses)) 
  message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}')) 
  keep <- which(kurtoses < threshold) 
  .res_qc_nam = list(NAM = NAM, keep = keep)
  
  res <- list()
  res[[paste0('NAM.T', suffix)]] <- t(.res_qc_nam[[1]])
  res[[paste0('keptcells', suffix)]] <- .res_qc_nam[[2]]
  res[[paste0('_batches', suffix)]] <- batches_vec
  
  if (verbose) message('Residualize NAM')
  N <- nrow(NAM)
  NAM_ <- scale(NAM, center = TRUE, scale = FALSE) 
  ncols_C <- 0
  if (!is.null(covs_mat)) { 
    covs_mat <- scale(covs_mat)
    ncols_C <- ncols_C + ncol(covs_mat)
  }
  if (is.null(covs_mat)) {
    M <- Matrix::Diagonal(n = N) 
  } else {
    M <- Matrix::Diagonal(n = N) - covs_mat %*% solve(t(covs_mat) %*% covs_mat, t(covs_mat))  
  }
  if (any(is.nan(M))) {
    M <- Matrix::Diagonal(n = N) 
    message('Warning: covs_mat is singular, using identity matrix') 
  }
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
  
  npcs <- min(npcs, nrow(data$samplem) - 1) 
  if (missing(npcs) | npcs > .5 * min(dim(NAM_))) {
    svd_res <- svd(scale(NAM_, center=FALSE, scale=TRUE))
  } else {
    svd_res <- RSpectra::svds(scale(NAM_, center=FALSE, scale=TRUE), k = npcs)
  }
  
  dim(svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v))
  
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
  
  M=res[[paste0('_M', suffix)]]
  r=res[[paste0('_r', suffix)]]
  
  yvals <- rcna_data$samplem[[test_var]]
  if (is(yvals, 'character') | is(yvals, 'factor') | is(yvals, 'integer') ) {
    stop('test_var is of class {class(yvals)}. It must be numeric variable for association testing.')
  }
  y = yvals 
  
  if (is.null(seed)) {
    set.seed(sample(1e6, 1))
  }
  if (force_permute_all) {
    batches_vec <- rep(1L, length(y))
  }
  
  U <- NAMsvd[[1]]
  sv <- NAMsvd[[2]]
  V <- NAMsvd[[3]]
  y <- scale(y)
  n <- length(y)
  
  if (is.null(ks)) {
    incr <- max(round(0.02*n), 1)
    maxnpcs <- min(4*incr, round(n/5))
    ks <- seq(incr, maxnpcs+1, incr)
  }
  
  
  .reg <- function(q, k) {
    Xpc <- U[, 1:k] 
    beta <- t(Xpc) %*% q 
    qhat <- Xpc %*% beta 
    return(list(qhat = qhat, beta = beta))
  }
  
  .stats <- function(yhat, ycond, k) {
    ssefull <- crossprod(yhat - ycond) 
    ssered <- crossprod(ycond) 
    deltasse <-  ssered - ssefull 
    f <- (deltasse / k) / (ssefull/n) 
    p <- -pf(f, k, n-(1+r+k), log.p = TRUE)   
    r2 <- 1 - ssefull/ssered 
    return(list(p=p, r2=r2))
  }
  
  .minp_stats <- function(z) {
    zcond <- scale(M %*% z, center = FALSE, scale = TRUE) 
    qhats <- purrr::map(ks, function(k) .reg(zcond, k)$qhat) 
    .tmp <- purrr::map2(qhats, ks, function(qhat, k) .stats(qhat, zcond, k)) 
    ps <- purrr::map_dbl(.tmp, 'p') 
    r2s <- purrr::map_dbl(.tmp, 'r2')
    k_ <- which.min(ps) 
    return(list(k=ks[k_], p=ps[k_], r2=r2s[k_]))
  }
  
  .tmp <- .minp_stats(y)
  k <- .tmp$k
  p <- .tmp$p
  r2 <- .tmp$r2
  if (k == max(ks)) {
    warning(glue::glue('data supported use of {k} NAM PCs, which is the maximum considered. Consider allowing more PCs by using the "ks" argument.'))      
  }
  ycond <- scale(M %*% y, center = FALSE, scale = TRUE)
  .tmp <- .reg(ycond, k) 
  yhat <- .tmp$qhat
  beta <- .tmp$beta
  r2_perpc <- (beta / as.numeric(sqrt(crossprod(ycond))))**2
  
  ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta/n)  
  rownames(ncorrs) <- rownames(V)
  
  set.seed(seed)
  y_ <- conditional_permutation(batches_vec, y, Nnull)
  .tmp <- apply(y_, 2, .minp_stats)
  nullminps <- purrr::map_dbl(.tmp, 'p')
  nullr2s <- purrr::map_dbl(.tmp, 'r2')
  
  pfinal <- (sum(nullminps <= p+1e-8) + 1)/(Nnull + 1)
  if (sum(nullminps <= p+1e-8) == 0) {
    warning('global association p-value attained minimal possible value. Consider increasing Nnull')
  }
  
  fdrs <- NULL
  fdr_5p_t <- NULL 
  fdr_10p_t <- NULL
  fdr_20p_t <- NULL
  fdr_30p_t <- NULL
  fdr_40p_t <- NULL
  fdr_50p_t <- NULL
  
  if (local_test) {   
    message('computing neighborhood-level FDRs')
    Nnull <- min(1000, Nnull)
    y_ <- y_[, 1:Nnull]
    ycond_ <- scale(M %*% y_, center = FALSE, scale = TRUE)
    gamma_ <- crossprod(U[, 1:k], ycond_)
    nullncorrs <- abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_ / n)))
    
    maxcorr <- max(abs(ncorrs))
    fdr_thresholds <- seq(maxcorr/4, maxcorr, maxcorr/400)
    fdr_vals <- empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)
    fdrs = data.frame(
      threshold = head(fdr_thresholds, -1),
      fdr = fdr_vals, 
      num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t)) 
    )
    if (min(fdrs$fdr) > 0.05) {      
      fdr_5p_t <- NULL
    } else {
      fdr_5p_t <- min(subset(fdrs, fdr < 0.05)$threshold)      
    }
    if (min(fdrs$fdr) > 0.10) {      
      fdr_10p_t <- NULL
    } else {
      fdr_10p_t <- min(subset(fdrs, fdr < 0.1)$threshold)
    }
    if (min(fdrs$fdr) > 0.20) {      
      fdr_20p_t <- NULL
    } else {
      fdr_20p_t <- min(subset(fdrs, fdr < 0.2)$threshold)
    }
    if (min(fdrs$fdr) > 0.30) {      
      fdr_30p_t <- NULL
    } else {
      fdr_30p_t <- min(subset(fdrs, fdr < 0.3)$threshold)
    }
    if (min(fdrs$fdr) > 0.40) {      
      fdr_40p_t <- NULL
    } else {
      fdr_40p_t <- min(subset(fdrs, fdr < 0.4)$threshold)
    }
    if (min(fdrs$fdr) > 0.50) {      
      fdr_50p_t <- NULL
    } else {
      fdr_50p_t <- min(subset(fdrs, fdr < 0.5)$threshold)
    }
  }
  
  res <- list(
    p = pfinal, 
    nullminps=nullminps,
    k=k,
    ncorrs=ncorrs, 
    fdrs=fdrs,
    fdr_5p_t=fdr_5p_t, 
    fdr_10p_t=fdr_10p_t,
    fdr_20p_t=fdr_20p_t,
    fdr_30p_t=fdr_30p_t,
    fdr_40p_t=fdr_40p_t,
    fdr_50p_t=fdr_50p_t,
    yhat=yhat, 
    ycond=ycond,
    ks=ks, 
    beta=beta,
    r2=r2, 
    r2_perpc=r2_perpc,
    nullr2_mean=mean(nullr2s), 
    nullr2_std=sd(nullr2s)
  )
  
  if (return_nam) {
    res[['NAM_embeddings']] <- nam_res$NAM_nbhdXpc
    res[['NAM_loadings']] <- nam_res$NAM_sampleXpc
    res[['NAM_svs']] <- nam_res$NAM_svs
    res[['NAM_varexp']] <- nam_res$NAM_varexp
  }
  
  seurat_object[['cna']] <- Seurat::CreateDimReducObject(
    embeddings = res$NAM_embeddings, 
    loadings = res$NAM_loadings, 
    stdev = res$NAM_svs, 
    assay = assay, 
    key = key,
    misc = res 
  )
  
  seurat_object@meta.data$cna_ncorrs <- ncorrs[colnames(seurat_object), , drop=TRUE]
  seurat_object@meta.data$cna_ncorrs_fdr05 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_5p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_5p_t)
    seurat_object@meta.data$cna_ncorrs_fdr05[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  
  seurat_object@meta.data$cna_ncorrs_fdr10 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_10p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_10p_t)
    seurat_object@meta.data$cna_ncorrs_fdr10[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  
  seurat_object@meta.data$cna_ncorrs_fdr20 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_20p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_20p_t)
    seurat_object@meta.data$cna_ncorrs_fdr20[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  
  seurat_object@meta.data$cna_ncorrs_fdr30 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_30p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_30p_t)
    seurat_object@meta.data$cna_ncorrs_fdr30[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  
  seurat_object@meta.data$cna_ncorrs_fdr40 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_40p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_40p_t)
    seurat_object@meta.data$cna_ncorrs_fdr40[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  
  seurat_object@meta.data$cna_ncorrs_fdr50 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_50p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_50p_t)
    seurat_object@meta.data$cna_ncorrs_fdr50[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  return(seurat_object)
}