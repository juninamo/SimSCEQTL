# FILE: R/generate_pseudo_features_fixedVar.R

#' Generate Pseudo-Features for Single Cells
#'
#' Creates a matrix of latent features (analogous to PCs) for single cells.
#' The variance of these features can be attributed to different metadata
#' variables like cell type, disease, batch, etc., according to specified ratios.
#'
#' @param data The cell-level metadata data frame.
#' @param n_features The number of latent features to generate.
#' @param cluster_features Integer vector specifying which features are associated with cell clusters.
#' @param disease_features Integer vector specifying which features are associated with disease status.
#' @param sex_features Integer vector specifying which features are associated with sex.
#' @param age_features Integer vector specifying which features are associated with age.
#' @param bmi_features Integer vector specifying which features are associated with BMI.
#' @param batch_features Integer vector specifying which features are associated with batch.
#' @param individual_features Integer vector specifying which features are associated with individual identity.
#' @param cluster_ratio Maximum proportion of variance a `cluster_feature` can explain.
#' @param disease_ratio Maximum proportion of variance a `disease_feature` can explain.
#' @param sex_ratio Maximum proportion of variance a `sex_feature` can explain.
#' @param age_ratio Maximum proportion of variance an `age_feature` can explain.
#' @param bmi_ratio Maximum proportion of variance a `bmi_feature` can explain.
#' @param batch_ratio Maximum proportion of variance a `batch_feature` can explain.
#' @param individual_ratio Maximum proportion of variance an `individual_feature` can explain.
#' @param ratio_variance Numeric, amount of random variation for variance ratios.
#' @param cluster_col Name of the cell type column in the metadata.
#' @param disease_col Name of the disease status column in the metadata.
#' @param sex_col Name of the sex column in the metadata.
#' @param age_col Name of the age column in the metadata.
#' @param bmi_col Name of the BMI column in the metadata.
#' @param batch_col Name of the batch column in the metadata.
#' @param individual_col Name of the individual ID column in the metadata.
#' @param common_var A numeric value to standardize the common variance across features.
#' @param seed An integer for reproducibility.
#' @param n_thread The number of parallel threads to use.
#' @param verbose Logical. If TRUE, progress messages will be printed.
#'
#' @return A matrix of pseudo-features where rows are cells and columns are features.
#'
#' @importFrom doParallel registerDoParallel stopCluster makeCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom dplyr %>% case_when group_by summarize
#' @importFrom stats rnorm runif
#'
#' @export
generate_pseudo_features_fixedVar = function(data,
                                             # Number of features (pseudo principal components)
                                             n_features = 20,
                                             # Features associated with each attribute
                                             cluster_features = 1:20,
                                             disease_features = 0,
                                             sex_features = 0,
                                             age_features = 0,
                                             bmi_features = 0,
                                             batch_features = 0,
                                             individual_features = 0,
                                             
                                             # Define the maximum proportion of variance for each attribute
                                             cluster_ratio = 0.25,
                                             disease_ratio = 0,
                                             sex_ratio = 0,
                                             age_ratio = 0,
                                             bmi_ratio = 0,
                                             batch_ratio = 0,
                                             individual_ratio = 0.1,
                                             ratio_variance = 0.5,
                                             
                                             cluster_col = "cell_type",
                                             disease_col = "disease",
                                             sex_col = "sex",
                                             age_col = "age",
                                             bmi_col = "bmi",
                                             batch_col = "batch",
                                             individual_col = "subject_id",
                                             
                                             common_var = 1, 
                                             
                                             seed = 1234,
                                             n_thread = 1,
                                             verbose = TRUE
) {
  
  set.seed(seed)
  disease_features = sample(disease_features)
  set.seed(seed*2)
  sex_features = sample(sex_features)
  set.seed(seed*3)
  age_features = sample(age_features)
  set.seed(seed*4)
  bmi_features = sample(bmi_features)
  set.seed(seed*5)
  batch_features = sample(batch_features)
  set.seed(seed*6)
  individual_features = sample(individual_features)
  
  # Register parallel backend to use multiple cores
  cl <- makeCluster(n_thread)
  registerDoParallel(cl)
  
  pcs_list = foreach(x = 1:n_features, .packages = c("dplyr")) %dopar% {
    
    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluster_col]))
    
    cell_sex <- as.integer(factor(data[,sex_col]))
    cell_age <- as.integer(factor(data[,age_col]))
    cell_bmi <- as.integer(factor(data[,bmi_col]))
    cell_batch <- as.integer(factor(data[,batch_col]))
    cell_diseases <- as.integer(factor(data[,disease_col]))
    cell_individual <- as.integer(factor(data[,individual_col]))
    
    var_all <- c()
    
    # Generate dummy features reflecting cell clusters
    cell_clusters <- as.integer(factor(data[,cluster_col]))
    for (j in 1:2){
      cell_clusters_tmp = cell_clusters
      set.seed(x*j)
      cluster_mean = sample(unique(cell_clusters_tmp))*10
      for (i in 1:length(unique(cell_clusters))){
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      assign(paste0("cell_clusters_means",j),cell_clusters_tmp)
    }
    
    variance <- 1 / cell_clusters_means2 # cell type specific variance
    set.seed(seed*x)
    
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters_means1), 
                       sd = #sqrt(variance)
                         common_var)
    var_all = c(var_all, variance)
    
    # Similarly, generate dummy PC for other covariates with changing seeds
    cell_covariate = cell_diseases
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*2)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*2)
    pc_disease = rnorm(n_cells, mean = scale(cell_covariates_means1), 
                       sd = #sqrt(variance)
                         common_var)
    var_all = c(var_all, variance)
    
    cell_covariate = cell_sex
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*3)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*3)
    pc_sex = rnorm(n_cells, mean = scale(cell_covariates_means1), 
                   sd = #sqrt(variance)
                     common_var)
    var_all = c(var_all, variance)
    
    # Treat age as fixed effect for PC mean
    cell_covariate = cell_age
    for (j in 2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*4)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*4)
    pc_age = rnorm(n_cells, mean = scale(data[,age_col]), 
                   sd = #sqrt(variance)
                     common_var)
    var_all = c(var_all, variance)
    
    cell_covariate = cell_bmi
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*5)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*5)
    pc_bmi = rnorm(n_cells, mean = scale(data[,bmi_col]), 
                   sd = #sqrt(variance)
                     common_var)
    var_all = c(var_all, variance)
    
    cell_covariate = cell_batch
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*6)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*6)
    pc_batch = rnorm(n_cells, mean = scale(cell_covariates_means1), 
                     sd = #sqrt(variance)
                       common_var)
    var_all = c(var_all, variance)
    
    cell_covariate = cell_individual
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*7)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*7)
    pc_individual = rnorm(n_cells, mean = scale(cell_covariates_means1), 
                          sd = #sqrt(variance)
                            common_var)
    var_all = c(var_all, variance)
    
    # Generate dummy PC for noise term
    set.seed(seed*x*100)
    variance = sample(var_all,n_cells)
    set.seed(seed*x*100)
    pc = rnorm(n_cells, mean = 0, 
               sd = #sqrt(variance)
                 common_var)
    
    cluster_ratio_tmp = cluster_ratio + runif(1, min = -cluster_ratio*ratio_variance, max = cluster_ratio*ratio_variance)
    disease_ratio_tmp = disease_ratio + runif(1, min = -disease_ratio*ratio_variance, max = disease_ratio*ratio_variance)
    sex_ratio_tmp = sex_ratio + runif(1, min = -sex_ratio*ratio_variance, max = sex_ratio*ratio_variance)
    age_ratio_tmp = age_ratio + runif(1, min = -age_ratio*ratio_variance, max = age_ratio*ratio_variance)
    bmi_ratio_tmp = bmi_ratio + runif(1, min = -bmi_ratio*ratio_variance, max = bmi_ratio*ratio_variance)
    batch_ratio_tmp = batch_ratio  + runif(1, min = -batch_ratio*ratio_variance, max = batch_ratio*ratio_variance)
    individual_ratio_tmp = individual_ratio  + runif(1, min = -individual_ratio*ratio_variance, max = individual_ratio*ratio_variance)
    if (!(x %in% cluster_features)) cluster_ratio_tmp = 0
    if (!(x %in% disease_features)) disease_ratio_tmp = 0
    if (!(x %in% sex_features)) sex_ratio_tmp = 0
    if (!(x %in% age_features)) age_ratio_tmp = 0
    if (!(x %in% bmi_features)) bmi_ratio_tmp = 0
    if (!(x %in% batch_features)) batch_ratio_tmp = 0
    if (!(x %in% individual_features)) individual_ratio_tmp = 0
    
    noise_ratio_tmp = (1
                       -cluster_ratio_tmp
                       -disease_ratio_tmp
                       -sex_ratio_tmp
                       -age_ratio_tmp
                       -bmi_ratio_tmp
                       -batch_ratio_tmp
                       -individual_ratio_tmp
    )
    if (noise_ratio_tmp < 0) noise_ratio_tmp = 0
    
    if(verbose){
      message(paste0("Feature",x,";\ncluster ratio = ",cluster_ratio_tmp,
                     "\ndisease ratio = ", disease_ratio_tmp,
                     "\nsex ratio = ", sex_ratio_tmp,
                     "\nage ratio = ", age_ratio_tmp,
                     "\nbmi ratio = ", bmi_ratio_tmp,
                     "\nbatch ratio = ", batch_ratio_tmp,
                     "\nindividual ratio = ", individual_ratio_tmp,
                     "\nnoise ratio = ", noise_ratio_tmp))
    }
    return(scale(pc) * noise_ratio_tmp +
             scale(pc_cluster) * cluster_ratio_tmp  +
             scale(pc_disease) * disease_ratio_tmp  +
             scale(pc_sex) * sex_ratio_tmp  +
             scale(pc_age) * age_ratio_tmp  +
             scale(pc_bmi) * bmi_ratio_tmp  +
             scale(pc_batch) * batch_ratio_tmp +
             scale(pc_individual) * individual_ratio_tmp
    )
  }
  
  stopCluster(cl)
  
  pcs <- do.call(cbind, pcs_list)
  return(pcs)
}
