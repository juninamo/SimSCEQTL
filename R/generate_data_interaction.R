# FILE: R/generate_data_interaction.R

#' Generate Simulated Single-Cell Data with eQTLs
#'
#' This is the main wrapper function to generate a complex simulated dataset.
#' It simulates cell metadata, gene expression counts based on latent features,
#' and pseudo-genotypes correlated with the expression of a specific gene in a
#' specific cell type to simulate eQTL effects.
#'
#' @param n_cells Base number of cells for major cell types per individual.
#' @param sd_celltypes Standard deviation for the number of cells.
#' @param n_major_cell_types Number of major cell types to simulate.
#' @param n_minor_cell_types Number of minor (rare) cell types to simulate.
#' @param relative_abundance The abundance ratio of minor to major cell types.
#' @param n_major_diff_celltypes Number of major cell types to set as differentially abundant in the 'disease' group.
#' @param n_minor_diff_celltypes Number of minor cell types to set as differentially abundant in the 'disease' group.
#' @param n_individuals Total number of individuals to simulate.
#' @param n_batchs Number of batches to distribute individuals into.
#' @param prop_sex Proportion of individuals assigned to sex '1'.
#' @param prop_disease Proportion of individuals assigned to disease group '1'.
#' @param fc_increase The fold change increase for differentially abundant cell types in the disease group.
#' @param n_features The number of latent features (pseudo principal components) to generate.
#' @param cluster_features Integer vector specifying which features should be associated with cell type clusters.
#' @param disease_features Integer vector specifying which features should be associated with disease status.
#' @param sex_features Integer vector specifying which features should be associated with sex.
#' @param age_features Integer vector specifying which features should be associated with age.
#' @param bmi_features Integer vector specifying which features should be associated with BMI.
#' @param batch_features Integer vector specifying which features should be associated with batch.
#' @param individual_features Integer vector specifying which features should be associated with individual identity.
#' @param cluster_ratio Maximum proportion of variance a `cluster_feature` can explain.
#' @param disease_ratio Maximum proportion of variance a `disease_feature` can explain.
#' @param sex_ratio Maximum proportion of variance a `sex_feature` can explain.
#' @param age_ratio Maximum proportion of variance an `age_feature` can explain.
#' @param bmi_ratio Maximum proportion of variance a `bmi_feature` can explain.
#' @param batch_ratio Maximum proportion of variance a `batch_feature` can explain.
#' @param individual_ratio Maximum proportion of variance an `individual_feature` can explain.
#' @param ratio_variance Numeric, amount of random variation to add to the variance ratios.
#' @param cluster_col Name of the cell type column in the metadata.
#' @param disease_col Name of the disease status column in the metadata.
#' @param sex_col Name of the sex column in the metadata.
#' @param age_col Name of the age column in the metadata.
#' @param bmi_col Name of the BMI column in the metadata.
#' @param batch_col Name of the batch column in the metadata.
#' @param individual_col Name of the individual ID column in the metadata.
#' @param dist_type The distribution to use for generating counts from features ("nb" for negative binomial or "poisson").
#' @param size_nb Size parameter for the negative binomial distribution (only used if `dist_type` is "nb").
#' @param sparsity The proportion of zero counts to introduce into the simulated count data.
#' @param eqtl_celltype The cell type in which to simulate the eQTL effect.
#' @param eqtl_gene The gene for which to simulate the eQTL effect.
#' @param cor_vals A numeric vector of target correlations for the eQTL simulation.
#' @param allele_freqs A numeric vector of allele frequencies for the eQTL simulation.
#' @param n_snps The number of simulated SNPs to generate for each condition.
#' @param seed An integer for reproducibility.
#' @param n_thread The number of parallel threads to use.
#' @param verbose Logical. If TRUE, progress messages will be printed.
#'
#' @return A list containing several data frames: `simulation_counts` (genes x cells count matrix),
#'   `dummy_data` (cell-level metadata), `bulk_df` (pseudobulk data),
#'   `genotype_df_list` (a list of simulated genotype matrices), and
#'   `cor_summary_df_list` (a summary of the simulated eQTL correlations).
#'
#' @importFrom dplyr %>% mutate left_join distinct bind_rows
#' @importFrom edgeR cpm
#' @importFrom stats qnorm rank rpois rnbinom cor cor.test quantile set.seed
#' @importFrom tibble rownames_to_column
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopCluster makeCluster
#' @importFrom parallel detectCores
#' @importFrom utils expand.grid
#' @importFrom iterators iter
#' @importFrom MASS mvrnorm
#'
#' @export
generate_data_interaction <- function(n_cells = 3000,
                                      sd_celltypes = 0.1,
                                      n_major_cell_types = 7,
                                      n_minor_cell_types = 3,
                                      relative_abundance = 0.1,
                                      n_major_interact_celltypes = 1,
                                      n_minor_interact_celltypes = 1,
                                      n_individuals = 30,
                                      n_batchs = 4,
                                      prop_sex = 0.5,
                                      prop_disease = 0.5,
                                      interaction_feature_1 = c("sex","age","disease","bmi"),
                                      interaction_feature_2 = c("sex","age","disease","bmi"),
                                      fc_interact = 0.1, 
                                      interaction_type = NULL,
                                      # Number of latent features
                                      n_features = 20,
                                      # Features associated with each attribute, ordered by variance explained
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
                                      dist_type = "nb",
                                      size_nb = 1,
                                      sparsity = 0,
                                      eqtl_celltype = "A",
                                      eqtl_gene = "Gene1",
                                      cor_vals = c(0, 0.15, 0.3),
                                      allele_freqs = c(0.05, 0.3),
                                      n_snps = 100,
                                      seed = 1234,
                                      n_thread = 1,
                                      verbose = TRUE
) {
  
  dummy_data <- generate_dummy_data_interaction(n_cells = n_cells,
                                                sd_celltypes = sd_celltypes,
                                                n_major_cell_types = n_major_cell_types,
                                                n_minor_cell_types = n_minor_cell_types,
                                                relative_abundance = relative_abundance,
                                                n_major_interact_celltypes = n_major_interact_celltypes,
                                                n_minor_interact_celltypes = n_minor_interact_celltypes,
                                                interaction_feature_1 = interaction_feature_1,
                                                interaction_feature_2 = interaction_feature_2,
                                                fc_interact = fc_interact,
                                                interaction_type = interaction_type,
                                                n_individuals = n_individuals,
                                                n_batchs = n_batchs,
                                                prop_sex = prop_sex,
                                                prop_disease = prop_disease,
                                                seed = seed
  )
  diff_cell_types = dummy_data[[2]]
  dummy_data = dummy_data[[1]]
  
  # Generate pseudo principal components for the data
  ## with the assumption that covariates are removed after harmonization
  pseudo_feature_matrix <- generate_pseudo_features(dummy_data,
                                                    n_features = n_features,
                                                    cluster_features = cluster_features,
                                                    disease_features = disease_features,
                                                    sex_features = sex_features,
                                                    age_features = age_features,
                                                    bmi_features = bmi_features,
                                                    batch_features = batch_features,
                                                    individual_features = individual_features,
                                                    cluster_ratio = cluster_ratio,
                                                    disease_ratio = disease_ratio,
                                                    sex_ratio = sex_ratio,
                                                    age_ratio = age_ratio,
                                                    bmi_ratio = bmi_ratio,
                                                    batch_ratio = batch_ratio,
                                                    individual_ratio = individual_ratio,
                                                    ratio_variance = ratio_variance,
                                                    cluster_col = cluster_col,
                                                    sex_col = sex_col,
                                                    age_col = age_col,
                                                    bmi_col = bmi_col,
                                                    batch_col = batch_col,
                                                    disease_col = disease_col,
                                                    individual_col = individual_col,
                                                    seed = seed,
                                                    n_thread = n_thread,
                                                    verbose = TRUE
  )
  colnames(pseudo_feature_matrix) <- paste0("Feature", 1:ncol(pseudo_feature_matrix))
  
  # Number of pseudo-genes to generate from the dummy features
  n_genes <- n_features
  
  # Start parallel backend
  cl <- makeCluster(n_thread)
  registerDoParallel(cl)
  
  # Initialize the result matrix
  simulation_results <- matrix(nrow = nrow(pseudo_feature_matrix), ncol = n_genes)
  
  # Parallel loop
  simulation_results <- foreach(i = 1:n_genes, .combine='cbind', .packages=c("stats")) %dopar% {
    if(dist_type == "poisson"){
      # Simulate gene expression
      set.seed(i)
      expression_data <- apply(scale(pseudo_feature_matrix[,i]), 1, function(x) rpois(1, lambda = exp(x - sparsity)))
    } else if (dist_type == "nb"){
      set.seed(i)
      expression_data <- apply(scale(pseudo_feature_matrix[,i]), 1, function(x) {
        mu <- exp(x - sparsity) # mean
        size <- size_nb    # this value should be adjusted based on analysis
        prob <- size / (size + mu)
        rnbinom(1, size, prob)
      })
    }
    return(expression_data)
  }
  
  # Convert results to a data frame
  simulation_counts <- as.data.frame(simulation_results)
  colnames(simulation_counts) <- paste0("Gene", 1:ncol(simulation_counts))
  simulation_counts[is.na(simulation_counts)] <- 0
  stopCluster(cl)
  
  dummy_data$nUMI = rowSums(simulation_counts)
  rownames(simulation_counts) = rownames(dummy_data)
  
  # Assume collapse_counts is a helper function available in the package
  all_collapse <- collapse_counts(t(simulation_counts), dummy_data, c("subject_id","cell_type"))
  colnames(all_collapse$counts_mat) <- paste0(all_collapse$meta_data$subject_id,
                                              "__",
                                              all_collapse$meta_data$cell_type)
  sum_exp = all_collapse$counts_mat
  
  # Log2 CPM normalization
  norm_exp <- log2(edgeR::cpm(sum_exp)+1)
  
  # Inverse normal transformation
  rn <- apply(norm_exp, 1, function(x){
    qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )
  })
  rn <- t(rn)
  
  all_collapse$meta_data = all_collapse$meta_data %>%
    dplyr::mutate(sample = paste0(subject_id,"__",cell_type)) %>%
    dplyr::left_join(.,
                     dummy_data[,c(individual_col,sex_col,age_col,bmi_col,batch_col,disease_col)] %>%
                       dplyr::distinct(),
                     by="subject_id")
  # all(all_collapse$meta_data$sample == colnames(all_collapse$counts_mat))
  # all(all_collapse$meta_data$sample == colnames(rn))
  
  bulk_df = merge(all_collapse$meta_data,
                  t(rn) %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("sample"),
                  by="sample")
  message("The simulated dataset contains ", nrow(bulk_df), " pseudobulk samples across ",
          length(unique(bulk_df$subject_id)), " individuals,\n",
          length(unique(bulk_df$cell_type)), " cell types,\n",
          ncol(rn), " genes,\n",
          "with ", nrow(simulation_counts), " single cells in total.")
  # print(head(bulk_df))
  
  # Set the cell type and gene pair for which to simulate a correlation (eQTL)
  Gene = eqtl_gene
  # Store the trait values that will be correlated with the pseudo-genotypes
  trait_values <- bulk_df[bulk_df$cell_type == eqtl_celltype, Gene]
  
  # Set the target correlation coefficients and allele frequencies
  tasks <- expand.grid(cor_val = cor_vals, allele_freq = allele_freqs)
  
  # Number of simulated SNPs to generate
  if(is.null(n_thread)) n_thread = detectCores()-1
  registerDoParallel(cores = n_thread)
  
  # List to store results
  results_list <- list()
  results_list <- foreach(task = iter(tasks, by = 'row') ) %dopar% {
    cor_val <- task$cor_val
    allele_freq <- task$allele_freq
    
    aggregate_trials <- function(cor_val, allele_freq, n_snps = 100, verbose = TRUE) {
      # Target correlation range
      target_cor_min <- cor_val - 0.05
      target_cor_max <- cor_val + 0.05
      
      # Prepare a list to store results
      correlations <- vector("list", n_snps)
      p_values <- vector("list", n_snps)
      
      # Parallel processing
      results <- foreach(i = 1:n_snps, .combine = 'rbind') %dopar% {
        generate_dosage_with_target_correlation <- function(trait_values, target_correlation, allele_freq, n_attempts = 5000, verbose = FALSE, n_trial) {
          n <- length(trait_values)
          # Calculate genotype assignment probabilities from allele frequency
          geno_probs <- c((1-allele_freq)^2, 2*allele_freq*(1-allele_freq), allele_freq^2)
          # Generate initial dosage_values
          set.seed(i) # For reproducibility
          dosage_values <- sample(0:2, n, replace = TRUE, prob = geno_probs)
          if(all(dosage_values == 0)) dosage_values[1] <- 1 # Ensure at least one non-zero value
          best_cor <- cor(trait_values, dosage_values)
          target_cor_min <- target_correlation - 0.015
          target_cor_max <- target_correlation + 0.015
          allele_freq_min <- allele_freq - 0.015
          allele_freq_max <- allele_freq + 0.015
          # Adjust to approach the target correlation and MAF
          for (attempt in 1:n_attempts) {
            set.seed(attempt * n_trial * i) # Change seed for each iteration
            idx <- sample(n, 1)
            current_dosage <- dosage_values[idx]
            # Randomly adjust the dosage
            possible_dosages <- setdiff(0:2, current_dosage)
            new_dosage <- sample(possible_dosages, 1)
            dosage_values[idx] <- new_dosage
            # Calculate the new correlation
            new_cor <- cor(trait_values, dosage_values)
            # Calculate the new MAF
            new_maf <- min(sum(dosage_values == 1) + 2*sum(dosage_values == 2), sum(dosage_values == 1) + 2*sum(dosage_values == 0)) / (2*n)
            # If the new correlation and MAF are closer to the target, update; otherwise, revert
            if (abs(new_cor - target_correlation) < abs(best_cor - target_correlation) && new_maf >= allele_freq_min && new_maf <= allele_freq_max) {
              best_cor <- new_cor
              # If a correlation and MAF within the target range are found, exit the loop
              if (new_cor >= target_cor_min && new_cor <= target_cor_max) {
                break
              }
            } else {
              dosage_values[idx] <- current_dosage # Revert the change
            }
          }
          return(dosage_values)
        }
        
        # Evaluate the result
        dosage_values <- generate_dosage_with_target_correlation(trait_values, cor_val, allele_freq, n_attempts = 1000, n_trial = i)
        cor_test_result <- cor.test(trait_values, dosage_values)
        cor_val_attempt <- cor_test_result$estimate
        p_val_attempt <- cor_test_result$p.value
        trait_values_cell <- simulation_counts[dummy_data$cell_type == eqtl_celltype, Gene]
        dosage_values_cell <- dplyr::left_join(dummy_data[dummy_data$cell_type == eqtl_celltype, ],
                                               data.frame(subject_id = bulk_df[bulk_df$cell_type == eqtl_celltype, individual_col],
                                                          dosage_values),
                                               by = "subject_id") %>%
          .$dosage_values
        cor_test_cell_result <- cor.test(trait_values_cell, dosage_values_cell)
        cor_val_attempt_cell <- cor_test_cell_result$estimate
        p_val_attempt_cell <- cor_test_cell_result$p.value
        
        success <- FALSE
        
        if(cor_val_attempt >= target_cor_min && cor_val_attempt <= target_cor_max && !is.na(as.numeric(cor_test_result$estimate))) {
          success <- TRUE
          return(c(dosage = dosage_values,
                   correlation = cor_val_attempt,
                   p_value = p_val_attempt,
                   correlation_cell = cor_val_attempt_cell,
                   p_value_cell = p_val_attempt_cell))
        }
        
        if (!success) {
          return(c(NA, NA, NA, NA, NA)) # If the correlation never reached the target range
        }
      }
      
      # Post-process the results
      genotypes <- results[,1:length(trait_values)]
      correlations <- as.numeric(results[, "correlation.cor"])
      p_values <- as.numeric(results[, "p_value"])
      correlations_cell <- as.numeric(results[, "correlation_cell.cor"])
      p_values_cell <- as.numeric(results[, "p_value_cell"])
      
      return(list(genotype_df = genotypes,
                  cor_summary = c(cor_val = cor_val,
                                  allele_freq = allele_freq,
                                  mean_cor = mean(correlations, na.rm = TRUE),
                                  ci_lower = quantile(correlations, probs = 0.025, na.rm = TRUE),
                                  ci_upper = quantile(correlations, probs = 0.975, na.rm = TRUE),
                                  power = mean(p_values < 0.05, na.rm = TRUE),
                                  mean_cor_cell = mean(correlations_cell, na.rm = TRUE),
                                  ci_lower_cell = quantile(correlations_cell, probs = 0.025, na.rm = TRUE),
                                  ci_upper_cell = quantile(correlations_cell, probs = 0.975, na.rm = TRUE),
                                  power_cell = mean(p_values_cell < 0.05, na.rm = TRUE),
                                  n_snps = n_snps,
                                  N = length(trait_values)
                  ))
      )
    }
    
    result <- aggregate_trials(cor_val, allele_freq, n_snps = n_snps, verbose = FALSE)
    return(result)
  }
  
  # Check if the ground truth correlations match the targets
  cor_summary_df_list <- lapply(results_list, function(x) x$cor_summary) %>% bind_rows()
  
  # For subsequent pseudobulk eQTL analysis, extract only the genotype_df elements
  genotype_df_list <- lapply(results_list, function(x) x$genotype_df)
  names(genotype_df_list) <- apply(tasks, 1, function(x) paste(x['cor_val'], x['allele_freq'], sep = "_"))
  
  return(list(simulation_counts = simulation_counts,
              dummy_data = dummy_data,
              bulk_df = bulk_df,
              genotype_df_list = genotype_df_list,
              cor_summary_df_list = cor_summary_df_list))
}
