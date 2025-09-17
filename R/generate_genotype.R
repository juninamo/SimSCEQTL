# FILE: R/generate_genotype.R

#' Generate Simulated Single-Cell Data with eQTLs
#'
#' This is the main wrapper function to generate a complex simulated dataset.
#' It simulates cell metadata, gene expression counts based on latent features,
#' and pseudo-genotypes correlated with the expression of a specific gene in a
#' specific cell type to simulate eQTL effects.
#'
#' @param bulk_df A data frame containing pseudobulk expression data with columns for cell type and gene expression.
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
generate_genotype <- function(bulk_df,
                              eqtl_celltype = "A",
                              eqtl_gene = "Gene1",
                              cor_vals = c(0, 0.15, 0.3),
                              allele_freqs = c(0.05, 0.3),
                              n_snps = 100,
                              seed = 1234,
                              n_thread = 1,
                              verbose = TRUE
) {
  
  message("Proceeding to simulate genotypes...")
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
}