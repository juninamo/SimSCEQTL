# FILE: R/generate_dummy_data.R

#' Generate Dummy Cell Metadata
#'
#' Simulates cell-level metadata including individual-level covariates
#' (sex, disease, age, batch, bmi) and cell type assignments. It can introduce
#' differential cell type abundance based on disease status.
#'
#' @param n_cells Base number of cells for major cell types per individual.
#' @param sd_celltypes Standard deviation for the number of cells.
#' @param n_major_cell_types Number of major cell types to simulate.
#' @param n_minor_cell_types Number of minor (rare) cell types to simulate.
#' @param relative_abundance The abundance ratio of minor to major cell types.
#' @param n_major_diff_celltypes Number of major cell types to make differentially abundant.
#' @param n_minor_diff_celltypes Number of minor cell types to make differentially abundant.
#' @param n_individuals Total number of individuals to simulate.
#' @param n_batchs Number of batches to distribute individuals into.
#' @param prop_sex Proportion of individuals assigned to sex '1'.
#' @param prop_disease Proportion of individuals assigned to disease group '1'.
#' @param fc_increase The fold change increase for differentially abundant cell types.
#' @param seed An integer for reproducibility.
#'
#' @return A list containing two elements: a data frame of the simulated cell metadata,
#'   and a character vector of the cell types designated as differentially abundant.
#'
#' @importFrom stats runif median
#' @importFrom dplyr %>% left_join group_by summarise filter
#'
#' @export
generate_dummy_data <- function(n_cells = 3000, # cells of major cell types per individual
                                sd_celltypes = 0.1,  # standard deviation of number of cells 
                                n_major_cell_types = 7,
                                n_minor_cell_types = 3,
                                relative_abundance = 0.1, # ratio between major and rare
                                n_major_diff_celltypes = 1,
                                n_minor_diff_celltypes = 1,
                                n_individuals = 30, # total individuals
                                n_batchs = 4,
                                prop_sex = 0.5, 
                                prop_disease = 0.5,  
                                fc_increase = 0.1, # additional FC of specified cell types which are from people with case groups compared to control groups (e.g., if you specify 0.5, FC comparing increased cell type between case and control groups will be 1.5 in total)
                                seed = 1234
) {
  n_cell_types = n_major_cell_types + n_minor_cell_types
  
  # Generate unique subject IDs
  subject_id <- paste0("SUB_", 1:n_individuals)
  set.seed(seed)
  sex = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_sex)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_sex)))
  )
  if (length(sex) != length(subject_id)) {
    sex = c(sex, rep(1,length(subject_id)-length(sex)))
  }
  set.seed(seed*2)
  disease = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_disease)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_disease))
    )
  )
  if (length(disease) != length(subject_id)) {
    disease = c(disease, rep(1,length(subject_id)-length(disease)))
  }
  set.seed(seed*3)
  age <- sample(18:60, n_individuals, replace = TRUE)
  set.seed(seed*4)
  bmi <- sample(15:35, n_individuals, replace = TRUE)
  batch <- rep(rep(1:n_batchs), length(subject_id))[1:length(subject_id)]
  dummy_data = data.frame(subject_id = subject_id,
                          sex = factor(sex, levels = c(0, 1)),
                          disease = factor(disease, levels = c(0, 1)),
                          age = age,
                          batch = factor(batch, levels = c(0:max(batch))),
                          bmi = bmi
  ) 
  
  # Major and rare cell type counts
  ## major_cell_types <- ceiling(n_cell_types / 2)
  major_cell_types <- n_major_cell_types
  ## rare_cell_types <- n_cell_types - major_cell_types
  rare_cell_types <- n_minor_cell_types
  
  celltype_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(celltype_df) <- c("cell_type", "subject_id")
  
  # Generate baseline of cell type data.frame 
  for (id in dummy_data$subject_id){
    set.seed(seed*5*grep(id,dummy_data$subject_id)*10)
    major_cell_counts <- round(runif(major_cell_types, n_cells-n_cells*sd_celltypes, n_cells+n_cells*sd_celltypes))
    set.seed(seed*6*grep(id,dummy_data$subject_id)*10)
    rare_cell_counts <- round(runif(rare_cell_types, n_cells*relative_abundance-n_cells*relative_abundance*sd_celltypes, n_cells*relative_abundance+n_cells*relative_abundance*sd_celltypes))
    cell_counts <- c(major_cell_counts, rare_cell_counts)
    for (i in 1:n_cell_types) {
      n <- cell_counts[i]
      celltype_df = rbind(celltype_df,
                          data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], n),
                                     subject_id = id)
      )
    }
  }
  
  # check
  ## celltype_df %>%
  ##   dplyr::group_by(subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  ## celltype_df %>%
  ##   dplyr::group_by(cell_type,subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  
  diff_clusters = c(1:n_major_diff_celltypes, (n_cell_types-n_minor_diff_celltypes+1):n_cell_types)
  diff_cell_types = LETTERS[diff_clusters]
  
  for (i in LETTERS[1:n_cell_types]) {
    if(i %in% diff_cell_types){
      
      abundance = dplyr::left_join(celltype_df,
                                   dummy_data,
                                   by="subject_id") %>%
        dplyr::filter(cell_type == i) %>%
        dplyr::group_by(subject_id) %>%
        dplyr::summarise(pro = dplyr::n()/n_cells,
                         count = dplyr::n())
      diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_increase)) - n_cells*median(abundance$pro))
      
      i# add diff cells only to disease == 1
      
      len = length(unique(dummy_data[dummy_data$disease == "1",]$subject_id)) * diff
      
      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "1",]$subject_id),diff)),
        celltype_df
      )
    }
  }
  dummy_data = merge(dummy_data,celltype_df,by="subject_id")
  # Shuffle rows
  set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(list(dummy_data,diff_cell_types))
}
