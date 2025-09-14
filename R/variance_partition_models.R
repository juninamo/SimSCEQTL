# FILE: R/variance_partition_models.R

#' Apply a Mixed Model and Calculate Variance Partition
#'
#' A wrapper to apply a specific lme4 mixed model to a column of data
#' and calculate the variance explained by each random effect using variancePartition.
#'
#' @param column_data A numeric vector to be used as the response variable.
#' @return A `varPart` object from the variancePartition package.
#' @importFrom lme4 lmer lmerControl
#' @importFrom variancePartition calcVarPart
#' @export
apply_model_to_column <- function(column_data) {
  # This function has a dependency on a global variable 'dummy_data'
  # which is not best practice. It should be passed as an argument.
  # For now, assuming it exists in the environment where this is called.
  form_test <- as.formula(paste0("column_data ~ (1|cell_type) + (1|sex) + (1|disease) + (1|batch) + (1|subject_id) + age + bmi"))
  fit <- lme4::lmer(form_test, data = dummy_data, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  var_stats <- variancePartition::calcVarPart(fit)
  return(var_stats)
}

#' Apply a GLMM Count Model and Calculate Variance Partition
#'
#' A wrapper to apply a specific lme4 generalized linear mixed model (Poisson)
#' and calculate the variance explained by each random effect.
#'
#' @param column_data A numeric vector of counts (response variable).
#' @param meta_data A data frame containing the covariates for the model.
#' @return A `varPart` object from the variancePartition package.
#' @importFrom lme4 glmer glmerControl
#' @importFrom variancePartition calcVarPart
#' @export
apply_countmodel_to_column <- function(column_data, meta_data) {
  form_test <- as.formula(paste0("column_data ~ (1|cell_type) + (1|sex) + (1|disease) + (1|batch) + (1|subject_id) + age + bmi"))
  fit <- lme4::glmer(form_test, family = "poisson", nAGQ = 0, data = meta_data, control = glmerControl(optimizer = "nloptwrap"))
  var_stats <- variancePartition::calcVarPart(fit)
  return(var_stats)
}