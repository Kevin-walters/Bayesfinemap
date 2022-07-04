#' Perform univariate logistic regression on case-control genotype data.
#' @param geno_array A (\code{num_cases} + \code{num_controls}) by \code{num_snps} by
#'   \code{num_datasets} array of genotypes.
#' @param num_cases The number of cases in \code{geno_array}.
#' @param num_controls The number of controls in \code{geno_array}.
#' @param num_snps The number of SNPs you are expecting \code{geno_array} to contain.
#' @param num_datasets The number of datasets you are expecting \code{geno_array} to contain.
#' @return A \code{num_snps} by three by \code{num_datasets} array. The three
#'   values of the second element of the array are the mle, variance of the mle
#'   and the p-value from a univariate logistic regression analysis (with null
#'   hypothesis that the effect size is zero).
#' @export

###########################################################################################
############################## logistic regression function ###############################
###########################################################################################

run_logistic_regression <- function(geno_array, num_cases, num_controls, num_snps, num_datasets){
  regression_output_by_dataset <- array(data = NA, dim = c(num_snps, 3 , num_datasets))
  if((num_cases + num_controls) != dim(geno_array)[1]){
    cat("sum of declared cases and controls does not match first dimension of genotype array")
    stop()
  }
  if(num_snps != dim(geno_array)[2]){
    cat("number of declared SNPs does not match second dimension of genotype array")
    stop()
  }
  if(num_datasets != dim(geno_array)[3]){
    cat("number of declared datasets does not match third dimension of genotype array")
    stop()
  }
  for(dataset in 1:num_datasets){
    pheno=c(rep(1, num_cases), rep(0, num_controls))
    geno_and_pheno_data <- cbind(pheno, geno_array[,, dataset])
    for(j in 1:num_snps) {
      model <- glm(geno_and_pheno_data[, 1] ~ geno_and_pheno_data[, j + 1], family=binomial)
      beta_hat <- coefficients(model)[2]
      var_beta_hat <- vcov(model)[2, 2]
      if(is.na(coefficients(model)[2]) == F) {
        regression_output_by_dataset[j, 1, dataset] <- beta_hat
        regression_output_by_dataset[j, 2, dataset] <- var_beta_hat
        regression_output_by_dataset[j, 3, dataset] <- 2 * (1 - pnorm(abs(beta_hat / sqrt(var_beta_hat))))
      }
    }
  } # end of 'for dataset in' loop
  return(regression_output_by_dataset)
}
