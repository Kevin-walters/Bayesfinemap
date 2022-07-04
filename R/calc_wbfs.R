#' Calculate Wakefield Bayes factors.
#' @param beta_hat The effect size estimate
#' @param v The variance of the effect size estimate.
#' @param w The variance of the Gaussian effect size prior
#' @examples calc_wbf(beta_hat = 0.1, v = 0.6, w = 0.04)
#' @export

calc_wbf <- function(beta_hat, v, w){
  sqrt(v / (v + w)) * exp(beta_hat ^ 2 * w / (2 * v * (v + w)))
}


