#' Calculate Laplace Bayes factors.
#' @param beta_hat The effect size estimate
#' @param v The variance of the effect size estimate.
#' @param lambda The rate parameter of the Laplace prior.
#' @examples calc_lbf(beta_hat = 0.1, v = 0.6, lambda = 5)
#' @export

calc_lbf <- function(beta_hat, v, lambda){
  Q_neg <- beta_hat + v * lambda
  Q_pos <- beta_hat - v * lambda
  D_neg <- exp(Q_neg ^ 2 / (2 * v)) * pnorm(0, Q_neg, sqrt(v))
  D_pos <- exp(Q_pos ^ 2 / (2 * v)) * (1 - pnorm(0, Q_pos, sqrt(v)))
  lambda * sqrt(pi * v / 2 ) * (D_neg + D_pos)
}



