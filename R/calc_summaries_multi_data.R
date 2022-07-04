#' Calculate Wakefield and Laplace Bayes factors.
#' @param log_reg_results_array a 3 dimensional array of dimension num_snps by
#' 2 by num_datasets. Column 1 contains the estimates of the effect sizes for
#' each SNP, column 2 contains the variances of the effect size estimates for
#' each SNP in each dataset.
#' @param num_snps The number of SNPs to be fine mapped.
#' @param num_datasets The number of datasets being supplied (designed for use
#' with simulated data from, for example, Hapgen).
#' @param calc_name A character string representing the statistic to be used.
#' Possibilities are lbf and wbf for Laplace and Wakefield Bayes factor; mean,
#' median and mode for posterior mean, median and mode; ci and hdi for 95
#' percent posterior credible or high density interval.
#' @inheritParams calc_lbf
#' @inheritParams calc_wbf
#' @export
calc_summary <- function(log_reg_results_array, lambda, num_snps,
                         num_datasets, calc_name, w){

  post_sum <- matrix(NA, nrow = num_snps, ncol = num_datasets)

  for(dataset in 1:num_datasets){
    beta_hat <- log_reg_results_array[, 1, dataset]
    v <- log_reg_results_array[, 2, dataset]
    max_ci_size <- rep(NA, num_snps)
    modal_value <- rep(NA, num_snps)

    Q_neg <- beta_hat + (v * lambda)
    Q_pos <- beta_hat - (v * lambda)

    Q_neg2 <- Q_neg ^ 2
    Q_pos2 <- Q_pos ^ 2

    Z_neg <- pnorm(0, Q_neg, sqrt(v))
    Z_pos <- pnorm(0, Q_pos, sqrt(v))

    C_neg <- exp(((beta_hat ^ 2) - (Q_neg ^ 2)) / (-2 * v))
    C_pos <- exp(((beta_hat ^ 2) - (Q_pos ^ 2)) / (-2 * v))

    D_neg <- C_neg * Z_neg
    D_pos <- C_pos *(1 - Z_pos)
    D <- D_neg + D_pos

    E_neg <- C_neg / D
    E_pos <- C_pos / D

    post_sum[, dataset] <- switch(calc_name,
     lbf = calc_lbf(beta_hat, v, lambda),
     wbf <- calc_lbf(beta_hat, v, w),
     mean = {
       cst <- sqrt(0.5 * v / pi)
       EV_pos <- cst * exp(-0.5 * Q_pos2 / v) + (Q_pos * (1 - Z_pos))
       EV_neg <- cst * exp(-0.5 * Q_neg2 / v) - (Q_neg * Z_neg)
       E_pos * EV_pos - E_neg * EV_neg
     },
     median = {
       m_neg <- sqrt(v) * qnorm(0.5 / E_neg) + Q_neg
       m_pos <- sqrt(v) * qnorm(1 - 0.5 / E_pos) + Q_pos
       ifelse(E_pos * (1 - Z_pos) >= 0.5, m_pos, m_neg)
     },
     mode = {
       for(c in 1:num_snps){
         if(Q_neg[c] > 0 & Q_pos[c] < 0) modal_value[c] <- 0
         if(Q_neg[c] < 0 & Q_pos[c] < 0) modal_value[c] <- Q_neg[c]
         if(Q_neg[c] > 0 & Q_pos[c] > 0) modal_value[c] <- Q_pos[c]
       }
       modal_value
     },
     ci = {
       for(c in 1:num_snps){
         low_tail <- E_neg[c] * Z_neg[c] #lower tail probability
         high_tail <- E_pos[c] * (1 - Z_pos[c]) # upper tail probability
         if(Q_neg[c] > 0 & Q_pos[c] < 0) max_ci_size[c]  <- 1 - 2 *
           min(low_tail, high_tail)
         if(Q_neg[c] < 0 & Q_pos[c] < 0) max_ci_size[c]  <- 1 - 2 * high_tail
         if(Q_neg[c] > 0 & Q_pos[c] > 0) max_ci_size[c]  <- 1 - 2 * low_tail
       }
       max_ci_size
     },
     hdi = {
       for(c in 1:num_snps){
         low_tail <- E_neg[c] * Z_neg[c] #lower boundary
         high_tail <- E_pos[c] * (1 - Z_pos[c]) #upper boundary
         if(Q_neg[c] > 0 & Q_pos[c] < 0) max_ci_size[c]  <- 0
         if(Q_neg[c] < 0 & Q_pos[c] < 0) max_ci_size[c]  <- 1 - 2 * high_tail
         if(Q_neg[c] > 0 & Q_pos[c] > 0) max_ci_size[c]  <- 1 - 2 * low_tail
       }
       max_ci_size
     }) # end of switch function
  } # end of dataset loop
  return(post_sum)
} # end of function definition
