#' Determine which SNPS are noteworthy.
#' @param bf A vector of Bayes factors.
#' @param r The ratio of the costs of false non-discovery to that of false discovery.
#' @param pi_0 The prior probability of no association.
#' @export

###########################################################################################
########################### Declare snps as noteworthy or not #############################
###########################################################################################

find_noteworthy <- function(bf, pi_0, r){
  as.numeric(pi_0 / (pi_0 + bf * (1 - pi_0)) < r / (1 + r))
}

calc_ppa <- function(bf, pi_0){
  bf * (1 - pi_0) / (pi_0 + bf * (1 - pi_0))
}
