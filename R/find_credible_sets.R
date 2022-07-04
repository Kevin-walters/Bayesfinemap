#' Calculate posterior probabilities of causality and form credible sets
#' assuming there is a single causal SNP.
#' @title Find credible sets
#' @param ppc A vector of posterior probabilities of causality
#' (defined as Bayes factor/sum of Bayes factors).
#' @param causal_snp The causal SNP number where SNPs are numbered from one to
#'  the number of SNPs.
#' @param thresh The cumulative probability threshold used to determine the
#' credible set.
#' @param bf vector of Bayes factors.
#' @examples find_cred_set(ppc = c(0.10, 0.20, 0.35, 0.40), causal_snp = 3, thresh = 0.9)
#' @return find_ppc returns the vector of posterior probabilities of causality
#' for each SNP; find_cred_set returns a length 2 vector where the first element
#' is the size of the credible set and the second element is T or F depending on
#' whether the credible set contains the causal SNP.
#' @export

###########################################################################################
#################################### Find credible sets ###################################
###########################################################################################

find_cred_set <- function(ppc, causal_snp, thresh){
  index_sort_ppc <- order(ppc, decreasing = T)
  ordered_ppc <- ppc[index_sort_ppc]
  ppc_cum <- cumsum(ordered_ppc)
  set_size <- min(which(ppc_cum > thresh))
  contains_causal <- as.numeric(is.element(causal_snp, index_sort_ppc[1:set_size]))
  return(c(set_size, contains_causal))
}

#' @describeIn find_cred_set Calculate posterior probabilities of causality

find_ppc <- function(bf) bf / sum(bf)
