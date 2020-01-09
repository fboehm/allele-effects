#' Calculate a p-value from a permutation test
#' 
#' @param observed_lod the test statistic, on the lod scale
#' @param maxlods a vector of max lods from a collection of permutations
#' @return a permutation test p-value
#' @details If the observed_lod is greater than all values in the vector maxlods, then a p-value of zero is returned. In practice, it means that you need more permutations.
#' @export

calc_pvalue <- function(observed_lod, maxlods){
  mean(observed_lod <= maxlods)
}
