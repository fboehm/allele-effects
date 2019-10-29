#' Collapse a genotype probabilities array according to allelic series
#' @param probs a genotype probabilities array (typically for one chromosome)
#' @param allele_map a numeric vector of integers mapping the founder alleles to QTL alleles
#' @return collapsed genotypes array
#' @export
collapse_geno <- function(probs, allele_map){
  map_matrix <- allele_map_vec_to_matrix(allele_map)
  list_of_matrices <- lapply(seq(dim(probs)[3]), function(x) probs[ , , x])
  sapply(X = list_of_matrices, FUN = function(x) {x %*% map_matrix}, simplify = "array")
}
