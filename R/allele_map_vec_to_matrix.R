#' Convert allele map vector to a binary matrix
#' @param allele_map a vector of integers with length f equal to the number of founders
#' @return a binary matrix, with dimension f by a, where a is the number of unique QTL alleles
#' @export
allele_map_vec_to_matrix <- function(allele_map){
  mat <- matrix(data = NA, nrow = length(allele_map), ncol = length(unique(allele_map)))
  for (i in 1:ncol(mat)){
    mat[ , i] <- as.numeric(allele_map == i - 1)
  }
  rownames(mat) <- LETTERS[1:nrow(mat)]
  colnames(mat) <- paste0("allele", 1:ncol(mat))
  return(mat)
}
