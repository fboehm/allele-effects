#' Call timbr for use on cluster 
#' 
#' @param input_vec a vector containing two named components, probe and trait
#' @param traits_df a data.frame with one trait per column and column name corresponding to trait name
#' @param addcovar additive covariates matrix to be appended to a column of 1's.
#' @param prior_M prior_M input for TIMBR::TIMBR()
#' @param genoprobs_array a n by 36 by m array of genotype probabilities, where n is the number of subjects and m the number of markers
#' @return output (a named list) of a single call to TIMBR
#' @export

call_timbr <- function(input_vec, 
                       traits_df, 
                       addcovar,
                       prior_M, 
                       genoprobs_array, 
                       A_matrix = mcv.data$prior.D$A){
  yy <- traits_df[ , colnames(traits_df) == input_vec$trait]
  gp <- genoprobs_array[ , , dimnames(genoprobs_array)[[3]] == input_vec$probe] 
  prior_d <- list(P = gp,
                             A = A_matrix, # Describes the mapping from full genoprobs to additive dosages
                             fixed.diplo = FALSE)
  
  out <- TIMBR::TIMBR(y = yy,
        Z = as.matrix(cbind(1, addcovar)), 
        prior.D = prior_d,
        prior.M = prior_M)
  return(out)
}
