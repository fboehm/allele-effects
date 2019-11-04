#' Call timbr for use on cluster 
#' 
#' @param mylist a list containing two named components, probe and trait
#' @param traits_df a data.frame with one trait per column and column name corresponding to trait name
#' @param addcovar additive covariates matrix to be appended to a column of 1's.
#' @param prior_M prior_M input for TIMBR::TIMBR()
#' @param genoprobs a n-long list containing, for each subject, a matrix of 36-state genotype probabilities at genome-wide collection of markers
#' @return output (a named list) of a single call to TIMBR
#' @export

call_timbr <- function(mylist, traits_df, prior_M, genoprobs, A_matrix = mcv.data$prior.D$A){
  yy <- traits_df[ , colnames(traits_df) == mylist$trait]
  gp <- genoprobs %>% # a 52-long list of matrices, one per mouse line
    purrr::map(.f = function(probe) {
      lapply(X = geno, 
             FUN = function(gg) gg[rownames(gg) %in% probe, drop = FALSE]
    )
  }
  ) %>%
    bind_cols() %>%
    as.matrix() %>%
    t()
  indices_gp <- match(rownames(traits_df), rownames(gp))
  qtl_full_genoprobs <- gp[indices_gp, ]
  prior_d <- list(P = qtl_full_genoprobs,
                             A = A_matrix, # Describes the mapping from full genoprobs to additive dosages
                             fixed.diplo = FALSE)
  
  out <- TIMBR(y = yy,
        Z = cbind(1, addcovar), 
        prior.D = prior_d,
        prior.M = prior_M)
  return(out)
}
