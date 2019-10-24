## Caveats
## - TIMBR can't currently include a kinship matrix
##   It's a Bayesian procedure and adding an additional
##   variance component with a dense covariance matrix
<<<<<<< HEAD
##   is challenging, if not infeasible computationally. 
##   My experience is that the CC and DO don't really need 
##   the kinship adjustment a majority of the time.
## - TIMBR needs the full genoprobs (36 state). For the 
##   CC, qtl2 forces the genoprobs to be homozygous, 
##   potentially a strong assumption. Nevertheless, here
##   I go along with that assumption and put zeros in at 
##   all the heterozygous states. It is particularly 
=======
##   is challenging, if not infeasible computationally.
##   My experience is that the CC and DO don't really need
##   the kinship adjustment a majority of the time.
## - TIMBR needs the full genoprobs (36 state). For the
##   CC, qtl2 forces the genoprobs to be homozygous,
##   potentially a strong assumption. Nevertheless, here
##   I go along with that assumption and put zeros in at
##   all the heterozygous states. It is particularly
>>>>>>> b50dac7b2bdef54879a4b3e3b38fa2d047ff27d7
##   important for the DO. If you generate the qtl2 data, that
##   shouldn't be a problem, but it is a major pain at JAX
##   because everyone tends to collapse to additive dosages
##   and throw away the full genoprobs.

### Look at using TIMBR to map founder-specific variants under multi-allelic QTL
#devtools::install_github("wesleycrouse/TIMBR)
library(TIMBR)

# CC data included in TIMBR
data(mcv.data)

# Specify allelic series prior
# Suggested by Wes
# Influences how much prior weight it places on more or less complicated allelic series
prior_M <- list(model.type = "crp", # crp - Chinese Restaurant Process
<<<<<<< HEAD
                prior.alpha.type = "gamma", 
                prior.alpha.shape = 1, 
=======
                prior.alpha.type = "gamma",
                prior.alpha.shape = 1,
>>>>>>> b50dac7b2bdef54879a4b3e3b38fa2d047ff27d7
                prior.alpha.rate = 2.333415)

## Grab QTL probability matrix
# qtl_chr - chromosome QTL is located on
# qtl_locus - name of QTL locus, needed to pull it from the genoprobs
n <- nrow(cc_genoprobs[[qtl_chr]][,, qtl_locus])
qtl_full_genoprobs <- cbind(cc_genoprobs[[qtl_chr]][,, qtl_locus],
                            matrix(0, ncol = 28, nrow = n)

prior_d <- list(P = qtl_full_genoprobs,
                A = mcv.data$prior.D$A, # Describes the mapping from full genoprobs to additive dosages
                fixed.diplo = FALSE)
#call TIMBR
results <- TIMBR(y = cc_phenotype[rownames(qtl_full_genoprobs)], # Order needs to match genoprobs
                 Z = cbind(rep(1, n), cc_addcovar),
<<<<<<< HEAD
                 prior.D = prior_d, 
=======
                 prior.D = prior_d,
>>>>>>> b50dac7b2bdef54879a4b3e3b38fa2d047ff27d7
                 prior.M = prior_M)
#expected number of alleles
mean(results$post.K)
#distribution of number of alleles
hist(results$post.K)
#posterior densities of allele effects
TIMBR.plot(results)


