TnSeq traits TIMBR analysis
================
Frederick J. Boehm
11/4/2019

Last modified: 2019-11-05 21:02:00.

``` r
library(TIMBR)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
tt <- read.csv("../data/neto_traits_by_probe3_annotated.csv")
neto <- tt %>%
  tidyr::pivot_longer(cols = V2:V35, values_to = "trait", names_to = "trait_name") %>%
  dplyr::filter(!is.na(trait)) %>%
  dplyr::select(- trait_name, - neto.n, - row, - n.traits, - cM, -hs, - chr)
```

``` r
if (!file.exists("../data/genotypes_array.rds")){
  fns <- dir("../data/36-state-genotypes")
  geno <- list()
  for (i in seq_along(fns)){
    load(file.path("../data/36-state-genotypes", fns[i]))
    geno[[i]] <- prsmth
  }
  names(geno) <- stringr::str_split_fixed(fns, ".genotype.probs.Rdata", 2)[, 1]
  g_transposed <- lapply(X = geno, FUN = t)
  ga <- do.call(abind::abind,c(g_transposed,list(along=0))) # 
  saveRDS(ga, "../data/genotypes_array.rds")
}
ga <- readRDS("../data/genotypes_array.rds")
```

``` r
##### From GK example code
# Specify allelic series prior
# Suggested by Wes
# Influences how much prior weight it places on more or less complicated allelic series
prior_M <- list(model.type = "crp", # crp - Chinese Restaurant Process
                prior.alpha.type = "gamma",
                prior.alpha.shape = 1,
                prior.alpha.rate = 2.333415)
```

``` r
if (!file.exists("../data/tnseq-traits.rds")) {
  load("../data/combined_hotspot_traits.Rdata")
  saveRDS(pheno, "../data/tnseq-traits.rds")
  pheno -> tnseq_traits
} else {
  tnseq_traits <- readRDS("../data/tnseq-traits.rds")
}
```

``` r
load("../data/reduced_map_qtl2_mapping_objects.Rdata")
```

``` r
data(mcv.data) # get A matrix
```

``` r
neto_list <- apply(FUN = as.list, X = neto, MARGIN = 1)
neto_small <- neto_list[1:3]
```

``` r
outfn <- "../data/timbr-tnseq-results-neto.rds"
# ensure that inputs to call_timbr all have subjects in same order!
subject_ids <- rownames(tnseq_traits)
indices_addcovar <- match(subject_ids, rownames(addcovar))
addcovar <- addcovar[indices_addcovar, ]
##
indices_gp <- match(subject_ids, rownames(ga))
gp <- ga[indices_gp, , ]
##
if (!file.exists(outfn)){
  timbr_out <- parallel::mclapply(neto_list, 
#  timbr_out <- lapply(neto_small, 
                                  FUN = qtl2effects::call_timbr, 
                                  mc.cores = parallel::detectCores(),
                                  traits_df = tnseq_traits,
                                  prior_M = prior_M, 
                                  genoprobs_array = gp,
                                  addcovar = addcovar
                                  )
  saveRDS(timbr_out, outfn)
}
```

``` r
timbr_out <- readRDS(outfn)
```
