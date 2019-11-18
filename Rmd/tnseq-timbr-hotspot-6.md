TIMBR: Allele effects in TnSeq Hotspots: One marker per trait
================
Frederick J. Boehm
11/15/2019

Last modified: 2019-11-17 19:00:31.

## Overview

We now consider only one marker per trait-hotspot pair. So, if a Neto
trait appears at more than one hotspot, it will be present more than
once below. However, if a trait is specific to a single hotspot, I
consider it at only one marker, its LOD peak marker within the hotspot.

``` r
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
library(TIMBR)
```

``` r
hots <- 1:10
```

``` r
probs <- readRDS("../data/genotypes_array.rds")
traits <- readRDS("../data/tnseq-traits.rds")
```

``` r
tt <- read.csv("../data/neto_traits_by_probe3_annotated.csv")
neto <- tt %>%
  tidyr::pivot_longer(cols = V2:V35, values_to = "trait", names_to = "trait_name") %>%
  dplyr::filter(!is.na(trait)) %>%
  dplyr::select(- trait_name, - neto.n, - row, - n.traits, - cM, -hs, - chr)
neto_plus <- tt %>%
  tidyr::pivot_longer(cols = V2:V35, values_to = "trait", names_to = "trait_name") %>%
  dplyr::filter(!is.na(trait)) %>%
  dplyr::select(- trait_name)

traits_timbr_annotated <- neto_plus %>%
  dplyr::select(probe, trait) %>%
  purrr::pmap( 
           .f = function(probe, trait){
             pheno <- traits[ , colnames(traits) == trait, drop = FALSE]
             geno <- probs[ , , dimnames(probs)[[3]] == probe]
             qtl2::fit1(genoprobs = geno, 
                        pheno = pheno, 
                        )
           }
             ) %>%
  purrr::map(.f = function(x){
    tibble::tibble(lod = x$lod)
  }) %>%
  bind_rows() %>%
  bind_cols(neto_plus) %>%
  dplyr::select(probe, trait, lod, hs) %>%
  dplyr::group_by(hs) %>%
  dplyr::group_by(trait) %>%
  dplyr::filter(lod == max(lod)) %>%
  dplyr::ungroup() %>%
  dplyr::ungroup()
```

``` r
(hot_indices <- neto_plus %>%
  dplyr::group_by(hs) %>%
  dplyr::tally() %>%
  dplyr::mutate(end = cumsum(n)) %>%
  dplyr::mutate(start = 1L + end - n)
)
```

    ## # A tibble: 10 x 4
    ##       hs     n   end start
    ##    <int> <int> <int> <int>
    ##  1     1     9     9     1
    ##  2     2    15    24    10
    ##  3     3    15    39    25
    ##  4     4  1905  1944    40
    ##  5     5    20  1964  1945
    ##  6     6  1265  3229  1965
    ##  7     7     8  3237  3230
    ##  8     8    80  3317  3238
    ##  9     9    66  3383  3318
    ## 10    10     2  3385  3384

## TIMBR setup

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
data(mcv.data) # get A matrix
```

``` r
tr_ann_sub <- traits_timbr_annotated %>%
  dplyr::filter(hs == params$hot) %>%
  dplyr::select(probe, trait)
neto_list <- apply(FUN = as.list, X = tr_ann_sub, MARGIN = 1)
neto_small <- neto_list[1:3]
```

``` r
outfn <- paste0("../data/timbr-tnseq-neto-traits-one-marker-per-trait-", params$hot, ".rds")
# ensure that inputs to call_timbr all have subjects in same order!
subject_ids <- rownames(traits)
##
indices_gp <- match(subject_ids, rownames(probs))
gp <- probs[indices_gp, , ]
##
if (!file.exists(outfn)){
  timbr_out <- parallel::mclapply(neto_list, 
  #timbr_out <- lapply(neto_small, 
                                  FUN = qtl2effects::call_timbr, 
                                  mc.cores = parallel::detectCores(),
                                  traits_df = traits,
                                  prior_M = prior_M, 
                                  genoprobs_array = gp,
                                  addcovar = NULL
                                  )
  saveRDS(timbr_out, outfn)
}
```
