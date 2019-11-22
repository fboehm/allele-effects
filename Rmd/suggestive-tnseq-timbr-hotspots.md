TIMBR: Allele effects in TnSeq Hotspots: One marker per trait for all 10 hotspots
================
Frederick J. Boehm
11/21/2019

Last modified: 2019-11-22 09:26:27.

## Overview

We now consider only one marker per trait-hotspot pair. The set of
traits differs from those in previous analyses. We now read a csv file
that I made by hand from Rick Baker’s html file.

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
aprobs <- readRDS("../data/aprobs.rds")
map <- readRDS("../data/map.rds")
probs <- readRDS("../data/genotypes_array.rds")
traits <- readRDS("../data/tnseq-traits.rds")
```

``` r
trait_ann <- readr::read_csv("../data/v2_byhand_hotspot_tnseq_traits.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   hs = col_double(),
    ##   trait = col_character()
    ## )

## Run `scan1` on each trait

``` r
pheno <- traits[ , match(trait_ann$trait, colnames(traits)),
                 drop = FALSE]
s1out <- qtl2::scan1(aprobs, pheno)
# check if every trait is represented in table, and where it should be!
hot_chr <- c("4", "8", "9", "10", "11", "11", "13", "14", "18", "X") 
hot_index <- 1:10
hs_ann <- tibble::tibble(hs = hot_index, chr = hot_chr) %>%
  dplyr::right_join(trait_ann, by = "hs")

ann_all <- qtl2::find_peaks(s1out, map = map, threshold = 4) %>%
  tibble::as_tibble() %>%
  dplyr::filter(as.character(chr) %in% hot_chr) %>%
  dplyr::full_join(hs_ann, by = c("lodcolumn" = "trait", "chr")) %>%
  dplyr::filter(!is.na(hs)) %>%
  dplyr::arrange(hs) %>%
  dplyr::rename(trait = lodcolumn) %>%
  dplyr::select( - lodindex) %>%
  dplyr::left_join(tibble(marker = names(unlist(map)), pos = unlist(map)), by = "pos" ) %>%
  dplyr::mutate(splitted = stringr::str_split(marker, 
                                                    pattern = "\\."
  )
                                                     %>%
                          sapply(FUN = function(x)x[2])
                ) 
```

    ## Warning: Column `chr` joining factor and character vector, coercing into
    ## character vector

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
 neto_list <- ann_all %>%
  dplyr::rename(probe = splitted) %>%
  dplyr::select(probe, trait) %>%
  apply(FUN = as.list, MARGIN = 1)
```

``` r
outfn <- paste0("../data/suggestive-timbr-tnseq-neto-traits-one-marker-per-trait-all.rds")
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
} else {
  timbr_out <- readRDS(outfn)
}
names(timbr_out) <- trait_ann$trait
```

``` r
annots <- hs_ann %>%
  dplyr::mutate(allele_series = purrr::map(timbr_out, .f = function(x){
    pp <- x$p.M.given.y[1]
    return(tibble::tibble(series = names(pp), 
                   probability = pp
                   )
    )
    }
    )
                )
print(annots %>% tidyr::unnest(allele_series))
```

    ## # A tibble: 69 x 5
    ##       hs chr   trait      series          probability
    ##    <dbl> <chr> <chr>      <chr>                 <dbl>
    ##  1     1 4     RVBD_0309  0,0,1,0,1,1,1,1      0.305 
    ##  2     1 4     RVBD_1204c 0,0,0,0,0,0,0,0      0.116 
    ##  3     1 4     RVBD_3127  0,1,1,0,0,0,0,1      0.282 
    ##  4     2 8     RVBD_0789c 0,1,0,0,0,1,1,1      0.178 
    ##  5     2 8     RVBD_0901  0,1,2,1,2,0,1,1      0.0881
    ##  6     2 8     RVBD_2161c 0,1,0,1,0,1,0,1      0.115 
    ##  7     2 8     RVBD_2466c 0,0,0,0,0,0,0,0      0.221 
    ##  8     2 8     RVBD_3208A 0,1,0,0,0,0,0,0      0.145 
    ##  9     2 8     RVBD_3236c 0,0,0,0,0,1,0,0      0.265 
    ## 10     2 8     RVBD_3394c 0,0,1,0,2,2,0,0      0.0716
    ## # … with 59 more rows

The annotations for each plot appear before the plot.

``` r
for (i in seq_along(timbr_out)){
  print(paste0(annots[i, ], collapse = " "))
  TIMBR::TIMBR.plot.haplotypes(timbr_out[[i]])
}
```

    ## [1] "1 4 RVBD_0309 list(RVBD_0309 = list(series = \"0,0,1,0,1,1,1,1\", probability = 0.3046))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-1.png)<!-- -->

    ## [1] "1 4 RVBD_1204c list(RVBD_1204c = list(series = \"0,0,0,0,0,0,0,0\", probability = 0.1159))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2.png)<!-- -->

    ## [1] "1 4 RVBD_3127 list(RVBD_3127 = list(series = \"0,1,1,0,0,0,0,1\", probability = 0.2815))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-3.png)<!-- -->

    ## [1] "2 8 RVBD_0789c list(RVBD_0789c = list(series = \"0,1,0,0,0,1,1,1\", probability = 0.178))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-4.png)<!-- -->

    ## [1] "2 8 RVBD_0901 list(RVBD_0901 = list(series = \"0,1,2,1,2,0,1,1\", probability = 0.0881))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-5.png)<!-- -->

    ## [1] "2 8 RVBD_2161c list(RVBD_2161c = list(series = \"0,1,0,1,0,1,0,1\", probability = 0.1151))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-6.png)<!-- -->

    ## [1] "2 8 RVBD_2466c list(RVBD_2466c = list(series = \"0,0,0,0,0,0,0,0\", probability = 0.2206))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-7.png)<!-- -->

    ## [1] "2 8 RVBD_3208A list(RVBD_3208A = list(series = \"0,1,0,0,0,0,0,0\", probability = 0.1454))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-8.png)<!-- -->

    ## [1] "2 8 RVBD_3236c list(RVBD_3236c = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.2647))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-9.png)<!-- -->

    ## [1] "2 8 RVBD_3394c list(RVBD_3394c = list(series = \"0,0,1,0,2,2,0,0\", probability = 0.0716))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-10.png)<!-- -->

    ## [1] "3 9 RVBD_0490 list(RVBD_0490 = list(series = \"0,1,0,0,0,0,0,1\", probability = 0.2454))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-11.png)<!-- -->

    ## [1] "3 9 RVBD_0491 list(RVBD_0491 = list(series = \"0,0,0,0,0,0,0,1\", probability = 0.1017))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-12.png)<!-- -->

    ## [1] "3 9 RVBD_0863 list(RVBD_0863 = list(series = \"0,0,0,0,0,0,0,0\", probability = 0.1721))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-13.png)<!-- -->

    ## [1] "3 9 RVBD_0961 list(RVBD_0961 = list(series = \"0,0,0,0,0,0,0,0\", probability = 0.3301))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-14.png)<!-- -->

    ## [1] "3 9 RVBD_1405c list(RVBD_1405c = list(series = \"0,1,0,0,0,0,0,1\", probability = 0.2474))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-15.png)<!-- -->

    ## [1] "3 9 RVBD_2257c list(RVBD_2257c = list(series = \"0,1,0,0,1,0,0,0\", probability = 0.1408))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-16.png)<!-- -->

    ## [1] "3 9 RVBD_2306A list(RVBD_2306A = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.3595))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-17.png)<!-- -->

    ## [1] "3 9 RVBD_3050c list(RVBD_3050c = list(series = \"0,0,0,1,0,1,0,0\", probability = 0.2881))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-18.png)<!-- -->

    ## [1] "3 9 RVBD_3813c list(RVBD_3813c = list(series = \"0,0,1,1,0,0,0,1\", probability = 0.108))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-19.png)<!-- -->

    ## [1] "3 9 RVBD_3881c list(RVBD_3881c = list(series = \"0,0,1,0,0,0,0,1\", probability = 0.2484))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-20.png)<!-- -->

    ## [1] "4 10 RVBD_0158 list(RVBD_0158 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.3725))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-21.png)<!-- -->

    ## [1] "4 10 RVBD_0437c list(RVBD_0437c = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.501))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-22.png)<!-- -->

    ## [1] "4 10 RVBD_0490 list(RVBD_0490 = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.3183))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-23.png)<!-- -->

    ## [1] "4 10 RVBD_1051c list(RVBD_1051c = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.2183))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-24.png)<!-- -->

    ## [1] "4 10 RVBD_1493 list(RVBD_1493 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.5663))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-25.png)<!-- -->

    ## [1] "4 10 RVBD_1692 list(RVBD_1692 = list(series = \"0,0,0,0,0,0,0,0\", probability = 0.2751))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-26.png)<!-- -->

    ## [1] "4 10 RVBD_2643A list(RVBD_2643A = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.114))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-27.png)<!-- -->

    ## [1] "4 10 RVBD_2854 list(RVBD_2854 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.6451))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-28.png)<!-- -->

    ## [1] "4 10 RVBD_3261 list(RVBD_3261 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.6531))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-29.png)<!-- -->

    ## [1] "4 10 RVBD_3262 list(RVBD_3262 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.3769))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-30.png)<!-- -->

    ## [1] "4 10 RVBD_3452 list(RVBD_3452 = list(series = \"0,1,1,1,2,1,1,2\", probability = 0.0949))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-31.png)<!-- -->

    ## [1] "5 11 RVBD_1103c list(RVBD_1103c = list(series = \"0,0,0,1,0,0,0,1\", probability = 0.1459))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-32.png)<!-- -->

    ## [1] "5 11 RVBD_2267c list(RVBD_2267c = list(series = \"0,1,2,1,1,1,1,1\", probability = 0.2302))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-33.png)<!-- -->

    ## [1] "5 11 RVBD_3658c list(RVBD_3658c = list(series = \"0,1,0,1,1,1,0,1\", probability = 0.1945))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-34.png)<!-- -->

    ## [1] "6 11 RVBD_0806c list(RVBD_0806c = list(series = \"0,0,1,1,0,1,0,1\", probability = 0.2016))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-35.png)<!-- -->

    ## [1] "6 11 RVBD_1782 list(RVBD_1782 = list(series = \"0,0,1,0,0,0,0,0\", probability = 0.3823))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-36.png)<!-- -->

    ## [1] "6 11 RVBD_2281 list(RVBD_2281 = list(series = \"0,0,1,0,0,0,0,1\", probability = 0.2244))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-37.png)<!-- -->

    ## [1] "7 13 RVBD_0740 list(RVBD_0740 = list(series = \"0,1,1,0,1,0,1,0\", probability = 0.3166))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-38.png)<!-- -->

    ## [1] "7 13 RVBD_1405c list(RVBD_1405c = list(series = \"0,1,1,1,0,1,1,1\", probability = 0.2396))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-39.png)<!-- -->

    ## [1] "7 13 RVBD_1793 list(RVBD_1793 = list(series = \"0,1,1,0,1,0,1,1\", probability = 0.3498))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-40.png)<!-- -->

    ## [1] "7 13 RVBD_1821 list(RVBD_1821 = list(series = \"0,0,1,0,0,0,0,0\", probability = 0.4276))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-41.png)<!-- -->

    ## [1] "7 13 RVBD_2039c list(RVBD_2039c = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.167))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-42.png)<!-- -->

    ## [1] "7 13 RVBD_2499c list(RVBD_2499c = list(series = \"0,0,1,1,1,1,1,1\", probability = 0.2194))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-43.png)<!-- -->

    ## [1] "7 13 RVBD_2515c list(RVBD_2515c = list(series = \"0,0,0,0,1,0,1,0\", probability = 0.0762))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-44.png)<!-- -->

    ## [1] "7 13 RVBD_3516 list(RVBD_3516 = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.8822))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-45.png)<!-- -->

    ## [1] "8 14 RVBD_0114 list(RVBD_0114 = list(series = \"0,0,1,1,1,1,1,0\", probability = 0.3475))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-46.png)<!-- -->

    ## [1] "8 14 RVBD_0590A list(RVBD_0590A = list(series = \"0,0,1,1,0,0,1,1\", probability = 0.2113))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-47.png)<!-- -->

    ## [1] "8 14 RVBD_1486c list(RVBD_1486c = list(series = \"0,0,1,1,0,1,0,1\", probability = 0.1761))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-48.png)<!-- -->

    ## [1] "8 14 RVBD_2284 list(RVBD_2284 = list(series = \"0,0,1,0,0,1,1,0\", probability = 0.4238))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-49.png)<!-- -->

    ## [1] "8 14 RVBD_2499c list(RVBD_2499c = list(series = \"0,0,1,0,0,0,0,0\", probability = 0.0645))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-50.png)<!-- -->

    ## [1] "8 14 RVBD_3369 list(RVBD_3369 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.4655))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-51.png)<!-- -->

    ## [1] "9 18 RVBD_0492c list(RVBD_0492c = list(series = \"0,1,1,1,1,1,1,1\", probability = 0.066))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-52.png)<!-- -->

    ## [1] "9 18 RVBD_0746 list(RVBD_0746 = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.5262))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-53.png)<!-- -->

    ## [1] "9 18 RVBD_0971c list(RVBD_0971c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.4842))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-54.png)<!-- -->

    ## [1] "9 18 RVBD_1705c list(RVBD_1705c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.4741))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-55.png)<!-- -->

    ## [1] "9 18 RVBD_2161c list(RVBD_2161c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.5739))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-56.png)<!-- -->

    ## [1] "9 18 RVBD_2318 list(RVBD_2318 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.2708))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-57.png)<!-- -->

    ## [1] "9 18 RVBD_2372c list(RVBD_2372c = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.6527))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-58.png)<!-- -->

    ## [1] "9 18 RVBD_2490c list(RVBD_2490c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.2619))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-59.png)<!-- -->

    ## [1] "9 18 RVBD_2949c list(RVBD_2949c = list(series = \"0,1,0,0,0,0,0,1\", probability = 0.2354))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-60.png)<!-- -->

    ## [1] "9 18 RVBD_2952 list(RVBD_2952 = list(series = \"0,0,0,1,0,1,0,0\", probability = 0.213))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-61.png)<!-- -->

    ## [1] "9 18 RVBD_3590c list(RVBD_3590c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.6815))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-62.png)<!-- -->

    ## [1] "9 18 RVBD_3653 list(RVBD_3653 = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.5678))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-63.png)<!-- -->

    ## [1] "9 18 RVBD_3848 list(RVBD_3848 = list(series = \"0,1,0,0,1,1,1,0\", probability = 0.1776))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-64.png)<!-- -->

    ## [1] "10 X RVBD_0120c list(RVBD_0120c = list(series = \"0,1,1,1,0,0,0,0\", probability = 0.1232))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-65.png)<!-- -->

    ## [1] "10 X RVBD_0501 list(RVBD_0501 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.1183))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-66.png)<!-- -->

    ## [1] "10 X RVBD_0502 list(RVBD_0502 = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.3994))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-67.png)<!-- -->

    ## [1] "10 X RVBD_1677 list(RVBD_1677 = list(series = \"0,0,0,1,1,0,1,0\", probability = 0.043))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-68.png)<!-- -->

    ## [1] "10 X RVBD_1850 list(RVBD_1850 = list(series = \"0,0,0,0,0,0,0,0\", probability = 0.4063))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-69.png)<!-- -->

``` r
readr::write_csv(annots %>% tidyr::unnest(allele_series), path = "../data/suggestive-tnseq-neto-traits-annotated-10-hotspots.csv")
```

## Use the other prior for allelic series

``` r
##### From GK example code
# Specify allelic series prior
# Suggested by Wes
# Influences how much prior weight it places on more or less complicated allelic series
prior_M <- list(model.type = "crp", # crp - Chinese Restaurant Process
                prior.alpha.type = "gamma",
                prior.alpha.shape = 2.3009322,
                prior.alpha.rate = 0.7488104)
```

``` r
outfn <- "../data/suggestive-timbr-tnseq-prior2.rds"
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
} else {
  timbr_out <- readRDS(outfn)
}
names(timbr_out) <- trait_ann$trait
```

``` r
annots <- hs_ann %>%
  dplyr::mutate(allele_series = purrr::map(timbr_out, .f = function(x){
    pp <- x$p.M.given.y[1]
    return(tibble::tibble(series = names(pp), 
                   probability = pp
                   )
    )
    }
    )
                )
print(annots %>% tidyr::unnest(allele_series))
```

    ## # A tibble: 69 x 5
    ##       hs chr   trait      series          probability
    ##    <dbl> <chr> <chr>      <chr>                 <dbl>
    ##  1     1 4     RVBD_0309  0,0,1,0,1,1,1,1      0.0675
    ##  2     1 4     RVBD_1204c 0,0,1,0,0,0,0,2      0.0467
    ##  3     1 4     RVBD_3127  0,1,1,0,0,2,0,1      0.0713
    ##  4     2 8     RVBD_0789c 0,1,0,0,0,1,1,1      0.0459
    ##  5     2 8     RVBD_0901  0,1,2,1,2,0,1,1      0.0227
    ##  6     2 8     RVBD_2161c 0,1,2,3,4,5,6,7      0.017 
    ##  7     2 8     RVBD_2466c 0,0,0,1,2,0,0,0      0.0409
    ##  8     2 8     RVBD_3208A 0,1,0,0,0,0,0,0      0.0238
    ##  9     2 8     RVBD_3236c 0,0,0,0,0,1,0,0      0.0352
    ## 10     2 8     RVBD_3394c 0,0,1,0,2,3,0,0      0.0243
    ## # … with 59 more rows

The annotations for each plot appear before the plot.

``` r
for (i in seq_along(timbr_out)){
  print(paste0(annots[i, ], collapse = " "))
  TIMBR::TIMBR.plot.haplotypes(timbr_out[[i]])
}
```

    ## [1] "1 4 RVBD_0309 list(RVBD_0309 = list(series = \"0,0,1,0,1,1,1,1\", probability = 0.0675))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-1.png)<!-- -->

    ## [1] "1 4 RVBD_1204c list(RVBD_1204c = list(series = \"0,0,1,0,0,0,0,2\", probability = 0.0467))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-2.png)<!-- -->

    ## [1] "1 4 RVBD_3127 list(RVBD_3127 = list(series = \"0,1,1,0,0,2,0,1\", probability = 0.0713))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-3.png)<!-- -->

    ## [1] "2 8 RVBD_0789c list(RVBD_0789c = list(series = \"0,1,0,0,0,1,1,1\", probability = 0.0459))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-4.png)<!-- -->

    ## [1] "2 8 RVBD_0901 list(RVBD_0901 = list(series = \"0,1,2,1,2,0,1,1\", probability = 0.0227))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-5.png)<!-- -->

    ## [1] "2 8 RVBD_2161c list(RVBD_2161c = list(series = \"0,1,2,3,4,5,6,7\", probability = 0.017))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-6.png)<!-- -->

    ## [1] "2 8 RVBD_2466c list(RVBD_2466c = list(series = \"0,0,0,1,2,0,0,0\", probability = 0.0409))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-7.png)<!-- -->

    ## [1] "2 8 RVBD_3208A list(RVBD_3208A = list(series = \"0,1,0,0,0,0,0,0\", probability = 0.0238))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-8.png)<!-- -->

    ## [1] "2 8 RVBD_3236c list(RVBD_3236c = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.0352))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-9.png)<!-- -->

    ## [1] "2 8 RVBD_3394c list(RVBD_3394c = list(series = \"0,0,1,0,2,3,0,0\", probability = 0.0243))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-10.png)<!-- -->

    ## [1] "3 9 RVBD_0490 list(RVBD_0490 = list(series = \"0,1,0,0,0,0,0,2\", probability = 0.0933))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-11.png)<!-- -->

    ## [1] "3 9 RVBD_0491 list(RVBD_0491 = list(series = \"0,1,0,0,1,1,0,2\", probability = 0.0358))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-12.png)<!-- -->

    ## [1] "3 9 RVBD_0863 list(RVBD_0863 = list(series = \"0,0,0,1,0,0,2,0\", probability = 0.0268))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-13.png)<!-- -->

    ## [1] "3 9 RVBD_0961 list(RVBD_0961 = list(series = \"0,0,0,0,1,0,2,0\", probability = 0.0203))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-14.png)<!-- -->

    ## [1] "3 9 RVBD_1405c list(RVBD_1405c = list(series = \"0,1,0,0,0,0,0,1\", probability = 0.0387))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-15.png)<!-- -->

    ## [1] "3 9 RVBD_2257c list(RVBD_2257c = list(series = \"0,1,0,0,2,0,0,0\", probability = 0.0229))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-16.png)<!-- -->

    ## [1] "3 9 RVBD_2306A list(RVBD_2306A = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.0584))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-17.png)<!-- -->

    ## [1] "3 9 RVBD_3050c list(RVBD_3050c = list(series = \"0,0,0,1,0,1,0,0\", probability = 0.058))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-18.png)<!-- -->

    ## [1] "3 9 RVBD_3813c list(RVBD_3813c = list(series = \"0,0,1,2,0,0,0,1\", probability = 0.0312))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-19.png)<!-- -->

    ## [1] "3 9 RVBD_3881c list(RVBD_3881c = list(series = \"0,0,1,0,0,0,0,1\", probability = 0.0676))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-20.png)<!-- -->

    ## [1] "4 10 RVBD_0158 list(RVBD_0158 = list(series = \"0,1,1,1,2,1,1,1\", probability = 0.1228))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-21.png)<!-- -->

    ## [1] "4 10 RVBD_0437c list(RVBD_0437c = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.1302))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-22.png)<!-- -->

    ## [1] "4 10 RVBD_0490 list(RVBD_0490 = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.0505))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-23.png)<!-- -->

    ## [1] "4 10 RVBD_1051c list(RVBD_1051c = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.0331))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-24.png)<!-- -->

    ## [1] "4 10 RVBD_1493 list(RVBD_1493 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.1538))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-25.png)<!-- -->

    ## [1] "4 10 RVBD_1692 list(RVBD_1692 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.089))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-26.png)<!-- -->

    ## [1] "4 10 RVBD_2643A list(RVBD_2643A = list(series = \"0,1,0,0,2,0,0,0\", probability = 0.0221))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-27.png)<!-- -->

    ## [1] "4 10 RVBD_2854 list(RVBD_2854 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.2318))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-28.png)<!-- -->

    ## [1] "4 10 RVBD_3261 list(RVBD_3261 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.2126))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-29.png)<!-- -->

    ## [1] "4 10 RVBD_3262 list(RVBD_3262 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.0678))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-30.png)<!-- -->

    ## [1] "4 10 RVBD_3452 list(RVBD_3452 = list(series = \"0,1,2,3,4,5,6,7\", probability = 0.033))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-31.png)<!-- -->

    ## [1] "5 11 RVBD_1103c list(RVBD_1103c = list(series = \"0,0,0,1,0,2,0,3\", probability = 0.0312))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-32.png)<!-- -->

    ## [1] "5 11 RVBD_2267c list(RVBD_2267c = list(series = \"0,1,2,3,1,1,1,1\", probability = 0.1225))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-33.png)<!-- -->

    ## [1] "5 11 RVBD_3658c list(RVBD_3658c = list(series = \"0,1,0,1,1,1,0,1\", probability = 0.0248))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-34.png)<!-- -->

    ## [1] "6 11 RVBD_0806c list(RVBD_0806c = list(series = \"0,0,1,1,0,1,0,1\", probability = 0.038))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-35.png)<!-- -->

    ## [1] "6 11 RVBD_1782 list(RVBD_1782 = list(series = \"0,0,1,0,0,0,0,0\", probability = 0.1946))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-36.png)<!-- -->

    ## [1] "6 11 RVBD_2281 list(RVBD_2281 = list(series = \"0,0,1,0,0,0,0,2\", probability = 0.0356))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-37.png)<!-- -->

    ## [1] "7 13 RVBD_0740 list(RVBD_0740 = list(series = \"0,1,1,0,1,0,1,0\", probability = 0.0481))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-38.png)<!-- -->

    ## [1] "7 13 RVBD_1405c list(RVBD_1405c = list(series = \"0,1,1,1,0,1,1,1\", probability = 0.0402))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-39.png)<!-- -->

    ## [1] "7 13 RVBD_1793 list(RVBD_1793 = list(series = \"0,1,1,0,1,0,1,1\", probability = 0.0654))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-40.png)<!-- -->

    ## [1] "7 13 RVBD_1821 list(RVBD_1821 = list(series = \"0,0,1,0,0,0,0,0\", probability = 0.0795))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-41.png)<!-- -->

    ## [1] "7 13 RVBD_2039c list(RVBD_2039c = list(series = \"0,0,1,2,0,0,0,0\", probability = 0.0316))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-42.png)<!-- -->

    ## [1] "7 13 RVBD_2499c list(RVBD_2499c = list(series = \"0,0,1,1,1,1,1,1\", probability = 0.0415))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-43.png)<!-- -->

    ## [1] "7 13 RVBD_2515c list(RVBD_2515c = list(series = \"0,1,2,3,4,5,6,7\", probability = 0.0264))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-44.png)<!-- -->

    ## [1] "7 13 RVBD_3516 list(RVBD_3516 = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.6452))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-45.png)<!-- -->

    ## [1] "8 14 RVBD_0114 list(RVBD_0114 = list(series = \"0,0,1,1,1,1,1,0\", probability = 0.0657))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-46.png)<!-- -->

    ## [1] "8 14 RVBD_0590A list(RVBD_0590A = list(series = \"0,0,1,1,0,0,1,1\", probability = 0.0391))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-47.png)<!-- -->

    ## [1] "8 14 RVBD_1486c list(RVBD_1486c = list(series = \"0,0,1,1,2,1,0,1\", probability = 0.0185))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-48.png)<!-- -->

    ## [1] "8 14 RVBD_2284 list(RVBD_2284 = list(series = \"0,0,1,0,0,1,1,0\", probability = 0.076))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-49.png)<!-- -->

    ## [1] "8 14 RVBD_2499c list(RVBD_2499c = list(series = \"0,0,1,0,0,0,0,2\", probability = 0.0237))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-50.png)<!-- -->

    ## [1] "8 14 RVBD_3369 list(RVBD_3369 = list(series = \"0,0,0,0,0,0,1,0\", probability = 0.1047))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-51.png)<!-- -->

    ## [1] "9 18 RVBD_0492c list(RVBD_0492c = list(series = \"0,1,2,3,4,5,6,7\", probability = 0.0232))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-52.png)<!-- -->

    ## [1] "9 18 RVBD_0746 list(RVBD_0746 = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.1526))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-53.png)<!-- -->

    ## [1] "9 18 RVBD_0971c list(RVBD_0971c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.1285))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-54.png)<!-- -->

    ## [1] "9 18 RVBD_1705c list(RVBD_1705c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.1231))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-55.png)<!-- -->

    ## [1] "9 18 RVBD_2161c list(RVBD_2161c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.2362))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-56.png)<!-- -->

    ## [1] "9 18 RVBD_2318 list(RVBD_2318 = list(series = \"0,0,0,0,1,0,0,0\", probability = 0.0574))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-57.png)<!-- -->

    ## [1] "9 18 RVBD_2372c list(RVBD_2372c = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.2441))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-58.png)<!-- -->

    ## [1] "9 18 RVBD_2490c list(RVBD_2490c = list(series = \"0,0,0,0,1,2,2,0\", probability = 0.048))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-59.png)<!-- -->

    ## [1] "9 18 RVBD_2949c list(RVBD_2949c = list(series = \"0,1,0,0,0,0,0,1\", probability = 0.0479))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-60.png)<!-- -->

    ## [1] "9 18 RVBD_2952 list(RVBD_2952 = list(series = \"0,0,0,1,0,1,0,0\", probability = 0.0336))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-61.png)<!-- -->

    ## [1] "9 18 RVBD_3590c list(RVBD_3590c = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.281))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-62.png)<!-- -->

    ## [1] "9 18 RVBD_3653 list(RVBD_3653 = list(series = \"0,0,0,0,0,1,0,0\", probability = 0.1904))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-63.png)<!-- -->

    ## [1] "9 18 RVBD_3848 list(RVBD_3848 = list(series = \"0,1,0,0,1,1,1,0\", probability = 0.0195))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-64.png)<!-- -->

    ## [1] "10 X RVBD_0120c list(RVBD_0120c = list(series = \"0,1,1,1,0,0,0,0\", probability = 0.0241))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-65.png)<!-- -->

    ## [1] "10 X RVBD_0501 list(RVBD_0501 = list(series = \"0,1,2,3,4,5,6,7\", probability = 0.0219))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-66.png)<!-- -->

    ## [1] "10 X RVBD_0502 list(RVBD_0502 = list(series = \"0,0,0,0,0,1,1,0\", probability = 0.1207))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-67.png)<!-- -->

    ## [1] "10 X RVBD_1677 list(RVBD_1677 = list(series = \"0,1,2,3,4,5,6,7\", probability = 0.0295))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-68.png)<!-- -->

    ## [1] "10 X RVBD_1850 list(RVBD_1850 = list(series = \"0,0,0,1,0,0,0,0\", probability = 0.0495))"

![](suggestive-tnseq-timbr-hotspots_files/figure-gfm/plots-2-69.png)<!-- -->

``` r
fn <- "../data/suggestive-tnseq-prior2.csv"
readr::write_csv(annots %>% tidyr::unnest(allele_series), path = fn)
```
