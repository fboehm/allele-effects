install.packages(c("devtools", "tidyverse"))

devtools::install_github("fboehm/qtl2effects")

library(qtl2effects)
library(tidyverse)

rmarkdown::render("../Rmd/tnseq-timbr.Rmd") # contains call to parallel::mclapply






