% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/call_timbr.R
\name{call_timbr}
\alias{call_timbr}
\title{Call timbr for use on cluster}
\usage{
call_timbr(
  input_vec,
  traits_df,
  addcovar = NULL,
  prior_M,
  genoprobs_array,
  A_matrix = mcv.data$prior.D$A
)
}
\arguments{
\item{input_vec}{a vector containing two named components, probe and trait}

\item{traits_df}{a data.frame with one trait per column and column name corresponding to trait name}

\item{addcovar}{additive covariates matrix to be appended to a column of 1's.}

\item{prior_M}{prior_M input for TIMBR::TIMBR()}

\item{genoprobs_array}{a n by 36 by m array of genotype probabilities, where n is the number of subjects and m the number of markers}
}
\value{
output (a named list) of a single call to TIMBR
}
\description{
Call timbr for use on cluster
}
