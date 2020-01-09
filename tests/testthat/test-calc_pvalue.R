library(testthat)
library(qtl2effects)

maxlods <- rep(c(0, 1), 10)

test_that("p-value is correctly calculated, considering discrete empirical distribution of max lods", {
  expect_equal(calc_pvalue(observed_lod = 0.9, maxlods = maxlods), 0.5)
  expect_equal(calc_pvalue(observed_lod = 0.1, maxlods = maxlods), 0.5)
  expect_equal(calc_pvalue(observed_lod = 2, maxlods = maxlods), 0)
  expect_equal(calc_pvalue(observed_lod = 1, maxlods = maxlods), 0.5)
  expect_equal(calc_pvalue(observed_lod = 0, maxlods = maxlods), 1)
})
