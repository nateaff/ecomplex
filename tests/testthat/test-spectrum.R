#============================================
# Test spectral features
#============================================

test_that("bp_spec returns correct data for vector input",{
  fs <- 1000
  freq <- c(1,20)
  freqs <- list(freq, c(20, 40))
  x <- rnorm(1000)
 
  expect_that(is.list(bp_pgram(x, fs, freq)), is_true())
  expect_that(length(bp_pgram(x, fs, freqs)), equals(length(freqs)))
})
