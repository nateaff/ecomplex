context('sample')

test_that("random downsampling samples the correct proportion", {
  n = 30 
  ds = 6
  rand_indices <- random_sample(n, ds)
  down_indices <- downsample_perm(n, ds)

  rand_lens <- unlist(lapply(rand_indices, length))
  down_lens <- unlist(lapply(down_indices, length))

  indices_sorted <- !unlist(lapply(rand_indices, is.unsorted))

  expect_that(length(rand_indices) == length(down_indices), is_true())
  expect_that(all(rand_lens == down_lens), is_true())
  expect_that(all(indices_sorted == TRUE), is_true())
})

