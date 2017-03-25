context('features')

test_that("feature output is the correct length", {
  f <- seq(0,1, length.out = 200)

  elift_feat <- clean_feature(ecomp_lift(f))
  expect_that(all.equal( dim(elift_feat),c(1,2) ), is_true())

  var_feat <- clean_feature(variance(f))
  expect_that(all.equal( dim(var_feat), c(1,1) ), is_true())

  hurst_feat <- clean_feature(hurst(f))
  expect_that(all.equal( dim(hurst_feat), c(1,1) ), is_true())

  perm_feat <- clean_feature(permutation_entropy(f))
  expect_that(all.equal( dim(perm_feat), c(1,1) ), is_true())

  band_feat <- clean_feature(bandpower(f))
  expect_that(all.equal( dim(band_feat), c(1,5) ), is_true())

  wvar_feat <- clean_feature(gmwm::wvar(f))
  expect_that(all.equal( dim(wvar_feat), c(1,4) ), is_true())

  sen_feat <- clean_feature(sample_entropy(f))
  expect_that(all.equal( dim(sen_feat), c(1,1) ), is_true())

  spec_feat <- clean_feature(spectral_entropy(f))
  expect_that(all.equal( dim(spec_feat), c(1,1) ), is_true())

})


test_that("feature output is numeric", {
  f <- seq(0,1, length.out = 200)

  elift_feat <- clean_feature(ecomp_lift(f))
  expect_that(all(apply(elift_feat, 2, is.numeric)), is_true())

  var_feat <- clean_feature(variance(f))
  expect_that(all(apply(var_feat, 2, is.numeric)), is_true())

  hurst_feat <- clean_feature(hurst(f))
  expect_that(all(apply(hurst_feat, 2, is.numeric)), is_true())

  perm_feat <- clean_feature(permutation_entropy(f))
  expect_that(all(apply(perm_feat, 2, is.numeric)), is_true())


  band_feat <- clean_feature(bandpower(f))
  expect_that(all(apply(band_feat, 2, is.numeric)), is_true())

  wvar_feat <- clean_feature(gmwm::wvar(f))
  expect_that(all(apply(wvar_feat, 2, is.numeric)), is_true())

  sen_feat <- clean_feature(sample_entropy(f))
  expect_that(all(apply(sen_feat, 2, is.numeric)), is_true())

  spec_feat <- clean_feature(spectral_entropy(f))
  expect_that(all(apply(spec_feat, 2, is.numeric)), is_true())

})





