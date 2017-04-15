context('ecomplex')

f <- function(){
  xx <- seq(0,15, 1/20)
  y1 <- sin(2*xx) + (1/3)*cos((.33*xx)) + exp(-xx) 
  y2 <- .8*cos(12*xx) + .5*cos(9*xx) + exp(-xx)
  ys <- normalize(y1 + y2)
}

test_that("ecomplex does not return NA's", {
  x <- rnorm(128)
  x <- f()
  fit1 <- ecomplex(x, ds = 6, max_degree = 5, method = "bspline")
  fit2 <- ecomplex(x, ds = 6, max_degree = 5, method = "cspline")
  fit3 <- ecomplex(x, ds = 6, max_degree = 5, method = "lift")
  
  expect_that(all(!is.na(fit1$fit)), is_true())
  expect_that(all(!is.na(fit2$fit)), is_true())
  expect_that(all(!is.na(fit3$fit)), is_true())

})


test_that("ecomplex throws appropriate errors or warnings for malformed data", {
  err1 <- c(1, 2, "c")
  err2 <- matrix(rep(1, 20))
  err3 <- data.frame(rep(1, 20))
  err4 <- c(1,2, NA)
  # errors  
  expect_error(ecomplex(err1))
  expect_error(ecomplex(err2))
  expect_error(ecomplex(err3))
  expect_error(ecomplex(err4))
  # warnings
  err5 <- c(1:101, Inf)
  # Length warning
  err6 <- c(1:20)
  expect_warning(ecomplex(err5))
  expect_warning(ecomplex(err6))
  # silent
  err7 <- c(1:101, NULL)
  expect_silent(ecomplex(err7))
})
