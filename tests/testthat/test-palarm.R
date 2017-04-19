context('palarm')

problem <- function(){
     c(0.82680105, 0.18305439, 0.00000000, 0.03675331, 0.92842346, 0.21954411,
       0.09442255, 0.26646925, 0.04393944, 0.32208523, 0.05664698, 0.13016508,
       0.06726348, 0.87707423, 0.34041029, 0.75410707, 0.88153270, 0.76655889,
       0.95112124, 0.76611597, 0.68941542, 0.73937119, 0.48055683, 0.20785725,
       0.11028071, 0.10255793, 0.06347670, 0.03794348)
}

test_that("palarm does not return nans/NAs", {

  test1 <- palarm(1:10)
  test2 <- palarm(problem())
  test3 <- palarm(runif(200))
  test4 <- palarm(rnorm(100))

  expect_that(any(is.na(test1$kout)), is_false())
  expect_that(any(is.na(test1$means)), is_false())
  expect_that(any(is.na(test2$kout)), is_false())
  expect_that(any(is.na(test2$means)), is_false())
  expect_that(any(is.na(test3$kout)), is_false())
  expect_that(any(is.na(test3$means)), is_false())
  expect_that(any(is.na(test4$kout)), is_false())
  expect_that(any(is.na(test4$means)), is_false())
})

test_that("palarm throws the correct errors",{
  expect_error(palarm(0))
})

