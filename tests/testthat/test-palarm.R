

test_that("palarm does not return nans/NAs", {

  test1 <- palarm(1:10)
  test2 <- palarm(problem())
  test3 <- palarm(runif(200))
  test4 <- palarm(rnorm(100))

  expect_that(any(is.na(test1$kout), is_false())
  expect_that(any(is.na(test1$means), is_false())
  expect_that(any(is.na(test2$kout), is_false())
  expect_that(any(is.na(test2$means), is_false())
  expect_that(any(is.na(test3$kout), is_false())
  expect_that(any(is.na(test3$means), is_false())
  expect_that(any(is.na(test4$kout), is_false())
  expect_that(any(is.na(test4$means), is_false())
})

test_that("palarm throws the correct erros"){

 expect_error(palarm(0))
}

problem <- function(){
     c(0.82680105, 0.18305439, 0.00000000, 0.03675331, 0.92842346, 0.21954411,
       0.09442255, 0.26646925, 0.04393944, 0.32208523, 0.05664698, 0.13016508,
       0.06726348, 0.87707423, 0.34041029, 0.75410707, 0.88153270, 0.76655889,
       0.95112124, 0.76611597, 0.68941542, 0.73937119, 0.48055683, 0.20785725,
       0.11028071, 0.10255793, 0.06347670, 0.03794348, 0.12344860, 0.34931837,
       0.75859703, 0.80512195, 0.99893581, 0.68643471, 0.84358726, 0.70004257,
       0.42722178, 0.25891216, 1.00000000, 0.88030578, 0.78629402, 0.59354777,
       0.81619840, 0.35353501, 0.36002230, 0.33113299, 0.19551102, 0.08413024,
       0.13117171, 0.09356542, 0.35028482, 0.16133553, 0.11510450, 0.13768223,
      0.05022654, 0.08964114, 0.11203881, 0.09971774, 0.10054290, 0.04440902)
}