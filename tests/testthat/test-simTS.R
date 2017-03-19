tsuite <- function(){
  # The time serires models : arma, arma, farima, logistic, quadratic, Mackey-Glass
  farima  <- farima(ar = c(0.2, -0.4), 
                    ma = c(0.4, 0.02), 
                    d  = 0.3)

  arma    <- arma(ar = c(0.4, 0.3, 0.2), 
                  ma = c(0.1, -0.5))
  
  log     <- logistic(r = 3.70)
  mg      <- mackeyglass(tau   = 17, 
                     beta  = 0.35, 
                     gamma = 0.2,
                     n     = 8, 
                     init  = rep(0.2),
                     noise = 0.5)

  weier   <- weierstrass(a = 0.8, b = 4)
  weiera  <- weierstrass_a(a = 0.4)
  weier_ran   <- weierstrass(a = 0.8, b = 4, random = FALSE)
  weiera_ran  <- weierstrass_a(a = 0.4, random = FALSE)

  test_fs <- list(arma = arma, 
                  log = log, 
                  weier = weier,
                  weiera = weiera,
                  weier_ran = weier_ran,
                  weiera_ran = weiera_ran, 
                  mg = mg, 
                  farima1 = farima)
  return(test_fs)
}

test_that("functions output the correct length", {
   len <- 500
   models <- tsuite()
   ys <- lapply(models, function(model) gen(model)(len))
   lens <- unlist(lapply(ys, length))
   names(lens) <- NULL
   expect_that(all.equal(lens, rep(500, length(lens))), is_true())
})








