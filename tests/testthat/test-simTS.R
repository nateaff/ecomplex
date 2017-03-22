tsuite <- function(){
  # The time serires models : arma, arma, farima, logistic, quadratic, Mackey-Glass
  farima  <- farima(ar = c(0.2, -0.4, 0.2), 
                    ma = c(0.4, 0.02), 
                    d  = 0.3)

  arma    <- arma(ar = c(0.4, 0.3), 
                  ma = c(0.1, -0.5))
  
  log     <- logistic(r = 3.70)
  mg      <- mackeyglass(tau   = 17, 
                     beta  = 0.35, 
                     gamma = 0.2,
                     n     = 8, 
                     init  = rep(0.2),
                     noise = 0.5)

  weiera   <- weierstrass(a = 0.8, b = 4)
  weierb  <- weierstrass(a = 0.8, b = 4, random = FALSE)
  weierc  <- weierstrass(a = 0.3)

  test_fs <- list(arma = arma, 
                  log = log, 
                  weiera = weiera,
                  weierb= weierb,
                  weierc = weierc,
                  mg = mg, 
                  farima = farima)
  return(test_fs)
}

test_that("functions output the correct length", {
   len    <- 500
   models <- tsuite()
   y      <- lapply(models, function(model) gen(model)(len))
   lens   <- unlist(lapply(y, length))
   names(lens) <- NULL
   expect_that(all(lens == len), is_true())
})


test_that("the one and two parameter weierstrass are the same", {
   len <- 500
   a <- 0.3; b <- 5
   aa <- -log(a)/log(b)   
   mod_a1 <- weierstrass(aa)
   mod_b1 <- weierstrass(aa,  random = FALSE)
   mod_c1 <- weierstrass(aa,  density = 100)
  
   mod_a2 <- weierstrass(a,b)
   mod_b2 <- weierstrass(a,b, random = FALSE)
   mod_c2 <- weierstrass(a,b, density = 100)
   
   mod3 <- weierstrass(a,b, density)
   set.seed(1)
   ya1 <- gen(mod_a1)(500)
   yb1 <- gen(mod_b1)(500)
   yc1 <- gen(mod_c1)(5)
   set.seed(1)
   ya2 <- gen(mod_a2)(500)
   yb2 <- gen(mod_b2)(500)
   yc2 <- gen(mod_c2)(5)
   
   expect_that(all(ya1 == ya2), is_true())
   expect_that(all(yb1 == yb2), is_true())
   expect_that(all(yc1 == yc2), is_true())
})







