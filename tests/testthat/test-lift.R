check_error <- function(ys, ts, omit, j){
    yrange <- (omit+j):(length(ys) - omit)
    trange <- (omit):(length(ys)-omit -j)
     abs(max(ts[trange] - ys[yrange]))
} 

f <- function(){
  xx <- seq(0,100, 1/20)
  y1 <- sin(2*xx) + (1/3)*cos((.33*xx)) + exp(-xx) 
  y2 <- .8*cos(12*xx) + .5*cos(9*xx) + exp(-xx)
  ys <- normalize(y1 + y2)
}


test_that("merge functions maintain correct order", {
   # test negative
   comp <- 1:100
   odd1 <- seq(1, 100, by = 2)
   even1 <- comp[-odd1]
   res1 <- merge1(odd1, even1)

   odd2 <- seq(1, 100, by = 3)
   even2 <- comp[-odd2]
   res2 <- merge2(odd2, even2)

   odd3 <- floor((5*(1:60)-1)/3)
   even3 <- comp[-odd3]
   res3 <- merge3(odd3, even3, 2)

   expect_that(all.equal(res1, comp), is_true())
   expect_that(all.equal(res2, comp), is_true())   
   expect_that(all.equal(res3, comp), is_true())   
   
})


test_that("inverseWT4 iterated interpolation error is below tolerance",{
  ys   <- f()
  ds   <- 4
  degs <- c(1,3, 5)
  tol  <- 0.3;   
  omit <- 10
  errs  <- double(ds*length(degs))
  ks   <- downsample_perm(length(ys), ds)
  
  for(j in 1:ds){
    for(k in seq_along(degs)){
      xs <- ys[ks[[j]]]
      res <- inverseWT1(xs, degs[k])
      res <- inverseWT1(res, degs[k])
      errs[(j-1)*3 + k] <- check_error(ys, res, omit, j-1)
      cat(sprintf("max err is %f \n", max(errs)))
    }
  }
      expect_that(max(errs) < tol ,is_true())
})




test_that("inverseWT5 interpolation error is below tolerance",{
  ys   <- f()
  ds   <- 6
  degs  <- c(1,3, 5)
  tol  <- 0.7;   
  omit <- 20
  err  <- double(ds*length(degs))
  ks <-downsample_perm(length(ys), ds)
  for(j in 1:length(ks)){
    for(k in seq_along(degs)){
      xs <- ys[ks[[j]]]
      res <- inverseWT5(xs, degs[k])
      err[(j-1)*3 + k] <- check_error(ys, res, omit, j-1)
      cat(sprintf("cur err is %f \n", err[j]))
    }
   }
      expect_that(max(err) < tol ,is_true())
})


# downsample 5
test_that("inverseWT4 error is below tolerance",{
  deg <- 3 
  tol  <- 0.3 
  omit <- 15
  ds   <- 5
  ys   <- f()
  err <- double(ds)
  ks  <- downsample_perm(length(ys), ds)
  for(j in 1:(length(ks))){
    ind <- ks[[j]]
    xs <- ys[ind]
    res <- inverseWT4(xs, deg)
    err[j] <- check_error(ys, res, omit, j-1)
    # cat(is.na(err[j]), "\n")
    cat(sprintf("cur err is %f02 \n", (err[j])))
  }
  expect_that(max(err) < tol, is_true())
})


test_that("inverseWT2 interpolation error is below tolerance",{
  ds   <- 3
  ys   <- f()
  tol  <- .05  
  omit <- 8
  degs <- c(1,3,5); 
  err  <- double(2)
  ks <-downsample_perm(length(ys), 3)
  for(j in 1:(length(ks))){
    for(k in seq_along(degs)){
      xs <- ys[ks[[j]]]
      res <- inverseWT2(xs, degs[k])
      err[j] <- check_error(ys, res, omit, j-1)
      cat(sprintf("cur err is %f02 \n", err[j]))
    }
  }
  expect_that(max(err) < tol ,is_true())
})


test_that("inverseWT1 interpolation error is below tolerance",{
  ys   <- f()
  ds   <- 2
  deg  <- c(1,3, 5)
  tol  <- 0.03;   
  omit <- 5
  err  <- double(ds*length(deg))
  ks <-downsample_perm(length(ys), ds)
  for(j in 1:ds){
    for(k in seq_along(deg)){
      xs <- ys[ks[[j]]]
      res <- inverseWT1(xs, deg[k])
      err[(j-1)*4 + k] <- check_error(ys, res, omit, j-1)
      cat(sprintf("max err is %f \n", max(err)))
    }
  }
      expect_that(max(err) < tol ,is_true())
})

