 #----------------------------------------------------------
# Three patterns are used for interpolation and all 
#  other patterns are combinations of these
#  1: o x o x o x o
#  2: o x x o x x o xx o
#  3: o x o o x o x o o x
#  
#  This number matches the numbering of the filters
#  and the merges. 
#  
#  The numbering on the inverse transforms (upsampling)
#  matches the initial downsampling patterns, which 
#  are formed by removing every other 1, 2, 3, 3, 5
#  points, respectively.
#----------------------------------------------------------

# Coefficient matrices for interpolating polynomials

pfilter <- function(degree, pos){
    N <- degree + 1
    eval <- seq(min(pos)-1, max(pos)+1, by = 2)
    cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
    ymat <- diag(N)
    for(j in 1:length(eval)){
        for(k in 1:N){
            ys <- ymat[k,]
            cmat[j,k] <- pracma::neville(pos, ys, eval[j])
        }
    }
    cmat
}

# could update to just evaluate the 0 position
pfilter2 <- function(degree, pos){
  xx <- (min(pos)-2):(max(pos)+2)
  eval <- xx[!(xx %in% pos)]
  N <- degree + 1
  cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
  ymat <- diag(N)
  for(j in 1:length(eval)){
      for(k in 1:N){
          ys <- ymat[k,]
          cmat[j,k] <- pracma::neville(pos, ys, eval[j])
      }
  }
  cmat
}

pfilter3 <- function(degree, pos){
  eval <- 0
  N <- degree + 1
  cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
  ymat <- diag(N)
  for(j in 1:length(eval)){
      for(k in 1:N){
          ys <- ymat[k,]
          cmat[j,k] <- pracma::neville(pos, ys, eval[j])
      }
  }
  cmat
}


# filter for single downsample: o x o 
filters1 <- function(){
     pts <- list( linear <- c(-1,1),
                  quad <- c(-1,1,3),
                  cubic <- c(-3,-1,1,3),
                  quart <- c(-3,-1,1,3,5), 
                  quintic <- c(-5,-3,-1,1,3,5), 
                  sextic <-c(1,1,1,1,1,1,1), # dummy
                  septic <- c(-7,-5,-3,-1,1,3,5,7)  )
  degrees <- 1:7
  ret <- lapply(degrees, function(x) pfilter(x, pts[[x]] ))
  ret
}


# filter for downsample = 2 : o xx o xx o 
filters2 <- function(){
  pts <- list( linear1 <- c(-1,2),
                linear1 <- c(-1,2),
                cubic1 <- c(-4,-1,2,5),
                cubic1 <- c(-4,-1,2,5),
                quintic1 <- c(-7, -4,-1,2,5, 8))
  degrees <- c(1,1,3,3, 5)
  ret <- lapply(1:length(degrees), function(k) pfilter2(degrees[k], pts[[k]] ))
  ret
}


# filter for downsample pattern: o x oo x o x  
filters3 <- function(){
  # Shifting the points gives the same filter
  pts <- list( linear1 <- c(-1,1),
                cubic1 <- c(-2, -1, 1, 3),
                cubic2 <- c(-3, -1, 1, 2),
                quintic1 <- c(-4,-2, -1, 1, 3, 4),
                quintic2 <- c(-4,-3, -1, 1, 2, 4))
  degrees <- c(1,3,3,5,5)
  ret <- lapply(1:length(degrees), function(k) pfilter3(degrees[k], pts[[k]] ))
  ret
}


# lshift <- function(x, t){
#     x[(1+t):length(x)]
# }



predict <- function(ys, filter){
    sum(ys*filter)
}


merge1 <- function(odd, even){
    len <- length(even) + length(odd)
    xx <- double(len)
    for(k in 1:length(even)){
        xx[2*k -1] <- odd[k]
        xx[2*k] <- even[k]
    }
    xx[1:len]
}


inverseWT1 <- function(xx, degree){

    filter <- filters1()[[degree]]
    offset <- ceiling(degree/2) 
    ps <- double(length(xx))
    #Note: handling of edge is not correct for differnt downsample patterns
    #TODO: Handle edges case by case or skip
    for(k in 1:offset){
      ps[k] <- predict(xx[1:(2*offset)], filter[k,])
      ps[length(ps) - offset + k] <- predict(xx[(length(ps)-(2*offset)+1):length(ps)], 
            filter[degree - k, ])
    }
    for(k in offset:(length(ps)- offset)){
        ps[k] <- predict(xx[(k - offset + 1):(k+offset)], filter[offset+1,])
    }
    res <- merge1(xx[1:length(xx)], ps[1:length(xx)])[1:(length(xx)*2)]
    res
}



# pattern:  0 x x 0 x x 0
inverseWT2 <- function(xx, degree){
    filter <- filters2()[[degree]]
    offset <- get_offset(degree)
    ps <- double(length(xx)*2)
    # xhead <- degree + 1
    xtail <- length(ps)- degree + 1 
    xseq <- seq(1, xtail, by = 2)

    for(j in (offset + 1):length(xseq)){
      k <- xseq[j]
      ps[k] <- predict(xx[(j - offset +1 ):(j + offset)], filter[degree + 2,])
      ps[k+1] <- predict(xx[(j - offset + 1):(j + offset)], filter[degree + 3,] )
    }
     ret <- merge2(xx, ps[1:length(ps)])
     # plot(ret)
     ret
}

inverseWT3 <- function(xx, degree){
  xx <- inverseWT1(xx, degree)
  yy <- inverseWT1(xx, degree)
  # ts.plot(yy); lines(xs, col = "blue")
  yy    
}

get_offset <- function(n){
  ceiling((n+1)/2)
}



inverseWT5 <- function(xx, degree){
  xx <- inverseWT1(xx, degree)
  yy <- inverseWT2(xx, degree)
  yy
}



# pattern :  0 x x x x 0 x x x x 0 
inverseWT4 <- function(xx, degree = 3){
  
  offset <- get_offset(degree)
  xx <- inverseWT2(xx, degree)
  
  filters <- filters3() 
  filter1 <- filters[[degree]]
  filter2 <- filters[[degree-1]]    
  ps <-  double(2*length(xx)/3)
 
  # last point used in prediction
  xtail <- length(xx)- degree + 1 
  xseq <- seq(offset, xtail, by = 3)
  # hack 
  # start <- ifelse(offset == 1, 1, 2)
  pseq <- seq(2, length(ps), by = 2)
  
  N <- min(length(xseq), length(pseq))
  for(i in offset:length(xseq)){
    j <- xseq[i]
    k <- pseq[i]
    ps[k] <- predict(xx[(j):(j + length(filter1) - 1)], filter1)
    ps[k+1] <- predict(xx[(j + 1):(j + length(filter1) )], filter2)
    # ps[k+1] <- predict(xx[(j - offset + 1):(j + offset)], filter[degree + 3,] )
  }
   ret <- merge3(xx, ps)
   # ts.plot(ret)
   ret
}

#
merge2 <- function(odd, even){
    xx <- double(length(even) + length(odd))
    xseq <- seq(1, length(xx), by = 3)
    len <- length(odd) + length(even)

    for(k in seq_along(odd)){
        # cat("insert", xseq[k], "\n")
        xx[xseq[k]] <- odd[k]
        xx[xseq[k]+1] <- even[2*k-1]
        xx[xseq[k]+2] <- even[2*k]
    }
    xx[1:len]
}

# 
merge3 <- function(odd, even, offset){

  n <- length(odd) + length(even)
  x <- double(n)
  oseq <- floor((5*(1:n)-1)/3)
  eseq <- (1:n)[-oseq]
  for(k in seq_along(odd)) x[oseq[k]] = odd[k]
  for(k in seq_along(even)) x[eseq[k]] = even[k] 
  
  x[1:n]
}


 
#----------------------------------------------------------
# error function for complexity measure
#
#----------------------------------------------------------


# generate <- function(mod) UseMethod("generate")

iwt_mod <- function(ds){
 # stopifnot(is.character(ds) || is.integer(ds))
 if(is.integer(ds)) type <- as.character(ds) 
 switch(type,
          "1" = { mod <- structure(list(ds = ds, degs = c(1,3,5), 
                                       omit = 0, inverse = inverseWT1), 
                                       class = "iwt0")},
          "2" = { mod <- structure(list(ds = ds, degs = c(1,3,5), 
                                       omit = 5, inverse = inverseWT1), 
                                       class = "iwt1")},
         "3" = { mod <- structure(list(ds = ds, degs = c(1,3,5), 
                                       omit = 8, inverse = inverseWT2), 
                                       class = "iwt2")},
         "4" = { mod <- structure(list(ds = ds, degs = c(1,3,5), 
                                       omit = 10, inverse = inverseWT3), 
                                       class = "iwt3")},
         "5" = { mod <- structure(list(ds = ds, degs = 3, 
                                       omit = 15, inverse = inverseWT4), 
                                       class = "iwt4")},
         "6" = { mod <- structure(list(ds = ds, degs = c(1,3,5),
                                       omit = 20, inverse = inverseWT5),
                                       class = "iwt5")}
         )
}


check_resid <- function(ys, ts, omit, j){
    yrange <- (omit+j):(length(ys) - omit)
    trange <- (omit):(length(ys)-omit -j)
    sum(abs(ts[trange] - ys[yrange]))/length(trange)
}



interp_err <- function(ys, iwt) {
  errs  <- matrix(0, nrow = iwt$ds, ncol = length(iwt$degs))
  ks <-downsample_perm(length(ys), iwt$ds)
  for(j in 1:iwt$ds){
    for(k in seq_along(iwt$degs)){
      xs <- ys[ks[[j]]]
      res <- iwt$inverse(xs, iwt$degs[k])      
      errs[j, k] <- check_resid(ys, res, iwt$omit, j-1)
      # cat(sprintf("cur err is %f \n", max(errs[j,k])))
    }
  }
  mean(unlist(apply(errs, 1, min)))
}

# 
interp_err.iwt3 <- function(iwt, ys) 
{
  errs  <- matrix(0, nrow = iwt$ds, ncol = length(iwt$degs))
  ks <-downsample_perm(length(ys), iwt$ds)
  for(j in 1:iwt$ds){
    for(k in seq_along(iwt$degs)){
      xs <- ys[ks[[j]]]
      res <- inverseWT1(xs, iwt$degs[k])
      res <- inverseWT1(res, iwt$degs[k])
      errs[j, k] <- check_resid(ys, res, iwt$omit, j-1)
      # cat(sprintf("cur err is %f \n", max(errs[j,k])))
    }
  }
  # average the best fit for each permutation
  mean(unlist(apply(errs, 1, min)))
}


#'  Compute epsilon-complexity using the lifting scheme
#'
#'  Computes the epsilon-complexity coefficient of a 
#'   time series using a wavelet lifting scheme.
#'
#' @param  xx A time series or vector 
#'
#' @return A \code{list} with :
#'
#' \tabular{ll}{
#' \code{coefficients} \tab The epsilon-complexity coefficients,\cr
#' \code{fit} \tab The complete linear model generated by lm(), \cr
#' \code{epsilons} \tab The mean sum of absolute errors at each level \cr
#' \code{S}     \tab The fraction of points removed at each level 
#'}
#' @export
lift_comp <- function(xx){

  S <- 1/(2:6)
  epsilons <- unlist(lapply((2:6), function(y) interp_err(xx, iwt_mod(y))))


  fit <- NA
   try({      
     fit <- lm(log(epsilons) ~ log(S))   
  }, silent=T)
  structure(list(fit = fit, epsilons = epsilons, 
                S = S), class = "ecomp_lift")
}

wavelet <- function(xx, degree){
  xx <- merge1(sin(seq(-1,0 , by = 0.1)), cos(seq(2, 3, by = 0.1)))
  ts.plot(xx)
  xx <- c(0,0,0,0,0,0,1,0,0,0,0,0,0)
  degree <- 1
  filter <- filters2()[[1]]
}

# cache filters  
pfilter <- memoise::memoise(pfilter)
pfilter2 <- memoise::memoise(pfilter2)
pfilter3 <- memoise::memoise(pfilter3)

filters1 <- memoise::memoise(filters1)
filters2 <- memoise::memoise(filters2)
filters3 <- memoise::memoise(filters3)

# profvis::profvis({
#   lift_comp(rnorm(5000))
#   })



