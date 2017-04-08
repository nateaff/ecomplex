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
#
# The general algorithm is : 
# 1. Filters-k hold the sampling pattern and are passed to
# 2. pfilters-k which create the matrix of time-domain 
#    filters used for interpolation.
# 3. The inverse functions interpolate points which are 
# 4. merged for depending on the downsampling pattern
#----------------------------------------------------------


#----------------------------------------------------------
# Pfilter functions build matrices of 
# time domain filters for various each downsampling pattern 
# and each position. 
#
# Note: Border filters aren't being used and these 
# functions can be simplified since only the centered 
# position (not the borders) are being interpolated.
#----------------------------------------------------------
pfilter1 <- function(degree, pos){
    N <- degree + 1
    eval <- seq(min(pos) - 1, max(pos) + 1, by = 2)
    cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
    ymat <- diag(N)
    for(j in 1:length(eval)){
      for(k in 1:N){
          y <- ymat[k, ]
          cmat[j, k] <- pracma::neville(pos, y, eval[j])
      }
    }
    cmat
}

# could update to just evaluate the 0 position
pfilter2 <- function(degree, pos){
  x <- (min(pos) - 2):(max(pos) + 2)
  eval <- x[!(x %in% pos)]
  N <- degree + 1
  cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
  ymat <- diag(N)
  for(j in 1:length(eval)){
      for(k in 1:N){
          y <- ymat[k, ]
          cmat[j, k] <- pracma::neville(pos, y, eval[j])
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
          y <- ymat[k,]
          cmat[j, k] <- pracma::neville(pos, y, eval[j])
      }
  }
  cmat
}

 
#----------------------------------------------------------
# Filters for each degree and downsampling pattern
#----------------------------------------------------------

# filter for single downsample: o x o 
filters1 <- function(){
     pts <- list( linear  = c(-1, 1),
                  quad    = c(-1, 1, 3),
                  cubic   = c(-3, -1, 1, 3),
                  quart   = c(-3, -1, 1, 3, 5), 
                  quintic = c(-5, -3, -1, 1, 3, 5), 
                  sextic  = c(1, 1, 1, 1, 1, 1, 1), # dummy
                  septic  = c(-7, -5, -3, -1, 1, 3, 5, 7))
  degrees <- 1:7
  ret <- lapply(degrees, function(x) pfilter1(x, pts[[x]] ))
  ret
}


# filter for downsample = 2 : o xx o xx o 
filters2 <- function(){
  pts <- list(  linear1  = c(-1, 2),
                linear2  = c(-1, 2),
                cubic1   = c(-4, -1, 2, 5),
                cubic2   = c(-4, -1, 2, 5),
                quintic1 = c(-7, -4, -1, 2, 5, 8))
  degrees <- c(1,1,3,3, 5)
  ret <- lapply(1:length(degrees), function(k) pfilter2(degrees[k], pts[[k]] ))
  ret
}


# filter for downsample pattern: o x oo x o x  
filters3 <- function(){
  # Shifting the points gives the same filter
  pts <- list( linear1   = c(-1,1),
                cubic1   = c(-2, -1, 1, 3),
                cubic2   = c(-3, -1, 1, 2),
                quintic1 = c(-4,-2, -1, 1, 3, 4),
                quintic2 = c(-4,-3, -1, 1, 2, 4))
  degrees <- c(1,3,3,5,5)
  ret <- lapply(1:length(degrees), function(k) pfilter3(degrees[k], pts[[k]] ))
  ret
}


predict <- function(y, filter){
    sum(y*filter)
}


get_offset <- function(n){
  ceiling((n + 1) / 2)
}


 
#----------------------------------------------------------
# The 'inverse' are the interpoloating functions for 
# each downsample pattern. 
#----------------------------------------------------------
inverseWT1 <- function(x, degree){

    filter <- filters1()[[degree]]
    offset <- ceiling(degree/2) 
    ps <- double(length(x))

    for(k in 1:offset){
      ps[k] <- predict(x[1:(2 * offset)], filter[k, ])
      ps[length(ps) - offset + k] <- predict(x[(length(ps) - (2 * offset) + 1):length(ps)], 
            filter[degree - k, ])
    }
    for(k in offset:(length(ps)- offset)){
        ps[k] <- predict(x[(k - offset + 1):(k+offset)], filter[offset+1, ])
    }
    res <- merge1(x[1:length(x)], ps[1:length(x)])[1:(length(x)*2)]
    res
}



# pattern:  0 x x 0 x x 0
inverseWT2 <- function(x, degree){
    filter <- filters2()[[degree]]
    offset <- get_offset(degree)
    ps <- double(length(x)*2)

    xtail <- length(ps)- degree + 1 
    xseq <- seq(1, xtail, by = 2)

    for(j in (offset + 1):length(xseq)){
      k <- xseq[j]
      ps[k] <- predict(x[(j - offset +1 ):(j + offset)], filter[degree + 2,])
      ps[k+1] <- predict(x[(j - offset + 1):(j + offset)], filter[degree + 3, ])
    }
     ret <- merge2(x, ps[1:length(ps)])
     ret
}

inverseWT3 <- function(x, degree){
  y <- inverseWT1(x, degree)
  inverseWT1(y, degree)
}


inverseWT5 <- function(x, degree){
  x <- inverseWT1(x, degree)
  inverseWT2(x, degree)
}



# pattern :  0 x x x x 0 x x x x 0 
inverseWT4 <- function(x, degree = 3){
  
  offset <- get_offset(degree)
  x <- inverseWT2(x, degree)
  
  filters <- filters3() 
  filter1 <- filters[[degree]]
  filter2 <- filters[[degree-1]]    
  ps <-  double(2*length(x)/3)
 
  xtail <- length(x)- degree + 1 
  xseq <- seq(offset, xtail, by = 3)
  pseq <- seq(2, length(ps), by = 2)
  
  N <- min(length(xseq), length(pseq))
  for(i in offset:length(xseq)){
    j <- xseq[i]
    k <- pseq[i]
    ps[k] <- predict(x[(j):(j + length(filter1) - 1)], filter1)
    ps[k+1] <- predict(x[(j + 1):(j + length(filter1) )], filter2)
  }
   ret <- merge3(x, ps)
   ret
}
 
#----------------------------------------------------------
# The merge functions merge interpolated points with 
# the previous points.
#----------------------------------------------------------

merge1 <- function(odd, even){
    len <- length(even) + length(odd)
    x <- double(len)
    for(k in 1:length(even)){
        x[2*k -1] <- odd[k]
        x[2*k] <- even[k]
    }
    x[1:len]
}

merge2 <- function(odd, even){
    x <- double(length(even) + length(odd))
    xseq <- seq(1, length(x), by = 3)
    len <- length(odd) + length(even)

    for(k in seq_along(odd)){
        # cat("insert", xseq[k], "\n")
        x[xseq[k]] <- odd[k]
        x[xseq[k]+1] <- even[2*k-1]
        x[xseq[k]+2] <- even[2*k]
    }
    x[1:len]
}

merge3 <- function(odd, even, offset){

  n <- length(odd) + length(even)
  x <- double(n)
  oseq <- floor((5*(1:n)-1)/3)
  eseq <- (1:n)[-oseq]
  for(k in seq_along(odd)) x[oseq[k]] = odd[k]
  for(k in seq_along(even)) x[eseq[k]] = even[k] 
  
  x[1:n]
}




# Returns a list holding parameters for the appropriate
# interpolation procedure based on the downsaple level "ds"
iwt_mod <- function(ds){
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


check_resid <- function(y, est, omit, j){
    yrange <- (omit + j):(length(y) - omit)
    erange <- (omit):(length(y) - omit - j)
    sum(abs(est[erange] - y[yrange])/length(erange))
}



interp_err <- function(y, iwt) {
  mem()
  big_num <- 1e6
  errs  <- matrix(big_num, nrow = iwt$ds, ncol = length(iwt$degs))
  ks <- downsample_perm(length(y), iwt$ds)
  for(j in 1:iwt$ds){
    for(k in seq_along(iwt$degs)){
      xs <- y[ks[[j]]]
      res <- iwt$inverse(xs, iwt$degs[k])      
      errs[j, k] <- check_resid(y, res, iwt$omit, j - 1)
    }
  }
  mean(unlist(apply(errs, 1, min)))
}


# cache filters  
# @importFrom memoise memoise
mem <- function(){
  pfilter1 <- memoise::memoise(pfilter1)
  pfilter2 <- memoise::memoise(pfilter2)
  pfilter3 <- memoise::memoise(pfilter3)

  filters1 <- memoise::memoise(filters1)
  filters2 <- memoise::memoise(filters2)
  filters3 <- memoise::memoise(filters3)
}