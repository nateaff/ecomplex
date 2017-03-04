#  Based on the function Palarm2 
#  Anatoly Zlotnik and Alexandra Piryatinska
#  See paper :Automated detection of neonate EEG sleep stages

#' Estimation of change points and mean on each
#' stable interval.
#'
#' @param x Input data.
#' @param delta1 A weighting parameter for phase I.
#' @param delta2 A weighting parmameter for phase II.
#' @param pf False detection probability.
#' @param m  Minimum interval size.
#' @param M  Number of terms in KS series.
#' @param epsilon Minimum relative distance of change point 
#'  from endpoints of local interval.
#'
#' @return 
#'
#' A \code{list} with :
#'
#' \tabular{ll}{
#' \code{kout} \tab Vector of indices of the detected change points,\cr
#' \code{means} \tab Vector of mean values on each stationary segment. \cr
#' }

palarm <- function(x, delta1 =1, 
                      delta2 = 0, 
                      pf = 0.1, 
                      m = 5, 
                      M = 1000,
                      epsilon = 0.02){
  xin <- x; x <- (x - mean(x))/sd(x) 
  N <- length(x)   

  thresh <- abs(ksinverse(pf, M))
  p0 <- 1; 
  kin <- double(0); 
  kone <- palarmf(x, kin, p0, delta1, thresh, m, epsilon)
  if(length(kone) == 0){
    kout <- kone;
    out <- mean(xin)*rep(N, 1)

  } else {
    kone <- sort(kone)
    thresh <- abs(ksinverse(pf/10, M))
    kout <- diagn(x, kone, delta2, m, thresh)
  }
  delta3=0.5;
  finalal(x, xin, kout, m, delta3, thresh) 
}


# Returns ksdistribution with parameters (pf,m)
ksdist <- function(pf, m){
  function(y){
    k <- 1:m
    s <- rep(-1, m)^(k+1)
    2*sum(s*exp(-(2*y^2)*(k^2))) -pf  
  }
}


ksinverse <- function(pf, M){
  ks <- ksdist(pf,M)    
  # Adding .005 to 1 to match behavior of 
  # Matlab fzero() function
  res <- pracma::fzero(ks, 1.005)
  res$x
}

#' Calculate statistics for change point and
#'  maximize to detect change point k
#'
#' @param X Input data
#' @param delta Detect parameter
#' @param m Minimum interval size
#'
#'  
#' @return 
#'
#' A \code{list} with :
#'
#' \tabular{ll}{
#' \code{ystat} \tab The test stastistic,\cr
#' \code{k} \tab The estimated change point \cr
#' }
ystat <- function(X, delta, m){
  N <- length(X)
  Y <- 1:(N-1)
  stopifnot(N > 1)
  stat <- rep(1, N)
  M1 <- cumsum(X[1:(N-1)])/Y
  M2 <- rev(cumsum(rev(X[2:N]))/Y)
  stat <- abs((((1- (Y/N))*(Y/N))^delta)*(M1-M2))
  
  k <- which.max(round(stat, 3))
  list(stat = stat, k = k)
}

# 2004 by Anatoly Zlotnik and Alexandra Piryatinska

#' Recursive algorithm for detection of change points by successive interval
#' bisection
#'
#' 
#' @param  X       local data
#' @param kin     vector containing estimated global change points
#' @param P       global index of first point of local data
#' @param delta   detect parameter
#' @param thresh  minimum statistic for change point characterization
#' @param m       minimum interval size
#' @param epsilon minimum relative distance of change point from endpoints
#'
#'  @return kout    updated vector of estimated global change points
palarmf <- function( X, 
                    kin = double(0), 
                    P = 1, 
                    delta = 1, 
                    thresh = 1.224, 
                    m = 5, 
                    epsilon = 0.02 ){
  # cat("palarmf called", "\n")
  res <- ystat(X, delta, m)
  k <- res$k 
  stat <- res$stat
  N <- length(X) 
  sigma <- sd(c(X[1:k] - mean(X[1:k]) , X[(k+1):N]- mean(X[(k+1):N])))/sqrt(N)
  
  level <- thresh*sigma
  d <- max(round(epsilon*N), 1)
  if(stat[k] < level){
    kout <- kin
  } else {
    X1 <- X[1:( k-d )] 
    P1 <- P
    X2 <- X[(k + d ):N]
    P2 <- P + k + d - 1;
    
    if( length(X1) <= m ) Ktemp1 <- kin
    else if (length(X1) > m) {
     Ktemp1 <- palarmf( X1, kin, P1, delta, thresh, m, epsilon)
    }
    if( length(X2) <= m )  Ktemp2 <- Ktemp1 
    else if( length(X2) > m ) {
       Ktemp2 <- palarmf( X2, Ktemp1, P2, delta, thresh, m, epsilon )
    }
    if( length(X1) > m && length(X2) > m ) {
      kout <- c( Ktemp2, P+k-1 )
    }
    else kout <- Ktemp2
  } 
  return(kout)
}

#' Checks estimated global change points for errors 
#'  and performs update
#'
#' @param x Input data
#' @param k A vector containing estimated global change points
#' @param delta Detect parameter
#' @param m Minimum interval size
#' @param thresh Minimum statistic for chnage point characterization
#' 
#' @return kout Estimated global change points
diagn <- function(x, kin, delta, m, thresh){
  N <- length(x); 
  Z <- length(kin); 
  kout <- double(0)
  # Some change points found
  if(length(kin) > 1) {
    b <- floor((kin[2] + kin[1])/2); a <- 0;
    kout <- check_pt(x[1:b], a, delta, m, thresh, kout)
    #cat(sprintf('a= %0.3f, b=%0.3f \n', a,b))
    for(i in 2:(length(kin) -1)){
      a <- floor((kin[i]+ kin[i-1])/2)+1
      b <- floor((kin[i+1] + kin[i])/2)
      kout <- check_pt(x[a:b], a, delta, m, thresh, kout)
    } # end for
    # final change point 
    a <- round((kin[Z]+ kin[Z-1])/2)
    b = N;
    kout <- check_pt(x[a:N], a, delta, m, thresh, kout)
    #cat(sprintf('a= %0.3f, b=%0.3f \n', a,b))
  } else {
    stats <- ystat(x, delta, m)
    if(stats$stat[stats$k] > (thresh*(sd(x)/sqrt(N))) ) {
      kout <- c(kout, stats$k)
    }
  }
  kout 
}


check_pt <- function(x, start, delta, m, thresh, kout){
  stats <- ystat(x, delta, m) 
  sigma <- sd(x)/sqrt(length(x)-1) 
  level <- thresh*sigma

  adjust <- 0
  if(start > 0) adjust <- 1
  pos <- start + stats$k - adjust

  if(stats$stat[stats$k] > level){
    #cat(sprintf('add k= %d \n', pos))
      kout <- c(kout, pos)
  }
  kout
}



#' Calculate final change point estimates and mean of each interval
#'  Based on Matlab function of Anatoly Zlotnik and Alexandra Piryatinska
#'    
#' @param X       Input data.
#' @param xin     The original time-series.
#' @param kin     A vector containing estimated global change points
#' @param m       minimum interval size.
#' @param delta   detect parameter
#' 
#' A \code{list} with :
#'
#' \tabular{ll}{
#' \code{meanK} \tab the mean values of the intervals,\cr
#' \code{kout} \tab the final estimated change points.\cr
#' }
finalal <- function(X, xin, kin, m, delta, thresh){

N <- length(X); 
Z <- length(kin);

kout <- double(0) 

if(length(kin)>1) {
    b <- floor((kin[2] + kin[1])/2);
    X1 <- X[1:b];
    stats <- ystat(X1,delta,m);
    
    if(stats$k > 1) {
        me <- mean(xin[1:(stats$k-1)])
        meanK <- me*rep(1, (stats$k - 1))
        } else {
        meanK <- double(0)
      }
      kout <- c(kout, stats$k)
    for(i in 2:(length(kin)-1)){
    # for i=2:length(kin)-1  
        a <- floor((kin[i] + kin[i-1] )/2) + 1
        b <- floor((kin[i+1] + kin[i])/2)
        # Xtemp <- X[a:b];
        stats <- ystat(X[a:b], delta,m)
        kout <- c(kout, a + stats$k - 1 )
        M1 <- mean(xin[(kout[i-1]+1):(kout[i]-1)])
        M2 <- M1*(rep(1, kout[i]-kout[i-1] ))
        meanK <- c(meanK, M2)
    }
    a <- round((kin[Z] + kin[Z-1])/2)
    X1 <- X[a:N]
    stats <- ystat(X1, delta, m);
    kout <- c(kout, a + stats$k - 1);
    M1 <- mean(xin[ (1 + kout[Z-1]): (kout[Z]-1)] )
    M2 <- M1*(rep(1 ,kout[Z]-kout[Z-1]))
    meanK <- c(meanK, M2)
    M3 <- mean(xin[(kout[Z] + 1):N])
    meanK <- c(meanK, M3*rep(1 , N - kout[Z]+1))
  } else {
     stats  <- ystat(X, delta, Fm); 
     kout <-  c(kout, stats$k)
     meanK <- c(rep(1, stats$k) * mean(xin[1:stats$k]), 
                rep(1, N - stats$k) * mean(xin[ stats$k:N] )) 
  }
  list(kout = kout, means = meanK)
}

