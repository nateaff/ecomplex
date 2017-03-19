
#' Compute the epsilon-complexity of a time series.
#'
#' @param xx A sequence of points. 
#' @param ds Number of times to downsample the input sequence.
#' @param method The interpolation or approximation method. One of
#'                c("bspline", "cspline")
#' @param max_degree The maximum order spline used in the approximation
#'              step
#'
#' @return A list of features or the slope coefficient of the fit.
#'
#'@export
ecomplexity <- function(xx, method = c("bspline", "cspline"), 
                             ds = 5, max_degree = 5){
  xx <- normalize(xx)
  epsilons <- double(ds-1)
  ds <- 2:ds
  method <- match.arg(method)
  for(k in ds){
    switch(method, 
    cspline = { epsilons[k-1] <- cspline_err(xx, k) },
    bspline = { epsilons[k-1] <- bspline_err(xx, k, max_degree) }, 
    ) # end switch
  }

  fit <- NA  #
  try({      
     fit <- lm(log(epsilons) ~ log(1/ds))   
  }, silent = T )
  structure(list(fit = fit, epsilons = epsilons, 
                ds = ds, max_degree = max_degree,
                method = method), 
                class = "ecomplexity")
}


#' Function returns result for a single downsample level.
#' 
#' @param ys           A vector or time series. 
#' @param sample_num   The amount the series is downsampled.
#' @param max_degree   The maximum degree spline polynomial to fit.
#' @export
bspline_err <- function(ys, sample_num, max_degree){
  xx <- 1:length(ys)
  df <- data.frame(x = xx,y = ys); 
  indices  <- downsample_perm(length(ys), sample_num);
  # minimum error for each permutation
  epsilons <- double(length(indices))
  for (k in 1:sample_num) {
    cur_knots = indices[[k]]
    temp      = 1:length(xx)
    hold_out  = temp[-cur_knots];

    # Mean absolute error for each degree spline
    errs <- matrix(1000, nrow = max_degree, 
                         ncol = length(hold_out) )
    for (d in 1:max_degree){
        basis  <- splines::bs(xx, knots = cur_knots, degree = d)
        yhat <- NA  
        try({      
           fit <- lm(y ~ basis, data = df);
           yhat <- stats::predict(fit)[hold_out]
           errs[d,] <- abs(ys[hold_out] - yhat)
        }, silent = T )
    }
    if(is.na(errs[d,])) epsilons[k] <- NA 
    else epsilons[k]  <- min(apply(errs, 1, sum))/length(ys) 
  }
  return(mean(epsilons))
}



# update with spline()
cspline_err <- function(ys, sample_num, max_degree = NULL){
  xx <- 1:length(ys)
  indices  <- downsample_perm(length(ys), sample_num);
  # errors for each permutation
  epsilons <- double(length(indices))
  for (k in 1:sample_num) {
    xs = indices[[k]]
    xout  = xx[-xs];
    yout <- spline(xs, ys[xs], xout = xout)
    epsilons[k]  <- sum(abs(yout$y - ys[-xs])) 
  }
  return(mean(epsilons))
}


#' Normalize sequence values to (0,1)
#'
#' @param  xx The sequence to normalize.
#' @return  The normalized sequence.
#' @export
normalize <- function(xx){
    if(!max(xx) == min(xx)){
      (xx - min(xx))/(max(xx) - min(xx)) 
    } else {
      rep(0, length(xx))
    }
  }


#' Create list of all possible index patterns 
#'  for a given downsampling amount.
#' 
#' Creates a list of indices that correspond to 
#'  possible downsampling of a vector of length n
#'  at the given downsample rate.
#'
#' @param  n Length of list.
#' @param  ds Downsample rate.
#'
#' @return  List of indices 
#' @export
downsample_perm <- function(n,ds){
  xx   <- 1:n
  ind  <- vector("list", ds)    
  for (k in 1:ds){
    ind[[k]] <- xx[(xx - (k))%% ds == 0] 
  }
  ind
}