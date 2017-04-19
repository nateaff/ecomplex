#' Function returns result for a single downsample level.
#' 
#' @param y           A vector or time series. 
#' @param sample_num  The amount the series is downsampled.
#' @param max_degree  The maximum degree spline polynomial to fit.
#' @param err_norm    The norm used in computing the approximation error.
#' @param sample_type Downsampling pattern to use. 
#' @return The mean errors for given sample_num
#' @importFrom stats spline
cspline_err <- function(y, sample_num, max_degree, err_norm, sample_type) {
  x <- 1:length(y)
  
  switch(sample_type, 
         "step"   = { indices <- downsample_perm(length(y), sample_num) },
         "random"  = { indices <- random_sample(length(y), sample_num) })

  # indices  <- downsample_perm(length(y), sample_num);
  epsilons <- double(length(indices))
  switch(err_norm, 
                "mse" = { err_norm <- mse }, 
                "mae" = { err_norm <- mae }, 
                "max" = { err_norm <- maxerr})

  for (k in 1:sample_num) {
    ind  <- indices[[k]]
    xout <- x[-ind];
    try({
    yout <- spline(ind, y[ind], xout = xout)
    # Average assuming sample points fit is exact
    epsilons[k]  <- err_norm(yout$y, y[-ind]) / length(y) 
    }, silent = TRUE)
  }
  return(mean(epsilons))
}
