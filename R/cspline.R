#' Function returns result for a single downsample level.
#' 
#' @param y           A vector or time series. 
#' @param sample_num   The amount the series is downsampled.
#' @param max_degree   The maximum degree spline polynomial to fit.
#' @return The mean errors for given sample_num
#' @importFrom stats spline
cspline_err <- function(y, sample_num, max_degree = NULL) {
  x <- 1:length(y)
  indices  <- downsample_perm(length(y), sample_num);

  epsilons <- double(length(indices))
  for (k in 1:sample_num) {
    ind  <- indices[[k]]
    xout <- x[-ind];
    try({
    yout <- spline(ind, y[ind], xout = xout)
    # Average assuming sample points fit is exact
    epsilons[k]  <- sum(abs(yout$y - y[-ind])) / length(y) 
    }, silent = TRUE)
  }
  return(mean(epsilons))
}