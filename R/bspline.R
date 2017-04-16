#' Function returns result for a single downsample level.
#' 
#' @param y           A vector or time series. 
#' @param sample_num   The amount the series is downsampled.
#' @param max_degree   The maximum degree spline polynomial to fit.
#' @importFrom splines bs
bspline_err <- function(y, sample_num, max_degree) {
  x <- 1:length(y)
  df <- data.frame(x = x, y = y); 
  indices  <- downsample_perm(length(y), sample_num)
  # minimum error for each permutation
  epsilons <- double(length(indices))
  for (k in 1:sample_num) {
    cur_knots <- indices[[k]]
    ind       <- 1:length(x)
    hold_out  <- ind[-cur_knots];
    # errs holds the absolue errors for each index set
    errs <- matrix(0, nrow = max_degree, ncol = length(hold_out))
    for (deg in 1:max_degree) {
        basis  <- splines::bs(x, knots = cur_knots, degree = deg)
        yhat <- NA  
        try({      
           fit      <- lm(y ~ basis, data = df);
           # Average on full prediction 
           yhat     <- stats::predict(fit)[hold_out]
           errs[deg,] <- abs(y[hold_out] - yhat) / length(y)
        }, silent = TRUE )
    }
    if (any(is.na(errs[deg, ]))) { 
      epsilons[k] <- NA 
    } else { 
      epsilons[k] <- min(apply(errs, 1, sum)) 
    }
  }
  return(mean(epsilons))
}
