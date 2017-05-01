#' Return ecomplex coefficients
#'
#' Returns a dataframe with a column for each
#' coefficient from an ecomplex object.
#'
#' @param fit An object of type 'ecomplex'.
#'
#' @return  A dataframe with the two coefficients.
#' @method coef ecomplex
coef.ecomplex <- function(fit){
  data.frame(A = fit$A, B = fit$B) 
}

#' Plot an ecomplex object
#'
#' Plot the linear fit log(epsilons) ~ log(S) 
#'  where epsilons are the errors at each downsample level
#'  and S is the fraction of points retained at a given
#'  downsample level.
#'
#' @param x The list returned by callling ecomplex()
#' @param ... Additional arguments to base plot function.
#' @method plot ecomplex
plot.ecomplex <- function(x, ...){
  plot(log(x$epsilons) ~ log(x$S), 
       col = "chocolate2",
       pch = 16, 
       cex = 1,
       xlab = expression(log(epsilon)), 
       ylab = expression(log(S)), 
       main = paste0('Log-log fit of the errors ',
                      'epsilon',  
                      ' v. percent of values used S'), 
       ...)

  abline(x$fit, col = "gray10", lwd = '1.6')

}