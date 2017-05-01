#' Return ecomplex coefficients
#'
#' Returns a dataframe with a column for each
#' coefficient from an ecomplex object.
#'
#' @param fit An object of type 'ecomplex'.
#'
#' @return  A dataframe with the two coefficients.
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
#' @param fit The list returned by callling ecomplex()
#'
#' @export
plot.ecomplex <- function(fit){
  plot(log(fit$epsilons) ~ log(fit$S), 
       col = "chocolate2",
       pch = 16, 
       cex = 1,
       xlab = expression(log(epsilon)), 
       ylab = expression(log(S)), 
       main = paste0('Log-log fit of the errors ',
                      'epsilon',  
                      ' v. percent of values used S'))

  abline(fit$fit, col = "gray10", lwd = '1.6')

}