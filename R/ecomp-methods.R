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
