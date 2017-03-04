#' Check the depth of a list.
#'
#' Recursively checks the maximum depth of a list.
#'
#' @param this  A list.
#' @return The depth of the list.
#' @export
depth <- function(x) {
  ifelse(is.list(x), 1L + max(sapply(x, depth)), 0L)
}

#' Return 2-D Size of Input.
#'
#' Matlab style size() function for vectors, lists, 
#' matrices or data.frames.
#'
#' @param x  Any one or two dimensional R data structure.
#' @return A numeric vector of length two: c(#rows, #cols).
#' @export
size <- function(x){
  if(is.null(dim(x))){
    m <- 1
    n <- length(x)
  } else {
    m <- dim(x)[1]
    n <- dim(n)[2]
  }
  c(m,n)
}
