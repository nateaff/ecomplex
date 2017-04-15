#  ne v i l l e . R  Neville and Newton Interpolation

#' Evaluate the Lagrange polynomials at position.
#'  
#' Compute the value of the Lagrange interpolating polynomial at a 
#'  position xs. The function is borrowed from the \code{pracma} package.
#'
#' @param xs The point at which to evaluate the lagrange polynomial.
#' @param y  The funtion values f(x) of the function to be interpolated.
#' @param x  The x locations of the function values y = f(x).
#' @return Interpolated value.
#' @export
neville <- function(x, y, xs) {
    stopifnot(is.numeric(x), is.numeric(y))
    if (!is.numeric(xs))
        stop("Argument 'xs' must be empty or a numeric vector.")
    x <- c(x); y <- c(y)
    n <- length(x)
    if (length(y) != n)
        stop("Vectors 'x' and 'y' must be of the same length.")

    ys <- y
    for (k in 1:(n-1)) {
        y[1:(n-k)] <- ((xs - x[(k+1):n]) * y[1:(n-k)] +
                       (x[1:(n-k)] - xs) * y[2:(n-k+1)]) /
                       (x[1:(n-k)] - x[(k+1):n])
    }
    ys <- y[1]
    return(ys)
}
