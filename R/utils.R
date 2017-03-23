#' Normalize sequence values to (0,1).
#'
#' Computes (x - min x)/(max x - min x).
#'
#' @param  x The sequence to normalize.
#' @return  The normalized sequence.
#' @export
normalize <- function(x){
    if(!max(x) == min(x)){
      (x - min(x))/(max(x) - min(x)) 
    } else {
      rep(0, length(x))
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
  x   <- 1:n
  ind  <- vector("list", ds)    
  for (k in 1:ds){
    ind[[k]] <- x[(x - (k))%% ds == 0] 
  }
  ind
}