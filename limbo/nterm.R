#' Find Nterm Complexity of a Time Series.
#'
#' Estimate of the complexity of a function based 
#' on the rate of error decay as fewer terms are used
#' to approximate the function. 
#'
#' @param ys A vector, time series. 
#' @param verbose Print debugging information.
#' @param slope_only Return only the slope.
#'
#' @return A list or the slope coefficient of the nterm approximation
nterm <- function(ys, verbose = FALSE, slope_only = FALSE, ...){
  trim = FALSE
  if(verbose) cat("verbose ...")
  ys <- normalize(ys)
  # TODO parameter bin
  start <- 30
  nbins <- 100
  ylen <- length(ys)
  # TODO: return info about modwt parameters
  wc <- waveslim::modwt(ys, ...)
  # remove smooth coefficients   
  # xsort <- sort(abs(x[1:(length(x)- ylen)])) 
  wc_sort <- sort(abs(unlist(wc))) 
  bin_len <- as.integer(floor(length(wc_sort)/nbins))
  thresholds <- errors <- fits <- vector("list", nbins-start+1)
  for(k in start:nbins){
    
    cutoff <- wc_sort[k*bin_len]
    j <- k- start + 1
    thresholds[[j]] <- cutoff 
    # if(verbose) sprintf("theta = %2f", theta)
    wc_thresh <- manual_thresh2(wc, value = cutoff, 
                                max.level = length(wc))
    fits[[j]] <- waveslim::imodwt(wc_thresh)
    errors[[j]] <- sum((ys - fits[[j]])^2)/length(ys)
  }
  #Note : magic number 3
  
  if(trim){
    trimmed_errors <- errors[3:18]
    coefficients <- fit_terms(trimmed_errors)
  } else {
    coefficients <- fit_terms(errors)
  }
  # coefficients <- fit_terms(errors)
  res <- structure(
              list(errors = errors, 
              thresholds =  thresholds, 
              fits = fits, 
              coefficients = coefficients), class = "nterm")
  res 
}

# 
fit_terms <- function(y, x = 1:length(y)){
  fit <- lm(unlist(y) ~ x)
  fit$coefficients
} 

# Adapted from waveslim::manual_thresh
manual_thresh2 <- function (wc, max.level = 4, value, hard = TRUE) {
    wc.shrink <- wc
    if (hard) {
        for (i in names(wc)[1:max.level]) {
            wci <- wc[[i]]
            # unithresh <- factor * value
            unithresh <- value
            wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
        }
    } else {
        for (i in names(wc)[1:max.level]) {
            wci <- wc[[i]]
            unithresh <- factor * value
            wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) * 
                (abs(wci) > unithresh)
        }
    }
    wc.shrink
}



#' Find the nterm complexity of a function.
#'
#' Estimate of the complexity of a function based 
#' on the rate of error decay as fewer terms are used
#' to approximate the function. 
#'
#' @param ys A vector, time series.
#' @param verbose Print debugging information.
#' @param slope_only Return only the slope.
#'
#' @return A list or the slope coefficient of the nterm approximation
nterm2 <- function(ys, verbose = FALSE, slope_only = FALSE, ...){
  trim = FALSE
  if(verbose) cat("verbose ...")
  ys <- normalize(ys)
  # param :
  # number of bins
  start <- 30
  nbins <- 100
  ylen <- length(ys)
  wc <- waveslim::modwt(ys, ...)
  wc_sort <- sort(abs(unlist(wc))) 
  bin_len <- as.integer(floor(length(wc_sort)/nbins))
  Ns <- thresholds <- errors <- fits <- vector("list", nbins-start+1)
  for(k in start:nbins){
    cutoff <- wc_sort[k*bin_len]
    j <- k- start + 1
    Ns[[j]] <- (length(wc_sort) - k*bin_len)
    thresholds[[j]] <- cutoff 
    # if(verbose) sprintf("theta = %2f", theta)
    wc_thresh <- manual_thresh2(wc, value = cutoff, 
                                max.level = length(wc))
    fits[[j]] <- waveslim::imodwt(wc_thresh)
    errors[[j]] <- sum((ys - fits[[j]])^2)/length(ys)   
  }
  index <- 1:length(errors)
  errors <- lapply(index, function(x) log(errors[[x]])/log(Ns[[x]]))
  coefficients <- fit_terms(errors[3:60], unlist(Ns[3:60]))
  res <- structure(
              list(errors = errors, 
              thresholds =  thresholds, 
              Ns = Ns,
              fits = fits, 
              coefficients = coefficients), class = "nterm2")
  res 
}

#' Find Nterm Complexity of a Time Series.
#'
#' Find the nterm complexity of a function
#'
#' @param ys A vector, time series (or matrix ... ?) 
#' @param verbose Print debugging information
#'
#' @return A list or the slope coefficient of the nterm approximation
nterm3 <- function(ys, verbose = FALSE, ...){
  trim = FALSE
  if(verbose) cat("verbose ...")
  ys <- normalize(ys)

  wc <- waveslim::modwt(ys, n.levels = 4 ,...)  
  wc <- wc[-length(wc)] 
  
  stat <- unlist(lapply(wc, function(x) sum(abs(x)^2)))
  coeffs <- fit_terms(log(stat))  
  
  res <- structure(list(coefficients = coeffs,  stat = stat,  
                        js = 1:length(stat)), 
                        class = "nterm_test")
  res 
}








