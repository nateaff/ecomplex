#' Extract set of features from a time series.
#'
#' @param data A time series, matrix of dataframe. 
#' @param features A list or vector of ts_features.
#' @param id Optional prefix for column headers.
#' @param multicore Parallelize feature extraction (Linux only).
#' @param verbose Print status to console.
#' 
#' @return A data frame with row of features for each time.
#'  series in data. (TODO: single feature, system 
#'  check for parallelization)
#'  
#' @export 
extract_features <- function(data, features, id = "1", multicore = TRUE, verbose = FALSE){

  # TODO: switch for system type
  if(multicore){
   cores <- parallel::detectCores()
   if(verbose) cat("Using", cores, " cores ... \n")
  }
  # TODO : single feature 
  if(is.vector(data) || class(data) == "ts"){
    cat("ts")
    ret <- extract_one_feature(as.vector(data), features)
    ret$id <- id
  } else {
    if(is.matrix(data)) {
      data <- as.data.frame(data) 
    }
    stopifnot(is.data.frame(data))

    f <- function(feature){
      if(multicore){
        ret <- parallel::mclapply(data, function(x) extract_one_feature(x, feature), 
                        mc.cores = cores)
      } else {
        ret <- lapply(data, function(x) extract_one_feature(x, feature))
      }
    do.call(rbind, ret)
    }

    ret <- do.call(cbind,(lapply(features, f))) 
    ret$id <- id
  } # end else
  ret 
} 


extract_one_feature <- function(tseries, feature, ...){
  clean_feature(feature(tseries, ...))
}

#'export
clean_feature <- function(x, ...) UseMethod("clean_feature")


clean_feature.FractalDim <- function(feature){
  ret <- plyr::unrowname(data.frame(fd = feature$fd)) 
  names(ret) <- paste0("fd_", feature$methods)
  ret
}


clean_feature.fd_variogram <- function(feature){
  ret <- plyr::unrowname(data.frame(fd = feature$fd)) 
  names(ret) <- paste0("fd_", feature$methods)
  class(ret) <- "fd_variogram"
  ret
}

clean_feature.bandpower <- function(feature){
  plyr::unrowname(data.frame((t(unlist(feature)))))
}


clean_feature.ecomp_bspline <- function(feature){
  ret <- plyr::unrowname(data.frame(t(feature$fit$coefficients))) 
  names(ret) <- paste0("bspline_", c("A", "B"))  
  ret
}

clean_feature.ecomp_cspline <- function(feature){
 ret <- plyr::unrowname(data.frame(t(feature$fit$coefficients)))
 names(ret) <- paste0("cspline_", c("A", "B"))  
 ret
}

clean_feature.ecomp_lift <- function(feature){
  ret <- plyr::unrowname(data.frame(t(feature$fit$coefficients))) 
  names(ret) <- paste0("lift_", c("A", "B"))  
  ret
}

clean_feature.sample_entropy <- function(feature){
  plyr::unrowname(data.frame(sample_entropy = feature[1])) 
}

clean_feature.hurst <- function(feature){
  plyr::unrowname(data.frame(hurst = feature$Hs))
}

clean_feature.variance <- function(feature){
   plyr::unrowname(data.frame(var = feature[1]))
}

clean_feature.permutation_entropy <- function(feature){
  plyr::unrowname(data.frame(p_entropy = feature[1]))
}

clean_feature.wvar <- function(feature){
  wav_var <- feature$variance[1:4]
  fnames <- paste0("wvar_", feature$scales[1:4])
  ret <- plyr::unrowname(data.frame(t(wav_var)))
  names(ret) <- fnames
  ret
}
 

clean_feature.default <- function(feature){
  plyr::unrowname(data.frame(feature = feature[1]))
}

#' Compute the permutation entropy of a time series.
#'
#' Calculates the permutation entropy of a time 
#'  series using the pdc package. The miminimum
#'  entropy over the embedding dimensions 3:7 is
#'  returned.
#'
#' @param xx A time series or vector.
#'
#' @return The permutation entropy of the time series.
#' @export
permutation_entropy <- function(xx){
  ret <- pdc::entropyHeuristic( xx )
  row <- which(ret$entropy.values[, 2] == ret$m)
  ret <- ret$entropy.values[row, 3]
  class(ret) <- "permutation_entropy"
  ret
}

#' Compute the sample entropy of the data
#'
#' @param xx The data 
#'
#' @return Sample entropy
#' @export
sample_entropy <- function(xx){
  cat("sample entropy \n" )
  ret <- pracma::sample_entropy(xx)
  # structure(list(sample_entropy = res), class = "sample_entropy")
  class(ret) <- "sample_entropy"
  ret
}

#' Corrected Hurst exponent
#'
#' The pracma implementation of the corrected Hurst
#'  exponent
#' 
#' @param xx The data 
#'
#' @return The features
#' @export
#
hurst <- function(xx){
  ret <- pracma::hurstexp(xx, d = 50, display = FALSE)
  class(ret) <- "hurst"
  ret
}

#' Variance wrapper
#'
#'
#' @param xx The data 
#'
#' @return The variance
#' @export
variance <- function(xx){
  cat("var \n")
  ret <- var(xx)
  class(ret) <- "variance"
  ret
}

#' Compute epsilon complexity using lifting 
#'
#'
#' @param xx The data 
#'
#' @return The features
#' @export
#' 
ecomp_lift <- function(xx){
  cat("ecomp lift \n")
  res <- lift_comp(xx)
  class(res) <-"ecomp_lift"
  res
}

#' Compute epsilon complexity using bsplines
#'
#'
#' @param xx The data 
#'
#' @return The features
#' @export
#' 
ecomp_bspline <- function(xx){
  cat("ecomp \n")
  res <- ecomplexity(xx, ds = 5, method = "bspline", max_degree = 4)
  class(res) <- "ecomp_bspline"
  res
}

ecomp_cspline <- function(xx){
  res <- ecomplexity(xx, ds = 5, method = "cspline")
  class(res) <- "ecomp_cspline"
  res
}


#' Compute bandpower
#'
#' Computes the bandpower on a default set 
#'   
#' @param xx The data
#' 
#' @return A data frame of power in each band
#' @export 
bandpower <- function(xx){
  cat("Bandpower \n")

  freqs <- list(
  delta = c(0.5,4),
  theta = c(4,8),
  alpha = c(8, 12),
  beta = c(12, 30),
  gamma = c(30, 100))

  if(!is.ts(xx)){
    fs <- 1220
  } else {
    fs <- frequency(xx)
  }
  res <- bp_pgram(xx, fs=fs, freqs=freqs)
  class(res) <- "bandpower"
  res
}


#' Variogram-based estimate of fractal dimension. 
#' 
#' This function is a wrapper for the fd.estimate 
#' function of the fractaldim package.
#'
#' @param xx The data
#'
#' @return The results from fd.estimate
#' @export
fd_variogram <- function(xx){
  cat("fd_variogram \n")
  fractaldim::fd.estimate(xx, methods = "variogram") 
}

#Temp 

#' Extract eeg features 
#'
#' @param x The number of time series generated for each group
#' @param features The features set to run 
#' @param multicore Parallelize computations (Linux only)
#'
#' @return A data frame with features as columns and an id
#'  column identifying group membership
#' @export
eeg_features <- function(x, features = NULL, multicore = TRUE){
  if (is.null(features)){
  features <- c(ecomp_bspline, bandpower, fd_variogram)
  }

  df <- extract_features(x, features, id = "1", multicore = multicore)
  df
} 

