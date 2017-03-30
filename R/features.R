#' Extract set of features from a time series.
#'
#' @param data A time series, matrix of dataframe. 
#' @param features A list or vector of ts_features.
#' @param id Optional prefix for column headers.
#' @param ncores Parallelize feature extraction (Linux only).
#' @param verbose Print status to console.
#' 
#' @return A data frame with row of features for each time.
#'  series in data. 
#'  
#' @export 
extract_features <- function(data, features, id = "1", 
                                             ncores = NULL, 
                                             verbose = FALSE){

  # TODO: switch for system type
  data <- as.data.frame(data)
  if(anyNA(data)) stop("Data contains NA values")
  
  if(!is.null(ncores)){
   cores <- parallel::detectCores()
   if(verbose) cat("Using", cores, " cores ... \n")
  }
  if(length(features) == 1) features <- list(features)
  
  f <- function(feature){
    if (!is.null(ncores)) {
      ret <- parallel::mclapply(data, 
                                function(x) extract_one_feature(x, feature), 
                                mc.cores = ncores)
    } else {
      ret <- lapply(data, function(x) extract_one_feature(x, feature))
    }
  do.call(rbind, ret)
  }

  ret <- do.call(cbind, (lapply(features, f))) 
  ret$id <- id
  ret 
} 


extract_one_feature <- function(tseries, feature, ...){
  clean_feature(feature(tseries, ...))
}

#' Returns dataframe with feature(s)
#'
#'@param feature The feature to clean
#'export
clean_feature <- function(feature) UseMethod("clean_feature")


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
  # change to Hal = (R/S -AL)
  # plyr::unrowname(data.frame(hurst = feature$Hs))
  plyr::unrowname(data.frame(hurst = feature$Hal))

}

clean_feature.variance <- function(feature){
   plyr::unrowname(data.frame(var = feature[1]))
}

clean_feature.permutation_entropy <- function(feature){
  plyr::unrowname(data.frame(p_entropy = feature[1]))
}

clean_feature.spectral_entropy <- function(feature){
    plyr::unrowname(data.frame(spec_entropy = feature[1]))
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
#' @param x A time series or vector.
#'
#' @return The permutation entropy of the time series.
#' @export
permutation_entropy <- function(x){
  cat("permutation entropy \n")
  ret <- pdc::entropyHeuristic( x )
  row <- which(ret$entropy.values[, 2] == ret$m)
  ret <- ret$entropy.values[row, 3]
  class(ret) <- "permutation_entropy"
  ret
}

#' Compute the sample entropy of the data
#'
#' @param x The data 
#'
#' @return Sample entropy
#' @export
sample_entropy <- function(x){
  cat("sample entropy \n" )
  ret <- pracma::sample_entropy(x)
  # structure(list(sample_entropy = res), class = "sample_entropy")
  class(ret) <- "sample_entropy"
  ret
}

#' Corrected Hurst exponent
#'
#' The pracma implementation of the corrected Hurst
#'  exponent
#' 
#' @param x The data 
#'
#' @return The features
#' @export
#
hurst <- function(x){
  cat("hurst \n")
  ret <- pracma::hurstexp(x, d = 50, display = FALSE)
  class(ret) <- "hurst"
  ret
}

#' Variance wrapper
#'
#'
#' @param x The data 
#'
#' @return The variance
#' @export
variance <- function(x){
  cat("var \n")
  ret <- var(x)
  class(ret) <- "variance"
  ret
}

#' Compute epsilon complexity using lifting 
#'
#'
#' @param x The data 
#'
#' @return The features
#' @export
ecomp_lift <- function(x){
  cat("ecomp lift \n")
  res <- ecomplex(x, ds = 5, method = "lift", max_degree = 4)
  class(res) <- "ecomp_lift"
  res
}

#' Compute epsilon complexity using bsplines
#'
#'
#' @param x The data 
#'
#' @return The features
#' @export
ecomp_bspline <- function(x){
  cat("ecomp bspline \n")
  res <- ecomplex(x, ds = 5, method = "bspline", max_degree = 4)
  class(res) <- "ecomp_bspline"
  res
}

#' Compute epsilon complexity using bsplines
#'
#'
#' @param x The data 
#'
#' @return The features
#' @export
ecomp_cspline <- function(x){
  cat("ecomp cspline \n")
  res <- ecomplex(x, ds = 5, method = "cspline")
  class(res) <- "ecomp_cspline"
  res
}


#' Compute spectral entropy
#'
#' Computes the entropy of the binned spectrogram.
#'  Wrapper of ForeCA package function.
#'
#' @param  x Time series
#'
#' @return  return 
#' @export
spectral_entropy <- function(x){
  cat("spectral entropy \n")
  ret <- ForeCA::spectral_entropy(x)
  class(ret) <- "spectral_entropy"
  ret
}

#' Compute bandpower
#'
#' Computes the bandpower on a default set 
#'   
#' @param x The data
#' 
#' @return A data frame of power in each band
#' @export 
bandpower <- function(x){
  cat("Bandpower \n")

  freqs <- list(
  delta = c(0.5,4),
  theta = c(4,8),
  alpha = c(8, 12),
  beta = c(12, 30),
  gamma = c(30, 100))

  if(!is.ts(x)){
    fs <- 1220
  } else {
    fs <- frequency(x)
  }
  res <- bp_pgram(x, fs=fs, freqs=freqs)
  res <- lapply(res, log)
  class(res) <- "bandpower"
  res
}


#' Variogram-based estimate of fractal dimension. 
#' 
#' This function is a wrapper for the fd.estimate 
#' function of the fractaldim package.
#'
#' @param x The data
#'
#' @return The results from fd.estimate
#' @export
fd_variogram <- function(x){
  cat("fd_variogram \n")
  fractaldim::fd.estimate(x, methods = "variogram") 
}
