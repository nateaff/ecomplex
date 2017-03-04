#' Extract set of features from a time series.
#'
#' @param data A time series, matrix of dataframe. 
#' @param features A list or vector of ts_features.
#' @param prefix Optional prefix for column headers.
#' @return A data frame with row of features for each time.
#'  series in data. (TODO: single feature, system 
#'  check for parallelization)
#'  
#' @export 
extract_features <- function(data, features, id = "1", para = TRUE, verbose = FALSE, ...){

  # TODO: switch for system type
  if(para){
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
      if(para){
        ret <- parallel::mclapply(data, function(x) extract_one_feature(x, feature), 
                        mc.cores = cores)
      } else {
        ret <- lapply(data, function(x) extract_one_feature(x, feature))
      }
    k <- k + 1
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


clean_feature <- function(x, ...) UseMethod("clean_feature")


clean_feature.FractalDim <- function(feature){
  res <- plyr::unrowname(data.frame(fd = feature$fd)) 
  names(res) <- paste0("fd_", feature$methods)
  res
}


clean_feature.fd_variogram <- function(feature){
  res <- plyr::unrowname(data.frame(fd = feature$fd)) 
  names(res) <- paste0("fd_", feature$methods)
  class(res) <- "fd_variogram"
  res
}

clean_feature.bandpower <- function(feature){
  res <- plyr::unrowname(data.frame((t(unlist(feature)))))
  names(res) <- unlist(lapply(1:length(res), function(x) paste0("bp", x)))
  res
}

clean_feature.nterm <- function(feature){ 
  plyr::unrowname(data.frame(nterm_coeff = feature$coefficients[2])) 
}

clean_feature.ecomplexity <- function(feature){
 res <-  plyr::unrowname(data.frame(ecomp_coeff = feature$fit$coefficients[2])) 
 names(res) <- paste0("fd_", feature$method)
 res
}

clean_feature.ecomp_bspline <- function(feature){
  plyr::unrowname(data.frame(ecomp_bspline = feature$fit$coefficients[2])) 
}

clean_feature.ecomp_cspline <- function(feature){
  plyr::unrowname(data.frame(ecomp_cspline = feature$fit$coefficients[2])) 
}


clean_feature.ecomp_adapt <- function(feature){
  plyr::unrowname(data.frame(ecomp_adapt = feature$fit$coefficients[2])) 
}

clean_feature.sample_entropy <- function(feature){
  plyr::unrowname(data.frame(sample_entropy = feature$sample_entropy)) 
}

clean_feature.nterm3 <- function(feature){
  plyr::unrowname(data.frame(nterm3 = feature$coefficients[2]))
}

clean_feature.nterm2 <- function(feature){
  plyr::unrowname(data.frame(nterm2 = feature$coefficients[2]))
}


clean_feature.default <- function(feature){
  plyr::unrowname(data.frame(feature = feature))
}


clean_feature.wvar <- function(feature){
  wav_var <- feature$variance[1:6]
  fnames <- paste0("wvar_", feature$scales[1:6])
  ret <- plyr::unrowname(data.frame(t(wav_var)))
  names(ret) <- fnames
  ret
}
 

# change class to differentiate features
sample_entropy2 <- function(xx){
  res <- pracma::sample_entropy(xx)
  structure(list(sample_entropy = res), class = "sample_entropy")
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

#' Compute bandpowep
#'
#' Computes the bandpower on a default set 
#'   
#' @param xx The data
#' 
#' @return A data frame of power in each band
#' @export
#' 
bandpower <- function(xx){
  cat("Bandpower \n")

  freqs <- list(
  delta = c(0.5,4),
  theta = c(4,8),
  alpha1 = c(8, 12),
  beta1 = c(12, 20),
  beta2 = c(20, 40),
  gamma1 = c(40, 60), 
  gamma2 = c(60, 100))

  if(!is.ts(xx)){
    fs <- 1220
  } else {
    fs <- frequency(xx)
  }
  res <- bp_pgram(xx, fs=fs, freqs=freqs)
  class(res) <- "bandpower"
  res
}

ecomp_adapt <- function(xx){
  res <- ecomplexity(xx, ds = 5, method = "adlift")
  class(res) <- "ecomp_adapt"
  res
}

#' fd variogram
#'
#' fd 
#'
#' @param xx The data
#'
#' @return The results from fd.estimate
#' @export
#' 
fd_variogram <- function(feature){
  cat("fd_variogram \nterm")
  fractaldim::fd.estimate(feature, methods = "variogram") 
}

#Temp 

#' Extract eeg features 
#'
#' @param x The number of time series generated for each group
#' @param features The features set to run 
#' @param parallel Parallelize computations (Linux only)
#'
#' @return A data frame with features as columns and an id
#'  column identifying group membership
#' @export
eeg_features <- function(x, features = NULL, parallel = TRUE){
  if (is.null(features)){
  features <- c(ecomp_bspline, bandpower, fd_variogram)
  }
  df <- extract_features(x, features, id = "1", para = para)
  df
} 


#' Test a default suite of time series features.
#'
#' Extracts default set of time series features from
#'  two groups of ARMA(2,2) processes. ARMA processes 
#'  are generated using arima.sim() using different initial
#'  parameters for each group.
#'
#' @param n The number of time series generated for each group.
#' @param seed The seed for the random number generator.
#'
#' @return A data frame with features as columns and an id
#'  column identifying group membership.
#' @export
test_features <- function(n = 10, seed = 2017, features =  NULL){
  set.seed(seed)
  if(is.null(features)){
  features <- c(nterm, sample_entropy2, ecomp_bspline, 
                bandpower, fd_variogram, gmwm::wvar)
  }

  ts_mat1 <- replicate(n, arima.sim(n = 200, list(ar = c(0.8897, -0.4858), 
                ma = c(-0.2279, 0.2488)), sd = sqrt(0.1796)))
  
  ts_mat2 <- replicate(n, arima.sim(n = 200,
   list(ar = c(-0.71, 0.18), ma = c(0.92, 0.14)), sd = sqrt(0.291)))


  
  df1 <- extract_features(ts_mat1, features, id = "1")
  df2 <- extract_features(ts_mat2, features, id = "2")
  comb <- rbind(df1, df2)
  comb 
}

