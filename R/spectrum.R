
#' Find the power within specified frequency bands.
#'
#' Estimates the power within the given frequency 
#'  bands using Welch's method. The periogram is 
#'  then smoothed using wavelet thresholding. The 
#'  sum of power in bins within the frequency band
#'  is used as the estimate of that band's power.
#'
#' @param x The input signal
#' @param fs The signal sampling rate
#' @param freqs A vector or list of vectors of the
#'  bandwidth's whose power should be computed
#' @param psd The output from spec.pgram(). 
#' @param plot_pgram Plots default periodogram if TRUE 
#' @return  A double or list of doubles representing the
#'  signal's frequency in the given bandwidths
#' 
#' @export
#' 
bp_pgram <- function(x, fs, freqs, 
                    psd = NULL, 
                    plot_pgram = FALSE){
  if(!is.list(freqs)){
    freqs <- list(freqs)
  }
  if(!is.ts(x)){ 
    x <- ts(x, frequency = fs)
  }
  if(is.null(psd)){ 
    pgram <- spec.pgram(x, plot = plot_pgram, n.used = 100)
  }
  lapply(freqs, function(x) bandpower_one(pgram, x)) 
}


bandpower_one <- function(pgram, bin){
  ran <- which(pgram$freq > bin[1] & pgram$freq < bin[2])
  sum(pgram$spec[ran])
}





