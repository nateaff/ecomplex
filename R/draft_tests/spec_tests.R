#=====================================
# Test functions
#=====================================

bp_fft <- function(){
  x <- ts(fex1()$x, frequency = 1000)
  P = abs(2*fft(x)/100)^2; 
  Fr = 0:((length(x)-1)/2)
  plot(Fr, P[1:length(P)], type = 'l', xlab = '  frq')
}

fex1 <- function(){
# example from shumway p. 178
  x1 <- 2*cos(2*pi*1:1000*6/1000) + 3*sin(2*pi*1:1000*6/1000)
  x2 <- 4*cos(2*pi*1:1000*10/1000) + 5*sin(2*pi*1:1000*10/1000)
  x3 <- 6*cos(2*pi*1:1000*40/1000) + 7*sin(2*pi*1:1000*40/1000)
  x <- x1 + x2 + x3
  return(list(x1 = x1, x2 = x2, x3 = x3, x = x))
}


test_pgram <- function(){
  set.seed(1)
  fs <- 200
  fxs <- fex1()
  x <- ts(fxs$x, frequency = fs)

  freqs <- list(
    delta <- c(1,4),
    theta <- c(4,8),
    alpha1 <- c(8, 12),
    beta1 <- c(13, 20),
    beta2 <- c(20, 40),
    gamma <- c(40, 100))

  bp_pgram(x, fs, freqs, plots = TRUE)   
}

