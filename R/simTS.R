#' Create an ARMA model.
#'
#' @param ar Autoregressive parameter.
#' @param ma Moving average parameter.
#' @param d Differencing parameter. (TODO: not used)
#' @return An ARMA model.
#' @export
arma_mod <- function(ar, ma, d ){
    structure(list(ar=ar, ma = ma, d = d), class = "arma")
}

#' Create a FARIMA model.
#'
#' @param ar The autoregressive parameters.
#' @param ma The moving average parameters.
#' @param d The long term memory parameter.
#' @return A function model of class "farima".
#' @export
farima_mod <- function(ar, ma, d){
  # check conditions on farima
  structure(list(ar=ar, ma = ma, d = d), class = "farima")
}


#' Create a logistic function model.
#'
#' @param r  Model parameter.
#' @return A logistic function model.
#' @export
logistic_mod <- function(r){
  structure(list(r = r), class = "logistic")  
}

#' Create a Mackey-Glass model
#'
#' @param tau The time lag 
#' @param init A vector used to initialize the Mackey-Glass function
#' @return A model for a Mackey-Glass equation
#' @export
mackeyglass_mod <- function(tau, init){
  structure(list(tau = tau, init = init), class = "mackeyglass")
}

#' Create a quadratic function model. 
#'
#'
#' @param x0 Initial value (TODO: ).
#' @return A quadratic function model.
#' @export
quadratic_mod <- function(x0){
  structure(list(x0 = x0), class = "quadratic")
}



#' Generate a function from a model
#'
#' @param mod A function model
#' @return Returns a function based on the model 
#'  parameters in mod 
#' @export
generate <- function(mod) UseMethod("generate")

#' @export
generate.farima <- function(mod){
  ar = mod$ar; ma = mod$ma; d = mod$d
  function(n){
    fracdiff::fracdiff.sim(n, ar = ar, ma = ma, d = d)$series
  }
}

#' @export
generate.arma <- function(mod){
  phi    <- mod$ar 
  theta  <- mod$ma
  burnin <- 200 
  sd <- 1
  function(n){
    p  <- size(phi)[1]
    q  <- size(theta)[1]
    n1 <- burnin + n;
    a  <- rnorm(n1+q, 0,sd);
    z  <- double(p)

    for(i in seq_along(1:n1)){
       zt = z[i:(i+p-1)]*phi[p:1] + a[i+q] - a[i:(i+q-1)]*theta[q:1];
       z <- c(z, zt)
    }
    z[(burnin + 1 + p):(n1 + p)];
  
  }
}


#' @export
generate.mackeyglass <- function(mod){
  stopifnot(class(x) == "mackeyglass")

  tau <- mod$tau
  # check appropriate range of yinit?
  init <- mod$init
  function(n){
    # add burnin 
    burnin <- 50
    end <- burnin+n
    
    mg <- function(t, y, parms, tau) {
    lag <- t - tau

    if (lag <= 0)
      lag <- 0.5
    else 
    lag <- deSolve::lagvalue(lag)

    dy <- 0.2 * lag * 1/(1+lag^10) - 0.1 * y
    list(c(dy))
    }
    
    times <- seq(from = 0, to = end*100, by = 100)
    y <- deSolve::dede(y = init, times = times, 
              func = mg, parms = NULL, tau = tau)
    y[(burnin+1):end, 2]
  }
}


#' @export
generate.logistic <- function(mod){
  r <-mod$r
  function(n, x0 = 0.1){
    res <- vector("double", n)
    res[1] <- x0
    
    for(k in seq_along(res)){
      if(k != 1){
        res[k] <- r*res[k-1]*(1-res[k-1])
      }
    }
    return(res)
  }
}

#' @export
generate.quadratic <- function(mod){
  x0 <- mod$x0
  function(n){
    res <- vector("double", n)
    res[1] <- x0
    
    for(k in seq_along(res)){
      if(k != 1){
        res[k] <- res[k-1]^2 -2
      }
    }
    return(res)
  }
}


#' Generate A Default suite of function models.
#'
#' @return A list of two lists of function models with
#'          the same functions and varying default parameters.
#' @export
fsuite <- function(){
  # The time serires models : arma, arma, farima, logistic, quadratic, Mackey-Glass

  # Farima
  farima1 <- farima_mod(ar = c(0.1, -0.5), ma = c(0.6, 0.01), d = 0.35) 
  farima2 <- farima_mod(ar = c(0.2, -0.4), ma = c(0.4, 0.02), d = 0.35)  
  # Arma1
  arma1 <- arma_mod(ar = c(-0.1, 0.3, 0.1), ma = c(0.2, 0.1), d = 0)
  arma2 <- arma_mod(ar = c(0.4, 0.3, 0.2), ma = c(0.1, -0.5), d = 0)
  # Arma2
  arma3 <- arma_mod(ar  = c(0.5, -0.6, 0.9), ma = c(0.5, 0.6), d = 0)
  arma4 <- arma_mod(ar = c(-0.2, -0.3, -0.8), ma = c(0.2, 0.1), d = 0)
  # Logistic 
  log1 <- logistic_mod(r = 3.98)
  log2 <- logistic_mod(r = 3.87)
  # Mackey-Glass
  mg1 <- mackeyglass_mod(tau = 100, init = rep(0.5))
  mg2 <- mackeyglass_mod(tau = 200, init = rep(0.5))
  # Quadratic  
  quad <- quadratic_mod(x0 = 0.2)

  group1 <- list(arma1 = arma1, 
                 arma2 = arma2, 
                 log1 = log1, 
                 quad = quad, 
                 # mg1 = mg1, 
                 farima1 = farima1)
  group2 <- list(arma3 = arma3, 
                 arma4 = arma4, 
                 log2 = log2, 
                 quad = quad, 
                 # mg2 = mg2, 
                 farima2 = farima2)

  return(list(group1 = group1, group2 = group2))
}


