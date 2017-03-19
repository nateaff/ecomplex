#' Create an ARMA model.
#'
#' @param ar x parameter.
#' @param ma Moving average parameter.
#'
#' @return An ARMA model.
#' @export
arma <- function(ar, ma){
    structure(list(ar=ar, ma = ma), class = "arma")
}


#' Create an Weierstrass function model.
#'
#' Creates a model for a Weierstrass or random 
#'  phase Weierstrass function. If random = FALSE
#'  the function is evaluated on (0,1), otherwise
#'  on (1,n).  
#'
#' @param a The amplitudes parameter.
#' @param b The frequencies parameter.
#' @param random If true a random phase ~unif(0,1) is added.
#' @return The parametrized model.
#' @export
weierstrass <- function(a, b, random = TRUE){
    structure(list(a = a, b = b, random = random), class = "weierstrass")
}

#' Single parameter Weierstrass function model.
#'
#' Creates a model for a Weierstrass or random 
#'  phase Weierstrass function with single parameter
#'  a. The resulting function is a-Holder continuous.
#'  If random, theu function is evaluated on (1,n), 
#'  otherwise on (0,1), 
#'
#' @param a The holder constant of the function. Equal to 
#'           -log(a)/log(b) for the two parameter Weierstrass
#' @param random If true generates a random phase model.
#' @return The parametrized model.
#' @export
weierstrass_a <- function(a, random = TRUE){
    structure(list(a = a, random = random), class = "weierstrass2")
}


#' Create a FARIMA model.
#'
#' @param ar Autoregressive parameters.
#' @param ma Moving average parameters.
#' @param d Long term memory parameter.
#' @return A parametrized FARIMA model.
#' @export
farima <- function(ar, ma, d){
  # check conditions on farima
  structure(list(ar=ar, ma = ma, d = d), class = "farima")
}


#' Create a logistic function model.
#'
#' @param r  Model parameter.
#' @return A logistic function model.
#' @export
logistic <- function(r){
  structure(list(r = r), class = "logistic")  
}

#' Create a quadratic function model. 
#'
#'
#' @param x0 Initial value (TODO: ).
#' @return A quadratic function model.
#' @export
quadratic <- function(x0){
  structure(list(x0 = x0), class = "quadratic")
}



#' Generate a function from a model
#'
#' @param mod A function model
#' @return Returns a function based on the model 
#'  parameters in mod 
#' @export
gen <- function(mod) UseMethod("gen")

#' @export
gen.farima <- function(mod){
  ar = mod$ar; ma = mod$ma; d = mod$d
  function(n){
    fracdiff::fracdiff.sim(n, ar = ar, ma = ma, d = d)$series
  }
}

#' @export 
gen.weierstrass <- function(mod){ 
  a <- mod$a 
  b <- mod$b
  random <- mod$random
  if(!(0 < a) || !(a < 1)) warning("Parameter 'a' is not in (0,1)") 
  if( a*b < 1) warning("'a*b' is not greater than 1")  
  n = 0:50
  theta = 0
  w <- function(y){
    f <- function(x){
         if(random){ theta <- runif(length(n));  
          }
      sum((a^n)*cos(2*pi*(b^n)*x + theta))
      }
    unlist(Map(f,y))
  }
  if(random) cat("true \n")
  function(n){
    xx <- 1:n
    if(!random){ xx <- seq(0, 1, length.out = n)}
    w(xx)
  }
}

#' @export
gen.weierstrass2 <- function(mod){ 
  a <- mod$a 
  random <- mod$random
  if(!(0 < a) || !(a < 1)) warning("Parameter 'a' is not in (0,1)") 
  n = 0:100
  theta = 0
  w <- function(y){
    f <- function(x){
      if(random){ theta <- runif(length(n)) }
      sum((2^(-a*n)*cos((2^n)*2*pi*x + theta)))
    }
  unlist(Map(f,y))
  }
  function(n){
    xx <- 1:n
    if(!random){ xx <- seq(0,1, length.out = n) }
    w(xx)
  }
}

#' @export
gen.arma <- function(mod){
  # Does not check for unit roots
  phi <- mod$ar 
  theta <- mod$ma 
  sd <- mod$sd
  burnin <- 200
  function(n){
    p <- length(phi)
    q <- length(theta)
    n1 = burnin + n;
    a = rnorm(n1 + q, 0,1);
    z <- double(p)

    for(i in seq_along(1:n1)){
       zt = z[i:(i+p-1)]*phi[p:1] + a[i+q] - a[i:(i+q-1)]*theta[q:1];
       z <- c(z, zt)
    }
    z[(burnin + 1 + p):(n1 + p)];
  
  }
}


#' Create a Mackey-Glass model
#'
#' @param tau The time lag 
#' @param init A value used to initialize the discrete Mackeyglass
#'              equation
#' @param beta  The beta parameter
#' @param gamma The gamma parameter
#' @param n     The n parameter
#' @param noise The noise ratio to add to the function
#' @return A model for a Mackey-Glass equation
#' @export
mackeyglass <- function(tau = 17, 
                        beta = 0.25, 
                        gamma = 0.1, 
                        n = 10, 
                        init = 0.5, 
                        noise = 0){

  structure(list( tau = tau, 
                  beta = beta, 
                  gamma = gamma, 
                  n = n,
                  noise = noise,
                  init = init), 
                  class = "mackeyglass")
}


#' @export 
gen.mackeyglass <- function(mod){
  tau   <- mod$tau
  init  <- mod$init
  gamma <- mod$gamma
  beta  <- mod$beta 
  n     <- mod$n
  noise <- mod$noise
  burnin <- max(tau, 500)
   
  # return function
  function(N) {
  len   <- burnin + N
  # initialize to length of delay
  x0 <- rep(init, tau)
  xx <- double(len)
  xx[1] = x0[tau] + beta*x0[1]/(1 + x0[1]^(n)) - gamma*x0[tau]
  
  for(t in 2:tau){
    xx[t] <= xx[t-1]  + beta*x0[t]/(1 + x0[t]^n) - gamma*xx[t-1] 
  }
    for(t in (tau + 1):len){
      xx[t] <- xx[t-1] + beta * xx[t-tau]/(1 + xx[t-tau]^n) - gamma*xx[t-1]
    }
    xx <- xx + rnorm(len)*noise*sd(xx)
    xx[(burnin+1):len]
  }
}



#' @export
gen.logistic <- function(mod){
  r <-mod$r
  function(n, x0 = 0.1){
    xx <- vector("double", n)
    xx[1] <- x0
    
    for(k in seq_along(xx)){
      if(k != 1){
        xx[k] <- r*xx[k-1]*(1 - xx[k-1])
      }
    }
    return(xx)
  }
}

#' @export
gen.quadratic <- function(mod){
  x0 <- mod$x0
  function(n){
    xx <- vector("double", n)
    xx[1] <- x0
    
    for(k in seq_along(xx)){
      if(k != 1){
        xx[k] <- xx[k-1]^2 -2
      }
    }
    return(xx)
  }
}

#' @export
gen.default <- function(mod){
  warning(sprintf("No method for class %s found", class(mod)))
  numeric(0)
}



