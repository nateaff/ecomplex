#' Create an ARMA model.
#'
#' @param ar Autoregressive parameter.
#' @param ma Moving average parameter.
#' @param d Differencing parameter. (TODO: not used)
#' @return An ARMA model.
#' @export
arma <- function(ar, ma, d ){
    structure(list(ar=ar, ma = ma, d = d), class = "arma")
}


#' Create an Weierstrass function model.
#'
#' Creates a model for a Weierstrass or random 
#'  phase Weierstrass function 
#'
#' @param a The amplitudes parameter.
#' @param b The frequencies parameter.
#' @param random Create a random phase model.
#' @return The parametrized model.
#' @export
weierstrass <- function(a, b, random = TRUE){
    structure(list(a = a, b = b, random = random), class = "weierstrass")
}


#' Create a FARIMA model.
#'
#' @param ar The autoregressive parameters.
#' @param ma The moving average parameters.
#' @param d The long term memory parameter.
#' @return A function model of class "farima".
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

#' Create a Mackey-Glass model
#'
#' @param tau The time lag 
#' @param init A vector used to initialize the Mackey-Glass function
#' @return A model for a Mackey-Glass equation
#' @export
mackeyglass <- function(tau, init){
  structure(list(tau = tau, init = init), class = "mackeyglass")
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
  if( a*b < 1) warning("Parameter 'b' is not greater than 1")  
  n = 1:20 
  theta = 0
  
  function(y){
    f <- function(x){
      if(random){ theta <- runif(length(n)) }
      sum((a^n)*cos(2*pi*((b^n)*x + theta)))
    }
  unlist(Map(f,y))
  }
}

#' @export
gen.arma <- function(mod){
  # Does not check for unit roots
  phi <- mod$ar 
  theta <- mod$ma 
  function(n){
    p <- size(phi)[1]
    q <- size(theta)[1]
    n1 = 200+n;
    a = rnorm(n1+q, 0,1);
    z <- double(p)

    for(i in seq_along(1:n1)){
       zt = z[i:(i+p-1)]*phi[p:1] + a[i+q] - a[i:(i+q-1)]*theta[q:1];
       z <- c(z, zt)
    }
    z[(201+p):(n1+p)];
  
  }
}


#' @export
gen.mackeyglass <- function(mod){
  # check appropriate range of inits
  cat("This takes a second ... \n") 
  tau <- mod$tau
  init <- mod$init
  function(n){
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
gen.logistic <- function(mod){
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
gen.quadratic <- function(mod){
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

  farima1 <- farima(ar = c(0.1, -0.5), ma = c(0.6, 0.01), d = 0.3) 
  farima2 <- farima(ar = c(0.2, -0.4), ma = c(0.4, 0.02), d = 0.3)  

  arma1 <- arma(ar = c(-0.1, 0.3, 0.1), ma = c(0.2, 0.1), d = 0)
  arma2 <- arma(ar = c(0.4, 0.3, 0.2), ma = c(0.1, -0.5), d = 0)

  arma3 <- arma(ar  = c(0.5, -0.6, 0.9), ma = c(0.5, 0.6), d = 0)
  arma4 <- arma(ar = c(-0.2, -0.3, -0.8), ma = c(0.2, 0.1), d = 0)

  log1 <- logistic(r = 3.98)
  log2 <- logistic(r = 3.87)

  mg1 <- mackeyglass(tau = 100, init = rep(0.5))
  mg2 <- mackeyglass(tau = 200, init = rep(0.5))

  weier1 <- weierstrass(a = 0.4, b = 4)
  weier2 <- weierstrass(a = 0.8, b = 4)


  group1 <- list(arma1 = arma1, 
                 arma2 = arma2, 
                 log1 = log1, 
                 weier = weier1, 
                 mg1 = mg1, 
                 farima1 = farima1)
  group2 <- list(arma3 = arma3, 
                 arma4 = arma4, 
                 log2 = log2, 
                 weier2 = weier2, 
                 mg2 = mg2, 
                 farima2 = farima2)

  return(list(group1 = group1, group2 = group2))
}


