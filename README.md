
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/nateaff/ecomplex.svg?branch=master)](https://travis-ci.org/nateaff/ecomplex)

ecomplex
========

R package for efficiently computing the epsilon-complexity coefficients of a time series. The coefficients are estimates of the complexity of a time series. The complexity coefficients are computed by finding a trend in the best approximation error as the series is successively downsampled and approximated.

An example
==========

The complexity coefficients can be used a feature in classification and clustering tasks. Here we generate two groups of time series, each group drawn from an ARMA(2,2) model with one parameter changed. A plot of the coefficients for each series shows the two groups are fairly well separated in the coefficient space.

``` r
library(ecomplex)
set.seed(1)
reps   <- 100; n <- 500
group1 <- replicate(reps, arima.sim(n = n, 
                              list(ar = c(0.89, -0.49), 
                                   ma = c(-0.23, 0.25))))
group2 <- replicate(reps, arima.sim(n = n, 
                              list(ar = c(0.69, -0.49),
                                   ma = c(-0.23, 0.25))))
ecomp1 <- apply(group1, 2, ecomplex)
ecomp2 <- apply(group2, 2, ecomplex)
coeffs <- lapply(c(ecomp1, ecomp2), function(fit) c(fit$A, fit$B))

df     <- data.frame(do.call(rbind, coeffs))
df$id  <- factor(rep(c(1,2), each = reps))
names(df) <- c('A', 'B', 'id')
palette(c('gray20', 'chocolate3'))
with(df, plot(B, A, col = id, lwd = 2))
```

![](figures/README-arimasim-1.png)

An example from each time series group plotted on the same axis.

``` r
plot(c(group1[1:500,1]), xlim = c(0,1000), 
                         ylab = '', 
                         col = 'gray20', 
                         lwd = 1.2, 
                         type ='l')
lines(501:1000, group2[, 1], col = 'chocolate3', lwd = 1.2)
abline(v = c(500), lwd = 3, col = 'gray20')
```

![](figures/README-ts-1.png)

``` r
palette('default')
```

The basic algorithm
===================

The `ecomplex` function successively down samples and approximates a time series. The coefficients are the parameters of a log-log regression of the set of approximation errors on the fraction of sample points retained for each approximation. Roughly, the coefficients measure the amount of information (in terms of sample points) needed to approximate a function within some error epsilon.

For mathematical details see Darkhovsky and Piryatinska, [Binary classification of multi-channel EEG records based on the epsilon-complexity of continuous vector functions](https://arxiv.org/pdf/1610.01633.pdf).

Installation
============

``` r
# install.packages('devtools')
devtools::install_github('nwaff/ecomplex')
```

Future work
===========

The package is in development but the `ecomplex` interface should be fairly stable. The `ecomplex` function computes the error on single variable time series but an option might be added to compute the epsilon-complexity coefficients for a multivariate time series.

The `palarm` function included in the package is a change point detection algorithm. There will likely be minor changes to return type of the function in the near future.
