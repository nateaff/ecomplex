
<!-- Note: edits should be made in the READ.Rmd file -->
ecomplex
========

R package for efficiently computing the epsilon-complexity coefficients of a time series. The coefficients are intended to characterize the functional complexity or regularity of a time series.

An example
==========

The complexity coefficients have been used as a feature in classification and clustering tasks. Here we look at how well the coefficients discriminate between two sets of time serires drawn from ARMA(2,2) models with one AR parameter changed.

``` r
library(ecomplex)
set.seed(1)
reps   <- 30
group1 <- replicate(30, arima.sim(n = 500, list(ar = c(0.89, -0.49), 
                                                ma = c(-0.23, 0.25))))
         

group2 <- replicate(30, arima.sim(n = 500, list(ar = c(0.59, -0.49),
                                                ma = c(-0.23, 0.25))))
ecomp1 <- apply(group1, 2, ecomplex)
ecomp2 <- apply(group2, 2, ecomplex)
coeffs <- lapply(c(ecomp1, ecomp2), function(fit) c(fit$A, fit$B))

df     <- data.frame(do.call(rbind, coeffs))
df$id  <- factor(rep(c(1,2), each = 30))
names(df) <- c("A", "B", "id")
with(df, plot(B, A, col = id))
```

![](tools/README-arimasim-1.png)

Basic algorithm
===============

The `ecomplex` function successively down samples and approximates a time series. The minimum absolute error fit at each downsample is calculated. A log-log linear model is fit to the errors regressed against the percentage of data used in the approximation and the coefficients of the linear model are the two epsilon-complexity coefficients.

Denote our time series \(x_t\) then the \(\varepsilon\)-complexity coefficients are computed as follows:

-   The observations \(x_t\) are sampled on a grid \(\alpha \mathbb{Z}\), with \(\alpha \in \{ 2,3,4,5, 6 \}\)
-   For each \(\alpha\) find \(f_{\alpha}\) in a set of approximation methods \(F\) that has the minimum error \(|f_{\alpha} - x_t| = \varepsilon_{\alpha}\)
-   Let \(S_{\alpha} = 1/\alpha\). Then the \(\varepsilon-\)complexity coefficients \(A\) and \(B\) of the linear fit \[
      \log(\varepsilon_{\alpha}) = A + B \cdot \log(S_{\alpha}).
    \]

See [Binary classification of multi-channel EEG records based on the epsilon-complexity of continuous vector functions](https://arxiv.org/pdf/1610.01633.pdf)

Installation
============

``` r
# install.packages('devtools')
devtools::install_github('nwaff/ecomplex')
```

Future work
===========

The package is in development but the `ecomplex` interface should be stable. The package also includes the `palarm` function which is a change point algorithm. The interface and return types of that function will likely be updated in the near future. The simulation and feature extraction functions were used in a separate analysis and will be moved to a separate package to reduce the number of package dependencies.
