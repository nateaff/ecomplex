---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# ecomplex

R package for efficiently computing the epsilon-complexity coefficients of a time series. The coefficients are a 
measure of the complexity of a time series.

# An example 

The complexity coefficients have been used as a feature in classification and clustering tasks. Here we look at how well the coefficients discriminate between two sets of time serires drawn from ARMA(2,2) models with one AR parameter changed.


```r
library(ecomplex)
set.seed(1)
reps   <- 30
group1 <- replicate(30, arima.sim(n = 500, list(ar = c(0.89, -0.49), 
                                                ma = c(-0.23, 0.25))))
         

group2 <- replicate(30, arima.sim(n = 500, list(ar = c(0.69, -0.49),
                                                ma = c(-0.23, 0.25))))
ecomp1 <- apply(group1, 2, ecomplex)
ecomp2 <- apply(group2, 2, ecomplex)
coeffs <- lapply(c(ecomp1, ecomp2), function(fit) c(fit$A, fit$B))

df     <- data.frame(do.call(rbind, coeffs))
df$id  <- factor(rep(c(1,2), each = 30))
names(df) <- c("A", "B", "id")

palette(c("gray20", "chocolate3"))
with(df, plot(B, A, col = id, lwd = 2.5))
```

![plot of chunk arimasim](tools/README-arimasim-1.png)

And a plot of a sample from each ARMA(2,2) group on the same 
axis. 


```r
plot(c(group1[,1]), xlim = c(0,1000), ylab = "", 
                                      main = title, 
                                      col = "gray20", 
                                      lwd = 1.5, type ='l')
lines(501:1000, group2[,1], col = "chocolate3", lwd = 1.5)
abline(v = c(500), lwd = 4, col = "gray20")
```

![plot of chunk ts](tools/README-ts-1.png)

```r
palette("default")
```

# The basic algorithm

The `ecomplex` function successively down samples and approximates a time series. The coefficients are the parameters of a log-log regression of the set of approximation errors on the fraction of sample points retained for each approximation. Roughly, the coefficients measure the amount of information (the fraction of samples) needed to approximate a function within some error epsilon.

See Darkhovsky and Piryatinska, [Binary classification of multi-channel EEG records based on the epsilon-complexity of continuous vector functions](https://arxiv.org/pdf/1610.01633.pdf).

# Installation


```r
# install.packages('devtools')
devtools::install_github('nwaff/ecomplex')
```
# Future work 

The package is in development but the `ecomplex` interface should be fairly stable. The `ecomplex` function computes the error on single variable time series but an option might be added to compute the epsilon-complexity coefficients for a multivariate time series.  

The `palarm` function included in the package is a change point detection algorithm. There will likely be minor changes to return 
type of the function in the near future.
