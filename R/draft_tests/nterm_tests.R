# Check if running
run_test <- function(n){
  groups <- test_functions()
  fs <- lapply(groups$group1, generate) 
  fxs <- lapply(fs, function(x) x(n))
  res <- lapply(fxs, function(x) nterm_approx(x, 
                      slope_only = FALSE, 
                      wf = "la8", 
                      n.levels = 4, 
                      boundary = "periodic") ) 
 res2 <- lapply(res, function(x) x$coefficients[2] )
 errs <- lapply(res, function(x) x$errors)
 return(list(funcs = groups, errors = errs, res = res, coeffs = res2))
}

