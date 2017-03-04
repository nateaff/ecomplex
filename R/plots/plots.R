
library(ggplot2)
library(reshape2)
library(ggthemes)
    
#' Diagnostic plot of time series features  
#'
#' Diagnostic plot of time series fatures returned by 
#'  (todo: name) test_features()
#'
#' @param data A data frame. 
#'   Data frame should be formatted with features in columns and an
#'   id column for   
#'   | feature_name1 | feature_name2 | ... | id | 
#'  
#'
#'
#'@export
plot_features <- function(data, type = c("histogram", "density")){
  # color blind friendly from: 
  # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

  data_long <- melt(data, id.vars = "id")

  type <- match.arg(type)
  switch(type,
  histogram = { gg <- ggplot(data_long, aes(x=value, fill=id)) + 
      geom_histogram(alpha = 0.5, aes(y = ..density..), 
                      bins = 25, 
                      position = "identity") },
   density = {  gg <- ggplot(data_long, aes(x=value, fill=id)) + 
                      geom_density(alpha = 0.5, aes(y = ..density..), 
                                   color = NA, 
                                   position = "identity")}
  ) #end switch
  gg + facet_wrap(~variable, ncol = 3, scales = "free") + 
     scale_fill_manual(values = cbPalette) +
     theme_minimal()
    #theme_par()
    #theme_hc()
    #theme_base()
}

#=============
# Plotting
# ============

# plot uses result$epsilons

plot_coeffs <- function(results, func_name){
  
  plot(log(1/(2:6)), log(results$epsilons), 
    pch  = 16, cex = 1, 
    col  = "darkblue", 
    ylab = "Log(epsilons)", 
    xlab = "Log(S)", 
    main = func_name)
  
  fit <- lm(log(results$epsilons) ~ log(1/(2:6))) 
  abline(fit, col = "gray20", lwd = "2")
}

# plot using using epsilons
plot_coeffs2 <- function(epsilons, ds, func_name){
  
  plot(log(1/(2:ds)), log(epsilons), 
    pch  = 16, cex = 1, 
    col  = "darkblue", 
    ylab = "Log(epsilons)", 
    xlab = "Log(S)", 
    main = func_name)
  
  fit <- lm(log(epsilons) ~ log(1/(2:ds))) 
  abline(fit, col = "gray20", lwd = "2")
}


plot_approx <- function(ys, yfit){
  plot(ys, type = "l", col = "gray20", 
                      lwd = 1.2, lty = 3, 
                      xlab = "", ylab = "")
  par(new = TRUE)
  plot(yfit, type = "l", lwd = 1.6, col = "darkblue", 
            xaxt = "none", yaxt = "none", ylab = "", xlab ="")
}


#===========================
# nterm fits  (old)
#==========================


# Used to plot multiple versions of fit
# Needs update
plotNtermFits <- function(ys, fits, pal){

  plot(ys[1:100], type = "l", col = "black", lwd = 1.7, lty = 1.7, xlab = "", ylab = "")
  
  for (k in 1:length(fits)) {
    par(new = TRUE)
    plot(fits[[k]][1:100], type = "l", lwd = 1.6, col = pal[k], 
              xaxt = "none", yaxt = "none", ylab = "", xlab ="")
  }
  par(new = TRUE)
   plot(ys[1:100], type = "l", col = "black", 
                lwd = 1.8, lty = 2, xlab = "", ylab = "")
 
}



# Plot linear fit of all functions
plot_errs <- function(y, fname){
  plot(log(y), 
    pch  = 16, cex = 0.8, 
    col  = "gray40", 
    ylab = "log(epsilons)", 
    xlab = "n", 
    main = fname)
  
  x <- 1:length(y)
  fit <- lm(log(y) ~ x) 
  abline(fit, col = "darkblue", lwd = "2")
}

