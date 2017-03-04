# library(Rcpp)
# library(pracma)

#-------------------------------------
# Build coefficient matrices for
# interpolating polynomial subdivision
# ------------------------------------

pfilter <- function(order, pos){
    
    N <- order + 1
    eval <- seq(min(pos)-1, max(pos)+1, by = 2)
    cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
    ymat <- diag(N)
    for(j in 1:length(eval)){
        for(k in 1:N){
            ys <- ymat[k,]
            cmat[j,k] <- pracma::neville(pos, ys, eval[j])
        }
    }
    cmat
}



filters1 <- function(){
     pts <- list( linear <- c(-1,1),
                  quad <- c(-1,1,3),
                  cubic <- c(-3,-1,1,3),
                  quart <- c(-3,-1,1,3,5), 
                  quintic <- c(-5,-3,-1,1,3,5), 
                  sextic <-c(1,1,1,1,1,1,1), # dummy
                  septic <- c(-7,-5,-3,-1,1,3,5,7)  )
  orders <- 1:7
  ret <- lapply(orders, function(x) pfilter(x, pts[[x]] ))
  ret
}




pfilter2 <- function(order, pos){
  xx <- (min(pos)-2):(max(pos)+2)
  eval <- xx[!(xx %in% pos)]
  #=============
  N <- order + 1
  cmat <- matrix(0, nrow = length(eval), ncol = length(pos))
  ymat <- diag(N)
  for(j in 1:length(eval)){
      for(k in 1:N){
          ys <- ymat[k,]
          cmat[j,k] <- pracma::neville(pos, ys, eval[j])
      }
  }
  cmat
}


filters2 <- function(){

  # Shifting the points gives the same filter
  pts <- list( linear1 <- c(-1,2),
                # linear2 <- c(-2,1),
                cubic1 <- c(-4,-1,2,5),
                # cubic2 <- c(-5,-2,1,4),
                quintic1 <- c(-7, -4,-1,2,5, 8))
                # quintic2 <- c(-8, -5,-2,1,4, 7))
  orders <- c(1,3,5)
  ret <- lapply(1:length(orders), function(k) pfilter2(orders[k], pts[[k]] ))
  ret
}




lshift <- function(x, t){
    x[(1+t):length(x)]
}


predict <- function(ys, filter){
    sum(ys*filter)
}



iswt_err <- function(ys, ds, deg){
  ind  <- downsample_perm(length(ys), ds);
  # errors for each permutation
  epsilons <- double(length(ind))
  filter <- filters1()[[deg]]
  reps = log(ds,2)
  for (k in 1:ds) {
    tys <- ys[ind[[k]]]
    yout <- ISWT(tys, deg, reps, filter)
    # hack
    adjust <- k + deg^2
    ran <- adjust:(length(yout) - adjust)
    cat("NA's", which(is.na(yout)), "\n")
    # stopifnot(which(is.na(yout))==0)
    epsilons[k]  <- sum(abs(yout[ran] - ys[ran + k -1])) 
  }
  # include last fit
  return(list(mu = mean(epsilons), eps = epsilons, yout = yout))
}







ISWT <- function(xx, deg = 3, reps, filter){
    
    xx <- interp_one(xx, deg, filter)
    reps <- reps -1
    if(reps > 0){
       xx <- Recall(xx, deg, reps, filter)
    }
    xx
}



interleave <- function(odd, even){
    xx <- double(length(even) + length(odd))
    for(k in 1:length(even)){
        xx[2*k -1] <- odd[k]
        xx[2*k] <- even[k]
    }
    xx
}


interp_one <- function(xx, deg, filter){
    # if(length(xx) %% 2 == 0){
    #     xx <- xx[1:(length(xx))]
    # }
    cat("len", length(xx), "\n")
    offset <- ceiling(deg/2) 
    ps <- double(length(xx))
    #Note: handling of edge is not correct for differnt downsample patterns
    #TODO: Handle edges case by case or skip
    for(k in 1:offset){
      ps[k] <- predict(xx[1:(2*offset)], filter[k,])
      ps[length(ps) - offset + k] <- predict(xx[(length(ps)-(2*offset)+1):length(ps)], 
            filter[deg - k,])
    }
    for(k in offset:(length(ps)- offset)){
        ps[k] <- predict(xx[(k - offset + 1):(k+offset)], filter[offset+1,])
    }
    res <- interleave(xx[1:length(xx)], ps[1:length(xx)])[1:(length(xx)*2)]
    res
}



interp2 <- function(xx, deg, filter){
    # if(length(xx) %% 2 == 0){
    #     xx <- xx[1:(length(xx))]
    # }
    # xx <- tys 
    # deg = 1
    # filter <- filters2()[[1]]
    # #====================
    cat("len", length(xx), "\n")
    cat("len filter", dim(filter)[2], "\n")
    offset <- ceiling(deg/2) 
    ps <- double(length(xx)*2)
    xhead <- deg + 1
    xtail <- length(ps)- xhead  
    #NOTE: this is not correct for different downsample patterns
    # for(k in 1:xhead){
    #   ps[k] <- predict(xx[1:(2*offset)], filter[k,])
    #   ps[xtail + k] <- predict(xx[(length(xx)-(2*offset)+1):length(xx)], 
    #         filter[(deg + 2) - k,])
    # }
    xseq <- seq(1, xtail, by = 2)
    for(j in (offset + 1):length(xseq)){
      k <- xseq[j]
      ps[k] <- predict(xx[(j - offset +1 ):(j + offset)], filter[deg + 2,])
      ps[k+1] <- predict(xx[(j - offset + 1):(j + offset)], filter[deg + 3,] )
    }
     ret <- interleave2(xx, ps[1:length(ps)])
     plot(ret)
     ret
}


interleave2 <- function(odd, evens){
    # odd <- c(1:10)
    # evens <- -1:-21
    
    xx <- double(length(evens) + length(odd))

    xseq <- seq(1, length(xx), by = 3)
    for(k in seq_along(odd)){
        # cat("insert", xseq[k], "\n")
        xx[xseq[k]] <- odd[k]
        xx[xseq[k]+1] <- evens[2*k-1]
        xx[xseq[k]+2] <- evens[2*k]
    }
    xx
}



test_interp2 <- function(deg = 3, add_noise = TRUE){
  # deg <- 3
  xx <- seq(0,5, 1/100)
  ys1 <- sin(2*xx) + (1/3)*cos((.33*xx)) + exp(-xx) 
  ys2 <- .8*cos(12*xx) + .5*cos(9*xx) + exp(-xx)
  ys <- ys1 + ys2

  ind <-downsample_perm(length(ys), 3)[[1]]
  tys <- ys[ind]

  if(add_noise){

    tys <- tys + rnorm(length(tys), 0, 0.5)
    }   
  # tys
  deg = 1
  tol = .1
  filter <- filters2()[[1]]
  res <- interp2(tys, deg, filter)
  dif <- res - ys 
  bad <- dif[abs(dif)  > .1]
  plot(res, cex = 0.3, col = "gray10")
  lines(res, col = "gray40", lwd = 1.2)
  # points(res, col = "gray10", cex = 0.5, pch = 16)  
  points(ys, col = "chocolate", cex = 0.5, pch = 16)
  lines(ys, col = "chocolate", lwd = 1.2)

  # remove (offset -1)    
  # predict(tys[1:4], ret[[3]][3,])
}





plot_isw <- function(ys, ds, deg, noise = TRUE){

    res <- iwt_err(ys, ds, deg)
    ran <- (k + deg^2):(length(yout)- adjust)
    plot(ys[ran + k - 1], col = "pink4", pch = 16, cex = 0.4)
    lines(ys[ran + k -1], col = "chocolate", lwd = 1)
    lines(res$yout[ran], col = "gray10", lwd = 1.5)
}


test_interp <- function(deg = 3, add_noise = TRUE){
  # deg <- 3
  xx <- seq(0,10, 1/100)
  ys1 <- sin(2*xx) + (1/3)*cos((.33*xx)) + exp(-xx) 
  ys2 <- .8*cos(12*xx) + .5*cos(9*xx) + exp(-xx)
  ys <- ys1 + ys2
  ds <- downsample_perm(length(ys), 8)
  # txs <- seq(1, length(ys), by = 2)
  tys <- ys[ds[[1]]] 
  if(add_noise) {
  tys <- tys + rnorm(length(tys), 0, 0.5)  
  }
  # tys
  filters <- filters1()
  res <- interp(tys, deg, filters1()[[deg]])
  res <- iswt_err(tys, 8, 3)  
  # ps <- double(length(tys))
  # offset <- ceiling(deg/2) 
  # for(k in 1:offset){
  #   ps[k] <- predict(tys[1:(2*offset)], ret[[deg]][k,])
  #   ps[length(ps) - offset + k] <- predict(tys[(length(ps)-(2*offset)+1):length(ps)], 
  #         ret[[deg]][deg - k,])
  # }
  # for(k in offset:(length(ps)- offset)){
  #     ps[k] <- predict(tys[(k - offset + 1):(k+offset)], ret[[deg]][offset+1,])
  # }
  # res <- interleave( tys[1:length(tys)], ps[1:length(tys)])[1:length(ys)]
  # which(is.na(res))
  plot(res$yout, cex = 0.3, col = "gray10")
  lines(res$yout, col = "gray40", lwd = 1.2)
  # points(res, col = "gray10", cex = 0.5, pch = 16)  
  points(ys, col = "chocolate", cex = 0.5, pch = 16)
  lines(ys, col = "chocolate", lwd = 1.2)

  # remove (offset -1)    
  # predict(tys[1:4], ret[[3]][3,])
}


# test_interp <- function(){
#     xx <- seq(0,10, 1/100)
#     ys <- sin(2*xx) + (1/3)*cos((.33*xx)) + exp(-xx)
#     txs <- seq(1, length(ys), by = 2)
#     tys <- ys[txs]  
#     # tys
#     ret <- filter_list()
#     ps <- double(length(txs))

#     for(k in 2:(length(ps)- 2)){
#         ps[k+1] <- predict(tys[k:(k+deg)], ret[[deg]][ceil(deg/2)+1,])
#     }
#     res <- interleave( tys[1:501], ps[1:501])
#     plot(res, cex = 0.3)
#     lines(res, type = 'l', col = "blue", cex = 0.3)
#     # predict(tys[1:4], ret[[3]][3,])
# }

test_pfilter <- function(){
    order <- 3 
    pos <- c(-3,-1,1,3)
    cmat <- matrix(c( ))
}

# test_interleave <- function(){
#     n <- 12 
#     x <- 1:n 
#     stopifnot(all(x == interleave(seq(1,11, by = 2), 
#         seq(2,12, by = 2))
# }

