#################################################
# Building custom optimization functions
#################################################

#' Beta optimization
#' @export
betaOptimize <- function(range = c(0,1)){
  # Return a function that takes expcounts and return 
  return(function(expCounts){.betaOptimize(expCounts, range = range)})
}

#' Linear Regression optimization
#' @export
linregOptimize <- function(range1 = 1:100, range2 = 1:100){
  return(function(expCounts){.linregOptimize(expCounts, range1 = range1, range2 = range2)})
}

#################################################
# Built-in optimization functions
# Take expectation counts and return a potential
#################################################

.linregOptimize <- function(expCounts, range1 = 1:nrow(expCounts), range2 = 1:ncol(expCounts)){
  # Ranges are only used for output
  SP_xy <- sum(1:nrow(expCounts) %*% t(1:ncol(expCounts)) * expCounts)
  S_x   <- sum(1:nrow(expCounts) %*% t(rep(1,ncol(expCounts))) * expCounts)
  S_y   <- sum(rep(1,nrow(expCounts)) %*% t(1:ncol(expCounts)) * expCounts)
  USS_x <- sum( (1:nrow(expCounts))**2 %*% t(rep(1,ncol(expCounts))) * expCounts)
  USS_y <- sum(rep(1,nrow(expCounts)) %*% t((1:ncol(expCounts))**2) * expCounts)
  N <- sum(expCounts)
  
  SSD_x  <- USS_x-S_x**2/N
  SSD_y  <- USS_y-S_y**2/N
  SPD_xy <- SP_xy - S_x*S_y/N
  
  alpha <- SPD_xy / SSD_x
  beta  <- (S_y - S_x* alpha) / N
  var   <- (SSD_y - SPD_xy * SPD_xy/SSD_x) / N
  if(var  <= 0)
    stop("Variance 0 or less")
  
  pot <- t(sapply(1:nrow(expCounts), FUN=function(x){
    diff( pnorm(c(-Inf, 2:ncol(expCounts), Inf) ,x*alpha+beta+0.5, sqrt(var)) )
  }))
  
  str <- "Linear Regression Potential\n"
  str <- paste0(str, "alpha:\t", signif(alpha, 5), "\nbeta:\t", signif(beta, 5), '\nvar:\t', signif(var, 5), '\n')
  return(list(pot = pot, str = str))    
}

.rowOptimize <- function(expCounts){
  pot <- sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')
  str <- "Row normalized multinomial potential\n"
  str <- paste0(str, '\t', paste(1:ncol(pot), collapse = '\t'), '\n')
  for(i in 1:nrow(pot))
    str <- paste0(str, i, '\t', paste(signif(pot[i,],5), collapse = '\t'), '\n')
  
  return(list(pot = pot, str = str))
}

.normOptimize <- function(expCounts){
  str <- "norm-potential update\n"
  i <- 0
  pot <- t(apply(expCounts, 1, FUN=function(x){
    xs <- sum(x*seq_along(x) )
    xm <- xs/sum(x)+1/2
    xv <- (sum(x*seq_along(x)**2)-xs**2/sum(x))/sum(x)
    i <<- i + 1
    str <<- paste0(str, i, "\tmean:\t", signif(xm, 5), "\tvar:\t", signif(xv, 5), '\n')
    diff(pnorm( c(-Inf, 2:length(x), Inf), xm, sqrt(xv)))
  }))
  return(list(pot = pot, str = str))
}

.betaOptimize <- function(expCounts, range = c(0,1)){
  # Method of moments
  str <- "beta-potential update\n"
  i <- 0
  pot <- t(apply(expCounts, 1, FUN=function(x){
    val <- (seq_along(x)-0.5)/length(x)
    xs  <- sum(x*val )
    xm  <- xs/sum(x)
    xv  <- (sum(x*val**2)-xs**2/sum(x) )/sum(x)
    a   <- xm*(xm*(1-xm)/xv-1)
    b   <- a*(1/xm-1)
    
    i <<- i + 1
    str <<- paste0(str, i, "\talpha:\t", signif(a, 5), "\tbeta:\t", signif(b, 5), '\n')
    diff(pbeta( 0:length(x)/length(x), a, b))
  }))
  return(list(pot = pot, str = str))
}