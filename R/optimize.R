#################################################
# Building custom optimization functions
#################################################

#' Beta optimization
#' @param range     range of variable
#' @export
betaOptimize <- function(range = c(0,1)){
  # Return a function that takes expcounts and return 
  return(function(expCounts){.betaOptimize(expCounts, range = range)})
}

#' Linear Regression optimization
#' Optimize a normal linear regression
#' @param range1    range of independent variable
#' @param range2    range of dependent variable
#' @export
linregOptimize <- function(range1 = c(0,100), range2 = c(0,100) ){
  return(function(expCounts){.linregOptimize(expCounts, range1 = range1, range2 = range2)})
}

#' Normal optimization
#' @param range     range of variable
#' @export
normOptimize <- function(range = c(0,100)){
  return(function(expCounts){.normOptimize(expCounts, range = range)})
}

#' Fixed link optimization
#' Known linear relationship between two variables but unknown variance
#' @param range1    range of independent variable
#' @param range2    range of dependent variable
#' @param alpha     known slope
#' @param beta      known intercept
#' @export
fixedlinkOptimize <- function(range1 = c(0,100), range2 = c(0,100), alpha = 1, beta = 0){
  return(function(expCounts){.fixedlinkOptimize(expCounts, range1 = range1, range2 = range2, alpha = alpha, beta = beta)})
}

#################################################
# Built-in optimization functions
# Take expectation counts and return a potential
#################################################

.linregOptimize <- function(expCounts, range1 = c(0,nrow(expCounts)), range2 = c(0,ncol(expCounts))){
  SP_xy <- sum(row(expCounts) * col(expCounts) * expCounts)
  S_x   <- sum(row(expCounts) * expCounts) # x runs along rows
  S_y   <- sum(col(expCounts)* expCounts) # y runs along columns
  USS_x <- sum(row(expCounts)**2 * expCounts)
  USS_y <- sum(col(expCounts)**2 * expCounts)
  N <- sum(expCounts)
  
  SSD_x  <- USS_x-S_x**2/N
  SSD_y  <- USS_y-S_y**2/N
  SPD_xy <- SP_xy - S_x*S_y/N
  
  alpha <- SPD_xy / SSD_x
  beta  <- (S_y - S_x* alpha) / N
  var   <- (SSD_y - SPD_xy * SPD_xy/SSD_x) / N
  if(var  <= 0)
    stop("Variance 0 or less")
  
  # Scale to range
  alphaScaled <- alpha / (diff(range1)/nrow(expCounts)) * (diff(range2)/ncol(expCounts))
  betaScaled <- (beta-0.5) * (diff(range2)/ncol(expCounts)) + range2[1]
  varScaled <- var * (diff(range2)/ncol(expCounts))**2
  
  # Use potential generator
  pot <- linregPotential(dim = dim(expCounts), 
                         range1 = range1,
                         range2 = range2,
                         alpha = alphaScaled,
                         beta = betaScaled,
                         var = varScaled)
  
  str <- "Linear Regression Potential\n"
  str <- paste0(str, "alpha:\t", signif(alphaScaled, 5), "\nbeta:\t", signif(betaScaled, 5), '\nvar:\t', signif(varScaled, 5), '\n')
  return(list(pot = pot, str = str, alpha = alphaScaled, beta = betaScaled, var = varScaled))    
}

.logregOptimize <- function(expCounts, range1 = c(0,nrow(expCounts)), range2 = c(0,ncol(expCounts))){
  #values <- .midpoints()
  n <- nrow(expCounts)
  m <- ncol(expCounts)
  values_x <- matrix(log(.midpoints(range1[1], range1[2], n)), n, m, byrow = FALSE)
  values_y <- matrix(.midpoints(range2[1], range2[2], m), n, m, byrow = TRUE)
  SP_xy <- sum(values_x * values_y * expCounts)
  S_x   <- sum(values_x * expCounts) # x runs along rows
  S_y   <- sum(values_y * expCounts) # y runs along columns
  USS_x <- sum(values_x**2 * expCounts)
  USS_y <- sum(values_y**2 * expCounts)
  N <- sum(expCounts)
  
  SSD_x  <- USS_x-S_x**2/N
  SSD_y  <- USS_y-S_y**2/N
  SPD_xy <- SP_xy - S_x*S_y/N
  
  alpha <- SPD_xy / SSD_x
  beta  <- (S_y - S_x* alpha) / N
  var   <- (SSD_y - SPD_xy * SPD_xy/SSD_x) / N
  if(var  <= 0)
    stop("Variance 0 or less")
  
  # Use potential generator
  pot <- logregPotential(dim = dim(expCounts), 
                         range1 = range1,
                         range2 = range2,
                         alpha = alpha,
                         beta = beta,
                         var = var)
  
  str <- "Log Regression Potential\n"
  str <- paste0(str, "alpha:\t", signif(alpha, 5), "\nbeta:\t", signif(beta, 5), '\nvar:\t', signif(var, 5), '\n')
  return(list(pot = pot, str = str, alpha = alpha, beta = beta, var = var))    
}

.fixedlinkOptimize <- function(expCounts, range1 = c(1,nrow(expCounts)), range2 = c(1,ncol(expCounts)), alpha = 1, beta = 0){
  # Find means
  m <- matrix(.midpoints(range1[1], range1[2], nrow(expCounts))*alpha+beta, nrow(expCounts), ncol(expCounts))
  v <- matrix(.midpoints(range2[1], range2[2], ncol(expCounts)), nrow(expCounts), ncol(expCounts), byrow = T)
  
  # Var
  var <- sum((m-v)**2*expCounts) / sum(expCounts)
  
  # Return
  pot <- linregPotential(dim = dim(expCounts), 
                         range1 = range1,
                         range2 = range2,
                         alpha = alpha,
                         beta = beta,
                         var = var)
  str <- "Fixed link regression potential\n"
  str <- paste0(str, "alpha:\t", signif(alpha, 5), "\nbeta:\t", signif(beta, 5), '\nvar:\t', signif(var, 5), '\n')
  return(list(pot = pot, str = str, range1 = range1, range2 = range2, alpha = alpha, beta = beta, var = var))
}

.rowOptimize <- function(expCounts){
  pot <- sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')
  str <- "Row normalized multinomial potential\n"
  str <- paste0(str, '\t', paste(1:ncol(pot), collapse = '\t'), '\n')
  for(i in 1:nrow(pot))
    str <- paste0(str, i, '\t', paste(signif(pot[i,],5), collapse = '\t'), '\n')
  
  return(list(pot = pot, str = str))
}

.independentOptimize <- function(expCounts){
  pot <- matrix(colSums(expCounts)/sum(expCounts),
                nrow = nrow(expCounts),
                ncol = ncol(expCounts),
                byrow = T)
  str <- "Row normalize independent of parent\n"
  str <- paste0(str, paste(1:ncol(pot), collapse = '\t'), '\n')
  str <- paste0(str, paste(signif(pot[1,],5), collapse = '\t'), '\n')
  
  return(list(pot = pot, str = str))
}

.normOptimize <- function(expCounts, range = c(0,ncol(expCounts))){
  str <- "norm-potential update\n"
  meanScaled <- rep(0, nrow(expCounts))
  varScaled <- rep(0, nrow(expCounts))
  
  for(i in 1:nrow(expCounts)){
    x <- expCounts[i,]
    xs <- sum(x*seq_along(x) )
    xm <- xs/sum(x)
    xv <- (sum(x*seq_along(x)**2)-xs**2/sum(x))/sum(x)
    
    meanScaled[i] <- (xm-0.5) * (diff(range)/ncol(expCounts)) + range[1]
    varScaled[i]  <- xv * (diff(range)/ncol(expCounts))**2
    str <- paste0(str, i, "\tmean:\t", signif(meanScaled[i], 5), "\tvar:\t", signif(varScaled[i], 5), '\n')
  }
  
  # Use potential generator
  pot <- normalPotential(dim = dim(expCounts),
                         range = range,
                         means = meanScaled, 
                         vars = varScaled)
  
  return(list(pot = pot, str = str, means = meanScaled, vars = varScaled))
}

.betaOptimize <- function(expCounts, range = c(0,1)){
  # Method of moments
  str <- "beta-potential update\n"
  val <- ((1:ncol(expCounts)-0.5)/ncol(expCounts))*diff(range)+range[1]
  breaks <- seq(range[1], range[2], length.out = ncol(expCounts)+1 )
  alphas <- rep(0, nrow(expCounts))
  betas <- rep(0, nrow(expCounts))
  
  for(i in 1:nrow(expCounts)){
    x   <- expCounts[i,]
    xs  <- sum(x*val )
    xm  <- xs/sum(x)
    xv  <- (sum(x*val**2)-xs**2/sum(x) )/sum(x)
    a   <- xm*(xm*(1-xm)/xv-1)
    b   <- a*(1/xm-1)
    
    alphas[[i]] <- a
    betas[[i]]  <- b
    str <- paste0(str, i, "\talpha:\t", signif(a, 5), "\tbeta:\t", signif(b, 5), '\n')
  }
  
  pot <- betaPotential(dim = dim(expCounts),
                       range = range,
                       alphas = alphas,
                       betas = betas)
  
  return(list(pot = pot, str = str, alphas = alphas, betas = betas))
}