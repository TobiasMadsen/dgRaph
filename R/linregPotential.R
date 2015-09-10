#' Linear Regression Potential
#' 
#' Initialize a linear regression to reasonable (random) defaults. If alpha, beta and var is provided they will be used
#' @param dim     A vector with dimensions of potential
#' @param range1  Range of independent variable
#' @param range2  Range of dependent variable
#' @param alpha   Slope of linear regression
#' @param beta    Intercept of linear regression
#' @param var     Variance
#' @export
linregPotential <- function(dim = c(100,100), range1 = c(0,100), range2 = c(0,100), alpha = NULL, beta = NULL, var = NULL){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  
  if(is.null(alpha) & is.null(beta) & is.null(var)){
    # Draw alpha beta and var
    # Increasing or decreasing
    monotonicity <- ifelse(runif(1) > .5, 1, -1)
    alpha <- monotonicity*runif(1, 0, diff(range2)/3/diff(range1) )
    beta <- runif(1, range2[1]+diff(range2)/3, range2[2]-diff(range2)/3)
    var <- diff(range2)**2/25
  } else {
    # Check mean, var and alpha
    if( is.null(alpha) | is.null(beta) | is.null(var))
      stop("Provide either all of alpha, beta and var or none of them")
    if( ! var > 0 )
      stop("var must be positive")
  }
  
  # Means
  means <- .midpoints(range1[1], range1[2], dim[1])*alpha+beta
  
  mat <- t(sapply(means, FUN=function(x){
    (head(dnorm(seq(range2[1], range2[2], length.out = dim[2] + 1),x,sqrt(var)),-1) +
      tail(dnorm(seq(range2[1], range2[2], length.out = dim[2] + 1),x,sqrt(var)),-1)) /
      dim[2]*diff(range2)/2
  }))
  pot <- list(mat = mat,
              param = list(alpha = alpha, beta = beta, var = var, range1 = range1, range2 = range2, dim = dim))
  class(pot) <- c("linregPotential", "potential")
  return(pot)
}

update.linregPotential <- function(pot, expCounts){
  if(! all(dim(expCounts) == dim(pot)))
    stop("dimensions of potential and expectation counts does not match")

  range1 <- pot$param$range1
  range2 <- pot$param$range2
  
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
}
