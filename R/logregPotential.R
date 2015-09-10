#' Log Regression Potential
#' 
#' Initialize a log regression to reasonable (random) defaults. 
#' \deqn{y = \alpha \log(x) + \beta + \epsilon}
#' If alpha, beta and var is provided they will be used
#' @param dim     A vector with dimensions of potential
#' @param range1  Range of independent variable
#' @param range2  Range of dependent variable
#' @param alpha   Slope of linear regression
#' @param beta    Intercept of linear regression
#' @param var     Variance
#' @export
logregPotential <- function(dim = c(100,100), range1 = c(0,100), range2 = c(0,100), alpha = NULL, beta = NULL, var = NULL){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  if(min(range1) < 0)
    stop("for log regression range1 should be positive")
  
  if(is.null(alpha) & is.null(beta) & is.null(var)){
    # Draw alpha beta and var
    # Increasing or decreasing
    monotonicity <- ifelse(runif(1) > .5, 1, -1)
    
    drange1 <- diff(range(log(.midpoints(range1[1], range1[2], dim[1]))))
    alpha <- monotonicity*runif(1, 0, diff(range2)/3/drange1 )
    beta <- runif(1, range2[1]+diff(range2)/3, range2[2]-diff(range2)/3)
    var <- diff(range2)**2/100
  } else {
    # Check mean, var and alpha
    if( is.null(alpha) | is.null(beta) | is.null(var))
      stop("Provide either all of alpha, beta and var or none of them")
    if( ! var > 0 )
      stop("var must be positive")
  }
  
  # Means
  means <- log(.midpoints(range1[1], range1[2], dim[1]))*alpha+beta
  
  mat <- t(sapply(means, FUN=function(x){
    (head(dnorm(seq(range2[1], range2[2], length.out = dim[2] + 1),x,sqrt(var)),-1) +
       tail(dnorm(seq(range2[1], range2[2], length.out = dim[2] + 1),x,sqrt(var)),-1)) /
      dim[2]*diff(range2)/2
  }))
  pot <- list(mat = mat,
              param = list(dim = dim,
                range1 = range1,
                range2 = range2,
                alpha = alpha,
                beta = beta,
                var = var))
  class(pot) <- c("logregPotential", "potential")
  return(pot)
}

update.logregPotential <- function(pot, expCounts){
  range1 <- pot$param$range1
  range2 <- pot$param$range2

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
}
