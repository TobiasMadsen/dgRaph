#' Potentials
#' @param dfg     discrete factor graph object
#' @return A list of current factors
#' @export
potentials <- function(dfg){
    stopifnot(is.dfg(dfg))

    dfg$facPot
}

#' Potentials <- 
#' @param dfg     discrete factor graph object
#' @export
"potentials<-" <- function(dfg, value){
  # Perform input check
  if( ! is.list(value) | ! all(sapply(value, is.matrix)))
    stop("Potentials must be a list of matrices")
  if(length(value) != length(dfg$facPot))
    stop("Potentials did not have correct length")
  
  sapply(seq_along(dfg$facPot), FUN=function(i){
      if( !all( dim(dfg$facPot[[i]]) == dim(value[[i]]) )){
        if( !all( dim(value[[i]]) == c(0,0) )){
          stop(i,"th potential doesn't match")
        }
      }
    })
  
  for(i in seq_along(dfg$facPot)){
    if( all(dim(value[[i]]) == c(0,0))){
      value[[i]] <- dfg$facPot[[i]]
    }
  }
  
  # Set potentials
  dfg$facPot <- value
  dfg
}

#' Linear Regression Potential
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
  
  t(sapply(means, FUN=function(x){
    (head(dnorm(seq(range2[1], range2[2], length.out = dim[2] + 1),x,sqrt(var)),-1) +
      tail(dnorm(seq(range2[1], range2[2], length.out = dim[2] + 1),x,sqrt(var)),-1)) /
      dim[2]*diff(range2)/2
  }))
}

#' Normal Potential
#' Initialize a norm potential if means and vars are not provided they will be initialized at random such that the whole range is covered.
#' @param dim     A vector with dimensions of potential
#' @param means   A vector of means. Provide a mean for each class.
#' @param vars    A vector of vars. Provide a variance for each class.
#' @export
normalPotential <- function(dim = c(1,100), range = c(0,100), means = NULL, vars = NULL){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  if(length(range) != 2)
    stop("range should be a vector of length 2")
  
  if(is.null(means) & is.null(vars)){
    # Draw means
    means <- runif(dim[1], range[1]+0.3*diff(range), range[2]-0.3*diff(range))
    vars <- rep( diff(range)**2/25, dim[1])
  } else {
    # Check means and vars
    if( ! dim[1] == length(means))
      stop("Provide as many means as dim[1]")
    if( ! dim[1] == length(vars))
      stop("Provide as many vars as dim[1]")
    if( ! all(vars > 0))
      stop("vars must be positive")
  }
  
  t(sapply(seq_along(means), FUN=function(i){
    (head(dnorm(seq(range[1], range[2], length.out = dim[2] + 1),means[i],sqrt(vars[i])),-1) +
       tail(dnorm(seq(range[1], range[2], length.out = dim[2] + 1),means[i],sqrt(vars[i])),-1)) /
    dim[2]*diff(range)/2
  }))
}

#' Multinomial Potential
#' Initialize a multinomial potential with random starting point
#' @param dim     A vector with dimensions of potential
#' @export
multinomialPotential <- function(dim = c(1,5)){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  m <- matrix(runif(prod(dim)), dim[1], dim[2])
  sweep(m, 1, STATS = rowSums(m), FUN = '/')
}

#' Beta Potential
#' Initialize a beta distributed potential
#' @param dim     A vector with dimensions of potential
#' @param range   The limits of the binning scheme.
#' @param alphas  Provide an alpha for each class. If alphas and betas are not provided they will be selected randomly and such that the whole range is covered.
#' @param betas   Provide a beta for each class. See alphas
#' @export
betaPotential <- function(dim = c(1, 100), range = c(0,1), alphas = NULL, betas = NULL){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  if(length(range) != 2)
    stop("range should be a vector of length 2")
  
  # Draw means
  if(is.null(alphas) & is.null(betas)){
    m <- runif(dim[1], min = range[1]+0.2*diff(range), max = range[1]+0.8*diff(range))
    v <- diff(range)**2 * 0.02
    alphas <- m*(m*(1-m)/v-1)
    betas <- alphas*(1/m-1)
  } else {
    # Check alphas and betas
    if( ! dim[1] == length(alphas))
      stop("Provide as many alphas as dim[1]")
    if( ! dim[1] == length(betas))
      stop("Provide as many betas as dim[1]")
    if( ! all(betas > 0))
      stop("Betas must be positive")
    if( ! all(alphas > 0))
      stop("Alphas must be positive")
  }
  
  t(sapply(seq_along(alphas), FUN=function(i){
    dbeta(.midpoints(range[1], range[2], length.out = dim[2]), alphas[i], betas[i]) / dim[2]*diff(range)
  }))
}
