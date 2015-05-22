#' Potentials
#' @param dfg     discrete factor graph object
#' @return A list of current factors
potentials <- function(dfg){
    stopifnot(is.dfg(dfg))

    dfg$dfgmodule$getPotentials()
}

#' Linear Regression Potential
#' Initialize a linear regression to reasonable defaults
#' @param dim     A vector with dimensions of potential
#' @export
linregPotential <- function(dim = c(100,100)){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  
  # Increasing or decreasing
  monotonicity <- ifelse(runif(1) > .5, 1, -1)
  
  # Standard deviation
  sd <- dim[2]/8
  m  <- dim[2]/3*monotonicity*(1:dim[1]-.5)/dim[1] +
    ifelse(monotonicity == 1, dim[2]/3, dim[2]*2/3)
  
  t(sapply(m, FUN=function(x){
    diff(pnorm( c(-Inf, 1:(dim[2]-1)+0.5, Inf), x, sd) )
  }))
}

#' Normal Potential
#' Initialize a norm potential with random starting points
#' @param dim     A vector with dimensions of potential
#' @export
normalPotential <- function(dim = c(1,100)){
  if(length(dim) != 2)
      stop("dim should be a vector of length 2")
  
  # Draw some means
  # TODO: Avoid means close to boundary
  m <- (sample(dim[2], dim[1])/2)+dim[2]/4
  t(sapply(m, FUN=function(x){
    diff(pnorm( c(-Inf, 1:(dim[2]-1)+0.5, Inf), x, dim[2]/5) )
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
#' Initialize a beta distributed potential with random starting point
#' @param dim     A vector with dimensions of potential
#' @param range   The limits of the binning scheme.
#' @export
betaPotential <- function(dim = c(1, 100), range = c(0,1)){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  if(length(range) != 2)
    stop("range should be a vector of length 2")
  
  # Draw means
  m <- runif(dim[1], min = range[1]+0.2*diff(range), max = range[1]+0.8*diff(range))
  v <- diff(range)**2 * 0.02
  t(sapply(m, FUN=function(x){
    a <- x*(x*(1-x)/v-1)
    b <- a*(1/x-1)
    diff(pbeta( range[1]+diff(range)*(0:dim[2]/dim[2]), a, b ) )
  }))
}

#' Fixed Normal Potential
#' @param dim     A vector with dimensions of potential
#' @param range1  The limits of binning scheme 1st variable
#' @param range2  The limits of binning scheme 2nd variable
#' @param a       Slope
#' @param b       Intercept
#' @param sd      standard deviation
#' @export
fixedNormalPotential <- function(dim = c(100, 100), range1 = c(0, 100), range2 = c(0, 100), a = 1, b = 0, sd = diff(range2)/5){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  if(length(range1) != 2)
    stop("range1 should be a vector of length 2")
  if(length(range2) != 2)
    stop("range1 should be a vector of length 2")
  
  # Means
  a_ <- a * dim[2] / diff(range2) / dim[1] * diff( range1 ) # slope in "bin" space
  b_ <- (b-range1[1])/(range1[2]-range1[1])*dim[1] # intercept in "bin" space
  m <- b_ + 1:dim[1] * a_
  sd_ <- sd * dim[2] / diff(range2) # standard deviation in "bin" space
  
  t(sapply(m, FUN=function(x){
    diff(pnorm( c(-Inf, 1:(dim[2]-1)+0.5, Inf), x, sd_) )
  }))
}