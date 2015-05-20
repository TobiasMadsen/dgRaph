#' Potentials
#' @param dfg     discrete factor graph object
#' @return A list of current factors
potentials <- function(dfg){
    stopifnot(is.dfg(dfg))

    dfg$dfgmodule$getFactorPotentials()
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
#' @export
betaPotential <- function(dim = c(1, 100)){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  # Draw some Means
  m <- runif(dim[1])*0.8+0.1 # range [0.1, 0.9]
  t(sapply(m, FUN=function(x){
    a <- x*(x*(1-x)/(0.02)-1)
    b <- a*(1/x-1)
    diff(pbeta( 0:dim[2]/dim[2], x, a, b ) )
  }))
}