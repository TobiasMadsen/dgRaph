#' Normal Potential
#' 
#' Initialize a norm potential if means and vars are not provided they will be initialized at random such that the whole range is covered.
#' @param dim     A vector with dimensions of potential
#' @param range   A vector with two entries providing the range of the variable.
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
  
  mat <- t(sapply(seq_along(means), FUN=function(i){
    (head(dnorm(seq(range[1], range[2], length.out = dim[2] + 1),means[i],sqrt(vars[i])),-1) +
       tail(dnorm(seq(range[1], range[2], length.out = dim[2] + 1),means[i],sqrt(vars[i])),-1)) /
      dim[2]*diff(range)/2
  }))
  pot <- list(mat = mat, 
              param = list(means = means, vars = vars, range = range, dim = dim))
  class(pot) <- c("normalPotential", "potential") 
  return(pot)
}


# Optimization
update.normalPotential <- function(pot, expCounts){
  if(! all(dim(expCounts) == dim(pot)))
    stop("dimensions of potential and expectation counts does not match")
  
  range <- pot$param$range
  
  meanScaled <- rep(0, nrow(expCounts))
  varScaled <- rep(0, nrow(expCounts))
  
  for(i in 1:nrow(expCounts)){
    x <- expCounts[i,]
    xs <- sum(x*seq_along(x) )
    xm <- xs/sum(x)
    xv <- (sum(x*seq_along(x)**2)-xs**2/sum(x))/sum(x)
    
    meanScaled[i] <- (xm-0.5) * (diff(range)/ncol(expCounts)) + range[1]
    varScaled[i]  <- xv * (diff(range)/ncol(expCounts))**2
  }
  
  # Use potential generator
  pot <- normalPotential(dim = dim(expCounts),
                         range = range,
                         means = meanScaled, 
                         vars = varScaled)
}
