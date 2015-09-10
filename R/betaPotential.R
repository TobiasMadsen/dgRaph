#' Beta Potential
#' 
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
  
  pot <- list(
    mat = t(sapply(seq_along(alphas), FUN=function(i){
      dbeta(.midpoints(range[1], range[2], length.out = dim[2]), alphas[i], betas[i]) / dim[2]*diff(range)})),
    param = list(alphas = alphas, betas = betas, range = range, dim = dim)    
  )
  class(pot) <- c("betaPotential", "potential") 
  return(pot)
}

#' Beta optimization
#' @param range     range of variable
#' @export
betaOptimize <- function(range = c(0,1)){
  # Return a function that takes expcounts and return 
  return(function(expCounts){.betaOptimize(expCounts, range = range)})
}

update.betaPotential <- function(pot, expCounts){
  if(! all(dim(expCounts) == dim(pot)))
    stop("dimensions of potential and expectation counts does not match")
  
  range <- pot$param$range
  
  # Method of moments
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
  }
  
  pot <- betaPotential(dim = dim(expCounts),
                       range = range,
                       alphas = alphas,
                       betas = betas)
  pot
}
