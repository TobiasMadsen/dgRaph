#' Mean link potential
#' 
#' Initialize a meanlink potential. The first neighbour encodes the conditional mean in the second neighbour. The variance can still be learned from data. Internally a linregPotential where alpha and beta are fixed trough when training.
#' @param dim     A vector with dimensions of potential
#' @param range1  Range of independent variable
#' @param range2  Range of dependent variable
#' @param alpha   Slope of linear regression
#' @param beta    Intercept of linear regression
#' @param var     Variance
#' @export
meanlinkPotential <- function(dim = c(100,100), range1 = c(0,100), range2 = c(0,100), alpha = 1, beta = 0, var = NULL){
  if(is.null(var))
    var <- diff(range2)**2/25
  pot <- linregPotential(dim = dim,
                         range1 = range1,
                         range2 = range2,
                         alpha = alpha,
                         beta = beta,
                         var = var)
  class(pot) <- c("meanlinkPotential","linregPotential","potential")
  return(pot)
}

update.meanlinkPotential <- function(pot, expCounts){
  range1 <- pot$param$range1
  range2 <- pot$param$range2
  alpha <- pot$param$alpha
  beta <- pot$param$beta
  var <- pot$param$var
  
  # Find means
  m <- matrix(.midpoints(range1[1], range1[2], nrow(expCounts))*alpha+beta, nrow(expCounts), ncol(expCounts))
  v <- matrix(.midpoints(range2[1], range2[2], ncol(expCounts)), nrow(expCounts), ncol(expCounts), byrow = T)
  
  # Var
  var <- sum((m-v)**2*expCounts) / sum(expCounts)
  
  # Return
  pot <- meanlinkPotential(dim = dim(expCounts), 
                           range1 = range1,
                           range2 = range2,
                           alpha = alpha,
                           beta = beta,
                           var = var)
}
