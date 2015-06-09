simulate <- function(x, ...){
  UseMethod("simulate")
}

simulate.default <- function(x, ...){
  stop("simulate not implemented for class: ", class(x))
}

#' Simulate from DFG
#' @param dfg     Discrete factor graph object
#' @param n       Number of samples
#' @export
simulate.dfg <- function(dfg, n){
  # Input check
  if(!is.numeric(n))
    stop("n should be numeric")
  if(n < 0)
    stop("n < 0")
  
  # Sample
  ret <- data.frame(dfg$dfgmodule$simulate(n))
  colnames(ret) <- dfg$varNames
  return(ret)
}