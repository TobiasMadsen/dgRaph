#' Simulate
#' Simulate from a distribution
#' @param x   Distribution from which to simulate.
#' @param ... Arguments to be passed to methods. Typically \code{n} number of simulations will be provided.
simulate <- function(x, ...){
  UseMethod("simulate")
}

#' Simulate from DFG
#' @param x       Discrete factor graph object
#' @param n       Number of samples
#' @param ...     Arguments to be passed to methods.
#' @export
simulate.dfg <- function(x, n, ...){
  .simulate.dfg(x, n)
}

.simulate.dfg <- function(dfg, n, module = NULL){
  # Input check
  if(!is.numeric(n))
    stop("n should be numeric")
  if(n < 0)
    stop("n < 0")

  if(is.null(module))
    module <- .build(dfg)
  
  # Sample
  ret <- data.frame(module$simulate(n))
  colnames(ret) <- dfg$varNames
  return(ret)
}
