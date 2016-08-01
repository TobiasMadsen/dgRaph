#' Expectation
#' Generic function for finding expectancy of function under a distribution x
#' @param x   distribution under which the expecation should be evaluated
#' @param ... Arguments to be passed to methods.
expect <- function(x, ...){
  UseMethod("expect")
}

#' Expectations in DFG
#' @param x   dfg object 
#' @param facScores   A list of matrices similar to factor potentials indicating a function which expectation should be evaluated
#' @param data dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param ... Arguments to be passed to methods.
expect.dfg <- function(x, facScores, data = NULL, ...){
  .expect.dfg(x, facScores, data = data)
}

.expect.dfg <- function(dfg, facScores, data = NULL, module = NULL){
  # Check dimensions
  stopifnot( is.list(facScores), all(sapply(facScores, is.matrix)))
  stopifnot( length(facScores) == length(dfg$facPot))
  
  facPotBg <- potentials(dfg)
  stopifnot( all(sapply(seq_along(facPotBg), FUN=function(i){all(dim(facPotBg[[i]])==dim(facScores[[i]]))})) )
  
  if(is.null(module))
    module <- .build(dfg)

  if(is.null(data)){
    ret <- module$expect(facScores)
    names(ret) <- c("Likelihood", "Expect")
  } else {
    .checkInputData(dfg, data)
    ret <- module$expectCondData(facScores, as.matrix(data))
    colnames(ret) <- c("Likelihood", "Expect")
  }
  
  ret  
}

#' Expectation 2nd order moment
#' Generic function for finding likelihood, expectancy and 2nd order moment of function under a distribution x
#' @param x   distribution under which the expecation should be evaluated
#' @param ... Arguments to be passed to methods.
expect2 <- function(x, ...){
  UseMethod("expect2")
}

#' Expectations and 2nd order moments in DFG
#' @param x   dfg object 
#' @param facScores   A list of matrices similar to factor potentials indicating a function which expectation should be evaluated
#' @param data        dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param ... Arguments to be passed to methods.
expect2.dfg <- function(x, facScores, data = NULL,...){
  .expect2.dfg(x, facScores, data = data)
}

.expect2.dfg <- function(dfg, facScores, data = NULL, module = NULL){
  # Check dimensions
  stopifnot( is.list(facScores), all(sapply(facScores, is.matrix)))
  stopifnot( length(facScores) == length(dfg$facPot))
  
  facPotBg <- potentials(dfg)
  stopifnot( all(sapply(seq_along(facPotBg), FUN=function(i){all(dim(facPotBg[[i]])==dim(facScores[[i]]))})) )
  
  if(is.null(module))
    module <- .build(dfg)
  
  ret <- module$gamma(facScores)
  names(ret) <- c("Likelihood", "Expect", "Gamma")
  ret  
}
