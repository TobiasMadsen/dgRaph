expect <- function(x, ...){
  UseMethod("expect")
}

expect.default <- function(x){
  stop("expect not implemented for class: ", class(x))
}

#' Expectations in DFG
#' 
expect.dfg <- function(dfg, facScores){
  .expect.dfg(dfg, facScores)
}

.expect.dfg <- function(dfg, facScores, module = NULL){
  # Check dimensions
  stopifnot( is.list(facScores), all(sapply(facScores, is.matrix)))
  stopifnot( length(facScores) == length(dfg$facPot))
  
  facPotBg <- potentials(dfg)
  stopifnot( all(sapply(seq_along(facPotBg), FUN=function(i){all(dim(facPotBg[[i]])==dim(facScores[[i]]))})) )
  
  if(is.null(module))
    module <- .build(dfg)
  
  ret <- module$expect(facScores)
  names(ret) <- c("Likelihood", "Expect")
  ret  
}
