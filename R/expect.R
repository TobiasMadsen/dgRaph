expect <- function(x, ...){
  UseMethod("expect")
}

expect.default <- function(x){
  stop("expect not implemented for class: ", class(x))
}

#' Expections in DFG
#' 
expect.dfg <- function(x, facScores){
  # Check dimensions
  stopifnot( is.list(facScores), all(sapply(facScores, is.matrix)))
  stopifnot( length(facScores) == length(x$facPot))
  facPotBg <- potentials(x)
  stopifnot( all(sapply(seq_along(facPotBg), FUN=function(i){all(dim(facPotBg[[i]])==dim(facScores[[i]]))})) )
  ret <- x$dfgmodule$expect(facScores)
  names(ret) <- c("Likelihood", "Expect")
  ret
}

