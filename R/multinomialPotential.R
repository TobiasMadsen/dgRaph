#' Multinomial Potential
#' 
#' Initialize a multinomial potential with random starting point
#' @param dim     A vector with dimensions of potential
#' @param mat     Optionally give fixed matrix potential
#' @param independent The two nodes are independent, i.e. rows are identical. This is mainly useful for encoding additional indepedence assumptions in a fixed graph structure.
#' @param pseudocount Add pseudocount to each cell. Default is 0.
#' @export
multinomialPotential <- function(dim = c(1,5), mat = NULL, independent = FALSE, pseudocount = 0L){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  if(is.null(mat)){
    if( ! independent ){
      mat <- matrix(runif(prod(dim)), dim[1], dim[2])
      mat <- sweep(mat, 1, STATS = rowSums(mat), FUN = '/')
    }
    if( independent ){
      r <- runif(dim[2])
      mat <- matrix( r/sum(r), dim[1], dim[2], byrow = T)
    }
  }
  pot <- list(mat = mat,
              param = list(independent = independent, pc = pseudocount))
  class(pot) <- c("multinomialPotential", "potential")
  pot
}

# Optimization
update.multinomialPotential <- function(pot, expCounts){
  if(pot$param$pc != 0)
    expCounts <- expCounts + pot$param$pc
    
  if(! pot$param$independent)
    return( multinomialPotential(dim = dim(expCounts),
                                 mat = sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')) )

  # Else: Doesn't depend on parent i.e. all rows are identical
  return( multinomialPotential(dim = dim(expCounts),
                               mat = matrix(colSums(expCounts)/sum(expCounts),
                                 nrow(expCounts),
                                 ncol(expCounts),
                                 byrow = TRUE),
                               independent = TRUE) )
}
