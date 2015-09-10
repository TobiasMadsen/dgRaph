#' Multinomial Potential
#' 
#' Initialize a multinomial potential with random starting point
#' @param dim     A vector with dimensions of potential
#' @export
multinomialPotential <- function(dim = c(1,5)){
  if(length(dim) != 2)
    stop("dim should be a vector of length 2")
  m <- matrix(runif(prod(dim)), dim[1], dim[2])
  sweep(m, 1, STATS = rowSums(m), FUN = '/')
}
