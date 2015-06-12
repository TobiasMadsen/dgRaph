#' Kullback-Leibler divergence
#' 
#' Calculates Kullback-leibler divergence KL(dfg1, dfg2) 
#' between two factor grahs with identical structure.
#' @param dfg1    Discrete factor graph object describing background distribution
#' @param dfg2    Discrete factor graph object describing foreground distribution
kl <- function(dfg1, dfg2){
  if(! is.dfg(dfg1))
    stop("dfg1 must be a dfg object")
  if(! is.dfg(dfg2))
    stop("dfg2 must be a dfg object")
  
  # Check same structure
  if( length(dfg1$facPot) != length(dfg2$facPot))
    stop("The two factor graphs must have same structure: facPot not identical")
  if( ! all(sapply(seq_along(dfg1$facPot), 
                 FUN=function(i){ all(dim(dfg1$facPot[[i]]) == dim(dfg2$facPot[[i]])) }
                 )))
    stop("The two factor graphs must have same structure: facPot not identical")
  if( ! all(dfg1$varDim == dfg2$varDim))
    stop("The two factor graphs must have same structure: varDim not identical")
  
  # Scores
  scores <- lapply(seq_along(dfg1$facPot), FUN=function(i){log(dfg1$facPot[[i]]/dfg2$facPot[[i]])})
  res <- expect(dfg1, scores)[2]
  names(res) <- "kl"
  res
}