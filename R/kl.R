#' Kullback-Leibler divergence
#' 
#' Calculates Kullback-leibler divergence KL(dfg1, dfg2) 
#' between two factor grahs with identical structure. Beware that if for some potential f, and some
#' variable configuration, x. f(x) > 0 in dfg1 and f(x) = 0 in dfg2, the KL divergence is infinite, 
#' here however, we ignore such entries and put f(x) = 0 in dfg1.
#' @param dfg1    Discrete factor graph object describing background distribution
#' @param dfg2    Discrete factor graph object describing foreground distribution
kl <- function(dfg1, dfg2){
  .kl(dfg1, dfg2)
}

.kl <- function(dfg1, dfg2, module = NULL){
  if(! is.dfg(dfg1))
    stop("dfg1 must be a dfg object")
  if(! is.dfg(dfg2))
    stop("dfg2 must be a dfg object")

  if(is.null(module))
    module <- .build(dfg1)
  
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
  scores <- lapply(seq_along(dfg1$facPot), FUN=function(i){
    m <- log(dfg1$facPot[[i]]/dfg2$facPot[[i]])
    
    m[dfg1$facPot[[i]] == 0 | dfg2$facPot[[i]] == 0] <- 0
    m
  })
  res <- .expect.dfg(dfg1, scores, module = module)[2]
  names(res) <- "kl"
  res
}
