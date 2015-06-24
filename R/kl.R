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
  .compareDfgs(dfg1, dfg2)

  # If not same nb structure remap
  if( ! all( sapply(seq_along(dfg1$facNbs), function(i){all(dfg1$facNbs[[i]] == dfg2$facNbs[[i]])})))
    dfg2 <- .remapFacNbsDfg(dfg2, match(dfg2$facNbs, dfg1$facNbs))
  
  # If not same potential map remap
  if( ! length(dfg1$potMap) == length(dfg2$potMap) |
      ! all(dfg1$potMap == dfg2$potMap) ){
    cm <- .commonMap(dfg1$potMap, dfg2$potMap)
    dfg1 <- .remapPotMapDfg(dfg1, cm)
    dfg2 <- .remapPotMapDfg(dfg2, cm)
  }
  
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
