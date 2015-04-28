#' Most Probable State
#' @param x         points to evaluate tail probabilities in
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @return A dataframe with a column for each variable and the most probable state for each observation
mps <- function(data, dfg){
  # Do some input checking
  
  # Calculate mps
  res <- apply(data, 1, FUN=function(x){
    obs <- x
    obs[is.na(obs)] <- 0
    dfg$dfgmodule$maxProbState(obs ,!is.na(x))
  })
  
  res <- as.data.frame(t(res))
  names(res) <- dfg$varNames
  res
}