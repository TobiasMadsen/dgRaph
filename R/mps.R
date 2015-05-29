#' Most Probable State
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @return A dataframe with a column for each variable and the most probable state for each observation
mps <- function(data, dfg){
  # Correct number of columns
  stopifnot(ncol(data) == length(dfg$varNames))
  
  # Correct data type
  if(is.matrix(data)){
      stopifnot(is.numeric(data))
  } else if(is.data.frame(data)){
      # TODO NA's has type logical
  } else {
      stop("data must be either a data.frame or a matrix")
  }

  # Check range of variables

  
  # Calculate mps
  res <- apply(data, 1, FUN=function(x){
    obs <- x
    obs[is.na(obs)] <- 1
    dfg$dfgmodule$maxProbState(obs ,!is.na(x))
  })
  
  res <- as.data.frame(t(res))
  names(res) <- dfg$varNames
  return(res)
}
