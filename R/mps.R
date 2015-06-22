#' Most Probable State
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param dataList  provides support for partially observed data. See online documentation.
#' @return A dataframe with a column for each variable and the most probable state for each observation
mps <- function(data, dfg, dataList = list()){
  .mps(data, dfg, dataList)
}

.mps <- function(data, dfg, dataList = list(), module = NULL){
  .checkInputData(dfg, data, dataList)
  
  if(is.null(module))
    module <- .build(dfg)
  
  # Calculate mps
  ret <- data.frame(module$mps(as.matrix(data), dataList))
  colnames(ret) <- dfg$varNames
  return(ret)
}
