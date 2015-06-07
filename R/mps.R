#' Most Probable State
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @return A dataframe with a column for each variable and the most probable state for each observation
mps <- function(data, dfg, dataList = list()){
  .checkInputData(dfg, data, dataList)

  # Calculate mps
  ret <- data.frame(dfg$dfgmodule$mps(as.matrix(data), dataList))
  colnames(ret) <- dfg$varNames
  return(ret)
}
