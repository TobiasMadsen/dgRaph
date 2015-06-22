#' Calculate expection counts
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param dataList  provides support for partially observed data. See online documentation.
#' @return A list of matrices containing the expectation counts for each factor
#' @export
facExpectations <- function(data, dfg, dataList = list()){
  .facExpectations(data, dfg, dataList)
}

.facExpectations <- function(data, dfg, dataList = list(), module = NULL){
  .checkInputData(dfg, data, dataList)
  
  if(is.null(module))
    module <- .build(dfg)
  
  # Calculate expectation counts
  res <- module$facExpCounts(as.matrix(data), dataList)
  return(res)
}
