#' Calculate expection counts
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @return A list of matrices containing the expectation counts for each factor

facExpectations <- function(data, dfg, dataList = list()){
    .checkInputData(dfg, data, dataList)

    # Calculate expectation counts
    res <- dfg$dfgmodule$facExpCounts(as.matrix(data), dataList)
    #names(res) <- dfg$facNames
    return(res)
}
