#' Calculate likelihood of observation
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param log       calculate loglikelihood
#' @return A vector of likelihoods for each observation
#' @export
likelihood <- function(data, dfg, log = FALSE, dataList = list()){
  # Correct number of columns
  stopifnot(ncol(data) == length(dfg$varNames))
  
  # Correct data type
  #if(is.matrix(data))
  #    stopifnot(is.numeric(data))
  # TODO NA's has type logical
  #if(is.data.frame(data))
  #    stopifnot(all(lapply(data, is.numeric)))
  
  # Calculate likelihood
  if( log ){
    dfg$dfgmodule$calcLogLikelihood(as.matrix(data), dataList)
  } else{
    dfg$dfgmodule$calcLikelihood(as.matrix(data), dataList)
  }
}
