#' Calculate likelihood of observation
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param log       calculate loglikelihood
#' @return A vector of likelihoods for each observation
#' @export
likelihood <- function(data, dfg, log = FALSE, dataList = list()){
  .likelihood(data, dfg, log, dataList)
}

# Hidden option: provide module to avoid rebuilding module for repeated evaluation
.likelihood <- function(data, dfg, log = FALSE, dataList = list(), module = NULL){
  .checkInputData(dfg, data, dataList)
  
  if(is.null(module))
    module <- .build(dfg)
  
  # Calculate likelihood
  if( log ){
    module$calcLogLikelihood(as.matrix(data), dataList)
  } else{
    module$calcLikelihood(as.matrix(data), dataList)
  }
}


