#' Calculate likelihood of observation
#' This method calculates the likelihood of an observation using the sum-product algorithm. 
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param log       calculate loglikelihood
#' @param dataList  provides support for partially observed data. See online documentation.
#' @return A vector of likelihoods for each observation
#' @examples
#' varDim <- c(2,2)
#' facPot <- list(matrix(c(0.7,0.3),1,2),
#'                matrix(c(0.8,0.4,0.2,0.6),2,2))
#' facNbs <- list(c(1),c(1,2))
#' mydfg <- dfg(varDim, facPot, facNbs)
#' 
#' dat <- matrix(c(1,2),1,2)
#' likelihood( dat, dfg = mydfg)
#' dat <- matrix(c(NA,2),1,2)
#' likelihood( dat, dfg = mydfg)
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


