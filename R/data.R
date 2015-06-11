.checkInputData <- function(dfg, data, dataList = list()){
  # Matrix or data frame
  if( !is.matrix(data) && !is.data.frame(data))
    stop("Data should either be matrix or data frame")
  
  # Number of columns matches number of variables
  if( ! ncol(data) == length(dfg$varDim) )
    stop("number of columns doesn't match number of variables")
  
  # Columns are numeric or logical
  if( is.matrix(data) ){
    if( ! (is.logical(data) || is.numeric(data)) )
      stop("Data should be either numeric or logical")
  } 
  if( is.data.frame(data)){
    if( ! all( sapply(data, is.numeric) | sapply(data, is.logical) ) )
      stop("Data should be either numeric or logical")
  }
  
  # Within range
  if( ! suppressWarnings( all( sapply(data, max, na.rm = T) <= dfg$varDim) ) )
    stop("Data outside range")
  if( ! suppressWarnings( all( sapply(data, min, na.rm = T) > 0 )))
    stop("Data outside range")
  
  # DataList is a list of lists of numeric vectors
  if( !is.list(dataList)){
    stop("dataList must be a list")
  }
  
  sapply(dataList, FUN = function(x){
    if( is.list(x) ){
      if(!all(sapply(x, is.numeric)))
        stop("dataList vectors must be numeric")
    } else if(! is.null(x)){
      stop("dataList elements must be either list or null")
    }
  })
}