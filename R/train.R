#' EM training
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. 
#'                  NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param optim     character vector with an optimization function for each potential. The built-in functions are
#'                  "row" optimizing a multinomial conditional distribution 
#'                  "norm" optimizing a discretized normal conditional distribution
#'                  "linreg" optizing a normal linear regression of x2 on x1
#' @param optimFun A named list with additional optimization functions. 
#'                  Refer to the optimization function by entry name in "optim".
#' @return A dataframe with a column for each variable and the most probable state for each observation
#' @export
train <- function(data, dfg, optim = NULL, optimFun = NULL, threshold = 1e-9, iter.max = 200){
  # Output
  sprintf("Training...\n")
    
  # Correct number of columns
  stopifnot(ncol(data) == length(dfg$varNames))
  
  # Correct data type
  if(is.matrix(data))
      stopifnot(is.numeric(data))
  if(is.data.frame(data)){
      # TODO: All NA column has type logical 
      #stopifnot(all(lapply(data, is.numeric)))
  }

  # Iterate till convergence
  curLik <- -Inf
  iter <- 0
  cat("Iterations:")
  while( TRUE && iter < iter.max){
      iter <- iter + 1
      cat(".")
      
      # Calculate likelihood of data
      oldLik <- curLik
      curLik <- sum(log(likelihood(data, dfg)))
      if( !(curLik-oldLik > threshold))
          break
      
      # Obtain expectation counts
      expCounts <- facExpectations(data, dfg)

      # Update potentials
      newPotentials <- lapply( seq_along(expCounts), FUN=function(i){
        if(is.null(optim) || optim[i] == 'row') # Default optimization
          return( .rowOptimize(expCounts[[i]]) )
        if(optim[i] == 'noopt')
          return( matrix(0,0,0) )
        if(optim[i] == 'norm')
          return( .normOptimize(expCounts[[i]]))
        if(optim[i] == 'beta')
          return( .betaOptimize(expCounts[[i]]))
        if(!is.null(optimFun) && optim[i] %in% names(optimFun)){
          if(!is.function(optimFun[[optim[i]]]))
            stop("optimFun should contain functions that take expectation counts and return potentials")
          return( optimFun[[optim[i]]](expCounts[[i]]) )
        }
        stop(paste0("No matching optimization function found for: ", optim[i]))
      })

      dfg$dfgmodule$resetPotentials( newPotentials);
  }
  cat("\n")
  
  # TODO summarize
  
}