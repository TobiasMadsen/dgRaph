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
  while( TRUE && iter < iter.max){
      iter <- iter + 1
      cat(sprintf("EM algorithm iteration %i current likelihood: %f\n", iter, curLik))
      
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
  
}

#################################################
# Built-in optimization functions
# Take expectation counts and return a potential
#################################################

.rowOptimize <- function(expCounts){
  sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')
}

.normOptimize <- function(expCounts){
  t(apply(expCounts, 1, FUN=function(x){
    xs <- sum(x*seq_along(x) )
    xm <- xs/sum(x)+1/2
    xv <- (sum(x*seq_along(x)**2)-xs**2/sum(x))/sum(x)
    diff(pnorm( c(-Inf, 2:length(x), Inf), xm, sqrt(xv)))
  }))
}

.betaOptimize <- function(expCounts){
  # Method of moments
  t(apply(expCounts, 1, FUN=function(x){
    val <- (seq_along(x)-0.5)/length(x)
    xs  <- sum(x*val )
    xm  <- xs/sum(x)
    xv  <- (sum(x*val**2)-xs**2/sum(x) )/sum(x)
    a   <- xm*(xm*(1-xm)/xv-1)
    b   <- a*(1/xm-1)
    diff(pbeta( 0:length(x)/length(x), a, b))
  }))
}