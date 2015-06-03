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
  
  # Info from potential updates
  strPotential <- list()

  # Iterate till convergence
  likVec <- rep(0, iter.max)
  curLik <- -Inf
  iter <- 0
  cat("Iterations:")
  while( TRUE && iter < iter.max){
      iter <- iter + 1
      cat(".")
      
      # Calculate likelihood of data
      oldLik <- curLik
      curLik <- sum(likelihood(data, dfg, log = T))
      likVec[iter] <- curLik
      lastIteration <- !(curLik-oldLik > threshold)
      
      # Obtain expectation counts
      expCounts <- facExpectations(data, dfg)
      
      # Update potentials
      newPotentials <- lapply( seq_along(expCounts), FUN=function(i){
        if(is.null(optim) || optim[i] == 'row'){ # Default optimization
          opt <- .rowOptimize(expCounts[[i]])
          if(lastIteration)
            strPotential[[i]] <<- opt[['str']]
          return(opt[['pot']])
        }
        if(optim[i] == 'noopt'){
          if(lastIteration)
            strPotential[[i]] <<- "No optimization performed\n"
          return( matrix(0,0,0) )
        }
        if(optim[i] == 'norm'){
          opt <- .normOptimize(expCounts[[i]])
          if(lastIteration)
            strPotential[[i]] <<- opt[['str']]
          return(opt[['pot']])
        }
        if(optim[i] == 'beta'){
          opt <- .betaOptimize(expCounts[[i]])
          if(lastIteration)
            strPotential[[i]] <<- opt[['str']]
          return(opt[['pot']])
        }
        if(!is.null(optimFun) && optim[i] %in% names(optimFun)){
          if(!is.function(optimFun[[optim[i]]]))
            stop("optimFun should contain functions that take expectation counts and return potentials")
          opt <- optimFun[[optim[i]]](expCounts[[i]])
          if(lastIteration)
            strPotential[[i]] <<- opt[['str']]
          return(opt[['pot']])
        }
        stop(paste0("No matching optimization function found for: ", optim[i]))
      })

      potentials(dfg) <- newPotentials
      
      if(lastIteration)
        break;
  }
  cat("\n")
  
  # Summary
  # Iterations / Convergence
  if(iter == iter.max){
    cat("EM-algorithm did not converge in", iter.max, "iterations\n")
  } else{
    cat("EM-algorithm converged after", iter, "iterations\n")
  }
  cat("Likelihood:", curLik, "\n\n")
  plot(likVec[1:iter], type = 'l', xlab = "Iteration", ylab = "likelihood", main = "EM-convergence") 
  
  # Output from optimization functions(i.e. parameters)
  for(i in seq_along(strPotential)){
    cat(i,'th potential\n',sep = '')
    cat(strPotential[[i]])
    cat('\n')
  }
  
}