#' EM training
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. 
#'                  NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @param threshold Stop training when difference in likelihood between two iterations is below threshold
#' @param iter.max  Maximal number of iterations of the EM-algorithm.
#' @param dataList  provides support for partially observed data. See online documentation.
#' @param verbose  A TRUE/FALSE statement enabling console output of information about the training process
#' @return A discrete factor graph object with updated potentials
#' @export
train <- function(data, dfg, threshold = 1e-9, iter.max = 200, dataList = list(), verbose = FALSE){
  #Rprof(filename = "Rprof_dgRaph_train.out", append = TRUE, line.profiling = T)
  .checkInputData(dfg, data, dataList)
  module <- .build(dfg)
  
  # Output
  if (verbose) cat("Training...\n")
  
  # Info from potential updates
  strPotential <- list()
  
  # Iterate till convergence
  likVec <- rep(0, iter.max)
  curLik <- -Inf
  iter <- 0
  if (verbose) cat("Iterations:")
  while( iter < iter.max){
    iter <- iter + 1
    if (verbose) cat(".")
    
    # Calculate likelihood of data
    oldLik <- curLik
    curLik <- sum(.likelihood(data, dfg, log = T, dataList = dataList, module = module))
    likVec[iter] <- curLik
    lastIteration <- !(curLik-oldLik > threshold)
    
    # Obtain expectation counts
    expCounts <- .facExpectations(data, dfg, dataList = dataList, module = module)
    
    # Update potentials
    for(i in seq_along(expCounts)){
      potentials(dfg)[[i]] <- update(potentials(dfg)[[i]], expCounts[[i]])
    }
    
    # Expose to C++ side
    newPotentials <- lapply(potentials(dfg), as.matrix)
    module$resetPotentials( newPotentials )
    
    if(lastIteration)
      break;
  }
  if (verbose) cat("\n")
  
  # Summary
  # Iterations / Convergence
  if (verbose) if(iter == iter.max){
    cat("EM-algorithm did not converge in", iter.max, "iterations\n")
  } else{
    cat("EM-algorithm converged after", iter, "iterations\n")
  }
  if (verbose) cat("Likelihood:", curLik, "\n\n")
  if (verbose) plot(likVec[1:iter], type = 'l', xlab = "Iteration", ylab = "likelihood", main = "EM-convergence") 
  
  # Output from optimization functions(i.e. parameters)
  if (verbose) for(i in seq_along(strPotential)){
    cat(i,'th potential\n',sep = '')
    cat(strPotential[[i]])
    cat('\n')
  }
  
  invisible(dfg)
}