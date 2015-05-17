#' EM training
#' @param data      dataframe or matrix with observed data. The mapping between columns and rows will be performed automatically. NA will be interpreted as a missing variable
#' @param dfg       discrete factor graph object
#' @return A dataframe with a column for each variable and the most probable state for each observation
train <- function(data, dfg){
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
  while( TRUE ){
      iter <- iter + 1
      cat(sprintf("EM algorithm iteration %i current likelihood: %f\n", iter, curLik))
      
      # Calculate likelihood of data
      oldLik <- curLik
      curLik <- sum(log(likelihood(data, dfg)))
      if( !(curLik-oldLik > 0.0001) )
          break
      
      # Obtain expectation counts
      expCounts <- facExpectations(data, dfg)

      # TODO: Map

      # Update potentials
      newPotentials <- lapply( expCounts, FUN=function(x){
          sweep( x, 1, STATS = rowSums(x), FUN = '/')
      })

      dfg$dfgmodule$resetFactorPotentials( newPotentials);
  }
  
  # Update facPot
  
}
