# Copy dfg object

.copy <- function(dfg){
  if(! is.dfg(dfg))
    stop("Object must be a dfg object")
  
  dfgCopy <- dfg(varDim = dfg$varDim, 
                 facPot = dfg$facPot,
                 facNbs = dfg$facNbs,
                 potMap = dfg$potMap, 
                 varNames = dfg$varNames, 
                 facNames = dfg$facNames)
}