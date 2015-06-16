# Build RDFG module corresponding to dfg object
.build <- function(dfg){
  if(! is.dfg(dfg))
    stop("dfg should be a dfg object")
  
  new("RDFG", dfg$varDim, dfg$facPot, dfg$facNbs, dfg$potMap-1)
}