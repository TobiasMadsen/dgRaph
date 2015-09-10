# Build RDFG module corresponding to dfg object
.build <- function(dfg){
  if(! is.dfg(dfg))
    stop("dfg should be a dfg object")
  
  new("RDFG", 
      dfg$varDim, 
      lapply(dfg$facPot,as.matrix), 
      dfg$facNbs, 
      dfg$potMap-1)
}