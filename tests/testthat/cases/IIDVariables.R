IIDVariables <- function(N){
  # IID variables
  varDim <- rep(2,N)
  facPot <- list(matrix(c(0.7,0.3),1,2))
  facNbs <- lapply(1:N, FUN=function(i){i})
  potMap <- rep(1, N)
  mydfg <- dfg(varDim = varDim, facPot = facPot, facNbs = facNbs, potMap = potMap)
}