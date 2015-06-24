nNeighbors <- function(N){
  varDim <- rep(100, N+1)
  facPot <- list(matrix(0.01, 100, 100))
  facNbs <- lapply(1:N, function(i){c(1,i+1)})
  potMap <- rep(1, N)
  dfg(varDim, facPot, facNbs, potMap)
}