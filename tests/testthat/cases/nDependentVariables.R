nDependentVariables <- function(N){
  varDim <- rep(2,N)
  facPot <- list(matrix(c(0.5,0.5),1,2),
                 matrix(c(0.75,0.25,0.25,0.75),2,2))
  facNbs <- c(list(c(1L)),
              lapply(2L:N, FUN=function(x){c(x-1,x)}))
  potMap <- c(1, rep(2, N-1))
  dfg(varDim, facPot, facNbs, potMap)
}