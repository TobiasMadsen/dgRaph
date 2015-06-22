twoDependentVariables <- function(){
    varDim <- rep(2,2)
    facPot <- list(matrix(c(0.7,0.3),1,2),
                   matrix(c(0.75,0.25,0.25,0.75),2,2))
    facNbs <- list(c(1L),
                   c(1L,2L))
    dfg(varDim, facPot, facNbs, varNames = c('x', 'y'))  
}