twoDependentVariables <- function(){
    varDim <- rep(2,2)
    facPot <- c(list(matrix(c(0.7,0.3),1,2)),
                list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
    facNbs <- c(list(c(1L)),
                list(c(1L,2L)))
    dfg(varDim, facPot, facNbs, varNames = c('x', 'y'))  
}