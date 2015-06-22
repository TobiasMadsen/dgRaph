fourIndependentVariables <- function(){
    varDim <- rep(4,4)
    facPot <- list(matrix(c(0.05,0.05,0.20,0.70),1,4),
                   matrix(c(0.05,0.05,0.70,0.20),1,4),
                   matrix(c(0.70,0.20,0.05,0.05),1,4),
                   matrix(c(0.05,0.70,0.20,0.05),1,4))
    facNbs <- list(c(1L),
                   c(2L),
                   c(3L),
                   c(4L))
    
    mydfg <- dfg(varDim, facPot, facNbs, varNames = c('x', 'y', 'z', 'w'))
}