fourIndependentVariables <- function(){
    varDim <- rep(4,4)
    facPot <- c(list(matrix(c(0.05,0.05,0.20,0.70),1,4)),
                list(matrix(c(0.05,0.05,0.70,0.20),1,4)),
                list(matrix(c(0.70,0.20,0.05,0.05),1,4)),
                list(matrix(c(0.05,0.70,0.20,0.05),1,4)))
    facNbs <- c(list(c(1L)),
                list(c(2L)),
                list(c(3L)),
                list(c(4L)))
    
    mydfg <- dfg(varDim, facPot, facNbs, varNames = c('x', 'y', 'z', 'w'))
}