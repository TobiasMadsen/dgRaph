#Functions for calculating Moment Generating Function
#Random variable defined as log likelihoodratio between two factorgraphs
#foreground and null, under null.

#Convert to sets of factor potentials to fun_a see "sumProduct.pdf"
facPotToFunA <- function(facPot1, facPot2, theta){
  fun_a <- lapply(seq_along(facPot1), function(i){
    facPot1[[i]]*exp(theta*(log(facPot2[[i]])-log(facPot1[[i]])))
  })
  return(rapply( fun_a, f=function(x) ifelse(is.nan(x),0,x), how="replace" ))
}

#Convert to sets of factor potentials to fun_b see "sumProduct.pdf"
facPotToFunB <- function(facPot1, facPot2){
  fun_b <- lapply(seq_along(facPot1), function(i){
    log(facPot2[[i]])-log(facPot1[[i]])
  })
  return(rapply( fun_b, f=function(x) ifelse(is.nan(x),0,x), how="replace" ))
}

#' Saddlepoint estimation of tail probabilities
#' @param  t        value at which to evaluate tail
#' @param  dfg      Probalistic graphical model specifying null model
#' @param  facPotFg Foreground model 
saddlePointTails <- function(t, dfg, facPotFg){
  stopifnot(is.dfg(dfg))
  require(numDeriv, quietly = TRUE)
  
  #Check compatible dimensions
  stopifnot( is.list(facPotFg), all(sapply(facPotFg, is.matrix)))
  stopifnot( length(facPotFg) == length(dfg$facPot))
  stopifnot( all(sapply(seq_along(facPot), FUN=function(i){all(dim(facPot[[i]])==dim(facPotFg[[i]]))})) )
  
  #Solve d/dt k(theta) = t
  #where k(theta) is the cumulant transform
  theta <- uniroot(function(x){
    res <- PGMExpectCpp(dfg$varDim, 
                        facPotToFunA(dfg$facPot, facPotFg, x),
                        facPotToFunB(dfg$facPot, facPotFg),
                        dfg$facNbs)
    return(res[2]/res[1]-t)
    }, c(-5,5))$root
  
  #Numerically find d/dt^2 k(t)
  kud2 <- grad(function(x){
    res <- PGMExpectCpp(varDim, facPotToFunA(dfg$facPot, facPotFg, x), facPotToFunB(dfg$facPot, facPotFg),facNbs)
    return(res[2]/res[1])
  }, theta)
  
  
  #Find MGF in theta
  phi <- PGMExpectCpp(varDim, facPotToFunA(dfg$facPot, facPotFg, theta), facPotToFunB(dfg$facPot, facPotFg),facNbs)[1]
  
  #Calculate 
  la <- theta*sqrt(kud2)
  return( phi*exp(-theta*t)*exp(la**2/2)*(1-pnorm(la)) )
}
