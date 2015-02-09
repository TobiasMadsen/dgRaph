#' ImportanceSampling for tail estimation
#' @param x         points to evaluate tail probabilities in
#' @param n         number of samples
#' @param q         find quantiles instead of probabilities see details
#' @param alpha     Tuning parameter, alpha=0 is naive sampling the higher alpha the more extreme observations
#' @param pgm       Probabilistic graphical model specifying null model
#' @param facPotFg  Foreground model
#' @return A dataframe with columns, x, tail estimate and confidence intervals
tailIS <- function(x=NULL, n, q=NULL, alpha=0.5, dfg, facPotFg){
  stopifnot( (!is.null(x)) | (!is.null(q)) )
  if( !is.null(x) )
    stopifnot(is.numeric(x))
  if( !is.null(q) )
    stopifnot(is.numeric(q))
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(alpha))
  stopifnot(is.dfg(dfg))
  
  #Check if compatible dimensions
  stopifnot( is.list(facPotFg), all(sapply(facPotFg, is.matrix)))
  stopifnot( length(facPotFg) == length(dfg$facPot))
  stopifnot( all(sapply(seq_along(dfg$facPot), FUN=function(i){all(dim(dfg$facPot[[i]])==dim(facPotFg[[i]]))})) )
  
  
  if(!is.null(x)){
    #Find expected score for each alpha and find nearest mean
    expScore <- sapply( alpha, FUN=function(a){ dfg$dfgmodule$calculateExpectedScoreIS(a, facPotFg)})
    x_sample_assign <- sapply( x, FUN=function(x){ which.min(abs(x-expScore)) })#the estimate at each x will be calculated from 
    
    ret <- data.frame(x = numeric(0), p = numeric(0), low = numeric(0), high = numeric(0), p_lower = numeric(0), alpha = numeric(0))
    
    for(i in seq_along(alpha)){
      x_alpha <- x[ which(x_sample_assign == i) ]
      if( length(x_alpha) == 0) #check if any samples are assigned
        next
      a <- alpha[i]
      
      #Make samples
      samples <- dfg$dfgmodule$makeImportanceSamples(n, a, facPotFg)
      
      #Tail probabilities P(S > x)
      cdf_upper_tail <- sapply(x_alpha, function(s){
        if(a < 0){ #estimate lower tail
          obs <- samples$weights*(samples$scores < s)
          error <- qnorm(0.975)*sd(obs)/sqrt(n)
          meanobs <- mean(obs)
          c( 1-(meanobs+error), 1-meanobs, 1-(meanobs-error), meanobs )
        } else{ #estimate upper tail
          obs <- samples$weights*(samples$scores > s)
          error <- qnorm(0.975)*sd(obs)/sqrt(n)
          meanobs <- mean(obs)
          return(c(meanobs-error, meanobs, meanobs+error, 1-meanobs))
        }
      })
      
      ret <- rbind( ret, 
                    data.frame(x=x_alpha, p=cdf_upper_tail[2,], low=cdf_upper_tail[1,], high=cdf_upper_tail[3,], 
                               p_lower=cdf_upper_tail[4,], alpha=a) )
    }
   return(ret) 
  } else if(!is.null(q)){ #Find quantiles
    #Make samples
    samples <- dfg$dfgmodule$makeImportanceSamples(n, alpha, facPotFg)
    
    #Compute quantiles
    samples <- samples[with(samples, order(-scores)), ]
    cs_weights <- cumsum(samples$weights)
    quant <- sapply(q, function(q){
      if(tail(cs_weights,n=1) < q*n)
        return(NA)
      samples$scores[ which.max(cs_weights>=q*n) ]
    })
    return( data.frame(x=quant, p=q))
  }
  
}

#####################################################################
# Helper functions for saddlepoint approximations
#####################################################################

#Convert to sets of factor potentials to fun_a see "sumProduct.pdf"
.facPotToFunA <- function(facPot1, facPot2, theta){
  fun_a <- lapply(seq_along(facPot1), function(i){
    facPot1[[i]]*exp(theta*(log(facPot2[[i]])-log(facPot1[[i]])))
  })
  return(rapply( fun_a, f=function(x) ifelse(is.nan(x),0,x), how="replace" ))
}

#Convert to sets of factor potentials to fun_b see "sumProduct.pdf"
.facPotToFunB <- function(facPot1, facPot2){
  fun_b <- lapply(seq_along(facPot1), function(i){
    log(facPot2[[i]])-log(facPot1[[i]])
  })
  return(rapply( fun_b, f=function(x) ifelse(is.nan(x),0,x), how="replace" ))
}

#' Saddlepoint approximation for tail estimation
#' @param x         points to evaluate tail probabilities in
#' @param pgm       Probabilistic graphical model specifying null model
#' @param facPotFg  Foreground model
#' @return A dataframe with columns, x, tail estimate and confidence intervals
tailSaddle <- function(x, dfg, facPotFg){
  stopifnot(is.numeric(x))
  stopifnot(is.dfg(dfg))
  
  #Check if compatible dimensions
  stopifnot( is.list(facPotFg), all(sapply(facPotFg, is.matrix)))
  stopifnot( length(facPotFg) == length(dfg$facPot))
  stopifnot( all(sapply(seq_along(dfg$facPot), FUN=function(i){all(dim(dfg$facPot[[i]])==dim(facPotFg[[i]]))})) )
  
  #For numerical derivative
  require(numDeriv, quietly = TRUE)
  
  #Setup data structures
  cdf_upper_tail <- rep(NA, length(x))
  
  for(i in seq_along(x)){
    t <- x[i]
    
    #Solve d/dt k(theta) = t
    #where k(theta) is the cumulant transform
    theta <- tryCatch(
    {
      uniroot(function(z)
      {
        res <- PGMExpectCpp(dfg$varDim, 
                            .facPotToFunA(dfg$facPot, facPotFg, z),
                            .facPotToFunB(dfg$facPot, facPotFg),
                            dfg$facNbs)
        return(res[2]/res[1]-t)
      }, c(-5,5))$root
    },
    error = function(e){
      NULL
    }
    )
    
    #If theta not determinable(t is not within the range of the score)
    if(is.null(theta))
      next
    
    #Numerically find d/dt^2 k(t)
    kud2 <- grad(function(x){
      res <- PGMExpectCpp(dfg$varDim, .facPotToFunA(dfg$facPot, facPotFg, x), .facPotToFunB(dfg$facPot, facPotFg),dfg$facNbs)
      return(res[2]/res[1])
    }, theta)
    
    
    #Find MGF in theta
    phi <- PGMExpectCpp(dfg$varDim, .facPotToFunA(dfg$facPot, facPotFg, theta), .facPotToFunB(dfg$facPot, facPotFg), dfg$facNbs)[1]
    
    #Calculate 
    la <- theta*sqrt(kud2)
    cdf_upper_tail[i] <- phi*exp(-theta*t)*exp(la**2/2)*(1-pnorm(la))
  }
  
  data.frame(x=x, p=cdf_upper_tail)
}