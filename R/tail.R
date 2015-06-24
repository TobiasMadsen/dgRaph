#' ImportanceSampling for tail estimation
#' @param x         points to evaluate tail probabilities in
#' @param n         number of samples
#' @param alpha     Tuning parameter, alpha=0 is naive sampling the higher alpha the more extreme observations
#' @param dfg1      dfg object specifying null model
#' @param dfg2      dfg object specifying foreground model
#' @param observed  A boolean vector indicating which variables are observed
#' @return A dataframe with columns, x, tail estimate and confidence intervals
tailIS <- function(x=NULL, n = 1000, alpha=0.5, dfg1, dfg2, observed = NULL){
  if( is.null(x))
    stop("Specify scores x ")
  if( !is.numeric(x) )
    stop("x must be a numeric vector")
  if( !is.numeric(n) | length(n) != 1) 
    stop("n must be a single integer")
  if( !is.numeric(alpha))
    stop("alpha must be a numeric vector")
  if( !is.dfg(dfg1))
    stop("dfg1 must be a dfg object")
  if( !is.dfg(dfg2))
    stop("dfg1 must be a dfg object")
  if( !is.null(observed)){
    if( !is.logical(observed) | !is.vector(observed))
      stop("observed must be a logical vector")
    if( length(observed) != length(dfg$varDim))
      stop("observed must have same length equal to number of variables")
  }
  
  # Restructure dfgs
  .compareDfgs(dfg1, dfg2)
  
  # If not same nb structure remap
  if( ! all( sapply(seq_along(dfg1$facNbs), function(i){all(dfg1$facNbs[[i]] == dfg2$facNbs[[i]])})))
    dfg2 <- .remapFacNbsDfg(dfg2, match(dfg2$facNbs, dfg1$facNbs))
  
  # If not same potential map remap
  if( ! length(dfg1$potMap) == length(dfg2$potMap) |
        ! all(dfg1$potMap == dfg2$potMap) ){
    cm <- .commonMap(dfg1$potMap, dfg2$potMap)
    dfg1 <- .remapPotMapDfg(dfg1, cm)
    dfg2 <- .remapPotMapDfg(dfg2, cm)
  }
  
  # Data frame to return
  ret <- data.frame(x = numeric(0), p = numeric(0), low = numeric(0), high = numeric(0), p_lower = numeric(0), alpha = numeric(0))

  # Build bg module
  facPotBg <- potentials(dfg1)
  facPotFg <- potentials(dfg2)
  moduleBg <- .build(dfg1)
  
  for(i in seq_along(alpha)){
    if( length(x) == 0) #check if any samples are assigned
      next
    a <- alpha[i]
    
    #Make samples
    dfgIS <- .copy(dfg1)
    facPotIS <- lapply(seq_along(facPotBg), FUN=function(i){
      m <- facPotBg[[i]]**(1-a)*facPotFg[[i]]**a
      
      # Score Matrix
      m[facPotBg[[i]] == 0 | facPotFg[[i]] == 0] <- 0
      m
      })
    potentials(dfgIS) <- facPotIS
    moduleIS <- .build(dfgIS)
    samples <- .simulate.dfg(dfgIS, n, module = moduleIS)
    if(!is.null(observed)){
      # Change unobserved columns to NA
      samples[,which(observed)] <- NA
    }
    
    # Calculate weights and scores
    dfgFg <- .copy(dfg1)
    potentials(dfgFg) <- facPotFg
    moduleFg <- .build(dfgFg)
    ncIS <- .likelihood(matrix(NA, 1, length(dfg1$varDim)), dfgIS, log = T, module = moduleIS) # dfgIS is not normalized
    weights <- exp(.likelihood(samples, dfg1, log = T, module = moduleBg) - .likelihood(samples, dfgIS, log = T, module = moduleIS) + ncIS)
    scores <- .likelihood(samples, dfgFg, log = T, module = moduleFg) - .likelihood(samples, dfg1, log = T, module = moduleBg)
    
    #Tail probabilities P(S > x)
    cdf_upper_tail <- sapply(x, function(s){
      if(a < 0){ #estimate lower tail
        obs <- weights*(scores < s)
        error <- qnorm(0.975)*sd(obs)/sqrt(n)
        meanobs <- mean(obs)
        c( max(0,1-(meanobs+error)), 1-meanobs, min(1,1-(meanobs-error)), meanobs )
      } else{ #estimate upper tail
        obs <- weights*(scores > s)
        error <- qnorm(0.975)*sd(obs)/sqrt(n)
        meanobs <- mean(obs)
        return(c(max(0,meanobs-error), meanobs, min(1,meanobs+error), 1-meanobs))
      }
    })
    
    ret <- rbind( ret, 
                  data.frame(x=x, p=cdf_upper_tail[2,], low=cdf_upper_tail[1,], high=cdf_upper_tail[3,], 
                             p_lower=cdf_upper_tail[4,], alpha=a) )
  }
  
  ret
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
#' @param dfg1      dfg object specifying null model
#' @param dfg2      dfg object specifying foreground model
#' @return A dataframe with columns, x, tail estimate and confidence intervals
tailSaddle <- function(x, dfg1, dfg2){
  stopifnot(is.numeric(x))
  
  #Check if compatible dimensions
  # Restructure dfgs
  .compareDfgs(dfg1, dfg2)
  
  # If not same nb structure remap
  if( ! all( sapply(seq_along(dfg1$facNbs), function(i){all(dfg1$facNbs[[i]] == dfg2$facNbs[[i]])})))
    dfg2 <- .remapFacNbsDfg(dfg2, match(dfg2$facNbs, dfg1$facNbs))
  
  # If not same potential map remap
  if( ! length(dfg1$potMap) == length(dfg2$potMap) |
        ! all(dfg1$potMap == dfg2$potMap) ){
    cm <- .commonMap(dfg1$potMap, dfg2$potMap)
    dfg1 <- .remapPotMapDfg(dfg1, cm)
    dfg2 <- .remapPotMapDfg(dfg2, cm)
  }
  
  facPotBg <- potentials(dfg1)
  facPotFg <- potentials(dfg2)
  
  #Setup data structures
  cdf_upper_tail <- rep(NA, length(x))
  
  for(i in seq_along(x)){
    t <- x[i]
    
    moduleSaddle <- .build(dfg1)
    
    #Solve d/dt k(theta) = t
    #where k(theta) is the cumulant transform
    theta <- tryCatch(
    {
      uniroot(function(z)
      {
        moduleSaddle$resetPotentials( .facPotToFunA(facPotBg, facPotFg, z) )
        facScore <- .facPotToFunB(facPotBg, facPotFg)
        res <- .expect.dfg(dfg1, facScore, module = moduleSaddle)
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
      moduleSaddle$resetPotentials( .facPotToFunA(facPotBg, facPotFg, x) )
      facScore <- .facPotToFunB(facPotBg, facPotFg)
      res <- .expect.dfg(dfg1, facScore, module = moduleSaddle)
      return(res[2]/res[1])
    }, theta)
    
    
    #Find MGF in theta
    moduleSaddle$resetPotentials( .facPotToFunA(facPotBg, facPotFg, theta) )
    phi <- .likelihood(data = matrix(NA, 1, length(dfg1$varDim)), dfg = dfg1, module = moduleSaddle)[1]
    
    #Calculate 
    la <- theta*sqrt(kud2)
    cdf_upper_tail[i] <- phi*exp(-theta*t)*exp(la**2/2)*(1-pnorm(la))
  }
  
  data.frame(x=x, p=cdf_upper_tail)
}
