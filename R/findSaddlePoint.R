#################################################
# Find saddlepoint 
# Suprisingly useful helper function
#################################################

findSaddlePoint <- function(dfg1, dfg2, t){
  if( !is.dfg(dfg1))
    stop("dfg1 should be a dfg object")
  if( !is.dfg(dfg2))
    stop("dfg2 should be a dfg object")

  ## Find range of scores
  facPotBg <- lapply(potentials(dfg1), as.matrix)
  facPotFg <- lapply(potentials(dfg2), as.matrix)
  facScore <- lapply(seq_along(facPotBg), function(i){log(facPotFg[[i]]/facPotBg[[i]])})
  moduleSaddle <- .build(dfg1)
  expFacScore <- lapply(facScore, function(x){exp(x)})
  expFacScoreMinus <- lapply(facScore, function(x){exp(-x)})
  moduleSaddle$resetPotentials( expFacScore )
  mpsDat <- .mps(data = matrix(NA, 1,length(dfg1$varDim)), dfg = dfg1, module = moduleSaddle)
  maxValue <- log(.likelihood(data = mpsDat, dfg = dfg1, module = moduleSaddle))
  moduleSaddle$resetPotentials( expFacScoreMinus )
  mpsDat <- .mps(data = matrix(NA, 1,length(dfg1$varDim)), dfg = dfg1, module = moduleSaddle)
  minValue <- -log(.likelihood(data = mpsDat, dfg = dfg1, module = moduleSaddle))

  if(t > maxValue)
    stop("t is larger than the maximal attainable value")
  if(t < minValue)
    stop("t is smaller than the minimal attainable value")

  ## Do newton raphson
  theta <- 0
  facScore <- .facPotToFunB(facPotBg, facPotFg)

  iter <- 0
  cum <- rep(0,2)
  while(T){
    if(t >= maxValue | t <= minValue){
      cum <- c(NA, NA)
      break;
    }
    moduleSaddle$resetPotentials( .facPotToFunA(facPotBg, facPotFg, theta) )
    res <- unname(.expect2.dfg(dfg1, facScore, module = moduleSaddle))
    cum <- c(res[2]/res[1], (res[3]*res[1]-res[2]**2)/(res[1]**2) )
    
    theta <- theta - min(max((cum[1]-t)/ cum[2],-5),5)
    iter <- iter +1
    if(iter > 100)
      break;
    if(any(is.nan(cum)) | abs(cum[1] - t) < 0.001)
      break;
  }

  return(theta)
}
