#' Potentials
#' @param dfg     discrete factor graph object
#' @return A list of current factors
#' @examples
#' varDim <- rep(2,2)
#' facPot <- list(matrix(c(0.7,0.3),1,2),
#'                matrix(c(0.75,0.25,0.25,0.75),2,2))
#' facNbs <- c(list(c(1L)),
#'             list(c(1L,2L)))
#' mydfg <- dfg(varDim, facPot, facNbs, varNames = c('x', 'y'))  
#' 
#' potentials(mydfg)
#' @export
potentials <- function(dfg){
    stopifnot(is.dfg(dfg))

    dfg$facPot
}

#' Potentials <- 
#' @param dfg     discrete factor graph object
#' @param value   a list of new factor potentials
#' @examples
#' varDim <- rep(2,2)
#' facPot <- list(matrix(c(0.7,0.3),1,2),
#'                matrix(c(0.75,0.25,0.25,0.75),2,2))
#' facNbs <- c(list(c(1L)),
#'             list(c(1L,2L)))
#' mydfg <- dfg(varDim, facPot, facNbs, varNames = c('x', 'y'))  
#' 
#' newFacPot <- list(matrix(c(0.5,0.5),1,2),
#'                   matrix(c(0.4,0.6,0.7,0.3),2,2))
#' potentials(mydfg) <- newFacPot
#' @export
"potentials<-" <- function(dfg, value){
  # Perform input check
  if( ! is.list(value) | ! all(sapply(value, function(v){is.matrix(v) | is.potential(v)})))
    stop("Potentials must be a list of matrices")
  if(length(value) != length(dfg$facPot))
    stop("Potentials did not have correct length")
  
  sapply(seq_along(dfg$facPot), FUN=function(i){
      if( !all( dim(dfg$facPot[[i]]) == dim(value[[i]]) )){
        if( !all( dim(value[[i]]) == c(0,0) )){
          stop(i,"th potential doesn't match")
        }
      }
    })
  
  for(i in seq_along(dfg$facPot)){
    if( all(dim(value[[i]]) == c(0,0))){
      value[[i]] <- dfg$facPot[[i]]
    }
  }
  
  # Set potentials
  dfg$facPot <- value
  dfg
}

#' is.potential
#' 
#' Check if object is a potential
#' 
#' @param x object to be tested
is.potential <- function(x) inherits(x, "potential")

#' dim.potential
#' 
#' Get dimensions of a potential
#' @param x potential object
dim.potential <- function(obj) dim(obj$mat)
nrow.potential <- function(obj) nrow(obj$mat)
ncol.potential <- function(obj) ncol(obj$mat)

#' as.matrix.potential
#' 
#' Get potential matrix
#' @param obj potential obj
as.matrix.potential <- function(obj) obj$mat

#' update.potential
#' 
#' Default row normalization update given expectation counts
update.potential <- function(obj, expCounts){
  if(!all(dim(obj$matrix) == dim(data))){
    stop("dimensions should match")
  }
  
  if(is.matrix(obj)){
    obj <- sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')
  } else {
    obj$mat <- sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')
  }
  
  obj
}

summary.potential <- function(obj){
  cat("El potentialo looks like zis")
}

