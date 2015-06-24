#' dfg
#' Construct a probalistic graphical model
#' 
#' @param  varDim  Vector of integers containing variable dimensions
#' @param  facPot  List of matrices containing potentials
#' @param  facNbs  List of vector containing describing neighbouring variables for each factor
#' @param  potMap  A vector of integers mapping with corresponding potential for each factor
#' @param  varNames  A character vector of names for each variable
#' @param  facNames  A character vector of names for each potential
#' 
#' @examples
#' varDim <- c(2L,3L)
#' facPot <- list(matrix(c(1:6),2,3))
#' facNbs <- list( c(1L,2L))
#' mydfg <- dfg(varDim, facPot, facNbs)
dfg <- function(varDim,
                facPot,
                facNbs,
                potMap = 1:length(facPot),
                varNames = seq_along(varDim),
                facNames = (length(varDim)+seq_along(potMap))){
  #Check varDim is a vector of integers
  if( ! is.vector(varDim, mode="numeric") |  ! all(varDim %% 1 == 0) )
    stop("varDim must be a vector of integers")
  
  #Check facPot is a list of matrices
  if( ! is.list(facPot) | ! all(sapply(facPot, is.matrix)))
    stop("facPot must be a list of matrices")
    
  #Check facNbs is a list of integer vectors
  if( ! is.list(facNbs) | ! all(sapply(facNbs, is.vector)) )
    stop("facNbs must be a list of vectors")
  if( ! all(sapply(facNbs, function(x) all(x %% 1 ==0))) )
    stop("facNbs must be a list of integer vectors")
  
  #Check all facNbs are in 1:seq_along(varDim)
  if( max(unlist(facNbs)) > length(varDim) || min(unlist(facNbs)) <= 0 )
    stop("facNbs should be integers 1,2,...,length(varDim)")
  
  #Check that potential dimensions matches var dimensions
  if( ! length(facNbs) == length(potMap) )
    stop("facNbs must have same length as potMap")
  facNbsCheck <- sapply(seq_along(facNbs), function(i){
    if(length(facNbs[[i]]) == 1){
      return( 1 == nrow(facPot[[ potMap[i] ]]) & varDim[ facNbs[[i]][1] ] == ncol(facPot[[ potMap[i] ]]))
    } else if(length(facNbs[[i]]) == 2){
      return(varDim[ facNbs[[i]][1] ] == nrow(facPot[[ potMap[i] ]]) & varDim[ facNbs[[i]][2] ] == ncol(facPot[[ potMap[i] ]]) )
    } else{ #Too many neighbors
      return(FALSE)
    }
  })
  if(! all(facNbsCheck)){
    fac <- which.min(facNbsCheck)
    if(length(facNbs[[fac]]) == 1)
      stop("Potential ", fac, " has wrong dimensions.\n",
           "Neighboring variable has dimensions: 1 ",varDim[ facNbs[[fac]][1] ], '\n',
           "Potential has dimensions: ", nrow(facPot[[ potMap[fac] ]]), " ", ncol(facPot[[ potMap[fac] ]]), '\n')
    if(length(facNbs[[fac]]) == 2)
      stop("Potential ", fac, " has wrong dimensions.\n",
           "Neighboring variables has dimensions: ", varDim[ facNbs[[fac]][1]], " ", varDim[ facNbs[[fac]][2] ], '\n',
           "Potential has dimensions: ", nrow(facPot[[ potMap[fac] ]]), " ", ncol(facPot[[ potMap[fac] ]]), '\n')
    
    stop("Too many or too few neighbors at factor ", fac, ". Factors should have 1 or 2 neighbors\n")
  }
  
  # Check acyclic
  if(! checkAcyclic(facNbs))
    stop("Graph contains a cycle")

  structure(list(varDim=varDim,
                 facPot=facPot,
                 facNbs=facNbs,
                 potMap=potMap,
                 varNames=varNames,
                 facNames=facNames),
            class = "dfg")
}

#' is.dfg
#' 
#' Check if object is a dfg
#' 
#' @param x object to be tested
is.dfg <- function(x) inherits(x, "dfg")

#' plot.dfg
#' 
#' Plot the graphical structure of a dfg object
#' 
#' @param x   dfg object to plot
#' @param ... Arguments to be passed to plotting method. Typically \code{layout}, for connected graphs \code{layout = layout.reingold.tilford} gives good results.
#' @examples
#' library(dgRaph)
#' varDim <- rep(4,5)
#' facPot <- list(multinomialPotential(c(1,4)),
#'                multinomialPotential(c(4,4)))
#' facNbs <- list(c(1L),
#'                c(1L,2L),
#'                c(1L,3L),
#'                c(3L,4L),
#'                c(3L,5L))
#' potMap <- c(1,2,2,2,2)
#' facNames <- c("P",rep("I",4))
#' varNames <- c("Do","Re","Mi","Fa","Sol")
#' mydfg <- dfg(varDim, facPot, facNbs, potMap, varNames, facNames)
#' plot(mydfg)
#' plot(mydfg, layout = layout.reingold.tilford)
plot.dfg <- function(x, ...){
  if(require(igraph)){
    graph <- graph.empty(directed=F) + vertices(x$facNames, color="red", shape="rectangle", size2=18, size=12*nchar(x$facNames))
    graph <- graph + vertices(x$varNames, color="blue", shape="crectangle", size2=18, size=12*nchar(x$varNames))
    
    lapply(seq_along(x$facNbs),FUN=function(i){
      graph <<- graph + edges( c(i, length(x$facNbs)+x$facNbs[[i]][1]), mode="mutual")
      if(length(x$facNbs[[i]]) == 2)
        graph <<- graph + edges( c(i, length(x$facNbs)+x$facNbs[[i]][2]) )
    })
    
    plot(graph,
         main = "DFG",
         vertex.frame.color = "white",
         vertex.label.color = "white",
         vertex.label.family = "sans",
         edge.width = 3,
         edge.color = "black", ...)
  }
}

# Check if two dfgs has similar structure
# Only need to compare variable dimensions and neighbor structure
.compareDfgs <- function(dfg1, dfg2){
  if(! is.dfg(dfg1))
    stop("dfg1 must be a dfg object")
  if(! is.dfg(dfg2))
    stop("dfg2 must be a dfg object")
  
  if( ! all(dfg1$varDim == dfg2$varDim))
    stop("The two factor graphs must have same structure: varDim not identical")
  
  if( ! length(dfg1$facNbs) == length(dfg2$facNbs))
    stop("The two factor graphs must have the same number of neighbors")
  
  if(! all(dfg1$facNbs %in% dfg2$facNbs))
    stop("The two factor graphs must have the same neighbor structure. For now this also includes the order, that is c(1,2) is not the same as c(2,1)")
}

# Remap potentials dfg
# Generate a dfg with a new richer mapping of potentials
# the unique elements in potMap must be sorted
.remapPotMapDfg <- function(dfg, potMap){
  if(! is.dfg(dfg))
    stop("dfg must be a dfg object")
  
  if(length(dfg$potMap) != length(potMap))
    stop("length of potMap must be equal to the number of potentials in dfg")
  
  # Make new potentials
  dfg$facPot <- dfg$facPot[ dfg$potMap[ match(unique(potMap), potMap) ] ]
  
  # Set new potMap
  dfg$potMap <- potMap
    
  dfg # return modified object
}

# Change the order of factors
.remapFacNbsDfg <- function(dfg, facMap){
  # Check is a map
  if( ! all( facMap %in% seq_along(facMap)))
    stop("facMap must be a map")
  
  dfg$potMap[facMap] <- dfg$potMap
  dfg$facNbs[facMap] <- dfg$facNbs
  
  dfg # return modified object
}
