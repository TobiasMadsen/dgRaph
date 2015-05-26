#' dfg
#' Construct a probalistic graphical model
#' 
#' @param  varDim  Vector of integers containing variable dimensions
#' @param  facPot  List of matrices containing potentials
#' @param  facNbs  List of vector containing describing neighbouring variables for each factor
#' @param  potMap  A vector of integers mapping with corresponding potential for each factor
#' 
#' @examples
#' varDim <- c(2L,3L)
#' facPot <- list(matrix(c(1:6),2,3))
#' facNbs <- list( c(0L,1L))
#' mydfg <- dfg(varDim, facPot, facNbs)
dfg <- function(varDim,
                facPot,
                facNbs,
                potMap = 1:length(facPot),
                varNames = seq_along(varDim),
                facNames = (length(varDim)+seq_along(potMap)),
                optim    = NULL){
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
                 dfgmodule=new("RDFG", varDim, facPot, facNbs, potMap-1),
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

plot.dfg <- function(x){
  if(require(igraph)){
    graph <- graph.empty(directed=F) + vertices(x$facNames, color="red", shape="rectangle", size2=18, size=12*nchar(x$facNames))
    graph <- graph + vertices(x$varNames, color="blue", shape="crectangle", size2=18, size=12*nchar(x$varNames))
    
    lapply(seq_along(facNbs),FUN=function(i){
      graph <<- graph + edges( c(i, length(x$facNbs)+x$facNbs[[i]][1]), mode="mutual")
      if(length(x$facNbs[[i]]) == 2)
        graph <<- graph + edges( c(i, length(x$facNbs)+x$facNbs[[i]][2]) )
    })
    
    plot(graph, layout = layout.reingold.tilford,
         main = "DFG",
         vertex.frame.color = "white",
         vertex.label.color = "white",
         vertex.label.family = "sans",
         edge.width = 3,
         edge.color = "black")
  } else{
    warning("igraph not installed: Check manually that your graph is acyclic")
  }
}
