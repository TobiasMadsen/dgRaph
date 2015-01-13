#' PGM
#' Construct a probalistic graphical model
#' 
#' @param  varDim  Vector of integers containing variable dimensions
#' @param  facPot  List of matrices containing factor potentials for null model
#' @param  facNbs  List of vector containing describing neighbouring variables for each factor 
#' 
#' @examples
#' varDim <- c(2L,3L)
#' facPot <- list(matrix(c(1:6),2,3))
#' facNbs <- list( c(0L,1L))
#' my_pgm <- PGM(varDim, facPot, facNbs)
PGM <- function(varDim, facPot, facNbs){
  #Check varDim is a vector of integers
  stopifnot( is.vector(varDim, mode="numeric"),all(varDim %% 1 == 0))
  
  #Check facPot is a list of matrices
  stopifnot( is.list(facPot), sapply(facPot, is.matrix))
  
  #Check facNbs is a list of integer vectors
  stopifnot( is.list(facNbs), sapply(facNbs, is.vector), all(sapply(facNbs, function(x) all(x %% 1 ==0))) )
  
  #Check that potential dimensions matches var dimensions
  stopifnot( length(facNbs) == length(facPot) )
  stopifnot( all(sapply(seq_along(facNbs), function(i){
    if(length(facNbs[[i]]) == 1){
      return(varDim[ facNbs[[i]][1] ] == ncol(facPot[[i]]))
    } else if(length(facNbs[[i]]) == 2){
      return(varDim[ facNbs[[i]][1] ] == nrow(facPot[[i]]) & varDim[ facNbs[[i]][2] ] == ncol(facPot[[i]]) )
    } else{ #Too many neighbors
      return(FALSE)
    }
  } )) )
  
  graph <- NULL
  #Create graph
  if(require(igraph)){
    graph <- graph.empty(directed=F) + vertices(length(varDim)+seq_along(facPot), color="red")
    graph <- graph + vertices(seq_along(varDim), color="blue")
    V(graph)$shape <- c(rep("square",length(facNbs)),rep("circle",length(varDim)))
    lapply(seq_along(facNbs),FUN=function(i){
      graph <<- graph + edges( c(i, length(facNbs)+facNbs[[i]][1]), mode="mutual")
      if(length(facNbs[[i]]) == 2)
        graph <<- graph + edges( c(i, length(facNbs)+facNbs[[i]][2]) )
    })
  } else{
    warning("igraph not installed: Check manually that your graph is acyclic")
  }
  
  structure(list(varDim=varDim, facPot=facPot, facNbs=facNbs, graph=graph), class = "pgm")
}

#' PGM
#' 
#' Check if object is a PGM
#' 
#' @param x object to be tested
is.pgm <- function(x) inherits(x, "pgm")

plot.pgm <- function(x){
  require(igraph)
  if(!is.null(x$graph)){
    plot(x$graph, layout = layout.reingold.tilford,
         main = "PGM",
         vertex.size = 25,
         vertex.frame.color = "white",
         vertex.label.color = "white",
         vertex.label.family = "sans",
         edge.width = 3,
         edge.color = "black")
  }
}